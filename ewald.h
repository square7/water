#ifndef EWALD_H
#define EWALD_H
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
using namespace std;

const double a2b=1.889725989;
const int NDIM=3;

class dVec;

dVec operator*(double d, const dVec& r);
dVec operator+(const dVec& l, const dVec& r);
dVec operator-(const dVec& l, const dVec& r);
ostream& operator<<(ostream& os, const dVec& v);
double operator*(const dVec& l, const dVec& r);
dVec crossProd(const dVec& l, const dVec& r);

class dVec
{
 public:
  vector<double> vec;
  double& operator[](int i) {return vec[i];}
  const double& operator[](int i) const {return vec[i];}
  int size()const {return vec.size();}
  dVec(int dim=NDIM){for(int i = 0; i < dim; ++i){vec.push_back(0.0);}}
 dVec(const vector<double>& v):vec(v) {}
  dVec operator-() const
    {
      dVec ans;
      for(int i = 0; i < vec.size(); ++i)
	{ans.vec[i]=-vec[i];}
      return ans;
    }
};

void zeroVec(dVec& v);
void zeroVec(vector<dVec>& v);



class ParameterClass
{
 public:
  double getVol() { 
    double v1=a*crossProd(b, c);
    return v1;
  }
  void readCell(string cel) 
  {
    ifstream is(cel.c_str());
    double ele;
    vector<vector<double> > v(6);
    vector<vector<double> > transV(6,vector<double>(3));
    // first 3 lines are cells
    for(int i = 0; i < NDIM; ++i)
      {
	for(int j = 0; j < NDIM; ++j)
	  {
	    is >> ele;
	    v[i].push_back(ele);
	    transV[j][i]=ele;
	  }
      }
    a=dVec(v[0]);
    b=dVec(v[1]);
    c=dVec(v[2]);
    // second 3 lines are invCells
    for(int i = 3; i < 3+NDIM; ++i)
      {
	for(int j = 0; j < NDIM; ++j)
	  {
	    is >> ele;
	    v[i].push_back(ele);
	    transV[3+j][i-3]=ele;
	  }
      }
    is.close();
    /* ina=dVec(v[3]); */
    /* inb=dVec(v[4]); */
    /* inc=dVec(v[5]); */
    vol=getVol();
    ka=2.0*M_PI/vol*crossProd(b,c);
    kb=2.0*M_PI/vol*crossProd(c,a);
    kc=2.0*M_PI/vol*crossProd(a,b);
    // transpose abc and inv abc
    a=dVec(transV[0]);
    b=dVec(transV[1]);
    c=dVec(transV[2]);
    ina=dVec(transV[3]);
    inb=dVec(transV[4]);
    inc=dVec(transV[5]);
  }
  ParameterClass() {}
  double alpha, kcut, qO, qH, AOO, AHH, AOH, AHO, RO, RH, rho, rcut;
  dVec a,b,c, ka, kb, kc, ina, inb, inc;
  double vol;
};

class SystemClass
{
 public:
  SystemClass(){}
  ParameterClass params;
  vector<dVec> r, k, force;
  vector<double> q;
  vector<complex<double> > eikr;
  vector<char> type;
  vector<vector<dVec> > disp;
  vector<vector<double> > dist, dist2;
  SystemClass(string cel, string pos)
    {
      params.alpha=8.;
      params.kcut=8.0;
      params.rcut=5.0;
      params.qO=-2.0;
      params.qH=1.0;
      params.AOO=10.0;
      params.AHH=10.0;
      params.AOH=params.AHO=10.0;
      params.RO=.5;
      params.RH=.2;
      params.rho=.1;
      params.readCell(cel);
      readPos(pos);
      buildk();
      disp=vector<vector<dVec> >(r.size(), vector<dVec>(r.size()));
      dist=vector<vector<double> >(r.size(), vector<double>(r.size()));
      dist2=vector<vector<double> >(r.size(), vector<double>(r.size()));
      buildDisp();
    }
  void readPos(string pos)
  {
    ifstream is(pos.c_str());
    string buff;
    while(getline(is, buff))
      {
	istringstream iss(buff);
	char typ;
	vector<double> p(NDIM);
	iss >> typ >> p[0] >> p[1] >> p[2];
	type.push_back(typ);
	if(typ=='O')
	  {
	    q.push_back(params.qO);
	    r.push_back(dVec(p));
	  }
	else
	  {
	    q.push_back(params.qH);
	    r.push_back(dVec(p));
	  }
	force.push_back(dVec());
      }
  }
  bool valid(int x, int y, int z)
  {
    if (z > 0)
      return true;
    else if (z==0 && y>0)
      return true;
    else if (z==0 && y==0 && x>0)
      return true;
    else return false;

  }
  void buildk()
  {
    k.clear();
    dVec MaxkIndex;
    MaxkIndex[0]= (int) ceil(2.*params.kcut/params.ka[0]);
    MaxkIndex[1]= (int) ceil(2.*params.kcut/params.kb[1]);
    MaxkIndex[2]= (int) ceil(2.*params.kcut/params.kc[2]);
    cerr << MaxkIndex << endl;
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
	for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
          /* if (ix==0 && iy==0 && iz==0){ */
	  /*   continue; */
          /* } */
	  if(!valid(ix,iy,iz))
	    {continue;}
	  dVec currK = ix*params.ka;
	  currK = currK+iy*params.kb;
	  currK = currK+iz*params.kc;
	  if ((currK*currK<params.kcut*params.kcut)){
	    k.push_back(currK);
	  }
	}
      }
    }
    cerr << k.size() << endl;
  }
  void getDisp(int i, int j, double& d, double& d2, dVec& r12)
  {
    d=dist[i][j];
    d2=dist2[i][j];
    r12=disp[i][j];
  }
  void buildDisp() //update eikr at the same time
  {
    eikr.clear();
    for(int i = 0; i < r.size(); ++i)
      {
	disp[i][i]=dVec();
	dist[i][i]=dist2[i][i]=0;
	for( int j = i+1; j < r.size(); ++j)
	  {
	    GetDistDisp(i,j,disp[i][j], dist2[i][j]);
	    dist[i][j]=sqrt(dist2[i][j]);
	    disp[j][i]=-disp[i][j];
	    dist[j][i]=dist[i][j];
	    dist2[j][i]=dist2[i][j];
	  }
      }
    // calculate S=sum e^ikr
    for(int i = 0; i < k.size(); ++i)
      {
	complex<double> currVal = 0;
	dVec kps=k[i];
	for(int j = 0; j < r.size();++j)
	  {
	    currVal+=q[j]*exp(complex<double>(0,-1.0)*(kps*r[j]));
	  }
	eikr.push_back(currVal);
      }
  }
  void GetDistDisp(int i, int j, dVec& disp, double& dist2)
  {
    dVec dr=r[i]-r[j];
    dVec reducedR;
    reducedR[0]=params.ina*dr;
    reducedR[1]=params.inb*dr;
    reducedR[2]=params.inc*dr;
    for(int ii = 0; ii < NDIM; ++ii)
      {
	reducedR[ii]-=round(reducedR[ii]);
      }
    dr[0]=params.a*reducedR;
    dr[1]=params.b*reducedR;
    dr[2]=params.c*reducedR;
    // non-cubic cell, check nearest neighbors
    dist2 = dr*dr;
    disp = dr;
    //cerr << sqrt(dist2) << endl;
    //cerr << r[i]-r[j] << endl;
    //cerr << disp << dist2 << endl;
    /* for(int i = -1; i < 2; ++i) */
    /*   { */
    /* 	for(int j = -1; j < 2; ++j) */
    /* 	  { */
    /* 	    for(int k = -1; k < 2; ++k) */
    /* 	      { */
    /* 		dVec currR(dr); */
    /* 		currR=currR+i*params.a; */
    /* 		currR=currR+j*params.b; */
    /* 		currR=currR+k*params.c; */
    /* 		double temp = currR*currR; */
    /* 		if (temp < dist2) */
    /* 		  { */
    /* 		    dist2=temp; */
    /* 		    disp=currR; */
    /* 		  } */
    /* 	      } */
    /* 	  } */
    /*   } */
  }
};


double ewald(SystemClass& sys);
double screened(SystemClass& sys);


#endif
