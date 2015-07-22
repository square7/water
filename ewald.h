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
};



class ParameterClass
{
 public:
  double getVol() { 
    double v0=a[0]*b[1]*c[2]-a[0]*b[2]*c[1]
      +a[1]*b[2]*c[0]-a[1]*b[0]*c[2]
      +a[2]*b[0]*c[1]-a[2]*b[1]*c[0];
    double v1=a*crossProd(b, c);
    cerr << "v0 " << v0 << " v1 " << v1 << endl;
    return v1;
  }
  void readCell(string cel) 
  {
    ifstream is(cel.c_str());
    double ele;
    vector<vector<double> > v(6);
    // first 3 lines are cells
    for(int i = 0; i < NDIM; ++i)
      {
	for(int j = 0; i < NDIM; ++i)
	  {
	    is >> ele;
	    v[i].push_back(ele);
	  }
      }
    a=dVec(v[0]);
    b=dVec(v[1]);
    c=dVec(v[2]);
    // second 3 lines are invCells
    for(int i = 3; i < 3+NDIM; ++i)
      {
	for(int j = 0; i < NDIM; ++i)
	  {
	    is >> ele;
	    v[i].push_back(ele);
	  }
      }
    is.close();
    ina=dVec(v[3]);
    inb=dVec(v[4]);
    inc=dVec(v[5]);
    getVol();
    ka=2.0*M_PI/vol*crossProd(b,c);
    kb=2.0*M_PI/vol*crossProd(c,a);
    kc=2.0*M_PI/vol*crossProd(a,b);
  }
  ParameterClass() {}
  double alpha, kcut, qO, qH, AOO, AHH, AOH, AHO, RO, RH, rho;
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
  vector<char> type;
  vector<vector<dVec> > disp;
  vector<vector<double> > dist, dist2;
  SystemClass(string cel, string pos)
    {
      params.alpha=1.0;
      params.kcut=4.0;
      params.qO=-2.0;
      params.qH=1.0;
      params.AOO=10.0;
      params.AHH=10.0;
      params.AOH=params.AHO=10.0;
      params.RO=1.0;
      params.RH=1.0;
      params.rho=.1;
      params.readCell(cel);
      readPos(pos);
      buildk();
      disp=vector<vector<dVec> >(r.size(), vector<dVec>(r.size()));
      dist=vector<vector<double> >(r.size(), vector<double>(r.size()));
      dist2=vector<vector<double> >(r.size(), vector<double>(r.size()));
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
      }
  }
  void buildk()
  {
    dVec MaxkIndex;
    MaxkIndex[0]= (int) ceil(2.*params.kcut/params.ka[0]);
    MaxkIndex[1]= (int) ceil(2.*params.kcut/params.kb[1]);
    MaxkIndex[2]= (int) ceil(2.*params.kcut/params.kc[2]);
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      dVec currK = ix*params.ka;
      for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
	currK = currK+iy*params.kb;
	for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
	  currK = currK+iz*params.kc;
          if (currK[0]==0 && currK[1]==0 && currK[2]==0){
          }
	  if ((currK*currK<params.kcut*params.kcut)){
	    k.push_back(currK);
	  }
	}
      }
    }
  }
  void getDisp(int i, int j, double& d, double& d2, dVec& r12)
  {
    d=dist[i][j];
    d2=dist2[i][j];
    r12=disp[i][j];
  }
  void buildDisp()
  {
    for(int i = 0; i < r.size(); ++i)
      {
	for( int j = 0; j < r.size(); ++j)
	  {
	    GetDistDisp(i,j,disp[i][j], dist2[i][j]);
	    dist[i][j]=sqrt(dist2[i][j]);
	  }
      }
  }
  void GetDistDisp(int i, int j, dVec& disp, double dist2)
  {
    dVec dr=r[i]-r[j];
    dVec reducedR;
    reducedR[0]=params.ina*dr;
    reducedR[1]=params.inb*dr;
    reducedR[2]=params.inc*dr;
    for(int i = 0; i < NDIM; ++i)
      {
	reducedR[i]-=round(reducedR[i]);
      }
    dr[0]=params.a*reducedR;
    dr[1]=params.b*reducedR;
    dr[2]=params.c*reducedR;
    // non-cubic cell, check nearest neighbors
    dist2 = dr*dr;
    disp = dr;
    for(int i = -1; i < 2; ++i)
      {
	for(int j = -1; j < 2; ++j)
	  {
	    for(int k = -1; k < 2; ++k)
	      {
		dVec currR(dr);
		currR=currR+i*params.a;
		currR=currR+j*params.b;
		currR=currR+k*params.c;
		double temp = currR*currR;
		if (temp < dist2)
		  {
		    dist2=temp;
		    disp=currR;
		  }
	      }
	  }
      }
  }
};


double ewald(SystemClass& sys);
double screened(SystemClass& sys);


#endif
