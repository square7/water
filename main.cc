#include "ewald.h"
#include <cstdlib>

double err(SystemClass& sys, const vector<dVec>& f)
{
  double e=0.0;
  double sf=0;
  double sf2=0;
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec Fdiff=sys.force[i]-f[i];
      e+=Fdiff*Fdiff;
      sf+=f[i]*f[i];
      sf2+=sys.force[i]*sys.force[i];
      //cerr << i << " " << sys.force[i] << endl;
      //cerr << i << " " << f[i] << endl;
    }
  //cerr << sf << endl;
  //cerr << sf2 << endl;
  //cerr << e << endl;
  return e;
}

vector<double> optimizer(const vector<double>& p, const vector<dVec>& force ,double step=1e-2)
{  
  vector<double> params(p);
  vector<double> p0(p);
  vector<double> outp(p);
  SystemClass sys(string("10000.cel"), string("10000.pos"),params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10]);
  double energy0=screened(sys);
  double e0=err(sys,force);
  cerr << e0 << endl;
  for(int i = 3;i<11;++i) // tunable params
    {
      params=p0;
      params[i]+=step;
      SystemClass temp(string("10000.cel"), string("10000.pos"),params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10]);
      double energy1=screened(temp);
      double e1=err(temp,force);
      if(e1<e0)
	{outp[i]+=step*(e0-e1);}
      else
	{outp[i]-=step;}
    }
  return outp;
}

int main(int argc, char** argv)
{
  cout << setprecision(15) << fixed;
  ifstream is("params.txt");
  vector<double> params;
  vector<string> pNames;
  double ele;
  string pName;
  for(int i = 0; i < 11; ++i)
    {
      is>>pName>>ele;
      //cerr << pName << ele << endl;
      params.push_back(ele);
      pNames.push_back(pName);
    }
  cout << setprecision(15) << fixed;
  vector<dVec> force;
  ifstream isf("10000.for");
  string str;
  while(getline(isf,str))
    {
      istringstream iss(str);
      char typ;
      vector<double> f(3);
      iss>>typ>>f[0]>>f[1]>>f[2];
      force.push_back(dVec(f));
    }
  //SystemClass sys("10000.cel","10000.pos");
  //double e0=screened(sys);
  //cout << sys.force[10] << endl;
  //sys.r[10][2]+=1e-7;
  //sys.buildDisp();
  //double e1=screened(sys);
  //cout << (e1-e0)/1e-7 << endl;
  vector<double> newP=optimizer(params, force);
  for(int i = 0;i < pNames.size(); ++i)
    {
      cout << pNames[i] << " " << newP[i] << endl;
    }
}
