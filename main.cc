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
      cerr << i << " " << sys.force[i] << endl;
      cerr << i << " " << f[i] << endl;
    }
  //cerr << sf << endl;
  //cerr << sf2 << endl;
  //cerr << e << endl;
  return e;
}

vector<double> optimizer(const vector<double>& p, const vector<dVec>& force ,double step=5e-2)
{  
  vector<double> params(p);
  vector<double> p0(p);
  vector<double> outp(p);
  SystemClass sys(string("10000.cel"), string("10000.pos"),params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10], params[11],params[12],params[13],params[14],params[15],params[16],params[17],params[18],params[19]);
  double energy0=morse(sys);
  double e0=err(sys,force);
  cerr << e0 << endl;
  vector<int> tunable;
  //tunable.push_back(7);
  //tunable.push_back(9);
  //tunable.push_back(10);
  tunable.push_back(11);
  tunable.push_back(12);
  tunable.push_back(13);
  tunable.push_back(14);
  tunable.push_back(15);
  tunable.push_back(16);
  tunable.push_back(17);
  tunable.push_back(18);
  tunable.push_back(19);
  for(int pos = 0 ;pos<tunable.size();++pos) // tunable params
    {
      int i=tunable[pos];
      params=p0;
      params[i]+=step;
      SystemClass temp(string("10000.cel"), string("10000.pos"),params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10], params[11],params[12],params[13],params[14],params[15],params[16],params[17],params[18],params[19]);
      //SystemClass temp(string("10000.cel"), string("10000.pos"),params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10]);
      double energy1=morse(temp);
      double e1=err(temp,force);
      if(e1<e0)
      	{outp[i]*=1+step;}
      else
      	{outp[i]*=1-step;}
    }
  return outp;
}

int main()
{  
  cout << setprecision(15) << fixed;
  ifstream enef("10000.ene");
  double ene;
  enef>>ene;
  SystemClass sys("10000.cel","10000.pos");
  double e0=ele(sys);
  cout << e0 << " " << ene << endl;
  vector<int> OO(100),OH(100),HH(100);
  double dr=0.05;
  for (int i = 0; i < sys.r.size(); ++i)
    {
      for(int j = 0;j < sys.r.size(); ++j)
  	{
  	  int idx=int(sys.dist[i][j]/dr);
  	  if(sys.type[i]=='H' && sys.type[j]=='H')
  	    {
  	      if(idx<100)
  		{HH[idx]++;}
  	    }
  	  else if(sys.type[i]=='O' && sys.type[j]=='O')
  	    {
  	      if(idx<100)
  		{OO[idx]++;}
  	    }
  	  else
  	    {
  	      if(idx<100)
  		{OH[idx]++;}
  	    }
  	}
    }
  for(int i = 0; i < OO.size(); ++i)
    {
      cout << OO[i] << " " << HH[i] << " " << OH[i] << endl;
    }
}

int forcecheck_main(int argc, char** argv)
{
  cout << setprecision(15) << fixed;
  SystemClass sys("10000.cel","10000.pos");
  double e0=morse(sys);
  int pos=atoi(argv[1]);
  int dir=atoi(argv[2]);
  cout << sys.force[pos][dir] << endl;
  sys.r[pos][dir]+=1e-6;
  sys.buildDisp();
  double e1=morse(sys);
  cout << (e0-e1)/1e-6 << endl;
}

int old_main(int argc, char** argv)
{
  cout << setprecision(15) << fixed;
  ifstream is("params.txt");
  vector<double> params;
  vector<string> pNames;
  double ele;
  string pName;
  for(int i = 0; i < 20; ++i)
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
  vector<double> newP=optimizer(params, force);
  for(int i = 0;i < pNames.size(); ++i)
    {
      cout << pNames[i] << " " << newP[i] << endl;
    }
}
