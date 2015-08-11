#include "ewald.h"

dVec operator*(double d, const dVec& r)
{
  dVec out(r);
  for(int i = 0; i < out.size(); ++i)
    {
      out[i]*=d;
    }
  return out;
}

dVec operator+(const dVec& l, const dVec& r)
{
  dVec out;
  for(int i=0; i < l.size(); ++i)
    {
      out[i]=l[i]+r[i];
    }
  return out;
}
dVec operator-(const dVec& l, const dVec& r)
{
  dVec out;
  for(int i=0; i < l.size(); ++i)
    {
      out[i]=l[i]-r[i];
    }
  return out;
}

ostream& operator<<(ostream& os, const dVec& v)
{
  for(int i = 0; i < v.size(); ++i)
    {
      os << v[i] << " ";
    }
  return os;
}

double operator*(const dVec& l, const dVec& r)
{
  double ans=0;
  for(int i = 0; i < l.size(); ++i)
    {
      ans+=l[i]*r[i];
    }
  return ans;
}

void zeroVec(dVec& v){v=dVec();}
void zeroVec(vector<dVec>& v)
{
  for(int i = 0; i < v.size(); ++i)
    {
      v[i]=dVec();
    }
}




dVec crossProd(const dVec& l, const dVec& r)
{
  dVec ans;
  ans[0]=l[1]*r[2]-l[2]*r[1];
  ans[1]=l[2]*r[0]-l[0]*r[2];
  ans[2]=l[0]*r[1]-l[1]*r[0];
  return ans;
}

double ele(SystemClass& sys)
{
  double energy=0;
  // screened potential
  for (int k = 0; k < sys.k.size(); ++k)
    {
      dVec kps=sys.k[k];
      double k2=kps*kps;
      double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
      complex<double> S2=sys.eikr[k]*conj(sys.eikr[k])-64*6.;
      //cerr << "S2=" << S2 << " ";
      energy+=fac*S2.real();
      //cerr << "imag should be 0: " << S2.imag() << endl;
    }  
  // screened coulomb 
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec sk;
      for(int k = 0; k < sys.k.size(); ++k)
	{
	  dVec kps=sys.k[k];
	  double k2 = kps*kps;
	  double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
	  complex<double> temp = sys.eikr[k]*exp(complex<double>(0,1.0)*(kps*sys.r[i]));
	  complex<double> fac2 = (temp-conj(temp))/complex<double>(0, 1.0);
	  //cerr << fac2 <<endl;
	  sk=sk+fac*fac2.real()*kps;
	}
      sys.force[i]=sys.q[i]*sk;
    }
  return energy;
}

double morse(SystemClass& sys)
{
  double energy=0;
  // screened potential
  for (int k = 0; k < sys.k.size(); ++k)
    {
      dVec kps=sys.k[k];
      double k2=kps*kps;
      double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
      complex<double> S2=sys.eikr[k]*conj(sys.eikr[k])-64*6.;
      energy+=fac*S2.real();
    }  
  for (int i = 0; i < sys.r.size(); ++i)
    {
      for(int j = i+1;j < sys.r.size(); ++j)
	{
	  dVec r12;
	  double dist, dist2;
	  sys.getDisp(i,j,dist,dist2,r12);
	  if(dist>sys.params.rcut){continue;}
	  double fac=0;
	  double fac2=0;
	  if(sys.type[i]=='O' && sys.type[j]=='O')
	    {
	      fac=sys.params.DOO;
	      fac2=(dist-sys.params.rmoo)*sys.params.raoo;
	    }
	  else if(sys.type[i]=='H' && sys.type[j]=='H')
	    {
	      fac=sys.params.DHH;
	      fac2=(dist-sys.params.rmhh)*sys.params.rahh;
	    }
	  else
	    {
	      fac=sys.params.DOH;
	      fac2=(dist-sys.params.rmoh)*sys.params.raoh;
	    }
	  double fac3=1.-exp(-fac2);
	  energy+=fac*fac3*fac3;
	}
    }
  // Next calculate forces
  // screened coulomb 
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec sk;
      for(int k = 0; k < sys.k.size(); ++k)
	{
	  dVec kps=sys.k[k];
	  double k2 = kps*kps;
	  double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
	  complex<double> temp = sys.eikr[k]*exp(complex<double>(0,1.0)*(kps*sys.r[i]));
	  complex<double> fac2 = (temp-conj(temp))/complex<double>(0, 1.0);
	  //cerr << fac2 <<endl;
	  sk=sk+fac*fac2.real()*kps;
	}
      sys.force[i]=sys.q[i]*sk;
    }
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec sj;
      for(int j = 0; j < sys.r.size(); ++j)
	{
	  if(i==j){continue;}
	  dVec r12;
	  double dist, dist2;
	  sys.getDisp(i,j,dist,dist2,r12);
	  if(dist>sys.params.rcut){continue;}
	  double fac=0;
	  double fac2=0;
	  if(sys.type[i]=='O' && sys.type[j]=='O')
	    {
	      fac=sys.params.DOO*sys.params.raoo;
	      double temp=sys.params.raoo*(dist-sys.params.rmoo);
	      fac2=exp(-temp)*(1.-exp(-temp));
	    }
	  else if(sys.type[i]=='H' && sys.type[j]=='H')
	    {
	      fac=sys.params.DHH*sys.params.rahh;
	      double temp=sys.params.rahh*(dist-sys.params.rmhh);
	      fac2=exp(-temp)*(1.-exp(-temp));
	    }
	  else
	    {
	      fac=sys.params.DOH*sys.params.raoh;
	      double temp=sys.params.raoh*(dist-sys.params.rmoh);
	      fac2=exp(-temp)*(1.-exp(-temp));
	    }
	  double temp=-2*fac*fac2/dist;
	  sj=sj+temp*r12; // Need to check sign!!!
	}
      sys.force[i]=sys.force[i]+sj;
    }
  return energy;
}



double screened(SystemClass& sys)
{
  double energy=0;
  // screened potential
  for (int k = 0; k < sys.k.size(); ++k)
    {
      dVec kps=sys.k[k];
      double k2=kps*kps;
      double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
      complex<double> S2=sys.eikr[k]*conj(sys.eikr[k])-64*6.;
      //cerr << "S2=" << S2 << " ";
      energy+=fac*S2.real();
      //cerr << "imag should be 0: " << S2.imag() << endl;
    }  
      ///////////////////////////////////////
      //old
      // complex<double> fac2=0;
      // for(int i = 0; i < sys.r.size(); ++i)
      // 	{
      // 	  for(int j = 0; j < sys.r.size(); ++j)
      // 	    {
      // 	      dVec r12;
      // 	      double dist, dist2;
      // 	      sys.getDisp(i,j,dist,dist2,r12);
      // 	      fac2+=sys.q[i]*sys.q[j]*exp(complex<double>(0,-1.0)*(kps*r12));
      // 	    }
      // 	}
      // //cerr << "fac2=" << fac2 << endl;
      // energy+=fac*fac2.real();
      //}
      //////////////////////////////////////
      /////////////////////////////////////
  //1cerr << "kEnergy = " << energy << endl;
  // short range repulsive potential
  for (int i = 0; i < sys.r.size(); ++i)
    {
      for(int j = i+1;j < sys.r.size(); ++j)
	{
	  dVec r12;
	  double dist, dist2;
	  sys.getDisp(i,j,dist,dist2,r12);
	  if(dist>sys.params.rcut){continue;}
	  double fac=0;
	  double fac2=0;
	  if(sys.type[i]=='O' && sys.type[j]=='O')
	    {
	      fac=sys.params.AOO;
	      fac2=(sys.params.RO*2-dist)/sys.params.rho;
	    }
	  else if(sys.type[i]=='H' && sys.type[j]=='H')
	    {
	      fac=sys.params.AHH;
	      fac2=(sys.params.RH*2-dist)/sys.params.rho;
	    }
	  else
	    {
	      fac=sys.params.AOH;
	      fac2=(sys.params.RO+sys.params.RH-dist)/sys.params.rho;
	    }
	  energy+=fac*exp(fac2);
	}
    }
  // Next calculate forces
  // screened coulomb 
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec sk;
      for(int k = 0; k < sys.k.size(); ++k)
	{
	  dVec kps=sys.k[k];
	  double k2 = kps*kps;
	  double fac=2./sys.params.vol/k2*exp(-k2/2.0/sys.params.alpha);
	  complex<double> temp = sys.eikr[k]*exp(complex<double>(0,1.0)*(kps*sys.r[i]));
	  complex<double> fac2 = (temp-conj(temp))/complex<double>(0, 1.0);
	  //cerr << fac2 <<endl;
	  sk=sk+fac*fac2.real()*kps;
	}
      sys.force[i]=sys.q[i]*sk;
    }
  //////////////////////////////////////
  //old method
  // for(int i = 0; i < sys.r.size(); ++i)
  //   {
  //     dVec sk;
  //     for (int k = 0; k < sys.k.size(); ++k)
  // 	{
  // 	  dVec kps=sys.k[k];
  // 	  double k2=kps*kps;
  // 	  double fac=4.0/k2*exp(-k2/2.0/sys.params.alpha);
  // 	  double sj = 0;
  // 	  for(int j = 0; j < sys.r.size(); ++j)
  // 	    {
  // 	      dVec r12;
  // 	      double dist, dist2;
  // 	      sys.getDisp(i,j,dist,dist2,r12);
  // 	      sj+=sys.q[j]*sin(kps*r12);
  // 	    }
  // 	  fac*=sj;
  // 	  sk=sk+fac*kps;
  // 	}
  //     sys.force[i]=sys.q[i]*sk;
  //   }
  // short range force
  for(int i = 0; i < sys.r.size(); ++i)
    {
      dVec sj;
      for(int j = 0; j < sys.r.size(); ++j)
	{
	  if(i==j){continue;}
	  dVec r12;
	  double dist, dist2;
	  sys.getDisp(i,j,dist,dist2,r12);
	  if(dist>sys.params.rcut){continue;}
	  double fac=0;
	  double fac2=0;
	  if(sys.type[i]=='O' && sys.type[j]=='O')
	    {
	      fac=sys.params.AOO;
	      fac2=(sys.params.RO*2-dist)/sys.params.rho;
	    }
	  else if(sys.type[i]=='H' && sys.type[j]=='H')
	    {
	      fac=sys.params.AHH;
	      fac2=(sys.params.RH*2-dist)/sys.params.rho;
	    }
	  else
	    {
	      fac=sys.params.AOH;
	      fac2=(sys.params.RO+sys.params.RH-dist)/sys.params.rho;
	    }
	  double temp=fac*exp(fac2)/dist/sys.params.rho;
	  sj=sj+temp*r12; // Need to check sign!!!
	}
      sys.force[i]=sys.force[i]+sj;
    }
  return energy;
}

double ewald(SystemClass& sys)
{
  double energy = 0;
  // real space 
  double realEnergy=0.0;
  for(int i = 0; i < sys.r.size(); ++i)
    {
      for(int j = i+1; j < sys.r.size(); ++j)
	{
	  dVec r12;
	  double dist;
	  double dist2;
	  sys.getDisp(i,j,dist,dist2,r12);
	  double q1=sys.q[i];
	  double q2=sys.q[j];
	  double fac=q1*q2*erfc(dist*sqrt(sys.params.alpha))/dist;
	  realEnergy += fac;
	}
    }
  // k space
  double kEnergy=0;
  // calculate rho(k)
  for(int k = 0; k < sys.k.size(); ++k)
    {
      dVec kps=sys.k[k];
      double k2=kps*kps;
      complex<double> rho=0;
      for(int i = 0; i<sys.r.size(); ++i) 
	{
	  rho+=sys.q[i]*exp(complex<double>(0, 1.0)*(kps*sys.r[i]));
	}
      double fac=(4.0*M_PI/k2*(rho*conj(rho))*exp(-k2/4.0/sys.params.alpha)).real();
      kEnergy += fac / 2.0 / sys.params.vol;
    }
  return realEnergy+kEnergy;
}
