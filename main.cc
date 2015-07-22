#include "ewald.h"
#include "Timer.h"
#include <cstdlib>

int main(int argc, char** argv)
{
  cout << setprecision(15) << fixed;
  SystemClass sys(argv[1], argv[2]);
  // for(int i = 1; i < 100; ++i)
  //   {
  //     sys.params.kcut=.1+i*1.0;
  //     sys.buildk();sys.buildDisp();
  //     cout << i << " " << screened(sys) << " " << sys.force[0] << endl;
  //   }
  double e=screened(sys);
  cerr << e << endl;
  for(int i = 0; i < sys.r.size(); ++i)
    {
      cout << sys.type[i] << " " << sys.q[i] << " " << sys.force[i] << endl;
    }
  //  cout << sys.force[0] << endl;
  // sys.r[0][1]+=1e-6;
  // sys.buildDisp();
  // double e2=screened(sys);
  // cout << (e2-e)/1e-6 << endl;
  // cout << sys.force[0] << endl;
}
