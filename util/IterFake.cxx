
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

#include "SusyTest0/IterativeFakeCorrection.h"

using namespace std;

/*

Plotting Macro

*/

void help()
{
  cout << "  Options:"                            << endl;

  cout << "  -n specify number of iterations"     << endl;
  cout << "     default: 8"                       << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 10;
  int n   = 8;
  
  cout << "Iterative Fake Rate" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    if (strcmp(argv[i], "-n") == 0)
      n = atoi(argv[++i]);
    else
      {
	help();
	return 0;
      }
  }
  
  // Figure out which run option was chosen
  cout << "flags:" << endl;
  cout << "  dbg                      " << dbg         << endl;
  cout << endl;
  

  // Create instance of the class:
  IterativeFakeCorrection* iter = new IterativeFakeCorrection();
  iter->setDebug(dbg);
  iter->setNIter(n);
  iter->iterate();

  cout << endl;
  cout << "Plotting job done" << endl;

  return 0;
}
