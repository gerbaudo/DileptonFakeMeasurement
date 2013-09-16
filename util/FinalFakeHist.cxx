
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

//#include "SusyTest0/FinalFake.h"
#include "SusyTest0/FinalNewFake.h"
//#include "SusyTest0/FinalChannelSepFake.h"

using namespace std;

/*

Plotting Macro

*/

void help()
{
  cout << "  Options:"                            << endl;

  cout << "  -n number of events to process"      << endl;
  cout << "     defaults: -1 (all events)"        << endl;

  cout << "  -o Specify output file name"         << endl;
  cout << "     default is 'fakeout'"         << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 0;
  string output = "fakeout";
  
  cout << "Final Fake Rate" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "-o") == 0)
      output = argv[++i];
    else
      {
	help();
	return 0;
      }
  }
  
    cout << "flags:" << endl;
  cout << "  dbg                      " << dbg         << endl;
  cout << "  output file:             " << output      << endl;
  cout << endl;
  

  // Create instance of the class:
  //FinalFake* format = new FinalFake(output);
  FinalNewFake* format = new FinalNewFake(output);
  //FinalChannelSepFake* format = new FinalChannelSepFake(output);
  format->setDebug(dbg);
  format->buildRates();
  
  cout << endl;
  cout << "Plotting job done" << endl;

  delete format;

  return 0;
}
