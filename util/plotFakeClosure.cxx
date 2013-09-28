
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

#include "SusyTest0/FancyPlotting.h"
#include "SusyTest0/myHist.h"

using namespace std;

/*

Plotting Macro

*/

void help()
{
  cout << "  Options:"                            << endl;

  cout << "  -d sets the debug level"             << endl;
  cout << "     default is 0 (off)"               << endl;

  cout << "  -r determine what Region to plot"    << endl;
  cout << "     0 -- All (default)"               << endl;
  cout << "     1 -- SR1"                         << endl;
  cout << "     2 -- SR2"                         << endl;
  cout << "     3 -- SR3"                         << endl;
  cout << "     4 -- SR4"                         << endl;
  cout << "     5 -- SR5"                         << endl;
  cout << "     6 -- ZWindow"                     << endl;
  cout << "     7 -- VR1"                         << endl;
  cout << "     8 -- VR2"                         << endl;
  cout << "     9 -- OSInc"                       << endl;
  cout << "    10 -- SSInc"                       << endl;
  cout << "    11 -- VR3"                         << endl;
  cout << "    12 -- VR4"                         << endl;
  cout << "    13 -- VRTL"                        << endl;
  cout << "  --int display integral in legend"    << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 0;
  FPRunOption option = RO_ALL;
  bool integral = false;
  bool makeTable = false;
  
  cout << "FancyPlot" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "--int") == 0)
      integral = true;
    else if (strcmp(argv[i], "-r") == 0)
      option = (FPRunOption) atoi(argv[++i]);
    else if (strcmp(argv[i], "-tab") == 0)
      makeTable = true;
    else
      {
	help();
	return 0;
      }
  }
  
  // Figure out which run option was chosen
  string s_option = "ALL";
  if( option == RO_SR1 )          s_option = "SR1";
  else if( option == RO_SR2 )     s_option = "SR2";
  else if( option == RO_SR3 )     s_option = "SR3";
  else if( option == RO_SR4 )     s_option = "SR4";
  else if( option == RO_SR5 )     s_option = "SR5";
  else if( option == RO_ZWindow ) s_option = "Z Window";
  else if( option == RO_VR1 )     s_option = "VR1";
  else if( option == RO_VR2 )     s_option = "VR2";
  

  
  cout << "flags:" << endl;
  cout << "  dbg                      " << dbg      << endl;
  cout << "  Region to plot:          " << s_option << endl;
  cout << endl;
  

  // Create instance of the class:
  FancyPlotting* plot = new FancyPlotting();
  plot->init(option);
  plot->setDebug(dbg);
  if(integral) plot->addIntegral();
  if(makeTable) plot->makeTable();
  plot->DataMCAnaPlots();

  
  cout << endl;
  cout << "Plotting job done" << endl;

  return 0;
}
