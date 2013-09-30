
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

#include "SusyTest0/FakeClosurePlot.h"
#include "SusyTest0/myHist.h"

using namespace std;

void help()
{
  cout << "  Options:"                            << endl;
  cout << "  -d sets the debug level"             << endl;
  cout << "  -i inputdir"                         << endl;
  cout << "  -t inputtag"                         << endl;
  cout << "  -o outputdir"                        << endl;
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
  string inputdir, outputdir, tag;

  for(int i = 1; i < argc; i++) {
    const string opt(argv[i]);
    if      (opt=="-d"   )       dbg = atoi(argv[++i]);
    else if (opt=="-i"   )  inputdir = argv[++i];
    else if (opt=="-o"   ) outputdir = argv[++i];
    else if (opt=="-t"   )       tag = argv[++i];
    else if (opt=="--int")  integral = true;
    else if (opt=="-r"   )    option = (FPRunOption) atoi(argv[++i]);
    else { help(); return 0; }
  }
  bool missingReqOpt(inputdir.size()==0 || tag.size()==0 || outputdir.size()==0);
  if(missingReqOpt) { cout<<"missing required option"<<endl; help(); return 0; }
  string s_option = "ALL";
  if( option == RO_SR1 )          s_option = "SR1";
  else if( option == RO_SR2 )     s_option = "SR2";
  else if( option == RO_SR3 )     s_option = "SR3";
  else if( option == RO_SR4 )     s_option = "SR4";
  else if( option == RO_SR5 )     s_option = "SR5";
  else if( option == RO_ZWindow ) s_option = "Z Window";
  else if( option == RO_VR1 )     s_option = "VR1";
  else if( option == RO_VR2 )     s_option = "VR2";

  cout<<"options:" << endl;
  cout<<"  dbg                      "<<dbg      <<endl;
  cout<<"  tag                      "<<tag      <<endl;
  cout<<"  inputdir                 "<<inputdir <<endl;
  cout<<"  outputdir                "<<outputdir<<endl;
  cout<<"  Region to plot:          "<<s_option<<endl;
  cout<<endl;

  FakeClosurePlot plot;
  plot.setTag(tag);
  plot.setInputDir(inputdir);
  plot.setOuputDir(outputdir);
  plot.setDebug(dbg);
  plot.setIntegralOption(integral);
  plot.init(option);

  plot.DataMCAnaPlots();
  return 0;
}
