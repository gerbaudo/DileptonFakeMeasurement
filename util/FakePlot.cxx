
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

#include "SusyTest0/FakePlotting.h"
#include "SusyTest0/myHist.h"

using namespace std;

/*

Plotting Macro

*/

void help()
{
  cout << "  Options:"                            << endl;

  cout << "  -n number of events to process"      << endl;
  cout << "     defaults: -1 (all events)"        << endl;

  cout << "  -r run option number"                << endl;
  cout << "     0 Data (default)"                 << endl;
  cout << "     1 MC"                             << endl;
  cout << "     2 Data and MC"                    << endl;
  cout << "     3 Data - MC Rate"                 << endl;
  cout << "     4 Photon + jet"                   << endl;
  cout << "     5 Data/MC scale factors"          << endl;
  cout << "     6 MC Signal Regions rates"        << endl;
  cout << "     7 MC Composition in SR"           << endl;
  cout << "     8 TT/TL/LT/LL Matrix Pred Plots"  << endl;
  cout << "     9 dump Matrix Pred in SR"         << endl;
  cout << "     10 HF normalization"              << endl;
  cout << "     11 Dump Percentages"              << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 0;
  RunOption ro = RO_Data;
  
  cout << "FakePlot" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "-r") == 0)
      ro = (RunOption) atoi(argv[++i]);
    //if (strcmp(argv[i], "-h") == 0)
    else
      {
	help();
	return 0;
      }
  }
  
  // Figure out which run option was chosen

  string runOption = "Data";
  if(ro == RO_MC)     runOption = "MC";
  if(ro == RO_DataMC) runOption = "Data and MC";
  
  cout << "flags:" << endl;
  cout << "  dbg                      " << dbg         << endl;
  cout << "  Run Option               " << runOption   << endl;
  cout << endl;
  

  // Create instance of the class:
  FakePlotting* plot = new FakePlotting(ro);
  plot->init();
  plot->setDebug(dbg);

  if(ro == RO_MC)        plot->MCFakeRate();
  if(ro == RO_Data)      plot->DataFakeRate();
  if(ro == RO_DataMCSub) plot->DataRateMCSub();
  if(ro == RO_GJetCR)    plot->GammaJetCRRates();
  if(ro == RO_DataMCSF)  plot->DataMCSF(ro);
  if(ro == RO_DataMC)    plot->DataMCSF(ro);
  if(ro == RO_SRRates)   plot->MCSRRates();
  if(ro == RO_SRComp){   /*plot->Composition();*/ plot->dumpSRTable(); }
  if(ro == RO_TLInfo)    plot->TLPlots();
  if(ro == RO_SRDump)    plot->dumpSRFake();
  if(ro == RO_HFNorm)    plot->GetHFNorm();
  if(ro == RO_DumpPer)   plot->checkPercentages();

  cout << endl;
  cout << "Plotting job done" << endl;

  return 0;
}
