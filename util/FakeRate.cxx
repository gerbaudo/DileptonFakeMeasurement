
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyAna2012/MeasureFakeRate2.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;

/*

Measure Fake Rate

*/

void help()
{
  cout << "  Options:"                          << endl;
  cout << "  -n number of events to process"    << endl;
  cout << "     defaults: -1 (all events)"      << endl;

  cout << "  -k number of events to skip"       << endl;
  cout << "     defaults: 0"                    << endl;

  cout << "  -d debug printout level"           << endl;
  cout << "     defaults: 0 (quiet) "           << endl;

  cout << "  -f name of input filelist"         << endl;
  cout << "     defaults: ''"                   << endl;

  cout << "  -s sample name, for naming files"  << endl;
  cout << "     defaults: ntuple sample name"   << endl;

  cout << "  -mc measure MC rates"              << endl;
  cout << "      default is off"                << endl;

  cout << "  -mcTrig use MC triggers"           << endl;
  cout << "      default is off"                << endl;

  cout << "  -altIso to use 2011 isolation"    << endl;
  cout << "      default is off"                << endl;

  cout << "  --optCut to fill histos useful for"<<endl;
  cout << "           determining optimum cuts "<<endl;

  cout << "  -h print this help"                << endl;
}

int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  bool useAltIso = false;
  bool ismc       = false;
  bool useMCTrig  = false;
  string sample;
  string fileList;
  string output = "";
  bool optCuts = false;

  cout << "MeasureFakeRate2" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0)
      nEvt = atoi(argv[++i]);
    else if (strcmp(argv[i], "-k") == 0)
      nSkip = atoi(argv[++i]);
    else if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "-f") == 0)
      fileList = argv[++i];
    else if (strcmp(argv[i], "-o") == 0)
      output = argv[++i];
    else if (strcmp(argv[i], "-s") == 0)
      sample = argv[++i];
    else if (strcmp(argv[i], "-mc") == 0)
      ismc = true;
    else if (strcmp(argv[i], "-altIso") == 0)
      useAltIso = true;
    else if (strcmp(argv[i], "-mcTrig") == 0)
      useMCTrig = true;
    else if (strcmp(argv[i], "--optCut") == 0)
      optCuts = true;
    else
    {
      help();
      return 0;
    }
  }

  // If no input specified except sample name, use a standard fileList
  if(fileList.empty())
    return 0;

  // Save the file name
  string fname = output;
  if(output.empty()){
    int pos0 = fileList.find("/");  
    int pos1 = fileList.find(".txt");
    fname = fileList.substr(pos0+1,pos1-pos0-1);
  }
  sample = fname;

  cout << "flags:" << endl;
  cout << "  sample  " << sample   << endl;
  cout << "  nEvt    " << nEvt     << endl;
  cout << "  nSkip   " << nSkip    << endl;
  cout << "  dbg     " << dbg      << endl;
  cout << "  input   " << fileList << endl;
  cout << "  output  " << fname    << endl;
  cout << endl;

  // Build the input chain
  TChain* chain = new TChain("susyNt");
  ChainHelper::addFileList(chain, fileList);
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  cout << "Chain built" << endl;
  // Build the TSelector
  MeasureFakeRate2 *fakeRate = new MeasureFakeRate2();
  fakeRate->buildSumwMap(chain);
  fakeRate->setDebug(dbg);
  fakeRate->setSampleName(sample);
  fakeRate->setFileName(fname);
  fakeRate->setIsMC(ismc);
  fakeRate->setUseMCTrig(useMCTrig);
  if(useAltIso) fakeRate->setAltIso();
  fakeRate->setFindOptCut(optCuts);

  // Run the job
  if(nEvt<0) nEvt = nEntries;
  cout << endl;
  cout << "Total entries:   " << nEntries << endl;
  cout << "Process entries: " << nEvt << endl;
  if(nEvt>0) chain->Process(fakeRate, sample.c_str(), nEvt, nSkip);

  cout << endl;
  cout << "MeasureFakeRate job done" << endl;

  delete chain;
  return 0;
}
