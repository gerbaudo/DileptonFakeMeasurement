#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/utils.h"

/*

Measure Fake Rate

Imported from Matt's code.

*/


using std::cout;
using std::endl;
using std::string;

void usage(const char *exeName) {
  cout<<"Usage:"<<endl
      <<exeName<<" options"<<endl
      <<"\t"<<"-n [--num-event]   nEvt (default -1, all)"<<endl
      <<"\t"<<"-k [--num-skip]    nSkip (default 0)"     <<endl
      <<"\t"<<"-i [--input]       (file, list, or dir)"  <<endl
      <<"\t"<<"-o [--output]      output file"           <<endl
      <<"\t"<<"-s [--sample]      samplename"            <<endl
      <<"\t"<<"--mcTrig           use MC triggers"       <<endl
      <<"\t"<<"--altIso           use 2011 isolation"    <<endl
      <<"\t"<<"--optCut           determine optimum cuts"<<endl
      <<"\t"<<"-d [--debug]     : debug (>0 print stuff)"<<endl
      <<"\t"<<"-h [--help]      : print help"            <<endl
      <<endl;
}

int main(int argc, char** argv)
{

  ROOT::Cintex::Cintex::Enable();
  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  string sample;
  string input;
  string output;
  bool useAltIso = false;
  bool useMCTrig = false;
  bool optCuts   = false;

  int optind(1);
  while ((optind < argc)) {
    if(argv[optind][0]!='-') {
      if(dbg) cout<<"skip "<<argv[optind]<<endl;
      optind++;
      continue;
    }
    std::string sw = argv[optind];
    if     (sw=="-n"||sw=="--num-event"  ) { nEvt = atoi(argv[++optind]); }
    else if(sw=="-k"||sw=="--num-skip"   ) { nSkip = atoi(argv[++optind]); }
    else if(sw=="-d"||sw=="--debug"      ) { dbg = atoi(argv[++optind]); }
    else if(sw=="-i"||sw=="--input"      ) { input = argv[++optind]; }
    else if(sw=="-o"||sw=="--output"     ) { output = argv[++optind]; }
    else if(sw=="-s"||sw=="--sample"     ) { sample = argv[++optind]; }
    else if(sw=="--mcTrig"               ) { useMCTrig = true; ++optind; }
    else if(sw=="--altIso"               ) { useAltIso = true; ++optind; }
    else if(sw=="--optCut"               ) { optCuts = true; ++optind; }
    else if(sw=="-h"||sw=="--help"       ) { usage(argv[0]); return 0; }
    else cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end while(optind<argc)

  cout<<"flags:"                <<endl
      <<"  sample  "<<sample    <<endl
      <<"  nEvt    "<<nEvt      <<endl
      <<"  nSkip   "<<nSkip     <<endl
      <<"  dbg     "<<dbg       <<endl
      <<"  input   "<<input     <<endl
      <<"  output  "<<output    <<endl
      <<"  mcTrig  "<<useMCTrig <<endl
      <<"  altIso  "<<useAltIso <<endl
      <<"  optCut  "<<optCuts   <<endl
      <<endl;

  TChain* chain = new TChain("susyNt");
  bool inputIsFile(string::npos!=input.find(".root"));
  bool inputIsList(string::npos!=input.find(".txt"));
  bool inputIsDir (endswith(input, "/"));
  bool validInput(inputIsFile||inputIsList||inputIsDir);
  bool validOutput(endswith(output, ".root"));
  if(!validInput || !validOutput) {
    cout<<"invalid "
        <<(validOutput ? "input" : "output")<<" "
        <<"'"<<(validOutput ? input : output)<<"'"
        <<endl;
    usage(argv[0]); return 1;
  }
  if(inputIsFile) ChainHelper::addFile    (chain, input);
  if(inputIsList) ChainHelper::addFileList(chain, input);
  if(inputIsDir ) ChainHelper::addFileDir (chain, input);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();

  MeasureFakeRate2 mfr;
  mfr.setDebug(dbg);
  if(sample.size()) mfr.setSampleName(sample);
  if(output.size()) mfr.setFileName(output);
  mfr.setUseMCTrig(useMCTrig);
  if(useAltIso) mfr.setAltIso();
  if(optCuts)   mfr.setFindOptCut(optCuts);

  mfr.buildSumwMap(chain);
  chain->Process(&mfr, sample.c_str(), nEvt, nSkip);
  delete chain;
  return 0;
}
