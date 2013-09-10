
#include <cstdlib>
#include <iostream>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyTest0/MatrixPrediction.h"
#include "SusyTest0/TightProbability.h"
#include "SusyTest0/utils.h"

/*!
  Main executable to determine from MC the following two
  probabilities:
  - r = p(pass tight lepton requirement | real, prompt lepton passing loose lepton requirement)
  - f = p(pass tight lepton requirement | fake, non-prompt lepton passing loose lepton requirement)

  Both probabilities are determined from MC samples, separately for
  each physical process. Both probabilities are determined from events
  in the same signal region used for the search.

  The two main sources of fake elecrons are photon conversions and
  jets where the energy deposition happened mostly through EM showers.

  The main source of fake muons are semileptonic decays of hadrons
  within jets, mostly heavy flavored ones.

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
      <<"\t"<<"-o [--output]      samplename"            <<endl
      <<"\t"<<"-s [--sample]      output file"           <<endl
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
    else if(sw=="-h"||sw=="--help"       ) { usage(argv[0]); return 0; }
    else cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end while(optind<argc)

  cout<<"flags:"             <<endl
      <<"  sample  "<<sample <<endl
      <<"  nEvt    "<<nEvt   <<endl
      <<"  nSkip   "<<nSkip  <<endl
      <<"  dbg     "<<dbg    <<endl
      <<"  input   "<<input  <<endl
      <<"  output  "<<output <<endl
      <<endl;

  // Build the input chain
  TChain* chain = new TChain("susyNt");
  bool inputIsFile(string::npos!=input.find(".root"));
  bool inputIsList(string::npos!=input.find(".txt"));
  bool inputIsDir (endswith(input, "/"));
  if(inputIsFile) ChainHelper::addFile    (chain, input);
  if(inputIsList) ChainHelper::addFileList(chain, input);
  if(inputIsDir ) ChainHelper::addFileDir (chain, input);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();
  Susy::TightProbability tp;
  tp.setDebug(dbg);
  if(sample.size()) tp.setSampleName(sample);
  if(output.size()) tp.setOutputFilename(output);
  tp.buildSumwMap(chain);
  chain->Process(&tp, sample.c_str(), nEvt, nSkip);
  delete chain;
  return 0;
}
