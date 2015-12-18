#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TChain.h"

#include "SusyNtuple/ChainHelper.h"
#include "DileptonFakeMeasurement/MeasureFakeRate2.h"
#include "DileptonFakeMeasurement/utils.h"

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
      <<"\t"<<"--write-tuple      write tuple (HF CR)"   <<endl
      <<"\t"<<"-d [--debug]     : debug (>0 print stuff)"<<endl
      <<"\t"<<"-h [--help]      : print help"            <<endl
      <<endl;
}

int main(int argc, char** argv)
{

  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  bool writeTuple = false;
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
    else if(sw=="--write-tuple"          ) { writeTuple = true; }
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
      <<"  tuple   "<<writeTuple<<endl
      <<endl;

  TChain* chain = new TChain("susyNt");
  ChainHelper::addInput(chain, input);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();

  MeasureFakeRate2 mfr;
  mfr.setDebug(dbg);
  mfr.setWriteFakeNtuple(writeTuple);
  if(sample.size()) mfr.setSampleName(sample);
  if(output.size()) mfr.setFileName(output);

  chain->Process(&mfr, sample.c_str(), nEvt, nSkip);
  delete chain;
  return 0;
}
