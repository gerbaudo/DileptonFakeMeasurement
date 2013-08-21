
#include <cstdlib>
#include <iostream>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyTest0/MatrixPrediction.h"
#include "SusyTest0/TightProbability.h"


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
      <<"\t"<<"-n [--num-event]   nEvt"       <<endl
      <<"\t"<<"-k [--num-skip]    nSkip"      <<endl
      <<"\t"<<"-d [--debug]    :  debug"      <<endl
      <<"\t"<<"-F [--input-file]  inputFile"  <<endl
      <<"\t"<<"-f [--input-list]  inputList"  <<endl
      <<"\t"<<"-D [--input-dir]   inputDir"   <<endl
      <<"\t"<<"-o [--output-file] samplename" <<endl
      <<"\t"<<"-s [--sample]      output file"<<endl
      <<"\t"<<"-h [--help]      : print help" <<endl
      <<endl;
}

int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  string sample;
  string file;
  string fileList;
  string fileDir;
  string outFile;

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
    else if(sw=="-F"||sw=="--input-file" ) { file = argv[++optind]; }
    else if(sw=="-f"||sw=="--input-list" ) { fileList = argv[++optind]; }
    else if(sw=="-D"||sw=="--input-dir"  ) { fileDir = argv[++optind]; }
    else if(sw=="-o"||sw=="--output-file") { outFile = argv[++optind]; }
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
      <<"  input   "<<(file.size() ? file :
                       (fileList.size() ? fileList :
                        (fileDir.size() ? fileDir :
                         "None")))<<endl
      <<"  output  "<<outFile<<endl
      <<endl;

  // Build the input chain
  TChain* chain = new TChain("susyNt");
  ChainHelper::addFile(chain, file);
  ChainHelper::addFileList(chain, fileList);
  ChainHelper::addFileDir(chain, fileDir);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();
  Susy::TightProbability tp;
  tp.setDebug(dbg);
  if(sample.size()) tp.setSampleName(sample);
  if(outFile.size()) tp.setOutputFilename(outFile);
  tp.buildSumwMap(chain);
  chain->Process(&tp, sample.c_str(), nEvt, nSkip);
  delete chain;
  return 0;
}
