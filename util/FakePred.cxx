
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyTest0/MatrixPrediction.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyTest0/utils.h"

using namespace std;

/*

SusyPlotter - perform selection and dump cutflows

*/

void usage(const char *exeName, const char *defaultMatrixFile) {
  cout<<"Usage:"<<endl
      <<exeName<<" options"<<endl
      <<"\t"<<"-m [--matrix-file] input matrix file"     <<endl
      <<"\t"<<" (default "<<defaultMatrixFile<<")"       <<endl
      <<"\t"<<"-n [--num-event]   nEvt (default -1, all)"<<endl
      <<"\t"<<"-k [--num-skip]    nSkip (default 0)"     <<endl
      <<"\t"<<"-i [--input]       (file, list, or dir)"  <<endl
      <<"\t"<<"-o [--output]      output file"           <<endl
      <<"\t"<<"-s [--sample]      samplename"            <<endl
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
  string matrixFile(getRootCoreDir()+"/../SusyMatrixMethod/data/forDavide_Sep11_2013.root");
  int optind(1);
  while ((optind < argc)) {
    std::string sw = argv[optind];
    if(sw[0]!='-') { if(dbg) cout<<"skip "<<sw<<endl; optind++; continue; }
    if     (sw=="-m"||sw=="--matrix-file") { matrixFile = argv[++optind]; }
    else if(sw=="-n"||sw=="--num-event"  ) { nEvt = atoi(argv[++optind]); }
    else if(sw=="-k"||sw=="--num-skip"   ) { nSkip = atoi(argv[++optind]); }
    else if(sw=="-d"||sw=="--debug"      ) { dbg = atoi(argv[++optind]); }
    else if(sw=="-i"||sw=="--input"      ) { input = argv[++optind]; }
    else if(sw=="-o"||sw=="--output"     ) { output = argv[++optind]; }
    else if(sw=="-s"||sw=="--sample"     ) { sample = argv[++optind]; }
    else if(sw=="-h"||sw=="--help"       ) { usage(argv[0], matrixFile.c_str()); return 0; }
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
    usage(argv[0], matrixFile.c_str()); return 1;
  }
  if(inputIsFile) ChainHelper::addFile    (chain, input);
  if(inputIsList) ChainHelper::addFileList(chain, input);
  if(inputIsDir ) ChainHelper::addFileDir (chain, input);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();

  MatrixPrediction fakePred;
  fakePred.setMatrixFilename(matrixFile);
  fakePred.setDebug(dbg);
  fakePred.setSampleName(sample);
  fakePred.setOutputFilename(output);

  fakePred.buildSumwMap(chain);
  chain->Process(&fakePred, sample.c_str(), nEvt, nSkip);

  delete chain;
  return 0;
}
