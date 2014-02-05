
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyTest0/SusyPlotter.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;

/*

SusyPlotter - perform selection and plot histos

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
      <<"\t"<<"-S [--with-syst]   fill also syst histos" <<endl
      <<"\t"<<"-d [--debug]     : debug (>0 print stuff)"<<endl
      <<"\t"<<"-h [--help]      : print help"            <<endl
      <<"\t"<<"--WH-sample      : xsec from SusyXSReader"<<endl
      <<endl;
}

bool endswith(const string &s, const string &end) {
  //http://stackoverflow.com/questions/874134/find-if-string-endswith-another-string-in-c
  if(s.length()<end.length()) return false;
  else return (0==s.compare(s.length() - end.length(), end.length(), end));
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
  bool useSusyXSReader = false;
  bool doSyst = false;

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
    else if(sw=="-S"||sw=="--with-syst"  ) { doSyst = true; }
    else if(sw=="-h"||sw=="--help"       ) { usage(argv[0]); return 0; }
    else if(sw=="--WH-sample"            ) { useSusyXSReader = true; }
    else cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end while(optind<argc)

  cout<<"flags:"             <<endl
      <<"  sample  "<<sample <<endl
      <<"  nEvt    "<<nEvt   <<endl
      <<"  nSkip   "<<nSkip  <<endl
      <<"  dbg     "<<dbg    <<endl
      <<"  doSyst  "<<doSyst <<endl
      <<"  input   "<<input  <<endl
      <<"  output  "<<output <<endl
      <<endl;

  TChain* chain = new TChain("susyNt");
  bool inputIsFile(string::npos!=input.find(".root"));
  bool inputIsList(string::npos!=input.find(".txt"));
  bool inputIsDir (endswith(input, "/"));
  if(inputIsFile) ChainHelper::addFile    (chain, input);
  if(inputIsList) ChainHelper::addFileList(chain, input);
  if(inputIsDir ) ChainHelper::addFileDir (chain, input);
  if(dbg>0) chain->ls();
  Long64_t nEntries = chain->GetEntries();
  SusyPlotter* susyPlot = new SusyPlotter();
  susyPlot->setDebug(dbg);
  susyPlot->setSampleName(sample);
  susyPlot->setUseXsReader(useSusyXSReader);
  if(doSyst) susyPlot->toggleSystematics();
  if(output.length()) susyPlot->setOutputFilename(output);
  susyPlot->buildSumwMap(chain);
  nEvt = (nEvt>0 ? nEvt : nEntries);
  cout<<"Total entries:   "<<nEntries<<endl
      <<"Process entries: "<<nEvt    <<endl;
  chain->Process(susyPlot, sample.c_str(), nEvt, nSkip);
  cout<<endl
      <<"SusyPlotter job done"<<endl;

  delete chain;
  return 0;
}
