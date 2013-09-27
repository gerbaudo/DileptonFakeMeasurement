
#include <algorithm> // std::find
#include <cstdlib>
#include <iomanip> // std::setw
#include <iterator>  // std::distance
#include <string>
#include <vector>

#include "Cintex/Cintex.h"

#include "SusyTest0/FakePlotting.h"
#include "SusyTest0/myHist.h"
#include "SusyTest0/utils.h"

using namespace std;

/*

Plotting Macro

*/

void help(const vector<string> &opts)
{
  size_t i=0;
  cout << "  Options:"                             << endl;
  cout << "  -i inputdir (required)"               << endl;
  cout << "  -t tag (required, e.g. '_Sep_23')"    << endl;
  cout << "  -c itecorrfile(required for some mode)"<< endl;
  cout << "  -o outputdir (required)"              << endl;
  cout << "  -r run option number"                 << endl;
  size_t width(24);
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Data (default)"                 << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : MC"                             << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Data and MC"                    << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Data - MC Rate"                 << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Photon + jet"                   << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Data/MC scale factors"          << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : MC Signal Regions rates"        << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : MC Composition in SR"           << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : TT/TL/LT/LL Matrix Pred Plots"  << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : dump Matrix Pred in SR"         << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : HF normalization"               << endl; i++;
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Dump Percentages"               << endl; i++;
  cout << "  -h print this help"                << endl;
  cout<<endl
      <<"Example:"<<endl
      <<" FakePlot -r datamc-sf -t _Sep_23 "
      <<" -i out/fakerate/merged/ -o out/fakerate/merged/"
      <<" -c out/fakerate/merged/iterative_out_Sep23.root"
      <<endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 0;
  string inputdir, inputItecorrfile, outputdir, tag;
  RunOption ro = RO_Data;

  vector<string> txtOptions; // ugly duplication, but allows for mnemonic...fixme
  txtOptions.push_back("data"            );
  txtOptions.push_back("mc"              );
  txtOptions.push_back("data-and-mc"     );
  txtOptions.push_back("data-mc-rate"    );
  txtOptions.push_back("photon-jet"      );
  txtOptions.push_back("datamc-sf"       );
  txtOptions.push_back("mc-sr-rates"     );
  txtOptions.push_back("mc-sr-comp"      );
  txtOptions.push_back("matrix-pred-plot");
  txtOptions.push_back("matrix-pred-dump");
  txtOptions.push_back("hf-norm"         );
  txtOptions.push_back("percent-dump"    );

  cout << "FakePlot" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if      (strcmp(argv[i], "-d") == 0) { dbg = atoi(argv[++i]);        }
    else if (strcmp(argv[i], "-t") == 0) { tag              = argv[++i]; }
    else if (strcmp(argv[i], "-i") == 0) { inputdir         = argv[++i]; }
    else if (strcmp(argv[i], "-c") == 0) { inputItecorrfile = argv[++i]; }
    else if (strcmp(argv[i], "-o") == 0) { outputdir        = argv[++i]; }
    else if (strcmp(argv[i], "-r") == 0) {
      string sw(argv[++i]);
      if(isInt(sw)) ro = (RunOption) atoi(sw.c_str());
      else {
        vector<string>::iterator it = std::find(txtOptions.begin(), txtOptions.end(), sw);
        bool isValidOpt(it != txtOptions.end());
        if(isValidOpt) { ro = (RunOption)distance(txtOptions.begin(), it); }
        else { cout<<"invalid text option '"<<sw<<"'"<<endl; help(txtOptions); return 0; }
      }
    } // end if('-r')
    else { cout<<"invalid option '"<<argv[i]<<"'"<<endl; help(txtOptions); return 0; }
  } // end for(i)
  bool unspecifiedOutdir(outputdir.size()==0), unspecifiedCorr(inputItecorrfile.size()==0);
  bool unspecifiedInputdir(inputdir.size()==0), unspecifiedTag(tag.size()==0);
  if(unspecifiedTag)      { cout<<"tag required."      <<endl; help(txtOptions); return 0; }
  if(unspecifiedInputdir) { cout<<"inputdir required." <<endl; help(txtOptions); return 0; }
  if(unspecifiedOutdir)   { cout<<"outputdir required."<<endl; help(txtOptions); return 0; }

  bool corrFileIsRequired(ro==RO_DataMCSF || ro==RO_DataMC);
  if(corrFileIsRequired && unspecifiedCorr) {
    cout<<"inputItecorrfile required. (example 'corFake_Sep11_2013_forDavide.root')"<<endl;
    help(txtOptions);
    return 0;
  }

  cout <<"flags:" << endl;
  cout <<"  dbg                      " << dbg              << endl;
  cout <<"  tag                      " << tag              << endl;
  cout <<"  inputdir                 " << inputdir         << endl;
  cout <<"  outputdir                " << outputdir        << endl;
  cout <<"  inputItecorrfile         " << inputItecorrfile << endl;
  cout <<"  Run Option               "
       <<" "<<txtOptions[ro]<<" ("<<ro<<") "<< endl;
  cout << endl;

  // Create instance of the class:
  FakePlotting* plot = new FakePlotting(ro);
  plot->setDebug(dbg);
  plot->setTag(tag);
  plot->setInputDir(inputdir);
  plot->setOuputDir(outputdir);
  plot->setInputItercorrFile(inputItecorrfile);
  plot->init();

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
