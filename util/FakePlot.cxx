
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
  cout << "  Options:"                             << endl;
  cout << "  -i inputdir (required)"               << endl;
  cout << "  -t tag (required, e.g. '_Sep_23')"    << endl;
  cout << "  -c itecorrfile(required for some mode)"<< endl;
  cout << "  -o outputdir (required)"              << endl;
  cout << "  -r run option number"                 << endl;
  size_t i=0;
  size_t width(24);
  cout<<"  "<<setw(3)<<i<<" "<<setw(width)<<opts[i]<<" : Data/MC scale factors"          << endl; i++;
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
  RunOption ro = RO_DataMCSF;

  vector<string> txtOptions; // ugly duplication, but allows for mnemonic...fixme
  txtOptions.push_back("datamc-sf"       );
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

  bool corrFileIsRequired(ro==RO_DataMCSF);
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

  if(ro == RO_DataMCSF)  plot->DataMCSF();
  else cout<<"Obsolete or unsupported option"<<endl;
  cout << endl;
  cout << "Plotting job done" << endl;

  return 0;
}
