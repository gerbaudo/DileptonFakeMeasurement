
#include <cstdlib>
#include <string>

#include "Cintex/Cintex.h"

#include "SusyTest0/IterativeFakeCorrection.h"

using namespace std;

/*

Plotting Macro

*/

void help()
{
  cout << "  Options:"                            << endl;

  cout << "  -n specify number of iterations"     << endl;
  cout << "     default: 8"                       << endl;
  cout<<" --input-mc   filemc.root"  <<endl;
  cout<<" --input-data filedata.root"<<endl;
  cout<<" --output     fileout.root" <<endl;
  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int dbg = 0;
  int n   = 8;
  string inputMcFile = "out/fakerate/merged/allBkgButHf_Sep_14.root";
  string inputDataFile = "out/fakerate/merged/data_Sep_14.root";
  string outputFile = "corFake_Sep11_2013_forDavide.root";

  cout << "Iterative Fake Rate" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0) dbg = atoi(argv[++i]);
    if (strcmp(argv[i], "-n") == 0) n = atoi(argv[++i]);
    if(strcmp(argv[i],"--input-mc")  ==0) inputMcFile   = argv[++i];
    if(strcmp(argv[i],"--input-data")==0) inputDataFile = argv[++i];
    if(strcmp(argv[i],"--output")    ==0) outputFile    = argv[++i];

    else { help(); return 0; }
  }

  cout<<"options:"<<endl
      <<"\t inputMcFile   : "<<inputMcFile<<endl
      <<"\t inputDataFile : "<<inputDataFile<<endl
      <<"\t outputFile    : "<<outputFile<<endl
      <<"\t n-iterations  : "<<n<<endl
      <<"\t dbg-level     : "<<dbg<<endl
      <<endl;

  // Create instance of the class:
  IterativeFakeCorrection* iter = new IterativeFakeCorrection();
  iter->setDebug(dbg);
  iter->setNIter(n);
  iter->setInputMc(inputMcFile);
  iter->setInputData(inputDataFile);
  iter->setOutputFilename(outputFile);
  iter->iterate();

  cout << endl;
  cout << "Plotting job done" << endl;

  return 0;
}
