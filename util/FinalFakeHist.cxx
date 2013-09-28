#include "SusyTest0/FinalNewFake.h"

#include <iostream>
#include <string>

#include "Cintex/Cintex.h"

using namespace std;

void help()
{
  cout<<"  Options:"                         << endl;
  cout<<"  -i inputdir            (required)"<< endl;
  cout<<"  -t tag (e.g. '_Sep_23', required)"<< endl;
  cout<<"  -o outfile.root        (required)"<< endl;
  cout<<"  -d debuglevel"                    << endl;
  cout<<"  -h print this help"               << endl;
}

int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  int dbg = 0;
  string inputdir, outputfile, tag;
  for(int i = 1; i < argc; i++) {
    string opt(argv[i]);
    if      (opt=="-d")        dbg = atoi(argv[++i]);
    else if (opt=="-i")   inputdir = argv[++i];
    else if (opt=="-o") outputfile = argv[++i];
    else if (opt=="-t")        tag = argv[++i];
    else if (opt=="-h") { help(); return 0; }
    else { cout<<"invalid option '"<<opt<<"'"<<endl; help(); return 0; }
  }
  bool missingRequiredOption(inputdir.size()==0 || outputfile.size()==0 || tag.size()==0);
  if(missingRequiredOption) { cout<<"Missing required option."<<endl; help(); return 0; }

  cout<<"options:"<<endl;
  cout<<"\t inputdir     : "<<inputdir  <<endl;
  cout<<"\t tag          : "<<tag       <<endl;
  cout<<"\t output file  : "<<outputfile<<endl;
  cout<<"\t dbg          : "<<dbg       <<endl;
  cout<<endl;

  FinalNewFake format;
  format.setTag(tag);
  format.setInputDir(inputdir);
  format.setOuputFilename(outputfile);
  format.setDebug(dbg);
  format.initIoFiles();
  format.buildRates();

  return 0;
}
