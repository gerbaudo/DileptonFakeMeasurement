#include "SusyTest0/FinalNewFake.h"

#include <iostream>
#include <string>

#include "Cintex/Cintex.h"

using namespace std;

void help()
{
  cout<<"  Options:"            << endl;
  cout<<"  -d debuglevel"       << endl;
  cout<<"  -o outfile.root"     << endl;
  cout<<"  -h print this help"  << endl;
}

int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  int dbg = 0;
  string output = "fakeout";
  for(int i = 1; i < argc; i++) {
    string opt(argv[i]);
    if      (opt=="-d")    dbg = atoi(argv[++i]);
    else if (opt=="-o") output = argv[++i];
    else if (opt=="-h") { help(); return 0; }
    else { cout<<"invalid option '"<<opt<<"'"<<endl; help(); return 0; }
  }

  cout<<"options:"<<endl;
  cout<<"\t dbg          : "<<dbg   <<endl;
  cout<<"\t output file  : "<<output<<endl;
  cout<<endl;

  FinalNewFake format(output);
  format.setDebug(dbg);
  format.buildRates();

  return 0;
}
