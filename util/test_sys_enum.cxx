/*
  Unit test for Systematics.h (validate names and self-consistentcy)

  davide.gerbaudo@gmail.com
  Feb 2014
*/

#include <iostream>
#include <string>

#include "SusyTest0/Systematics.h"
#include "SusyNtuple/SusyDefs.h"

using std::cout;
using std::endl;
using std::string;
namespace swh = susy::wh;
using swh::Systematic;

int main(int argc, char** argv)
{
    size_t failing=0;
    cout<<"\tSusyNtSys\tswh::Systematic"<<endl;
    for(int i=0; i<NtSys_N; ++i){
        SusyNtSys sNt = static_cast<SusyNtSys>(i);
        Systematic sWh = swh::ntsys2sys(sNt);
        SusyNtSys sNt2 = swh::sys2ntsys(sWh);
        if(sNt!=sNt2) failing++;
        cout<<"["<<i<<"]"
            <<" : "<<swh::syst2str(sNt)
            <<"\t -> "<<swh::syst2str(sWh)
            <<" ; "<<sNt<<" : "<<sNt2
            <<(sNt==sNt2 ? "" : "\t <<<<<<<")
            <<endl;
    }
    if(failing) cout<<"\nFailing "<<failing<<" cases."<<endl;
    else        cout<<"All fine."<<endl;
    return 0;
}
