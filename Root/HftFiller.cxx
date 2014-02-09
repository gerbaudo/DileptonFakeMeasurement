#include "SusyTest0/HftFiller.h"

#include "SusyTest0/kinematic.h"
#include "HistFitterTree/HistFitterTree.h"

#include <cassert>

using std::cout;
using std::endl;
using std::string;
using susy::wh::HftFiller;

//----------------------------------------------------------
HftFiller::HftFiller()
{
}
//----------------------------------------------------------
HftFiller::~HftFiller()
{
    if(m_hftTrees.size()) {
        cout<<"HftFiller::~HftFiller: you forgot to call HftFiller::close"<<endl;
        assert(false); // throw instead
    }
}
//----------------------------------------------------------
bool HftFiller::fill(size_t systIndex, const susy::wh::kin::DilepVars &v)
{
    bool someBytesWritten(false);
    const float gev2mev = 1.0e3;
    if(HistFitterTree *t = m_hftTrees[systIndex]) {
        t->lept1Pt = v.pt0*gev2mev;
        t->lept2Pt = v.pt1*gev2mev;
        // todo: all the other variables.
//         = v.mll;
//         = v.detall;
//         = v.metrel;
//         = v.mlj;
//         = v.mljj;
//         = v.mtmax();
//         = v.mtllmet;
        t->WriteTree();
        someBytesWritten = true;
    }
    return someBytesWritten;
}
//----------------------------------------------------------
bool HftFiller::init(const std::string &mcid, const std::vector<std::string> &systematics)
{
    for(size_t i=0; i<systematics.size(); ++i)
        {
            cout<<"initializing "<<mcid<<" -> "<<systematics[i].c_str()<<endl;
        m_hftTrees.push_back(new HistFitterTree(mcid.c_str(), systematics[i].c_str()));
        }
    return m_hftTrees.size()>0;
}
//----------------------------------------------------------
bool HftFiller::close(float sumw)
{
    for(size_t i=0; i<m_hftTrees.size(); ++i){
        m_hftTrees[i]->setSumOfMcWeights(sumw);
        delete m_hftTrees[i];
    }
    m_hftTrees.clear();
    return true; // anything we should check?
}
//----------------------------------------------------------
