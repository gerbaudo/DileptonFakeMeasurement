#include "SusyTest0/HtfFiller.h"

#include "SusyTest0/kinematic.h"
#include "HistFitterTree/HistFitterTree.h"

using std::cout;
using std::endl;
using std::string;
using susy::wh::HtfFiller;

//----------------------------------------------------------
HtfFiller::HtfFiller(const std::string &mcid, float sumw):
    m_mcid(mcid),
    m_sumw(0.0)
{
}
//----------------------------------------------------------
HtfFiller::~HtfFiller()
{
    close();
}
//----------------------------------------------------------
bool HtfFiller::fill(size_t systIndex, const susy::wh::kin::DilepVars &v)
{
    bool someBytesWritten(false);
    const float GeV2MeV = 1.0e3;
    if(HistFitterTree *t = m_hftTrees[systIndex]) {
        t->lept1Pt = v.pt0*GeV2MeV;
        t->lept2Pt = v.pt1*GeV2MeV;
        // todo: all the other variables.
//         = v.mll;
//         = v.detall;
//         = v.metrel;
//         = v.mlj;
//         = v.mljj;
//         = v.mtmax();
//         = v.mtllmet;
        someBytesWritten = true;
    }
    return someBytesWritten;
}
//----------------------------------------------------------
bool HtfFiller::init(const std::vector<std::string> &systematics)
{
    for(size_t i=0; i<systematics.size(); ++i) m_hftTrees.push_back(new HistFitterTree(m_mcid.c_str(), systematics[i].c_str()));
    return m_hftTrees.size()>0;
}
//----------------------------------------------------------
bool HtfFiller::close()
{
    for(size_t i=0; i<m_hftTrees.size(); ++i){
        m_hftTrees[i]->setSumOfMcWeights(m_sumw);
        delete m_hftTrees[i];
    }
    m_hftTrees.clear();
    return true; // anything we should check?
}
//----------------------------------------------------------
