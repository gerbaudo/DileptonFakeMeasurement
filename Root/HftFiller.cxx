#include "SusyTest0/HftFiller.h"

#include "SusyTest0/kinematic.h"
#include "SusyTest0/utils.h"
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
void HftFiller::assignDilepVars(HistFitterTree* const tree, const susy::wh::kin::DilepVars &v)
{
    if(tree) {
        const float gev2mev = 1.0e3;
        tree->eventweight = v.weight;
        tree->L2qFlipWeight = v.qflipWeight;
        tree->isEE = v.isEe;
        tree->isEMU = v.isEm;
        tree->isMUMU = v.isMm;
        tree->isOS = !v.isSs;
        tree->L2nCentralLightJets = v.numCentralLightJets;
        tree->lept1Pt = v.pt0*gev2mev;
        tree->lept2Pt = v.pt1*gev2mev;
        tree->L2Mll = v.mll*gev2mev;
        tree->deltaEtaLl = v.detall;
        tree->L2METrel = v.metrel*gev2mev;
        tree->Ht = v.ht*gev2mev;
        tree->mlj = v.mlj*gev2mev;
        tree->mljj = v.mljj*gev2mev;
        tree->mtmax = v.mtmax()*gev2mev;
        tree->mtllmet = v.mtllmet*gev2mev;
    }
}
//----------------------------------------------------------
bool HftFiller::fill(size_t systIndex, const susy::wh::kin::DilepVars &v, unsigned int run, unsigned int event)
{
    bool someBytesWritten(false);
    if(HistFitterTree *t = m_hftTrees[systIndex]) {
        t->runNumber = run;
        t->eventNumber = event;
        assignDilepVars(t, v);
        t->WriteTree();
        someBytesWritten = true;
    }
    return someBytesWritten;
}
//----------------------------------------------------------
bool HftFiller::init(const std::string &mcid, const std::vector<std::string> &systematics)
{
    for(size_t i=0; i<systematics.size(); ++i) {
            const string &sysname = systematics[i];
            string filename = m_outdir+"/"+sysname+"_"+mcid+".root";
            cout<<"initializing "<<mcid<<" -> "<<sysname<<endl;
            m_hftTrees.push_back(new HistFitterTree(mcid.c_str(), sysname.c_str(), filename));
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
HftFiller& HftFiller::setOutputDir(const std::string dir)
{
    if(dirExists(dir)) m_outdir = dir;
    else {
        bool validDir(mkdirIfNeeded(dir).length()>0);
        if(validDir) m_outdir = dir;
        else {
            cout<<"HftFiller::setOutputDir: cannot set to '"<<dir<<"', using './'"<<endl;
            m_outdir = "./";
        }
    }
    return *this;
}
//----------------------------------------------------------
