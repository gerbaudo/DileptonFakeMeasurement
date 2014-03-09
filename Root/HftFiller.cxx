#include "SusyTest0/HftFiller.h"

#include "SusyTest0/kinematic.h"
#include "SusyTest0/utils.h"
#include "SusyTest0/XsecUncertainty.h"

#include "HistFitterTree/HistFitterTree.h"

#include <cassert>

using std::cout;
using std::endl;
using std::string;
using susy::wh::HftFiller;

//----------------------------------------------------------
HftFiller::HftFiller() :
    xsecRelativeUncertainty_(1.0)
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
void HftFiller::assignWeightVars(HistFitterTree* const tree, const susy::wh::HftFiller::WeightVariations &wv)
{
    if(tree) {
        tree->syst_BKGMETHODUP   = wv.qflipUp_;
        tree->syst_BKGMETHODDOWN = wv.qflipDo_;
        tree->syst_ETRIGREWUP   = wv.elTrigUp_;
        tree->syst_ETRIGREWDOWN = wv.elTrigDo_;
        tree->syst_MTRIGREWUP   = wv.muTrigUp_;
        tree->syst_MTRIGREWDOWN = wv.muTrigDo_;
        tree->syst_BJETUP   = wv.bTagUp_;
        tree->syst_BJETDOWN = wv.bTagDo_;
        tree->syst_XSUP   = wv.xsecUp_;
        tree->syst_XSDOWN = wv.xsecDo_;
    }
}
//----------------------------------------------------------
bool HftFiller::fill(size_t systIndex, const susy::wh::kin::DilepVars &v, unsigned int run, unsigned int event)
{
    WeightVariations wv;
    return fill(systIndex, v, run, event, wv);
}
//----------------------------------------------------------
bool HftFiller::fill(size_t systIndex, const susy::wh::kin::DilepVars &v, unsigned int run, unsigned int event,
                     const WeightVariations &wv)
{
    bool someBytesWritten(false);
    if(HistFitterTree *t = m_hftTrees[systIndex]) {
        t->runNumber = run;
        t->eventNumber = event;
        assignDilepVars(t, v);
        assignWeightVars(t, wv);
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
bool HftFiller::determineXsecUncertainty(const int dsid)
{
    XsecUncertainty xsecUnc;
    bool groupFound = xsecUnc.determineGroup(dsid);
    if(!groupFound)
        cout<<"HftFiller::determineXsecUncertainty : cannot determine group for dsid "<<dsid<<endl
            <<"  will use "<<xsecUnc.str()<<endl;
    cout<<"HftFiller::determineXsecUncertainty : using "<<xsecUnc.str()<<" for dsid "<<dsid<<endl;
    xsecRelativeUncertainty_ = xsecUnc.fractionalUncertainty();
    return groupFound;
}
//----------------------------------------------------------
