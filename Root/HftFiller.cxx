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
        susy::wh::HftFiller::WeightVariations w(wv); // Hft wants up above 1, do below; need to modify WeightVariations
        struct SwapFunc{ float tmp_; void operator()(float &a, float &b) {tmp_=a; a=b; b=tmp_;} } swap;
        if(w.qflipUp_  < w.qflipDo_ ) swap(w.qflipUp_,  w.qflipDo_);
        if(w.elTrigUp_ < w.elTrigDo_) swap(w.elTrigUp_, w.elTrigDo_);
        if(w.muTrigUp_ < w.muTrigDo_) swap(w.muTrigUp_, w.muTrigDo_);
        if(w.elEffUp_  < w.elEffDo_)  swap(w.elEffUp_,  w.elEffDo_);
        if(w.muEffUp_  < w.muEffDo_)  swap(w.muEffUp_,  w.muEffDo_);
        if(w.bTagUp_   < w.bTagDo_  ) swap(w.bTagUp_,   w.bTagDo_);
        if(w.cTagUp_   < w.cTagDo_  ) swap(w.cTagUp_,   w.cTagDo_);
        if(w.lTagUp_   < w.lTagDo_  ) swap(w.lTagUp_,   w.lTagDo_);
        if(w.xsecUp_   < w.xsecDo_  ) swap(w.xsecUp_,   w.xsecDo_);
        tree->syst_BKGMETHODUP   = w.qflipUp_;
        tree->syst_BKGMETHODDOWN = w.qflipDo_;
        tree->syst_ETRIGREWUP    = w.elTrigUp_;
        tree->syst_ETRIGREWDOWN  = w.elTrigDo_;
        tree->syst_MTRIGREWUP    = w.muTrigUp_;
        tree->syst_MTRIGREWDOWN  = w.muTrigDo_;
        tree->syst_ESFUP         = w.elEffUp_;
        tree->syst_ESFDOWN       = w.elEffDo_;
        tree->syst_MEFFUP        = w.muEffUp_;
        tree->syst_MEFFDOWN      = w.muEffDo_;
        tree->syst_BJETUP        = w.bTagUp_;
        tree->syst_BJETDOWN      = w.bTagDo_;
        tree->syst_CJETUP        = w.cTagUp_;
        tree->syst_CJETDOWN      = w.cTagDo_;
        tree->syst_BMISTAGUP     = w.lTagUp_;
        tree->syst_BMISTAGDOWN   = w.lTagDo_;
        tree->syst_XSUP          = w.xsecUp_;
        tree->syst_XSDOWN        = w.xsecDo_;
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
