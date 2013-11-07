#include "SusyTest0/TupleMaker.h"

#include "SusyNtuple/SusyNt.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

#include <algorithm> // transform
#include <functional> // unary_function
#include <iostream>
#include <iterator> // back_inserter

using susy::wh::TupleMaker;
using susy::wh::FourMom;
using Susy::Lepton;
using Susy::Jet;
using std::cout;
using std::endl;
using std::string;

//----------------------------------------------------------
TupleMaker::TupleMaker(const std::string &outFilename, const std::string &treename):
    file_(0),
    tree_(0)
{
    init(outFilename, treename);
}
//----------------------------------------------------------
TupleMaker::~TupleMaker()
{
    close();
}
//----------------------------------------------------------
// util functions to convert Lepton, Jet -> FourMom
FourMom lepton2FourMom (const Lepton *l)
{    
    return (  l && l->isMu()  ? FourMom().setMu(*l)
            : l && l->isEle() ? FourMom().setEl(*l)
            : FourMom());
}
FourMom jet2FourMom (const Jet *j) { return (j ? FourMom().set4mom(*j) : FourMom()); }
//----------------------------------------------------------
bool TupleMaker::fill(const Lepton &l0, const Lepton &l1,
                      const LeptonVector &otherLeptons, const JetVector &jets)
{
    bool someBytesWritten(false);
    if(tree_) {
        l0.isMu() ? l0_.setMu(l0) : l0_.setEl(l0);
        l1.isMu() ? l1_.setMu(l1) : l1_.setEl(l1);
        jets_.clear();
        lowptLepts_.clear();
        const LeptonVector &olps = otherLeptons;
        std::transform(jets.begin(), jets.end(), std::back_inserter(jets_),       jet2FourMom);
        std::transform(olps.begin(), olps.end(), std::back_inserter(lowptLepts_), lepton2FourMom);
        someBytesWritten = (tree_->Fill()>0);
    }
    return someBytesWritten;
}
//----------------------------------------------------------
bool TupleMaker::init(const std::string &outFilename, const std::string &treename)
{
    return (initFile(outFilename) && initTree(treename));
}
//----------------------------------------------------------
bool TupleMaker::initFile(const std::string &outFilename)
{
    file_ = TFile::Open(outFilename.c_str(), "recreate");
    if(!file_) cout<<"TupleMaker::initFile('"<<outFilename<<"') : failed to create file"<<endl;
    return (file_ && file_->IsOpen());
}
//----------------------------------------------------------
bool TupleMaker::initTree(const std::string &treename)
{
    bool initialized(false);
    TDirectory *startingDir = gDirectory;
    if(file_) {
        file_->cd();
        string title("TupleMaker tree");
        tree_ = new TTree(treename.c_str(), title.c_str());
        initTreeBranches();
        initialized = true;
    } else {
        cout<<"TupleMaker::initTree: invalid file, failed to create tree"<<endl;
    }
    if(startingDir) startingDir->cd(); // root is easily confused by pwd; cd back to where we were
    return initialized;
}
//----------------------------------------------------------
bool TupleMaker::initTreeBranches()
{
    bool initialized(false);
    if(tree_) {
        tree_->Branch("l0",    &l0_);
        tree_->Branch("l1",    &l1_);
        tree_->Branch("met",   &met_);
        tree_->Branch("jets",  &jets_);
        tree_->Branch("lepts", &lowptLepts_);
        tree_->Branch("pars",  &eventPars_);
    } else {
        cout<<"TupleMaker::initTreeBranches : invalid tree, failed to init branches"<<endl;
    }
    return initialized;
}
//----------------------------------------------------------
bool TupleMaker::close()
{
    bool closed(false);
    if(file_) {
        file_->Write();
        file_->Close();
        file_->Delete();
        file_ = 0;
        closed = true;
    }
    return closed;
}
//----------------------------------------------------------
