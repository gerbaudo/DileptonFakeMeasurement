#ifndef SusySelection_h
#define SusySelection_h

//////////////////////////////////////////////////////////
// General script to implement basic selection with all //
// signal region cut methods.                           //
//////////////////////////////////////////////////////////

// Common Packages

// Root Packages
#include "TTree.h"

// Susy Common
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyDefs.h"

#include "SUSYTools/SUSYObjDef.h"

#include "SusyXSReader/XSReader.h"

#include <fstream>

enum WeightType {
  WT_Raw = 0,   // raw counts
  WT_Evt,       // include gen weight (from ntuple or xsreader)
  WT_PU,        // include pileup weight
  WT_LSF,       // include lepton scale factor
  WT_Btag,      // include b-tag scale factor
  WT_Trig,      // include trigger weight
  WT_All,       // include all weights above
  WT_N
};

enum WH_SR {
  WH_SRSS1=0,
  WH_SRSS2,
  WH_SRSS3,
  WH_SRSS4,
  WH_SRN
};

// fw decl
class chargeFlip;

class SusySelection : public SusyNtAna
{
 public:
  typedef       LeptonVector  vl_t;  //!< just to make some decl shorter
  typedef const LeptonVector cvl_t;  //!< just to make some decl shorter
  typedef const TauVector    cvt_t;  //!< just to make some decl shorter
  typedef const JetVector    cvj_t;  //!< just to make some decl shorter
  typedef const Met          cmet_t; //!< just to make some decl shorter
  struct WeightComponents {
    WeightComponents(): susynt(1), gen(1), pileup(1), norm(1),lepSf(1), btag(1), trigger(1) {}
    float product() const { return susynt * lepSf * btag * trigger; }
    void reset() { susynt = gen = pileup = norm = lepSf = btag = trigger = 1.0; }
    float susynt; // from SusyNtTools::getEventWeight: includes gen, pu, xs, lumi, sumw
    float gen, pileup, norm; // breakdown of the above; norm is xs*lumi/sumw
    float lepSf, btag, trigger; // factors that we compute, not from upstream
  };
 public:
    SusySelection();
    virtual ~SusySelection(){};

    ofstream out;

    virtual void    Begin(TTree *tree); //!< called before looping on entries
    virtual void    Terminate(); //!< called after looping is finished
    virtual Bool_t  Process(Long64_t entry); //!< called at each event
    virtual void dumpEventCounters();
    bool selectEvent(); //!< event selection  based on event qtities (mostly...)
    // Signal regions
    bool passSR6base(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR7base(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR8base(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR9base(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR6(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR7(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR8(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    bool passSR9(cvl_t& leptons, cvj_t& jets, const Met* met, bool count=false);
    // std SR7 has at least 2jets + the requirements below
    // (but no counters, just so that the fit on one line)
    bool passSR7ge2j     (cvl_t& l, cvj_t& j, const Met* m) { return passSR7base(l,j,m) && passge2Jet(j); }
    bool passSR7ge3j     (cvl_t& l, cvj_t& j, const Met* m) { return passSR7base(l,j,m) && passge3Jet(j); }
    bool passSR7eq2jNfv  (cvl_t& l, cvj_t& j, const Met* m) { return passSR7base(l,j,m) && passeq2JetWoutFwVeto(j); }
    bool passSR7ge2jNfv  (cvl_t& l, cvj_t& j, const Met* m) { return passSR7base(l,j,m) && passge2JetWoutFwVeto(j); }
    bool passSR7Nj       (cvl_t& l, cvj_t& j, const Met* m) { return passSR7base(l,j,m) && passNj(j); }
    bool passSR7NjZttVeto(cvl_t& l, cvj_t& j, const Met* m) { return passSR7Nj(l,j,m)   && passZtautauVeto(l,j,m); }
    bool passSR7NjPtTot  (cvl_t& l, cvj_t& j, const Met* m) { return passSR7Nj(l,j,m)   && passPtTot(l,j,m); }
    bool passSR7NjMll    (cvl_t& l, cvj_t& j, const Met* m) { return passSR7Nj(l,j,m)   && passMllMax(l); }
    bool passSrSsBase();
    bool passSrSs(const WH_SR signalRegion,
                  vl_t &l, cvt_t &t, cvj_t &j, const Met* m);
    // Cut methods
    bool passHfor();
    bool passTrig2L(const LeptonVector& leptons);
    bool passTrig2LMatch(const LeptonVector& leptons);
    bool passTrig2LwithMatch(const LeptonVector& leptons);
    bool sameFlavor(const LeptonVector& leptons);
    bool oppositeFlavor(const LeptonVector& leptons);
    bool sameSign(const LeptonVector& leptons);
    //! for the SS selection we want to accept OS events with an el, for Qflip
    /*!  For MC events that are OS, but that we want to consider as SS
      with some charge-flip probability, the 4-mom of the q-flipped
      electron is modified if update4mom==true.
     */
    bool sameSignOrQflip(LeptonVector &leptons, Met &met,
                         const DiLepEvtType eventType,
                         bool update4mom, bool isMC);
    bool oppositeSign(const LeptonVector& leptons);
    bool passHtautauVeto(int hdecay);
    // Signal Region Cuts
    bool passJetVeto(const JetVector& jets);
    bool passMetRelMin(const Met *met, cvl_t &leptons, cvj_t &jets, float minVal);
    bool passbJetVeto(const JetVector& jets);
    bool passfJetVeto(const JetVector& jets);
    bool passge1Jet(const JetVector& jets);
    bool passge2Jet(const JetVector& jets);
    bool passge3Jet(const JetVector& jets);
    bool passeq2Jet(const JetVector& jets);
    bool passge2JetWoutFwVeto(const JetVector& jets);
    bool passeq2JetWoutFwVeto(const JetVector& jets);
    bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut);
    bool passMtLlMetMin(const LeptonVector& leptons, const Met* met, float minVal=50.0);
    bool passMtMinlmetMin(const LeptonVector& leptons, const Met* met, float minVal=50.0);
    bool passMT2(const LeptonVector& leptons, const Met* met, float cut);
    bool passHtMin(const cvl_t& l, cvj_t &j, const Met* met, float minVal);
    bool passNlepMin(const LeptonVector &leptons, size_t minVal);
    bool passNj(const JetVector& jets, int minNj=2, int maxNj=3);
    bool passZtautauVeto(cvl_t& l, cvj_t& j, const Met* m, float widthZpeak=40.0);
    bool passZllVeto(cvl_t& l, float mllLo, float mllHi);
    bool pass2LepPt(cvl_t& l, float minPt0, float minPt1); //!< assume pt0 > pt1
    bool passPtllMin(cvl_t& l, float minPt=50.0);
    bool passPtTot(cvl_t& l, cvj_t& j, const Met* m, float maxPtTot=50.0);
    bool passMllMax(const LeptonVector& leptons, float maxMll=80.0);
    bool passMllMin(const LeptonVector& leptons, float minVal);
    bool passDrllMax(const LeptonVector& leptons, float maxDr=2.0);

    // Idendification methods
    bool isRealLepton(const Lepton* lep);
    bool isFakeLepton(const Lepton* lep);
    bool isConvLepton(const Lepton* lep);
    bool isHFLepton(const Lepton* lep);
    bool isLFLepton(const Lepton* lep);
    bool isTrueDilepton(const LeptonVector &leptons);
    bool passMuonRelIso(const LeptonVector &leptons, float maxVal);
    bool passEleD0S(const LeptonVector &leptons, float maxVal);

    //! method that should be used to fill the histos
    float getEvtWeight(const LeptonVector &leptons, bool includeBTag=false, bool includeTrig=true);
    void setUseXsReader(bool val){ m_useXsReader = val; };
    void setUseMCTrig(bool useMCTrig){ m_useMCTrig = useMCTrig; };
    //! increment the counters for the all event weight types
    void increment(float counters[], bool includeLepSF=false, bool includeBtag=false);
    static float computeMt2(const TLorentzVector &l0, const TLorentzVector &l1,
                            const TLorentzVector &met);
    float computeChargeFlipProb(LeptonVector &leptons, Met &met,
                                uint systematic, bool update4mom);
    static int pdgIdFromLep(const Lepton *l);
 protected:
    //! call SusyNtAna::getEventWeight, replacing the ntuple xsec with the one from the reader
    float computeEventWeightXsFromReader(float lumi);
    float getXsFromReader();     //!< cache xsec from xsreader
    float getBTagWeight(const Event* evt);
    float getPythiaBbCcScaleFactor(uint datasetId, const LeptonVector &leptons) const;
    float getTriggerWeight2Lep(const LeptonVector &leptons);
    float getLeptonEff2Lep(const LeptonVector &leptons) const;
    void resetAllCounters();
    void initChargeFlipTool();
    ClassDef(SusySelection, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj
    XSReader* m_xsReader;

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger, i.e. toggle the matching in DilTrigLogic::passDil*()
    chargeFlip*         m_chargeFlip;   //!< tool providing the electron charge flip probability
    float               m_w;            // mc weight
    bool                m_useXsReader;  // use SusyXSReader to get the xsec for normalization
    float               m_xsFromReader; // cached xsec from reader
    DiLepEvtType        m_ET;           // Dilepton event type to store cf
    float               m_qflipProb;     //! charge flip probability
    TLorentzVector      m_unsmeared_lv0; //! cached lepton LV before charge-flip smearing
    TLorentzVector      m_unsmeared_lv1; //! see above
    Met                 m_unsmeared_met; //! cached met before charge-flip smearing
    WeightComponents    m_weightComponents;
    // Event counters
    float n_readin          [WT_N]; // [weight type]
    float n_pass_Grl        [WT_N];
    float n_pass_LarErr     [WT_N];
    float n_pass_TileErr    [WT_N];
    float n_pass_TTCVeto    [WT_N];
    float n_pass_GoodVtx    [WT_N];
    float n_pass_TileTrip   [WT_N];
    float n_pass_LAr        [WT_N];
    float n_pass_BadJet     [WT_N];
    float n_pass_FEBCut     [WT_N];
    float n_pass_BadMuon    [WT_N];
    float n_pass_Cosmic     [WT_N];
    float n_pass_hfor       [WT_N];
    float n_pass_HttVeto    [WT_N];
    float n_pass_ge2l       [WT_N];
    float n_pass_eq2l       [WT_N];
    float n_pass_mll        [WT_N];
    float n_pass_signalLep  [WT_N];

    // SR6 counts
    float                n_pass_SR6sign[ET_N][WT_N];
    float                n_pass_SR6flav[ET_N][WT_N];
    float                n_pass_SR6eq2j[ET_N][WT_N];
    float                n_pass_SR6eq2jNfv[ET_N][WT_N];
    float                n_pass_SR6ge1j[ET_N][WT_N];
    float                n_pass_SR6ge2j[ET_N][WT_N];
    float                n_pass_SR6ge2jNfv[ET_N][WT_N];
    float                n_pass_SR6metr[ET_N][WT_N];
    float                n_pass_SR6DrllMax     [ET_N][WT_N];
    float                n_pass_SR6PtllMin     [ET_N][WT_N];
    float                n_pass_SR6MllMax      [ET_N][WT_N];
    float                n_pass_SR6METRel      [ET_N][WT_N];
    float                n_pass_SR6MtLlmetMin  [ET_N][WT_N];
    float                n_pass_SR6MtMinlmetMin[ET_N][WT_N];
    float                n_pass_SR6ZtautauVeto [ET_N][WT_N];
    float                n_pass_SR6            [ET_N][WT_N];
    // SR7 counts
    float                n_pass_SR7sign[ET_N][WT_N];
    float                n_pass_SR7flav[ET_N][WT_N];
    float                n_pass_SR7eq2j[ET_N][WT_N];
    float                n_pass_SR7eq2jNfv[ET_N][WT_N];
    float                n_pass_SR7ge1j[ET_N][WT_N];
    float                n_pass_SR7ge2j[ET_N][WT_N];
    float                n_pass_SR7ge2jNfv[ET_N][WT_N];
    float                n_pass_SR7metr[ET_N][WT_N];
    float                n_pass_SR7DrllMax     [ET_N][WT_N];
    float                n_pass_SR7PtllMin     [ET_N][WT_N];
    float                n_pass_SR7MllMax      [ET_N][WT_N];
    float                n_pass_SR7METRel      [ET_N][WT_N];
    float                n_pass_SR7MtLlmetMin  [ET_N][WT_N];
    float                n_pass_SR7MtMinlmetMin[ET_N][WT_N];
    float                n_pass_SR7ZtautauVeto [ET_N][WT_N];
    float                n_pass_SR7[ET_N][WT_N];
    // SR8 counts
    float                n_pass_SR8sign[ET_N][WT_N];
    float                n_pass_SR8flav[ET_N][WT_N];
    float                n_pass_SR8eq2j[ET_N][WT_N];
    float                n_pass_SR8eq2jNfv[ET_N][WT_N];
    float                n_pass_SR8ge1j[ET_N][WT_N];
    float                n_pass_SR8ge2j[ET_N][WT_N];
    float                n_pass_SR8ge2jNfv[ET_N][WT_N];
    float                n_pass_SR8metr[ET_N][WT_N];
    float                n_pass_SR8[ET_N][WT_N];
    // SR9 counts
    float                n_pass_SR9sign[ET_N][WT_N];
    float                n_pass_SR9flav[ET_N][WT_N];
    float                n_pass_SR9eq2j[ET_N][WT_N];
    float                n_pass_SR9eq2jNfv[ET_N][WT_N];
    float                n_pass_SR9ge1j[ET_N][WT_N];
    float                n_pass_SR9ge2j[ET_N][WT_N];
    float                n_pass_SR9ge2jNfv[ET_N][WT_N];
    float                n_pass_SR9metr[ET_N][WT_N];
    float                n_pass_SR9[ET_N][WT_N];
    // SS counts
    float n_pass_flavor     [ET_N][WT_N]; // [event type][weight type]
    float n_pass_os         [ET_N][WT_N];
    float n_pass_ss         [ET_N][WT_N];
    float n_pass_tr2L       [ET_N][WT_N];
    float n_pass_tr2LMatch  [ET_N][WT_N];
    float n_pass_mcTrue2l   [ET_N][WT_N];
    float n_pass_category   [ET_N][WT_N];
    float n_pass_nSigLep    [ET_N][WT_N];
    float n_pass_tauVeto    [ET_N][WT_N];
    float n_pass_mllMin     [ET_N][WT_N];
    float n_pass_muIso      [ET_N][WT_N];
    float n_pass_elD0Sig    [ET_N][WT_N];
    float n_pass_fjVeto     [ET_N][WT_N];
    float n_pass_bjVeto     [ET_N][WT_N];
    float n_pass_ge1j       [ET_N][WT_N];
    float n_pass_lepPt      [ET_N][WT_N];
    float n_pass_mllZveto   [ET_N][WT_N];
    float n_pass_mWwt       [ET_N][WT_N];
    float n_pass_ht         [ET_N][WT_N];
    float n_pass_metRel     [ET_N][WT_N];
};

#endif // SusySelection_h
