// emacs -*- C++ -*-
#ifndef SusySelection_h
#define SusySelection_h

//////////////////////////////////////////////////////////
// General script to implement basic selection with all //
// signal region cut methods.                           //
//////////////////////////////////////////////////////////

// Common Packages

// Root Packages
#include "TTree.h"

#include "SusyTest0/TupleMaker.h"
// Susy Common
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyDefs.h"

#include "SusyXSReader/XSReader.h"
#include "SusyTest0/ProgressPrinter.h"
#include "SusyTest0/SsPassFlags.h"
#include "SusyTest0/DileptonChannel.h"

#include <fstream>

enum WeightTypes {
  kRaw = 0,   // raw counts
  kEvt,       // include gen weight (from ntuple or xsreader)
  kPU,        // include pileup weight
  kLSF,       // include lepton scale factor
  kBtag,      // include b-tag scale factor
  kTrig,      // include trigger weight
  kAll,       // include all weights above
  kWeightTypesN
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

namespace susy {
namespace wh {
namespace kin {
class DilepVars;
}
}
}

class SusySelection : public SusyNtAna
{
 public:
  typedef       LeptonVector  vl_t;  //!< just to make some decl shorter
  typedef const LeptonVector cvl_t;  //!< just to make some decl shorter
  typedef const TauVector    cvt_t;  //!< just to make some decl shorter
  typedef const JetVector    cvj_t;  //!< just to make some decl shorter
  typedef const Met          cmet_t; //!< just to make some decl shorter
  struct WeightComponents {
    WeightComponents() { reset(); }
    double product() const { return susynt * lepSf * btag * trigger * qflip * fake; }
    void reset() { susynt = gen = pileup = norm = lepSf = btag = trigger = qflip = fake = 1.0; }
    double susynt; // from SusyNtTools::getEventWeight: includes gen, pu, xs, lumi, sumw
    double gen, pileup, norm; // breakdown of the above; norm is xs*lumi/sumw
    double lepSf, btag, trigger, qflip, fake; // factors that we compute, not from upstream
    std::string str() const;
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
    bool passSrSsBase();
    SsPassFlags passSrSs(const WH_SR signalRegion,
                         vl_t &l, cvt_t &t, cvj_t &j, const Met* m, bool allowQflip);
    // Cut methods
    bool passHfor() { return passHfor(nt); }
    static bool passHfor(Susy::SusyNtObject &nto);
    static bool passTrig2L         (const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt);
    static bool passTrig2LMatch    (const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt);
    static bool passTrig2LwithMatch(const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt);
    bool passTrig2L(const LeptonVector& leptons) { return passTrig2L(leptons, m_trigObj, m_met->Et, nt.evt()); }
    bool passTrig2LMatch(const LeptonVector& leptons) { return passTrig2LMatch(leptons, m_trigObj, m_met->Et, nt.evt()); }
    bool passTrig2LwithMatch(const LeptonVector& leptons) { return passTrig2LwithMatch(leptons, m_trigObj, m_met->Et, nt.evt()); }
    //! for the SS selection we want to accept OS events with an el, for Qflip
    /*!  For MC events that are OS, but that we want to consider as SS
      with some charge-flip probability, the 4-mom of the q-flipped
      electron is modified if update4mom==true.
     */
    bool sameSignOrQflip(LeptonVector &leptons, Met &met,
                         const DiLepEvtType eventType,
                         bool update4mom, bool isMC);
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
    bool passNj(const JetVector& jets, int minNj=2, int maxNj=3);
    bool passMuonRelIso(const LeptonVector &leptons, float maxVal);
    static bool passEwkSs     (const LeptonVector& leptons, const JetVector& jets, const Met* met);
    static bool passEwkSsLoose(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    static bool passCrWhZVfakeEe(const susy::wh::kin::DilepVars &v);
    static bool passCrWhZVfakeEm(const susy::wh::kin::DilepVars &v);
    static bool passCrWhfakeEm  (const susy::wh::kin::DilepVars &v);
    static bool passCrWhZVMm    (const susy::wh::kin::DilepVars &v);
    static bool passCrWhfakeMm  (const susy::wh::kin::DilepVars &v);
    // just wrapping the channel-specific functions above
    static bool passCrWhZVfake(const susy::wh::kin::DilepVars &v);
    static bool passCrWhfake    (const susy::wh::kin::DilepVars &v);
    static bool passCrWhZV    (const susy::wh::kin::DilepVars &v);

    //
    static bool passSrWh1j(const susy::wh::kin::DilepVars &v);
    static bool passSrWh2j(const susy::wh::kin::DilepVars &v);

    void setUseXsReader(bool val){ m_useXsReader = val; };
    void setUseMCTrig(bool useMCTrig){ m_useMCTrig = useMCTrig; };
    void setWriteNtuple(bool val) { m_writeTuple = val; };
    void setTupleFile(const std::string &name) { m_writeTuple = true; m_outTupleFile = name; };
    //! increment the counters for the all event weight types
    void increment(float counters[], const WeightComponents &wc);
    float computeChargeFlipProb(LeptonVector &leptons, Met &met,
                                uint systematic, bool update4mom);
    //! any electron or muon, before pt cuts, before overlap removal
    static vl_t getAnyElOrMu(SusyNtObject &susyNt/*, SusyNtSys sys*/);
    static susy::wh::Chan getChan(const LeptonVector& leps); //!< compute lepton channel
    static SsPassFlags assignNjetFlags(const JetVector& jets, SsPassFlags f);
    //! determine whether a third lepton makes a Z candidate with a signal lepton
    static bool passThirdLeptonVeto(const Susy::Lepton* l0, const Susy::Lepton* l1, const LeptonVector& otherLeptons, bool verbose=false);
    //! ugly hack function : utils::filter seems not to work properly with SusyNtTools::isCentralLightJet...
    static JetVector filterClJets(const JetVector &jets);

 protected:
    //! call SusyNtAna::getEventWeight, replacing the ntuple xsec with the one from the reader
    float computeEventWeightXsFromReader(float lumi);
    float getXsFromReader();     //!< cache xsec from xsreader
    float getBTagWeight(cvj_t& jets, const Event* evt);
    float getTriggerWeight2Lep(const LeptonVector &leptons);
    void resetAllCounters();
    void initChargeFlipTool();
    void cacheStaticWeightComponents(); //! cache those weight components that do not depend on sel
    void computeNonStaticWeightComponents(cvl_t& leptons, cvj_t& jets);
    ClassDef(SusySelection, 2);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj
    XSReader* m_xsReader;
    susy::wh::TupleMaker m_tupleMaker;
    bool m_writeTuple;
    bool m_debugThisEvent;
    std::string m_outTupleFile;

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger, i.e. toggle the matching in DilTrigLogic::passDil*()
    chargeFlip*         m_chargeFlip;   //!< tool providing the electron charge flip probability
    float               m_w;            // mc weight
    bool                m_useXsReader;  // use SusyXSReader to get the xsec for normalization
    float               m_xsFromReader; // cached xsec from reader
    float               m_qflipProb;     //! charge flip probability
    TLorentzVector      m_unsmeared_lv0; //! cached lepton LV before charge-flip smearing
    TLorentzVector      m_unsmeared_lv1; //! see above
    Met                 m_unsmeared_met; //! cached met before charge-flip smearing
    WeightComponents    m_weightComponents;
    susy::ProgressPrinter m_printer;
    // Event counters
    float n_readin          [kWeightTypesN]; // [weight type]
    float n_pass_Grl        [kWeightTypesN];
    float n_pass_LarErr     [kWeightTypesN];
    float n_pass_TileErr    [kWeightTypesN];
    float n_pass_TTCVeto    [kWeightTypesN];
    float n_pass_GoodVtx    [kWeightTypesN];
    float n_pass_TileTrip   [kWeightTypesN];
    float n_pass_LAr        [kWeightTypesN];
    float n_pass_BadJet     [kWeightTypesN];
    float n_pass_FEBCut     [kWeightTypesN];
    float n_pass_BadMuon    [kWeightTypesN];
    float n_pass_Cosmic     [kWeightTypesN];
    float n_pass_hfor       [kWeightTypesN];
    float n_pass_HttVeto    [kWeightTypesN];
    float n_pass_ge2l       [kWeightTypesN];
    float n_pass_eq2l       [kWeightTypesN];
    float n_pass_mll        [kWeightTypesN];
    float n_pass_signalLep  [kWeightTypesN];

    // SR6 counts
    float                n_pass_SR6sign[ET_N][kWeightTypesN];
    float                n_pass_SR6flav[ET_N][kWeightTypesN];
    float                n_pass_SR6eq2j[ET_N][kWeightTypesN];
    float                n_pass_SR6eq2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR6ge1j[ET_N][kWeightTypesN];
    float                n_pass_SR6ge2j[ET_N][kWeightTypesN];
    float                n_pass_SR6ge2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR6metr[ET_N][kWeightTypesN];
    float                n_pass_SR6DrllMax     [ET_N][kWeightTypesN];
    float                n_pass_SR6PtllMin     [ET_N][kWeightTypesN];
    float                n_pass_SR6MllMax      [ET_N][kWeightTypesN];
    float                n_pass_SR6METRel      [ET_N][kWeightTypesN];
    float                n_pass_SR6MtLlmetMin  [ET_N][kWeightTypesN];
    float                n_pass_SR6MtMinlmetMin[ET_N][kWeightTypesN];
    float                n_pass_SR6ZtautauVeto [ET_N][kWeightTypesN];
    float                n_pass_SR6            [ET_N][kWeightTypesN];
    // SR7 counts
    float                n_pass_SR7sign[ET_N][kWeightTypesN];
    float                n_pass_SR7flav[ET_N][kWeightTypesN];
    float                n_pass_SR7eq2j[ET_N][kWeightTypesN];
    float                n_pass_SR7eq2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR7ge1j[ET_N][kWeightTypesN];
    float                n_pass_SR7ge2j[ET_N][kWeightTypesN];
    float                n_pass_SR7ge2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR7metr[ET_N][kWeightTypesN];
    float                n_pass_SR7DrllMax     [ET_N][kWeightTypesN];
    float                n_pass_SR7PtllMin     [ET_N][kWeightTypesN];
    float                n_pass_SR7MllMax      [ET_N][kWeightTypesN];
    float                n_pass_SR7METRel      [ET_N][kWeightTypesN];
    float                n_pass_SR7MtLlmetMin  [ET_N][kWeightTypesN];
    float                n_pass_SR7MtMinlmetMin[ET_N][kWeightTypesN];
    float                n_pass_SR7ZtautauVeto [ET_N][kWeightTypesN];
    float                n_pass_SR7[ET_N][kWeightTypesN];
    // SR8 counts
    float                n_pass_SR8sign[ET_N][kWeightTypesN];
    float                n_pass_SR8flav[ET_N][kWeightTypesN];
    float                n_pass_SR8eq2j[ET_N][kWeightTypesN];
    float                n_pass_SR8eq2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR8ge1j[ET_N][kWeightTypesN];
    float                n_pass_SR8ge2j[ET_N][kWeightTypesN];
    float                n_pass_SR8ge2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR8metr[ET_N][kWeightTypesN];
    float                n_pass_SR8[ET_N][kWeightTypesN];
    // SR9 counts
    float                n_pass_SR9sign[ET_N][kWeightTypesN];
    float                n_pass_SR9flav[ET_N][kWeightTypesN];
    float                n_pass_SR9eq2j[ET_N][kWeightTypesN];
    float                n_pass_SR9eq2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR9ge1j[ET_N][kWeightTypesN];
    float                n_pass_SR9ge2j[ET_N][kWeightTypesN];
    float                n_pass_SR9ge2jNfv[ET_N][kWeightTypesN];
    float                n_pass_SR9metr[ET_N][kWeightTypesN];
    float                n_pass_SR9[ET_N][kWeightTypesN];
    // SS counts
    float n_pass_flavor     [ET_N][kWeightTypesN]; // [event type][weight type]
    float n_pass_os         [ET_N][kWeightTypesN];
    float n_pass_ss         [ET_N][kWeightTypesN];
    float n_pass_tr2L       [ET_N][kWeightTypesN];
    float n_pass_tr2LMatch  [ET_N][kWeightTypesN];
    float n_pass_mcTrue2l   [ET_N][kWeightTypesN];
    float n_pass_category   [ET_N][kWeightTypesN];
    float n_pass_nSigLep    [ET_N][kWeightTypesN];
    float n_pass_tauVeto    [ET_N][kWeightTypesN];
    float n_pass_mllMin     [ET_N][kWeightTypesN];
    float n_pass_muIso      [ET_N][kWeightTypesN];
    float n_pass_elD0Sig    [ET_N][kWeightTypesN];
    float n_pass_fjVeto     [ET_N][kWeightTypesN];
    float n_pass_bjVeto     [ET_N][kWeightTypesN];
    float n_pass_ge1j       [ET_N][kWeightTypesN];
    float n_pass_eq1j       [ET_N][kWeightTypesN];
    float n_pass_ge2j       [ET_N][kWeightTypesN];
    float n_pass_lepPt      [ET_N][kWeightTypesN];
    float n_pass_mllZveto   [ET_N][kWeightTypesN];
    float n_pass_mWwt       [ET_N][kWeightTypesN];
    float n_pass_ht         [ET_N][kWeightTypesN];
    float n_pass_metRel     [ET_N][kWeightTypesN];
    float n_pass_3rdLep     [ET_N][kWeightTypesN];
    float n_pass_eq1jlepPt      [ET_N][kWeightTypesN];
    float n_pass_eq1jmllZveto   [ET_N][kWeightTypesN];
    float n_pass_eq1jmWwt       [ET_N][kWeightTypesN];
    float n_pass_eq1jht         [ET_N][kWeightTypesN];
    float n_pass_eq1jmetRel     [ET_N][kWeightTypesN];
    float n_pass_eq1j3rdLep     [ET_N][kWeightTypesN];
    float n_pass_ge2jlepPt      [ET_N][kWeightTypesN];
    float n_pass_ge2jmllZveto   [ET_N][kWeightTypesN];
    float n_pass_ge2jmWwt       [ET_N][kWeightTypesN];
    float n_pass_ge2jht         [ET_N][kWeightTypesN];
    float n_pass_ge2jmetRel     [ET_N][kWeightTypesN];
    float n_pass_ge2j3rdLep     [ET_N][kWeightTypesN];
};

#endif // SusySelection_h
