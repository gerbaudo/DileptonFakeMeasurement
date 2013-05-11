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

enum WeightType
{
  WT_Raw = 0,   // weight = 1;
  WT_Evt,       // weight = gen weight
  WT_PU,        // weight = pileup weight
  WT_PU1fb,     // weight = pileup weight for 1/fb
  WT_LSF,       // weight = lepton SF
  WT_Btag,      // weight = btag
  WT_Trig,      // Trigger weight
  WT_AllAB3,    // all weights for A-B3
  WT_AllAE,     // all weights for A-E
  WT_N
};

class SusySelection : public SusyNtAna
{

 public:
  typedef const LeptonVector cvl_t;  //!< just to make some decl shorter
  typedef const JetVector    cvj_t;  //!< just to make some decl shorter
  typedef const Met          cmet_t; //!< just to make some decl shorter

 public:
    SusySelection();
    virtual ~SusySelection(){};

    ofstream out;

    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);

    // Full event selection. Specify which leptons to use.
    bool selectEvent(bool doMll=true);
    bool selectAnaEvent(cvl_t& leptons, cvl_t& baseLeptons);

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

    // Cut methods
    bool passHfor();
    bool passNLepCut(const LeptonVector& leptons);
    bool passNBaseLepCut(const LeptonVector& baseLeptons);
    bool passTrigger(const LeptonVector& leptons);
    bool sameFlavor(const LeptonVector& leptons);
    bool oppositeFlavor(const LeptonVector& leptons);
    bool sameSign(const LeptonVector& leptons);
    bool oppositeSign(const LeptonVector& leptons);
    bool passMll(const LeptonVector& leptons, float mll = 20); // this one (by Matt) increments
    bool passHtautauVeto(int hdecay);
    // Signal Region Cuts
    bool passJetVeto(const JetVector& jets);
    bool passZVeto(const LeptonVector& leptons, float Zlow = 81.2, float Zhigh = 101.2);
    bool passMETRel(const Met *met, const LeptonVector& leptons, const JetVector& jets, float maxMet = 50.0);

    bool passbJetVeto(const JetVector& jets);
    bool passge1Jet(const JetVector& jets);
    bool passge2Jet(const JetVector& jets);
    bool passge3Jet(const JetVector& jets);
    bool passeq2Jet(const JetVector& jets);
    bool passge2JetWoutFwVeto(const JetVector& jets);
    bool passeq2JetWoutFwVeto(const JetVector& jets);
    bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut);
    bool passMtLlmetMin(const LeptonVector& leptons, const Met* met, float minVal=50.0);
    bool passMtMinlmetMin(const LeptonVector& leptons, const Met* met, float minVal=50.0);
    bool passMT2(const LeptonVector& leptons, const Met* met, float cut);
    bool passNj(const JetVector& jets, int minNj=2, int maxNj=3);
    bool passZtautauVeto(cvl_t& l, cvj_t& j, const Met* m, float widthZpeak=40.0);
    bool passPtllMin(cvl_t& l, float minPt=50.0);
    bool passPtTot(cvl_t& l, cvj_t& j, const Met* m, float maxPtTot=50.0);
    bool passMllMax(const LeptonVector& leptons, float maxMll=80.0);
    bool passDrllMax(const LeptonVector& leptons, float maxDr=2.0);

    // Idendification methods
    bool isRealLepton(const Lepton* lep);
    bool isFakeLepton(const Lepton* lep);
    bool isConvLepton(const Lepton* lep);
    bool isHFLepton(const Lepton* lep);
    bool isLFLepton(const Lepton* lep);
    bool isTrueDilepton(const LeptonVector &leptons);

    // Photon+jet MC methods
    bool passCheckMC(int mcRunNumber, float pt);
    float getPhotonXS(int mcRunNumber);

    // Dump cutflow - if derived class uses different cut ordering,
    // override this method
    virtual void dumpEventCounters();
    void dumpInterestingEvents(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    void dumpTrigFlag(uint flag);

    // debug check
    bool debugEvent();
    void dumpPreObjects();
    void dumpJets();
    void checkSys();

    // Save file name for easy writing
    void setFileName(string f){ m_fileName = f; };
    string getFileName(){ return m_fileName; };
    // cache xsec from xsreader
    float getXsFromReader();
    float computeEventWeightXsFromReader(float lumi);
    // Get Btag weight
    float getEvtWeight(const LeptonVector &leptons, bool includeBTag=false, bool includeTrig=true,
		       bool doMediumpp=false);
    float getBTagWeight(const Event* evt);

    // Some controls
    void setUse1fb(bool use1fb){ m_do1fb = use1fb; };
    bool is1fb(){ return isPeriodAB3(nt.evt()->run); };
    void setUseAD(bool useAD){ m_doAD = useAD; };
    void setUseXsReader(bool val){ m_useXsReader = val; };
    bool isB3(){
      uint run = nt.evt()->run;
      return 203169 <= run && run <= 203195;
    };
    void setUseMCTrig(bool useMCTrig){ m_useMCTrig = useMCTrig; };

    // Method to increment the counters for the event weight types
    void increment(float flag[], bool includeLepSF=false, bool includeBtag=false){
      flag[WT_Raw]   += 1.0;
      flag[WT_Evt]   += nt.evt()->w;
      flag[WT_PU]    += nt.evt()->w * nt.evt()->wPileup;
      flag[WT_PU1fb] += nt.evt()->w * nt.evt()->wPileupAB3;
      flag[WT_LSF]   += (includeLepSF ?
			 nt.evt()->w * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF :
			 nt.evt()->w);
      float btag = includeBtag ? getBTagWeight(nt.evt()) : 1.0;
      flag[WT_Btag]  += nt.evt()->w * btag;

      float trig = m_baseLeptons.size() == 2 && nt.evt()->isMC ?
	m_trigObj->getTriggerWeight(m_baseLeptons,
				    nt.evt()->isMC,
				    m_met->Et,
				    m_signalJets2Lep.size(),
				    nt.evt()->nVtx,
				    NtSys_NOM) : 1;

      //cout<<"\tTrigger weight: "<<trig<<endl;
      flag[WT_Trig] += trig * nt.evt()->w;

      float all = getEventWeightAB3() * btag * trig;
      all = includeLepSF ? all * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF : all;
      flag[WT_AllAB3] += all;

      float allAE = getEventWeightFixed(nt.evt()->mcChannel,LUMI_A_E) * btag * trig;
      allAE = includeLepSF ? allAE * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF : allAE;
      flag[WT_AllAE] += allAE;
    };

    // Miscellaneous methods
    void printLep(const Lepton* lep);
    void printJet(const Jet* jet);
    static float computeMt2(const TLorentzVector &l0, const TLorentzVector &l1,
			    const TLorentzVector &met);

    ClassDef(SusySelection, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj
    XSReader* m_xsReader;

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger

    string              m_fileName;     // File name
    float               m_w;            // mc weight

    bool                m_do1fb;        // For get weight method
    bool                m_doAD;         // do weights for A-B
    bool                m_useXsReader;  // use SusyXSReader to get the xsec for normalization
    float               m_xsFromReader; // cached xsec from reader

    bool                m_dumpCounts;   // Flag to dump counters

    // Cut variables
    float                m_nLepMin;      // min leptons
    float                m_nLepMax;      // max leptons
    bool                m_cutNBaseLep;  // apply nLep cuts to baseline leptons as well as signal

    DiLepEvtType        m_ET;           // Dilepton event type to store cf

    // Event counters
    float n_readin          [WT_N];
    float n_pass_Grl        [WT_N];
    float n_pass_LarErr     [WT_N];
    float n_pass_TileErr    [WT_N];
    float n_pass_TTCVeto    [WT_N];
    float n_pass_GoodVtx    [WT_N];
    float n_pass_TileTrip   [WT_N];
    float n_pass_LAr        [WT_N];
    float n_pass_BadJet     [WT_N];
    float n_pass_BadMuon    [WT_N];
    float n_pass_Cosmic     [WT_N];
    float n_pass_HttVeto    [WT_N];
    float n_pass_atleast2Lep[WT_N];
    float n_pass_exactly2Lep[WT_N];
    float n_pass_signalLep  [WT_N];
    float n_pass_flavor     [ET_N][WT_N];
    float n_pass_mll        [ET_N][WT_N];
    float n_pass_os         [ET_N][WT_N];
    float n_pass_ss         [ET_N][WT_N];
    float n_pass_evtTrig    [ET_N][WT_N];
    float n_pass_trigMatch  [ET_N][WT_N];

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
};

#endif
