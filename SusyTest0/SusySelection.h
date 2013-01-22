#ifndef SusySelection_h
#define SusySelection_h

//////////////////////////////////////////////////////////
// General script to implement basic selection with all //
// signal region cut methods.                           //
//////////////////////////////////////////////////////////

// Common Packages
#include "Mt2/mt2_bisect.h" // I don't want to recode this..
// #include "LeptonTruthTools/RecoTruthMatch.h" // DG not needed ?

// Root Packages
#include "TTree.h"

// Susy Common
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyDefs.h"

#include "SUSYTools/SUSYObjDef.h"

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
    bool selectAnaEvent(const LeptonVector& leptons, const LeptonVector& baseLeptons);

    // Signal regions
    bool passSR1(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR2(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR3(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR4(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR4b(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR5(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR6(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR7(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR8(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);
    bool passSR9(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false);

    // Some validation regions
    // These are under development!!!
    bool passVR1(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passVR2(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passVRTL(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passVR3(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passVR4(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passWWCR1(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passWWCR2(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passWWCR3(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passTOPCR(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // Bonus regions to look at same-sign muons
    bool passBR1(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passBR2(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passBR3(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passBR4(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // Cut methods
    bool passHfor();
    bool passNLepCut(const LeptonVector& leptons);
    bool passNBaseLepCut(const LeptonVector& baseLeptons);
    bool passTrigger(const LeptonVector& leptons);
    bool sameFlavor(const LeptonVector& leptons);
    bool oppositeFlavor(const LeptonVector& leptons);
    bool sameSign(const LeptonVector& leptons);
    bool oppositeSign(const LeptonVector& leptons);
    bool passMll(const LeptonVector& leptons, float mll = 20);

    // Signal Region Cuts
    bool passJetVeto(const JetVector& jets);
    bool passZVeto(const LeptonVector& leptons, float Zlow = 81.2, float Zhigh = 101.2);
    bool passMETRel(const Met *met, const LeptonVector& leptons,
		    const JetVector& jets, float maxMet = 100);
    bool passbJetVeto(const JetVector& jets);
    bool passge2Jet(const JetVector& jets);
    bool passeq2Jet(const JetVector& jets);
    bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut);
    bool passMT2(const LeptonVector& leptons, const Met* met, float cut);

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

    // Get Btag weight
    float getEvtWeight(const LeptonVector &leptons, bool includeBTag=false, bool includeTrig=true,
		       bool doMediumpp=false);
    float getBTagWeight(const Event* evt);

    // Some controls
    void setUse1fb(bool use1fb){ m_do1fb = use1fb; };
    bool is1fb(){ return isPeriodAB3(nt.evt()->run); };
    void setUseAD(bool useAD){ m_doAD = useAD; };
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
	m_trigObj->getTriggerWeight(m_baseLeptons,nt.evt()->isMC,NtSys_NOM) : 1;
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

    ClassDef(SusySelection, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger

    string              m_fileName;     // File name
    float               m_w;            // mc weight

    bool                m_do1fb;        // For get weight method
    bool                m_doAD;         // do weights for A-B

    bool                m_dumpCounts;   // Flag to dump counters

    // Cut variables
    float                m_nLepMin;      // min leptons
    float                m_nLepMax;      // max leptons
    bool                m_cutNBaseLep;  // apply nLep cuts to baseline leptons as well as signal

    DiLepEvtType        m_ET;           // Dilepton event type to store cf

    // Event counters
    float                n_readin[WT_N];
    float                n_pass_LAr[WT_N];
    float                n_pass_BadJet[WT_N];
    float                n_pass_BadMuon[WT_N];
    float                n_pass_Cosmic[WT_N];
    float                n_pass_atleast2Lep[WT_N];
    float                n_pass_exactly2Lep[WT_N];
    float                n_pass_signalLep[WT_N];
    float                n_pass_flavor[ET_N][WT_N];
    float                n_pass_mll[ET_N][WT_N];
    float                n_pass_os[ET_N][WT_N];
    float                n_pass_ss[ET_N][WT_N];
    float                n_pass_evtTrig[ET_N][WT_N];
    float                n_pass_trigMatch[ET_N][WT_N];

    // SR1 counts
    float                n_pass_SR1jv[ET_N][WT_N];
    float                n_pass_SR1Zv[ET_N][WT_N];
    float                n_pass_SR1MET[ET_N][WT_N];

    // SR2 counts
    float                n_pass_SR2jv[ET_N][WT_N];
    float                n_pass_SR2MET[ET_N][WT_N];

    // SR3 counts
    float                n_pass_SR3ge2j[ET_N][WT_N];
    float                n_pass_SR3Zv[ET_N][WT_N];
    float                n_pass_SR3bjv[ET_N][WT_N];
    float                n_pass_SR3mct[ET_N][WT_N];
    float                n_pass_SR3MET[ET_N][WT_N];

    // SR5 counts
    float                n_pass_SR5jv[ET_N][WT_N];
    float                n_pass_SR5MET[ET_N][WT_N];
    float                n_pass_SR5Zv[ET_N][WT_N];
    float                n_pass_SR5L0pt[ET_N][WT_N];
    float                n_pass_SR5SUMpt[ET_N][WT_N];
    float                n_pass_SR5dPhiMETLL[ET_N][WT_N];
    float                n_pass_SR5dPhiMETL1[ET_N][WT_N];

    // SR4 counts
    float                n_pass_SR4jv[ET_N][WT_N];
    float                n_pass_SR4Zv[ET_N][WT_N];
    float                n_pass_SR4MET[ET_N][WT_N];
    float                n_pass_SR4MT2[ET_N][WT_N];



};

#endif
