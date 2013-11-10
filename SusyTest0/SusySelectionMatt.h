#ifndef SUSYSELECTIONMATT_h
#define SUSYSELECTIONMATT_h

// Common Packages
#include "Mt2/mt2_bisect.h"
#include "LeptonTruthTools/RecoTruthMatch.h"

// Root Packages
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"

// Susy Common
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyXSReader/XSReader.h"
#include "SusyTest0/SusyAnaDefsMatt.h"
#include "SusyTest0/SsPassFlags.h"

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

class SusySelectionMatt : public SusyNtAna
{

  public:

    SusySelectionMatt();
    virtual ~SusySelectionMatt(){};
    virtual void    Begin(TTree *tree);
    virtual void    Terminate();
    virtual Bool_t  Process(Long64_t entry);
    virtual void dumpEventCounters();
    bool selectEvent(bool count=false);
    bool selectBaseEvent(bool doMll=true, bool count=false);
    bool selectAnaEvent(const LeptonVector& leptons, const LeptonVector& baseLeptons, bool count=false);
    // Signal regions
    // Completely Redefining this shit!!!

    // SRmT2
    bool passSRmT2a(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passSRmT2b(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passSRmT2c(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    
    // SRWW
    bool passSRWWa(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passSRWWb(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passSRWWc(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    
    // SR-SS
    //bool passSRSS(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);

    // SR-ZJets
    bool passSRZjets(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);

    // Control and validaiton regions
    bool passVRSS(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // CR for SR-WW and SR-mT2
    bool passCRWWMet(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRWWmT2(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRTopMet(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRTopmT2(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZVMet(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZVmT2a(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZVmT2b(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZVmT2c(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZVmT2d(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRPremT2(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // Cr for SR-ZJets
    bool passCRTopZjets(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passCRZXZjets(const LeptonVector& leptons, const JetVector& jets, const Met* met, bool count=false, bool loose=false);
    bool passPreSRZjets(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passCRZjets(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // Region to check Met sys
    bool passSimpleZ(const LeptonVector& leptons, const JetVector& jets, const Met* met);
    bool passSimpleZ2(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    bool passMuonRelIso(const LeptonVector &leptons, float maxVal);
    SsPassFlags passWhSS(const LeptonVector& leptons, const JetVector& jets, const Met* met);

    // Higgs jet counting..
    //int nHiggsSignalJets(const JetVector& jets);
    
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
    bool passBadMet(const Met* met, float cutval=0.8);

    // Signal Region Cuts
    bool passJetVeto(const JetVector& jets);
    int nL20Close(const JetVector& jets);
    int nB20Close(const JetVector& jets);
    int nF30Close(const JetVector& jets);

    bool passZVeto(const LeptonVector& leptons, float Zlow = 81.2, float Zhigh = 101.2);
    bool passMETRel(const Met *met, const LeptonVector& leptons, 
		    const JetVector& jets, float maxMet = 100);
    bool passbJetVeto(const JetVector& jets);
    bool passge2Jet(const JetVector& jets);
    bool passdPhi(TLorentzVector v0, TLorentzVector v1, float cut);
    bool passMT2(const LeptonVector& leptons, const Met* met, float cut);
    float getMt2(const LeptonVector& leptons, const Met* met);

    // Idendification methods
    bool isRealLepton(const Lepton* lep);
    bool isFakeLepton(const Lepton* lep);
    bool isConvLepton(const Lepton* lep);
    bool isHFLepton(const Lepton* lep);
    bool isLFLepton(const Lepton* lep);
    bool isQCDLepton(const Lepton* lep);
    bool isTrueDilepton(const LeptonVector &leptons);
    bool isFakeDilepton(const LeptonVector &leptons);
    void setFileName(string f){ m_fileName = f; };
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

    void setDoSusy(bool susy){ m_doSusy = susy; };

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
      
      float trig = m_baseLeptons.size() == 2 && nt.evt()->isMC && !m_useMCTrig ? 
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

      float allAE = getEventWeight(LUMI_A_L,true) * btag * trig;
      allAE = includeLepSF ? allAE * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF : allAE;
      flag[WT_AllAE] += allAE;

    };

    // Channels for plotting
    int getChan(const LeptonVector& leps);

    // Miscellaneous methods
    void printCounter(string cut, float counter[ET_N][WT_N], int weight);
    ClassDef(SusySelectionMatt, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger
    string              m_fileName;     // File name
    float               m_w;            // mc weight
    bool                m_do1fb;        // For get weight method
    bool                m_doAD;         // do weights for A-B

    bool                m_dumpCounts;   // Flag to dump counters

    TF1                 m_BoundLow;      // function for ZJets CR

    // Cut variables
    float                m_nLepMin;      // min leptons
    float                m_nLepMax;      // max leptons
    bool                m_cutNBaseLep;  // apply nLep cuts to baseline leptons as well as signal

    DiLepEvtType        m_ET;           // Dilepton event type to store cf

    bool                 m_doSusy;       // flag to turn on susy xs
    XSReader*            m_susyXS;       // susy xs object

    // Event counters
    float                n_readin[WT_N];
    float                n_pass_LAr[WT_N];
    float                n_pass_BadJet[WT_N];
    float                n_pass_BadMuon[WT_N];
    float                n_pass_Cosmic[WT_N];
    float                n_pass_atleast2Lep[WT_N];
    float                n_pass_exactly2Lep[WT_N];
    float                n_pass_mll20[WT_N];
    float                n_pass_signalLep[WT_N];
    float                n_pass_HFOR[WT_N];
    float                n_pass_HotSpot[WT_N];
    float                n_pass_TileError[WT_N];
    float                n_pass_FEBCut[WT_N];

    float                n_pass_flavor[ET_N][WT_N];
    float                n_pass_mll[ET_N][WT_N];    
    float                n_pass_signalTau[ET_N][WT_N];
    float                n_pass_os[ET_N][WT_N];
    float                n_pass_ss[ET_N][WT_N];
    float                n_pass_evtTrig[ET_N][WT_N];
    float                n_pass_trigMatch[ET_N][WT_N];
    float                n_pass_truth[ET_N][WT_N];


    // SRmT2a
    float                n_pass_SRmT2a_zv[ET_N][WT_N];
    float                n_pass_SRmT2a_jv[ET_N][WT_N];
    float                n_pass_SRmT2a_PtCut[ET_N][WT_N];
    float                n_pass_SRmT2a_mt2[ET_N][WT_N];
    float                n_pass_SRmT2b_mt2[ET_N][WT_N];
    float                n_pass_SRmT2c_mt2[ET_N][WT_N];


    // SR WW
    float                n_pass_SRWWa_zv[ET_N][WT_N];
    float                n_pass_SRWWa_jv[ET_N][WT_N];
    float                n_pass_SRWWa_PtCut[ET_N][WT_N];
    float                n_pass_SRWWa_ptll80[ET_N][WT_N];
    float                n_pass_SRWWa_metrel80[ET_N][WT_N];
    float                n_pass_SRWWa_mll120[ET_N][WT_N];

    // SR WW b
    float                n_pass_SRWWb_zv[ET_N][WT_N];
    float                n_pass_SRWWb_jv[ET_N][WT_N];
    float                n_pass_SRWWb_PtCut[ET_N][WT_N];
    float                n_pass_SRWWb_mll170[ET_N][WT_N];
    float                n_pass_SRWWb_mt2_90[ET_N][WT_N];

    // SR WW c
    float                n_pass_SRWWc_zv[ET_N][WT_N];
    float                n_pass_SRWWc_jv[ET_N][WT_N];
    float                n_pass_SRWWc_PtCut[ET_N][WT_N];
    float                n_pass_SRWWc_mt2_100[ET_N][WT_N];

    // SR ZJets
    float                n_pass_SRZjets_zw[ET_N][WT_N];
    float                n_pass_SRZjets_2ljets[ET_N][WT_N];
    float                n_pass_SRZjets_bfveto[ET_N][WT_N];
    float                n_pass_SRZjets_JetPt[ET_N][WT_N];
    float                n_pass_SRZjets_PtCut[ET_N][WT_N];
    float                n_pass_SRZjets_ptll80[ET_N][WT_N];
    float                n_pass_SRZjets_mjjw[ET_N][WT_N];
    float                n_pass_SRZjets_metrel80[ET_N][WT_N];
    float                n_pass_SRZjets_dRll[ET_N][WT_N];


    // CR WWMet
    float                n_pass_CRWWMet_OF[ET_N][WT_N];
    float                n_pass_CRWWMet_jv[ET_N][WT_N];
    float                n_pass_CRWWMet_PtCut[ET_N][WT_N];
    float                n_pass_CRWWMet_mll120[ET_N][WT_N];
    float                n_pass_CRWWMet_ptll40[ET_N][WT_N];
    float                n_pass_CRWWMet_metrel_60_80[ET_N][WT_N];

    // CR WWmT2 
    float                n_pass_CRWWmt2_OF[ET_N][WT_N];
    float                n_pass_CRWWmt2_jv[ET_N][WT_N];
    float                n_pass_CRWWmt2_PtCut[ET_N][WT_N];
    float                n_pass_CRWWmt2_mt2_50_90[ET_N][WT_N];

    // CR Top mT2
    float                n_pass_CRTopmt2_OF[ET_N][WT_N];
    float                n_pass_CRTopmt2_1bjet[ET_N][WT_N];
    float                n_pass_CRTopmt2_lfveto[ET_N][WT_N];
    float                n_pass_CRTopmt2_PtCut[ET_N][WT_N];
    float                n_pass_CRTopmt2_mt2_70[ET_N][WT_N];

    // CR Top Met
    float                n_pass_CRTopMet_OF[ET_N][WT_N];
    float                n_pass_CRTopMet_1bjet[ET_N][WT_N];
    float                n_pass_CRTopMet_lfveto[ET_N][WT_N];
    float                n_pass_CRTopMet_PtCut[ET_N][WT_N];
    float                n_pass_CRTopMet_mll120[ET_N][WT_N];
    float                n_pass_CRTopMet_ptll80[ET_N][WT_N];
    float                n_pass_CRTopMet_metrel80[ET_N][WT_N];


    // CR Top ZJets
    float                n_pass_CRTopZjets_SF[ET_N][WT_N];
    float                n_pass_CRTopZjets_zv[ET_N][WT_N];
    float                n_pass_CRTopZjets_2jets[ET_N][WT_N];
    float                n_pass_CRTopZjets_bjet[ET_N][WT_N];
    float                n_pass_CRTopZjets_fveto[ET_N][WT_N];
    float                n_pass_CRTopZjets_PtCut[ET_N][WT_N];
    float                n_pass_CRTopZjets_ptll80[ET_N][WT_N];
    float                n_pass_CRTopZjets_dRll[ET_N][WT_N];
    float                n_pass_CRTopZjets_metrel80[ET_N][WT_N];

    // CR ZV MET
    float                n_pass_CRZVMet_zw[ET_N][WT_N];
    float                n_pass_CRZVMet_jv[ET_N][WT_N];
    float                n_pass_CRZVMet_PtCut[ET_N][WT_N];
    float                n_pass_CRZVMet_ptll80[ET_N][WT_N];
    float                n_pass_CRZVMet_metrel80[ET_N][WT_N];

    // CR ZV mt2
    float                n_pass_CRZVmt2a_zw[ET_N][WT_N];
    float                n_pass_CRZVmt2a_jv[ET_N][WT_N];
    float                n_pass_CRZVmt2a_PtCut[ET_N][WT_N];
    float                n_pass_CRZVmt2a_mt2_90[ET_N][WT_N];
    float                n_pass_CRZVmt2b_mt2_120[ET_N][WT_N];
    float                n_pass_CRZVmt2c_mt2_150[ET_N][WT_N];
    float                n_pass_CRZVmt2d_mt2_100[ET_N][WT_N];

    // CR ZXZjets
    float                n_pass_CRZXZjets_SF[ET_N][WT_N];
    float                n_pass_CRZXZjets_zw[ET_N][WT_N];
    float                n_pass_CRZXZjets_2ljets[ET_N][WT_N];
    float                n_pass_CRZXZjets_bfveto[ET_N][WT_N];
    float                n_pass_CRZXZjets_JetPtCut[ET_N][WT_N];
    float                n_pass_CRZXZjets_PtCut[ET_N][WT_N];
    float                n_pass_CRZXZjets_dRll[ET_N][WT_N];
    float                n_pass_CRZXZjets_metrel[ET_N][WT_N];
    float                n_pass_CRZXZjets_lowerbound[ET_N][WT_N];

    float                n_pass_CRWHSS2lss  [ET_N][WT_N];
    float                n_pass_CRWHSStauv  [ET_N][WT_N];
    float                n_pass_CRWHSSmuiso [ET_N][WT_N];
    float                n_pass_CRWHSSeled0 [ET_N][WT_N];
    float                n_pass_CRWHSSnfj   [ET_N][WT_N];
    float                n_pass_CRWHSSnbj   [ET_N][WT_N];
    float                n_pass_CRWHSSnj    [ET_N][WT_N];
    float                n_pass_CRWHSS2lpt  [ET_N][WT_N];
    float                n_pass_CRWHSSzveto [ET_N][WT_N];
    float                n_pass_CRWHSSmwwt  [ET_N][WT_N];
    float                n_pass_CRWHSShtmin [ET_N][WT_N];
    float                n_pass_CRWHSSmetrel[ET_N][WT_N];
    float                n_pass_CRWHSS      [ET_N][WT_N];
};

#endif
