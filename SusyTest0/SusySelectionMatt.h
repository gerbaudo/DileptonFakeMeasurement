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
    // Method to increment the counters for the event weight types
    void increment(float flag[], bool includeLepSF=false, bool includeBtag=false){
      flag[WT_Raw]   += 1.0;
      flag[WT_Evt]   += nt.evt()->w;
      flag[WT_PU]    += nt.evt()->w * nt.evt()->wPileup;
      flag[WT_PU1fb] += nt.evt()->w * nt.evt()->wPileupAB3;
      flag[WT_LSF]   += (includeLepSF ? 
			 nt.evt()->w * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF :
			 nt.evt()->w);
      float btag = includeBtag ? 1.0 : 1.0;
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

      float allAE = 1.0 * btag * trig;
      allAE = includeLepSF ? allAE * m_baseLeptons[0]->effSF * m_baseLeptons[1]->effSF : allAE;
      flag[WT_AllAE] += allAE;

    };

    // Miscellaneous methods
    void printCounter(string cut, float counter[ET_N][WT_N], int weight);
    ClassDef(SusySelectionMatt, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger
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
