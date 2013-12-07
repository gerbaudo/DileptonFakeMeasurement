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


class SusySelectionMatt : public SusyNtAna
{

  public:

    SusySelectionMatt();
    virtual ~SusySelectionMatt(){};
    virtual void    Begin(TTree *tree);
    virtual void    Terminate();
    virtual Bool_t  Process(Long64_t entry);
    ClassDef(SusySelectionMatt, 1);

  protected:

    SUSYObjDef* m_susyObj;            // susy obj

    DilTrigLogic*       m_trigObj;      // My trigger logic class
    bool                m_useMCTrig;    // Use MC Trigger
    float               m_w;            // mc weight
    bool                m_do1fb;        // For get weight method
    bool                m_doAD;         // do weights for A-B

    bool                m_dumpCounts;   // Flag to dump counters

    // Cut variables
    float                m_nLepMin;      // min leptons
    float                m_nLepMax;      // max leptons
    bool                m_cutNBaseLep;  // apply nLep cuts to baseline leptons as well as signal

    DiLepEvtType        m_ET;           // Dilepton event type to store cf

    XSReader*            m_susyXS;       // susy xs object

};

#endif
