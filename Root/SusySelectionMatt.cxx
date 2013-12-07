#include "SusyTest0/SusySelectionMatt.h"

#include <iomanip>
#include "TCanvas.h"

#include "SusyTest0/criteria.h"
#include "SusyTest0/SusySelection.h" // passHfor

using namespace std;
using namespace Susy;

/*--------------------------------------------------------------------------------*/
// SusySelectionMatt Constructor
/*--------------------------------------------------------------------------------*/
SusySelectionMatt::SusySelectionMatt() :
  m_susyObj(NULL),
  m_trigObj(NULL),
  m_useMCTrig(false),
  m_w(1.0),
  m_do1fb(false),
  m_doAD(false),
  m_dumpCounts(true),
  m_nLepMin(2),
  m_nLepMax(2),
  m_cutNBaseLep(true),
  m_ET(ET_Unknown),
  m_susyXS(NULL)
{
  // Loop over weight types
  for(int w=0; w<WT_N; ++w){
    n_readin[w]       = 0;
    n_pass_LAr[w]     = 0;
    n_pass_BadJet[w]  = 0;
    n_pass_BadMuon[w] = 0;
    n_pass_Cosmic[w]  = 0;
    n_pass_atleast2Lep[w] = 0;
    n_pass_exactly2Lep[w] = 0;
    n_pass_mll20[w]       = 0;
    n_pass_signalLep[w]   = 0;
    n_pass_HFOR[w]        = 0;
    n_pass_HotSpot[w]     = 0;
    n_pass_TileError[w]   = 0;
    n_pass_FEBCut[w]      = 0;

    // The rest are channel specific.
    for(int i=0; i<ET_N; ++i){
      n_pass_flavor[i][w]   = 0;
      n_pass_signalTau[i][w]   = 0;
      n_pass_evtTrig[i][w]     = 0;
      n_pass_trigMatch[i][w]   = 0;
      n_pass_mll[i][w]         = 0;
      n_pass_ss[i][w]          = 0;
      n_pass_os[i][w]          = 0;
      n_pass_truth[i][w]       = 0;
      //
      n_pass_CRWHSS2lss  [i][w] = 0;
      n_pass_CRWHSStauv  [i][w] = 0;
      n_pass_CRWHSSmuiso [i][w] = 0;
      n_pass_CRWHSSeled0 [i][w] = 0;
      n_pass_CRWHSSnfj   [i][w] = 0;
      n_pass_CRWHSSnbj   [i][w] = 0;
      n_pass_CRWHSSnj    [i][w] = 0;
      n_pass_CRWHSS2lpt  [i][w] = 0;
      n_pass_CRWHSSzveto [i][w] = 0;
      n_pass_CRWHSSmwwt  [i][w] = 0;
      n_pass_CRWHSShtmin [i][w] = 0;
      n_pass_CRWHSSmetrel[i][w] = 0;
      n_pass_CRWHSS      [i][w] = 0;
    }
  }// end loop over weight types

  m_BoundLow = TF1("expo","expo",1.5,4.0);
  m_BoundLow.FixParameter(0, 6.11861);
  m_BoundLow.FixParameter(1,-1.118815);
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  if(m_dbg) cout << "SusySelectionMatt::Begin" << endl;

  // Specify 2-lep ana type
  setAnaType(Ana_2LepWH);
  setSelectTaus(true);
  string period("Moriond");
  m_trigObj = new DilTrigLogic(period, false/*No Reweight Utils!*/);
  if(m_useMCTrig) m_trigObj->useMCTrigger();
}

/*--------------------------------------------------------------------------------*/
Bool_t SusySelectionMatt::Process(Long64_t entry)
{
  return kTRUE;
}
/*--------------------------------------------------------------------------------*/
void SusySelectionMatt::Terminate()
{
  SusyNtAna::Terminate();
  if(m_dbg) cout << "SusySelectionMatt::Terminate" << endl;
}
