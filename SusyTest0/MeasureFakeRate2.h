#ifndef MeasureFakeRate2_h
#define MeasureFakeRate2_h


//////////////////////////////////////////////////////////
// Code to measure the fake rates and real efficiencies //
// from a variety of control regions. Output will be    //
// root files that can then be used to calc. weights.   //
//////////////////////////////////////////////////////////

// Root Packages
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TProfile.h"

// Susy Packages
#include "SusyTest0/SusySelectionMatt.h"
#include "SusyTest0/SusyAnaDefsMatt.h"
#include "SusyTest0/EffObject.h"
#include <fstream>
#include <vector>

using namespace std;
using namespace Susy;

typedef unsigned int uint;

enum CRPLOT
{
  CRP_met=0,
  CRP_mll,
  CRP_mt_tag,
  CRP_mt_probe,
  CRP_ht,
  CRP_N
};

class MeasureFakeRate2 : public SusySelectionMatt
{
 public:
  MeasureFakeRate2();
  virtual ~MeasureFakeRate2();
  virtual void    Begin(TTree *tree);
  virtual void    Terminate();
  virtual Bool_t  Process(Long64_t entry);
  // Initialize Histograms
  void initHistos(string outName);
  // Data Control Regions
  bool passRealCR(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                  ControlRegion CR);
  bool passHFCR(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                ControlRegion CR);
  bool passLFZjetCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  bool passLFWjetCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  bool passConvCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  bool passSignalRegion(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                        ControlRegion CR);
  bool passCFCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  // Monte Carlo Regions
  bool passMCReg(const LeptonVector &leptons, const JetVector& jets,
                 const Met* met, ControlRegion CR);
  void fillRatesHistos(const Lepton* lep, const JetVector& jets,
                       const Met* met, ControlRegion CR);
  void fillCrHistos(const Lepton* tag, const Lepton* probe, const JetVector& jets, const Met* met,
                    ControlRegion CR, CRPLOT CRP);
  // Miscellaneous
  bool passLFTrig(const LeptonVector &leps);
  LeptonSource getLeptonSource(const Lepton* l);
  bool isSignalWithEtcone(const Lepton* lep);
  bool isSignalWithoutEtcone(const Lepton* lep);
  bool isSignalWithoutPtcone(const Lepton* lep);
  bool isSignalWithoutd0Sig(const Lepton* lep);
  bool isSignalWithoutz0Sig(const Lepton* lep);
  bool isSignalWithoutIP(const Lepton* lep);
  DiLepPair getDilepPair(const Lepton* tag, const Lepton* probe);
  // Histograms for variables
  const int CR_N;
  static const int kNmaxControlRegions=64;
  //  static const int getNcontrolRegions()
  const std::vector<int> m_controlRegions; // where we compute SF and rates (pseudo t&p)
  const std::vector<int> m_signalRegions;  // where we compute fractions to make the weighted avg

  #define NEW(name)   name[LT_N][kNmaxControlRegions][Ch_N]
  #define NEWCR(name) name[LT_N][kNmaxControlRegions]

  EffObject* NEW(h_l_pt);
  EffObject* NEW(h_l_pt_coarse);
  EffObject* NEW(h_l_pt_heavy);
  EffObject* NEW(h_l_pt_others);
  EffObject* NEW(h_l_eta);
  EffObject* NEW(h_l_eta_coarse);
  EffObject* NEW(h_metrel);
  EffObject* NEW(h_metrel_fine);
  EffObject* NEW(h_metrel_coarse);
  EffObject* NEW(h_met);
  EffObject* NEW(h_met_fine);
  EffObject* NEW(h_met_coarse);

  EffObject* NEW(h_njets);
  EffObject* NEW(h_nlightjets);
  EffObject* NEW(h_nheavyjets);
  EffObject* NEW(h_nlightjetsNoB);

  EffObject* NEW(h_flavor);
  EffObject* NEW(h_l_type);
  EffObject* NEW(h_l_origin);

  EffObject* NEW(h_onebin);

  EffObject* NEW(h_heavy_d0sig);
  EffObject* NEW(h_light_d0sig);
  EffObject* NEW(h_conv_d0sig);

  // New params
  EffObject* NEW(h_ht);
  EffObject* NEW(h_ht_pt);
  EffObject* NEW(h_ht_wMet);
  EffObject* NEW(h_ht_pt_wMet);

  // 2D param
  EffObject2* NEW(h_l_pt_bjet);
  EffObject2* NEW(h_l_pt_eta);

  // CR Plots
  EffObject* NEWCR(h_met_cr);
  EffObject* NEWCR(h_mll_cr);
  EffObject* NEWCR(h_mt_tag_cr);
  EffObject* NEWCR(h_mt_probe_cr);
  EffObject* NEWCR(h_ht_cr);

  // Check for Conf
  EffObject* NEW(h_with_without_Etcone);


  #undef NEW
  #undef NEWCR

  // Histos to determine cuts
  #define NEWCUT(name) name[LT_N][kNmaxControlRegions]

  EffObject* NEWCUT(h_ptcone);
  EffObject* NEWCUT(h_etcone);
  EffObject* NEWCUT(h_d0Sig);
  EffObject* NEWCUT(h_z0Sig);

  EffObject* NEWCUT(h_dist_relptcone);
  EffObject* NEWCUT(h_dist_d0Sig);
  EffObject* NEWCUT(h_dist_z0Sig);
  EffObject* NEWCUT(h_dist_l_pt);
  EffObject* NEWCUT(h_dist_ptcone);

  EffObject2* NEWCUT(h_dist_mll_ptcone);

  TProfile* NEWCUT(p_npv_etcone);
  TProfile* NEWCUT(p_npv_ptcone);
  TProfile* NEWCUT(p_npv_ptconeElStyle);
  TProfile* NEWCUT(p_mu_etcone);
  TProfile* NEWCUT(p_mu_ptcone);
  TProfile* NEWCUT(p_mu_ptconeElStyle);

  #undef NEWCUT

 protected:
  TFile*       m_outFile;           // Output file
  LeptonVector m_probes;            // Probe lepton vector
  LeptonVector m_tags;              // Tag Lepton vector
  float        m_evtWeight;         // Event Weight
  float        m_metRel;            // Met Rel to be plotted
  int          m_ch;                // Set the channel

};
#endif
