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
#include "SusyTest0/FakeRegions.h"
#include "SusyTest0/FakeLeptonSources.h"
#include "SusyTest0/EffObject.h"
#include "SusyTest0/DileptonChannel.h"
#include <fstream>
#include <vector>

using namespace std;
using namespace Susy;
namespace sf = susy::fake;

class MeasureFakeRate2 : public SusySelectionMatt
{

 public:
  MeasureFakeRate2();
  virtual ~MeasureFakeRate2();
  virtual void    Begin(TTree *tree);
  virtual void    Terminate();
  virtual Bool_t  Process(Long64_t entry);
  void initHistos(string outName);
  // Data Control Regions
  bool passRealCR(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                  sf::ControlRegion CR);
  bool passHFCR(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                sf::ControlRegion CR);
  bool passConvCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  bool passSignalRegion(const LeptonVector &leptons, const JetVector& jets, const Met* met,
                        sf::ControlRegion CR);
  // Monte Carlo Regions
  bool passMCReg(const LeptonVector &leptons, const JetVector& jets,
                 const Met* met, sf::ControlRegion CR);
  void fillRatesHistos(const Lepton* lep, const JetVector& jets,
                       const Met* met, sf::ControlRegion CR);
  // Miscellaneous
  sf::LeptonSource getLeptonSource(const Lepton* l);
  const int CR_N;
  static const int kNmaxControlRegions=64;
  const std::vector<int> m_controlRegions; //!< where we compute SF and rates (pseudo t&p)
  const std::vector<int> m_signalRegions;  //!< where we compute fractions to make the weighted avg

  EffObject* h_l_pt         [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_pt_coarse  [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_eta        [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_eta_coarse [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_metrel       [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_met          [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_njets        [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_onebin       [LT_N][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_flavor       [LT_N][kNmaxControlRegions][susy::wh::Ch_N];

 protected:
  TFile*       m_outFile;           // Output file
  LeptonVector m_probes;            // Probe lepton vector
  LeptonVector m_tags;              // Tag Lepton vector
  float        m_evtWeight;         // Event Weight
  float        m_metRel;            // Met Rel to be plotted
  int          m_ch;                // Set the channel

};
#endif
