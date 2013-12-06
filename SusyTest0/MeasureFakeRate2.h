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
#include "SusyMatrixMethod/FakeRegions.h"
#include "SusyTest0/DileptonChannel.h"
#include "SusyTest0/FakeLeptonSources.h"
#include "SusyTest0/SusySelectionMatt.h"
#include "SusyTest0/EffObject.h"
#include <vector>

using namespace std;
using namespace Susy;

class MeasureFakeRate2 : public SusySelectionMatt
{

 public:
  MeasureFakeRate2();
  virtual ~MeasureFakeRate2();
  virtual void    Begin(TTree *tree);
  virtual void    Terminate();
  virtual Bool_t  Process(Long64_t entry);
  void initHistos(string outName);
  MeasureFakeRate2& setFileName(string f){ m_fileName = f; return *this; }
  // Data Control Regions
  bool passRealCR(const LeptonVector &leptons, const JetVector& jets, const Met* met, susy::fake::Region CR);
  bool passHFCR(const LeptonVector &leptons, const JetVector& jets, const Met* met, susy::fake::Region CR);
  bool passConvCR(const LeptonVector &leptons, const JetVector& jets, const Met* met);
  bool passSignalRegion(const LeptonVector &leptons, const JetVector& jets, const Met* met, susy::fake::Region CR);
  // Monte Carlo Regions
  bool passMCReg(const LeptonVector &leptons, const JetVector& jets, const Met* met, susy::fake::Region CR);
  void fillRatesHistos(const Lepton* lep, const JetVector& jets, const Met* met, size_t regionIndex);
  // Miscellaneous
  susy::fake::LeptonSource getLeptonSource(const Lepton* l);
  static const size_t kNmaxControlRegions=64, kNmaxLeptonTypes=2;
  const std::vector<susy::fake::Region> m_controlRegions; //!< where we compute SF and rates (pseudo t&p)
  const std::vector<susy::fake::Region> m_signalRegions;  //!< where we compute fractions to make the weighted avg
  const std::vector<susy::fake::Region> allRegions() const; //!< generate on the fly the sum of the two above
  enum LeptonType {kElectron, kMuon}; // DG Dec13: move this enum to a separate file when Fake is a separate package
  const std::vector<LeptonType> m_leptonTypes; //! types of leptons for which we will measure the fake probability
  static const std::string LeptonType2str(const LeptonType l);

  EffObject* h_l_pt         [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_pt_coarse  [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_eta        [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_l_eta_coarse [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_metrel       [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_met          [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_njets        [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_onebin       [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];
  EffObject* h_flavor       [kNmaxLeptonTypes][kNmaxControlRegions][susy::wh::Ch_N];

 protected:
  std::string  m_fileName;          // Outname file name
  TFile*       m_outFile;           // Output file
  LeptonVector m_probes;            // Probe lepton vector
  LeptonVector m_tags;              // Tag Lepton vector
  float        m_evtWeight;         // Event Weight
  float        m_metRel;            // Met Rel to be plotted
  int          m_ch;                // Set the channel

};
#endif
