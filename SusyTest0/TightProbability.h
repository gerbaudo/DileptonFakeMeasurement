#ifndef TIGHTPROBABILITY_H
#define TIGHTPROBABILITY_H

#include "SusyTest0/SusySelection.h"

#include "TH1F.h"

#include <string>
#include <vector>



class TDirectory;

namespace Susy {

//! A class to determine p(tight|real) and p(tight|fake)
/*!
  More details to come.
  Based on Matt's (mrelich@) procedure.

  davide.gerbaudo@gmail.com
  August 2013
 */

class TightProbability : public SusySelection
{
 public:
  enum LeptonOrigin {
    kReal=0,
    kHeavyFlavor,
    kLigthFlavor,
    kConversion,
    kMultijet,
    kUnknownOrigin,
    kLeptonOriginN
  };
  enum TightLoosePairType {
    kTightTight=0,
    kTightLoose,
    kLooseTight,
    kLooseLoose,
    kAnyAny,
    kTightLoosePairTypeN
  };
  //! container used instead of TEfficiency; also catches under/overflow
  struct NumDenHisto {
    // \todo implement automatic rebinning when computing ratio
    TH1F m_num, m_den;
    float m_min, m_max;
    float m_widthFirst, m_widthLast;
    NumDenHisto(string name, int nbins, float min, float max);
    NumDenHisto(string name, int nbins, const float* binEdges);
    void Fill(bool alsoFillNum, float weight, float value);
    void Sumw2() { m_num.Sumw2(); m_den.Sumw2(); }
    void SetDirectory(TDirectory* dir) {m_num.SetDirectory(dir); m_den.SetDirectory(dir);}
    void setMinMax();
  };

 public:
  TightProbability();
  virtual ~TightProbability();
  virtual void    Begin(TTree *tree);
  virtual void    Terminate();
  virtual Bool_t  Process(Long64_t entry);
  void initOutput(string outName);
  void finalizeOutput();
  TightProbability& setOutputFilename(const std::string &s);

  TightProbability::LeptonOrigin getLeptonOrigin(const Lepton* l);
  TightProbability::TightLoosePairType getTightLoosePairType(const Lepton* tag,
                                                             const Lepton* probe);
/*   bool isSignalWithEtcone(const Lepton* lep); */
/*   bool isSignalWithoutEtcone(const Lepton* lep); */
/*   bool isSignalWithoutPtcone(const Lepton* lep); */
/*   bool isSignalWithoutd0Sig(const Lepton* lep); */
/*   bool isSignalWithoutz0Sig(const Lepton* lep); */
/*   bool isSignalWithoutIP(const Lepton* lep); */

  float getPtcone(const Lepton* lep);
  float getEtcone(const Lepton* lep);
  void printLeptonOriginCounters() const;
 protected:
  NumDenHisto* histoForLepForOrigin(const Lepton *l, const LeptonOrigin &lo);
  void incrementLeptonOriginCounter(const LeptonOrigin &lo);
 protected:

  std::string  m_outFname;
  TFile*       m_outFile;           // Output file
  bool         m_isMC;              // is MC
  LeptonVector m_probes;            // Probe lepton vector
  LeptonVector m_tags;              // Tag Lepton vector
  float        m_evtWeight;         // Event Weight
  NumDenHisto  m_h_el_pt_any;  // electron probabilities
  NumDenHisto  m_h_el_pt_real;
  NumDenHisto  m_h_el_pt_hf;
  NumDenHisto  m_h_el_pt_lf;
  NumDenHisto  m_h_el_pt_conv;
  NumDenHisto  m_h_el_pt_mjet;
  NumDenHisto  m_h_el_pt_other;
  NumDenHisto  m_h_mu_pt_any;  // muon probabilities
  NumDenHisto  m_h_mu_pt_real;
  NumDenHisto  m_h_mu_pt_hf;
  NumDenHisto  m_h_mu_pt_lf;
  NumDenHisto  m_h_mu_pt_conv;
  NumDenHisto  m_h_mu_pt_mjet;
  NumDenHisto  m_h_mu_pt_other;
  std::vector< size_t > m_leptonOriginCounter;
};

} // end namespace Susy

#endif
