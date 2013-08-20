#ifndef TIGHTPROBABILITY_H
#define TIGHTPROBABILITY_H

#include <string>

#include "SusyTest0/SusySelection.h"


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

 protected:

  std::string  m_outFname;
  TFile*       m_outFile;           // Output file
  bool         m_isMC;              // is MC
  LeptonVector m_probes;            // Probe lepton vector
  LeptonVector m_tags;              // Tag Lepton vector
  float        m_evtWeight;         // Event Weight
  bool         m_AltIso;            // If true, use Alt isolation
  float        m_metRel;            // Met Rel to be plotted
  int          m_ch;                // Set the channel

};

} // end namespace Susy

#endif
