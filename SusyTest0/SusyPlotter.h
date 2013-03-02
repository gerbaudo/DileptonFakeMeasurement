#ifndef SusyAna_SusyPlotter_h
#define SusyAna_SusyPlotter_h

#include "SusyTest0/SusySelection.h"
#include "SusyTest0/SusyAnaDefs.h"

#include "TH1F.h"
#include "TH2F.h"


/*

    SusyPlotter - class for making analysis histograms

*/

// lepton channel
enum Chan {
  Ch_all = 0,
  Ch_ee,
  Ch_mm,
  Ch_em,
  Ch_N
};

// Lepton chan names
static string chanNames[] = {
  "all",
  "ee",
  "mm",
  "em"
};

// Plotting regions
enum PlotRegion{
  PR_SR1 = 0,
  PR_SR2,
  PR_SR3,
  PR_SR4,
  PR_SR4b,
  PR_NONE,
  PR_VR1,
  PR_VR2,
  PR_VR3,
  PR_VR4,
  PR_SSInc,
  PR_OSInc,
  PR_WWCR1,
  PR_WWCR2,
  PR_TOPCR,
  PR_BR1,
  PR_BR2,
  PR_BR3,
  PR_BR4,
  PR_SR6base,
  PR_SR7base,
  PR_SR8base,
  PR_SR9base,
  PR_SR6ge2j,
  PR_SR7ge2j,
  PR_SR8ge2j,
  PR_SR9ge2j,
  PR_SR6ge3j,
  PR_SR7ge3j,
  PR_SR8ge3j,
  PR_SR9ge3j,
  PR_SR6eq2jNfv,
  PR_SR7eq2jNfv,
  PR_SR8eq2jNfv,
  PR_SR9eq2jNfv,
  PR_SR6ge2jNfv,
  PR_SR7ge2jNfv,
  PR_SR8ge2jNfv,
  PR_SR9ge2jNfv,
  PR_SR6,
  PR_SR7,
  PR_SR8,
  PR_SR9,
  PR_N
};

static string PRNames[] = {
  "sr1",
  "sr2",
  "sr3",
  "sr4",
  "sr4b",
  "srnone",
  "vr1",
  "vr2",
  "vr3",
  "vr4",
  "ssinc",
  "osinc",
  "wwcr1",
  "wwcr2",
  "topcr",
  "br1",
  "br2",
  "br3",
  "br4",
  "sr6base",
  "sr7base",
  "sr8base",
  "sr9base",
  "sr6ge2j",
  "sr7ge2j",
  "sr8ge2j",
  "sr9ge2j",
  "sr6ge3j",
  "sr7ge3j",
  "sr8ge3j",
  "sr9ge3j",
  "sr6eq2jNfv",
  "sr7eq2jNfv",
  "sr8eq2jNfv",
  "sr9eq2jNfv",
  "sr6ge2jNfv",
  "sr7ge2jNfv",
  "sr8ge2jNfv",
  "sr9ge2jNfv",
  "sr6",
  "sr7",
  "sr8",
  "sr9"
};

class SusyPlotter : public SusySelection
{

  public:

    SusyPlotter();
    virtual ~SusyPlotter(){};



    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);

    // 
    // Method to fill histograms
    //
    void fillHistos(const LeptonVector& leps, const JetVector &jets, const Met* met, 
		    const float weight, PlotRegion PR = PR_NONE, uint sys = 0);

    //
    // Define some additional plot regions
    //
    
    bool passZwindow(const LeptonVector& leps);

    //
    // Miscellaneous methods
    //
    
    // compute Mt
    float Mt(TLorentzVector p1, TLorentzVector met){
      return sqrt(2*p1.Pt()*met.Et()*(1-cos(p1.DeltaPhi(met))));
    };

    // compute lepton channel
    int getChan(const LeptonVector& leps);

    // get list of systematics to consider
    // override in SusyMatrixMethod
    void setSysts();
 public:
    static float transverseMass(const TLorentzVector &lep, const TLorentzVector &met);

    ClassDef(SusyPlotter, 1);

  protected:

    std::vector<uint> m_systs;                // systematics to process
    std::vector<string> m_systNames;          // systematics to process

    std::string         m_histFileName;         // output histo file name
    TFile*              m_histFile;             // output histo file

    // control flags
    bool                m_doLepSF;              // use lepton efficiency SF
    bool                m_doTrigW;              // do trigger reweighting
    //bool                m_do1fb;                // use 1/fb weights
    bool                m_doFake;               // do Fake estimate
    bool                m_doCF;                 // do charge flip

    // preprocessor convenience - add more indices later
    #define DEFHIST( name ) h_ ## name[Ch_N][PR_N][40/*Guess for # of sys*/];

    // Pt
    TH1F* DEFHIST(l0_pt);
    TH1F* DEFHIST(l1_pt);
    TH1F* DEFHIST(e_pt);
    TH1F* DEFHIST(m_pt);

    TH1F* DEFHIST(j0_pt);
    TH1F* DEFHIST(j1_pt);

    // eta
    TH1F* DEFHIST(l0_eta);
    TH1F* DEFHIST(l1_eta);
    TH1F* DEFHIST(e_eta);
    TH1F* DEFHIST(m_eta);
    TH1F* DEFHIST(j0_eta);
    TH1F* DEFHIST(j1_eta);

    // Mass
    TH1F* DEFHIST(ll_M);
    TH1F* DEFHIST(jj_M);
    
    // MET
    TH1F* DEFHIST(met);
    TH1F* DEFHIST(metrel);
    TH1F* DEFHIST(met_refEle);
    TH1F* DEFHIST(met_refMuo);
    TH1F* DEFHIST(met_refJet);
    TH1F* DEFHIST(met_softJet);
    TH1F* DEFHIST(met_refGamma);
    TH1F* DEFHIST(met_refCell);

    // # of jets
    TH1F* DEFHIST(njets);
    TH1F* DEFHIST(nbjets);

    // Type and origin
    TH1F* DEFHIST(l_type);
    TH1F* DEFHIST(l_origin);

    // One bin
    TH1F* DEFHIST(onebin);

    // Charge plots
    TH1F* DEFHIST(sumQ);
    TH1F* DEFHIST(ll_M_pos);
    TH1F* DEFHIST(ll_M_neg);
    TH1F* DEFHIST(njets_pos);
    TH1F* DEFHIST(njets_neg);
    TH1F* DEFHIST(njets_mll_90_120);
    TH1F* DEFHIST(njets_mll_90_120_pos);
    TH1F* DEFHIST(njets_mll_90_120_neg);
    TH1F* DEFHIST(nbjets_mll_90_120);
    TH1F* DEFHIST(nbjets_mll_90_120_pos);
    TH1F* DEFHIST(nbjets_mll_90_120_neg);

    TH1F* DEFHIST(nfjets);
    TH1F* DEFHIST(nfjets_mll_90_120);

    TH1F* DEFHIST(llj_M);
    TH1F* DEFHIST(llj_M_pos);
    TH1F* DEFHIST(llj_M_neg);
    TH1F* DEFHIST(llj_M_mll_90_120);
    TH1F* DEFHIST(llj_M_mll_90_120_pos);
    TH1F* DEFHIST(llj_M_mll_90_120_neg);

    // Met mass plots
    TH1F* DEFHIST(met_j_M);
    TH1F* DEFHIST(met_ll_M);
    TH1F* DEFHIST(met_j_ll_M);
    TH1F* DEFHIST(met_j_Mt);
    TH1F* DEFHIST(met_ll_Mt);
    TH1F* DEFHIST(met_j_ll_Mt);

    TH1F* DEFHIST(met_l0_Mt);
    TH1F* DEFHIST(met_l1_Mt);

    TH1F* DEFHIST(met_ll_Mt_noj);
    TH1F* DEFHIST(met_ll_Mt_onej);
    TH1F* DEFHIST(met_ll_Mt_twoj);
    TH1F* DEFHIST(met_ll_Mt_ge3j);
    TH1F* DEFHIST(met_ll_Mt_oneOrtwoj);

    TH1F* DEFHIST(dPhi_llmet_j);
    TH1F* DEFHIST(dPhi_met_l0);
    TH1F* DEFHIST(dPhi_met_l1);
    TH1F* DEFHIST(dPhi_met_ll);
    TH1F* DEFHIST(dPhi_met_j);
    TH1F* DEFHIST(dPhi_ll_j);
    TH1F* DEFHIST(dPhi_l0_j);
    TH1F* DEFHIST(dPhi_l1_j);
    TH1F* DEFHIST(dPhi_l0_l1);

    TH1F* DEFHIST(dPhi_ll_jj);
    TH1F* DEFHIST(dPhi_l0_jj);
    TH1F* DEFHIST(dPhi_l1_jj);
    TH1F* DEFHIST(dR_l0_l1);
    TH1F* DEFHIST(dR_ll_jj);


    TH1F* DEFHIST(dPhi_woSig_llmet_j);
    TH1F* DEFHIST(dPhi_woSig_met_l0);
    TH1F* DEFHIST(dPhi_woSig_met_l1);
    TH1F* DEFHIST(dPhi_woSig_met_ll);
    TH1F* DEFHIST(dPhi_woSig_met_j);
    TH1F* DEFHIST(dPhi_woSig_ll_j);
    TH1F* DEFHIST(dPhi_woSig_l0_j);
    TH1F* DEFHIST(dPhi_woSig_l1_j);
    TH1F* DEFHIST(dPhi_woSig_l0_l1);

    TH1F* DEFHIST(dR_llmet_j);

    TH2F* DEFHIST(l0_l1_pt);

    TH1F* DEFHIST(ll_M_dPhiReg);

    // Finely binned plots for Daniel
    TH1F* DEFHIST(ll_M_fine);
    TH1F* DEFHIST(ll_M_finer);
    
    TH1F* DEFHIST(l0_qeta);
    TH1F* DEFHIST(l1_qeta);
    TH1F* DEFHIST(mt_l0_met);
    TH1F* DEFHIST(mt_l1_met);
    TH1F* DEFHIST(mt_l_met_min);
    TH1F* DEFHIST(mct_top_tag);

    // Test
    //TH1F* h_met_test[Ch_N][PR_N];
    //TH1F* h_met_test2[Ch_N][PR_N];


    #undef DEFHIST

};

#endif
