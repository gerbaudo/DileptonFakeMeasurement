#ifndef SusyAna_MatrixPrediction_h
#define SusyAna_MatrixPrediction_h

#include "DileptonFakeMeasurement/SusyPlotter.h"
#include "DileptonFakeMeasurement/PlotRegion.h"
#include "SusyMatrixMethod/FakeRegions.h"
#include "SusyMatrixMethod/MatrixLepton.h"
#include "SusyMatrixMethod/DiLeptonMatrixMethod.h"

#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

/*

    Matrix Prediction - Looper to run over data to make relevant histograms
    Histograms defined in SusyPlotter, maybe define extras in this class 
    relevant to the matrix method

*/

enum MatrixPair{ MP_TT = 0, MP_TL, MP_LT, MP_LL, MP_ALL, MP_N };
static string MPNames[] = {"TT","TL","LT","LL","ALL"};

enum WeightToggle { WT_ON = 0, WT_OFF, WTog_N };
static string WTNames[] = {"wOn","wOff"};

const int nFakeSys = 21;
static string FAKESYSNames[nFakeSys] =
{
  "NONE",
  "ALL_UP",
  "ALL_DN",
  "EL_RE_UP",
  "EL_RE_DOWN",
  "EL_FR_UP",
  "EL_FR_DOWN",
  "MU_RE_UP",
  "MU_RE_DOWN",
  "MU_FR_UP",
  "MU_FR_DOWN",
  "EL_SSOS_UP",
  "EL_SSOS_DOWN",
  "EL_WJETS",
  "EL_TTBAR",
  "EL_METREL",
  "MU_SSOS_UP",
  "MU_SSOS_DOWN",
  "MU_WJETS",
  "MU_TTBAR",
  "MU_METREL"
};

class MatrixPrediction : public SusyPlotter
{

  public:

    MatrixPrediction();
    virtual ~MatrixPrediction(){};
    virtual void    Begin(TTree *tree); //!< called before looping on entries
    virtual void    Terminate(); //!< called after looping is finished
    virtual Bool_t  Process(Long64_t entry);
    //! Extra histograms and plotting function specific to matrix method
    bool bookFakeHisto();
    void fillFakeHistos(const LeptonVector &baseLeps, const JetVector &jets, 
                        const Met* met,float weight, size_t regionIndex, uint sys);
    // Get the fake event weight given a signal region
    float getFakeWeight(const LeptonVector &baseLeps, susy::fake::Region region, float metRel, SusyMatrixMethod::SYSTEMATIC sys = SusyMatrixMethod::SYS_NOM);
    MatrixPair getMatrixPair(const LeptonVector &baseLeps); //!< Get the Matrix Pair type
    MatrixPrediction& setMatrixFilename(const std::string filename); // to be called before Begin
    MatrixPrediction& use2dParametrization() { m_use2dparametrization = true; return *this; }
    static std::string dilepDetails(const Susy::Event &event, const DiLepEvtType &ll,
                                    const LeptonVector &ls);
    std::string eventDetails(bool passSrSs, const Susy::Event &event, const DiLepEvtType &ll,
                             const LeptonVector &ls);
    ClassDef(MatrixPrediction, 3);

  protected:
    std::vector<uint> m_matrixSys;      // systematics to process
    //    SusyMatrixMethod::DiLeptonMatrixMethod* m_matrix;
    susy::fake::DileptonMatrixMethod* m_matrix;
    std::string m_matrixFilename;
    bool m_use2dparametrization;
    bool m_allconfigured;
    bool initMatrixTool();
    // Histograms
#define FAKEHIST( name ) hf_ ## name[susy::wh::Ch_N][susy::wh::kNumberOfPlotRegions][MP_N][WTog_N][SusyMatrixMethod::SYS_N_USER];
    TH1F* FAKEHIST( ll_M );
    TH1F* FAKEHIST( l0_pt );
    TH1F* FAKEHIST( l1_pt );

    TH1F* FAKEHIST( met );
    TH1F* FAKEHIST( metrel );
    
    TH1F* FAKEHIST( onebin );
    TH1F* FAKEHIST( evt_weight );

    TH2F* FAKEHIST( pt0vspt1 );

    TH1F* FAKEHIST( met_l0_Mt );
    TH1F* FAKEHIST( met_l1_Mt );

    TH1F* FAKEHIST( njets );
    TH1F* FAKEHIST( bjet_pt );
    TH1F* FAKEHIST( ljet_pt );
#undef FAKEHIST
};

#endif
