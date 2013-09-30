////////////////////////////////////////////////////
// Will house the control regions used to measure //
// the fake rates used in 2-lep analysis.         //
////////////////////////////////////////////////////

#include "SusyTest0/MeasureFakeRate2.h"
#include "SusyTest0/criteria.h"

float Htbins[] = {1,3,6,9,12};
int nHtbins = 4;


/*--------------------------------------------------------------------------------*/
// Fake Rate Constructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::MeasureFakeRate2() :
  m_outFile(NULL),
  m_evtWeight(1.),
  m_AltIso(false),
  m_metRel(0.),
  m_ch(0),
  m_findOptCut(false)
{

}

/*--------------------------------------------------------------------------------*/
// Fake Rate Destructor
/*--------------------------------------------------------------------------------*/
MeasureFakeRate2::~MeasureFakeRate2()
{


}

/*--------------------------------------------------------------------------------*/
// Begin -- Called at the beginning of the event loop
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::Begin(TTree* /*tree*/)
{

  if(m_dbg) cout << "MeasureFakeRate2::Begin" << endl;

  SusySelectionMatt::Begin(0);
  // Output file name
  cout<<"Creating: "<<out<<endl;
  initHistos(m_fileName);


}

/*--------------------------------------------------------------------------------*/
// Terminate
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::Terminate()
{

  if(m_dbg) cout << "MeasureFakeRate2::Terminate" << endl;

  cout<<"Writing file "<<m_outFile<<endl;
  m_outFile->Write();
  cout<<"Closing file"<<endl;
  m_outFile->Close();
  dumpEventCounters();
}

/*--------------------------------------------------------------------------------*/
// Initialize the histograms
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::initHistos(string outName)
{

  cout<<"Creating file: "<<outName<<endl;
  m_outFile = new TFile((outName).c_str(),"recreate");
  m_outFile->cd();
  cout<<"File created: "<<m_outFile<<endl;


  // All of the rates will be stored in efficiency objects.
  // These are much more convenient for handling the errors.

  #define EFFVAR(eff,name,nbins,bins)					\
    do{									\
      eff = new EffObject(name,nbins,bins);				\
    }while(0)
  #define EFFVAR2(eff,name,nbinsx,xbins,nbinsy,ybins)				\
    do{									\
      eff = new EffObject2(name,nbinsx,xbins,nbinsy,ybins);		\
    }while(0)
  #define EFF(eff,name,nbins,xmin,xmax)					\
    do{									\
      eff = new EffObject(name,nbins,xmin,xmax);			\
    }while(0)
  #define EFF2(eff,name,nbins,xmin,xmax,nbinsy,ymin,ymax)			\
    do{									\
      eff = new EffObject2(name,nbins,xmin,xmax,nbinsy,ymin,ymax);		\
    }while(0)

  #define PROFILE(prof,name,nbins,xmin,xmax)					\
    do{									\
      prof = new TProfile(name.c_str(),name.c_str(),nbins,xmin,xmax);	\
    }while(0)

  #define LABEL(eff, bin, label)			\
    do{							\
      eff->SetXLabel(bin, label);			\
    }while(0)



  // Loop over lepton type
  for(int il=0; il<LT_N; ++il){
    string lName = LTNames[il];

    // For the control regions we need the
    // plots from both Data and MC
    for(int icr=0; icr<CR_N; ++icr){
      string cName = CRNames[icr];

      if(!m_findOptCut){

	// Loop over the channels
	for(int ich=0; ich<Ch_N; ++ich){
	  string chan = chanNames[ich];

	  string base = lName + "_" + cName + "_" + chan + "_";

	  EFFVAR(h_l_pt[il][icr][ich], (base+"l_pt"),nFakePtbins,FakePtbins);
	  EFFVAR(h_l_pt_coarse[il][icr][ich], (base+"l_pt_coarse"),nCoarseFakePtbins,coarseFakePtbins);
	  EFFVAR(h_l_pt_heavy[il][icr][ich], (base+"l_pt_heavy"),nFakePtbins,FakePtbins);
	  EFFVAR(h_l_pt_others[il][icr][ich], (base+"l_pt_others"),nFakePtbins,FakePtbins);
	  EFFVAR(h_l_eta[il][icr][ich], (base+"l_eta"),nEtabins,Etabins);
	  EFFVAR(h_l_eta_coarse[il][icr][ich], (base+"l_eta_coarse"),nCoarseEtabins,CoarseEtabins);
	  EFFVAR(h_metrel[il][icr][ich], (base+"metrel"),nMetbins, Metbins);
	  EFFVAR(h_met[il][icr][ich], (base+"met"),nMetbins, Metbins);
	  EFF(h_metrel_fine[il][icr][ich], (base+"metrel_fine"),40,0,200);
	  EFF(h_met_fine[il][icr][ich], (base+"met_fine"),40,0,200);
	  EFFVAR(h_metrel_coarse[il][icr][ich], (base+"metrel_coarse"),nCoarseMetbins, coarseMetbins);
	  EFFVAR(h_met_coarse[il][icr][ich], (base+"met_coarse"),nCoarseMetbins, coarseMetbins);
	  EFFVAR(h_njets[il][icr][ich], (base+"njets"),nJetbins, Jetbins);
	  EFFVAR(h_nlightjets[il][icr][ich], (base+"nlightjets"),nJetbins, Jetbins);
	  EFFVAR(h_nheavyjets[il][icr][ich], (base+"nheavyjets"),nJetbins, Jetbins);
	  EFFVAR(h_nlightjetsNoB[il][icr][ich], (base+"nlightjetsNoB"),nJetbins, Jetbins);
	  EFF(h_onebin[il][icr][ich], (base+"onebin"), 1, -0.5, 0.5);

	  // d0sig
	  EFF(h_heavy_d0sig[il][icr][ich], (base+"heavy_d0sig"), 50, 0, 5);
	  EFF(h_light_d0sig[il][icr][ich], (base+"light_d0sig"), 50, 0, 5);
	  EFF(h_conv_d0sig[il][icr][ich], (base+"conv_d0sig"), 50, 0, 5);

	  EFF(h_l_type[il][icr][ich], (base+"l_type"), nType, Typemin, Typemax);
	  EFF(h_l_origin[il][icr][ich], (base+"l_origin"), nOrigin, Originmin, Originmax);

	  // Flavor counting
	  EFF(h_flavor[il][icr][ich], (base+"flavor"), LS_N, -0.5, LS_N-0.5);
	  for(int lbl=0; lbl<LS_N; ++lbl) LABEL(h_flavor[il][icr][ich], lbl+1, LSNames[lbl]);

	  // 2 D param
	  EFFVAR2(h_l_pt_bjet[il][icr][ich], (base+"l_pt_bjet"),
		  nCoarseFakePtbins,coarseFakePtbins,nJetbins,Jetbins);
	  EFFVAR2(h_l_pt_eta[il][icr][ich], (base+"l_pt_eta"),
		  nCoarseFakePtbins,coarseFakePtbins,nCoarseEtabins,CoarseEtabins);

	  EFF(h_ht[il][icr][ich], (base+"ht"), 20, 0, 400);
	  EFFVAR(h_ht_pt[il][icr][ich], (base+"ht_pt"), nHtbins, Htbins);
	  EFF(h_ht_wMet[il][icr][ich], (base+"ht_wMet"), 20, 0, 400);
	  EFFVAR(h_ht_pt_wMet[il][icr][ich], (base+"ht_pt_wMet"), nHtbins, Htbins);

	  EFF(h_with_without_Etcone[il][icr][ich], (base+"with_without_Etcone"),3, -0.5, 2.5);

	}// end loop over channels

      }// end if !find opt cuts

      // Find Opt cuts
      else{

	string base = lName + "_" + cName + "_";

	EFF(h_ptcone[il][icr], (base+"ptcone"), nOptIsobins, OptIsomin, OptIsomax);
	EFF(h_etcone[il][icr], (base+"etcone"), nOptIsobins, OptIsomin, OptIsomax);
	EFF(h_d0Sig[il][icr], (base+"d0sig"), nOptSigbins, OptSigmin, OptSigmax);
	EFF(h_z0Sig[il][icr], (base+"z0sig"), nOptZ0bins, OptZ0min, OptZ0max);

	EFF(h_dist_d0Sig[il][icr], (base+"dist_d0sig"), nOptSigbins, OptSigmin, OptSigmax);
	EFF(h_dist_z0Sig[il][icr], (base+"dist_z0sig"), nOptZ0bins, OptZ0min, OptZ0max);
	EFF(h_dist_relptcone[il][icr], (base+"dist_relptcone"), 60, 0, 0.3);
	EFF(h_dist_ptcone[il][icr], (base+"dist_ptcone"), 100, 0, 10);
	EFF(h_dist_l_pt[il][icr], (base+"dist_l_pt"), 100, 0, 200);

	// Make 2-D plot of m(ll) vs ptcone(30)
	EFF2(h_dist_mll_ptcone[il][icr], (base+"mll_ptcone"), 20, 0, 200, 100, 0, 10);

	PROFILE(p_npv_etcone[il][icr], (base+"npv_etcone"), nNPVbins, NPVmin, NPVmax);
	PROFILE(p_npv_ptcone[il][icr], (base+"npv_ptcone"), nNPVbins, NPVmin, NPVmax);
	PROFILE(p_npv_ptconeElStyle[il][icr], (base+"npv_ptconeElStyle"), nNPVbins, NPVmin, NPVmax);

	PROFILE(p_mu_etcone[il][icr], (base+"mu_etcone"), nMUbins, MUmin, MUmax);
	PROFILE(p_mu_ptcone[il][icr], (base+"mu_ptcone"), nMUbins, MUmin, MUmax);
	PROFILE(p_mu_ptconeElStyle[il][icr], (base+"mu_ptconeElStyle"), nMUbins, MUmin, MUmax);

      }// end if find opt cuts


      // Saving CR plots to look at how
      // cuts affect distribution
      string base = lName + "_" + cName + "_all_";
      EFF(h_met_cr[il][icr], (base+"met_cr"), 50, 0, 200);
      EFF(h_mt_tag_cr[il][icr], (base+"mt_tag_cr"), 50, 0, 200);
      EFF(h_mt_probe_cr[il][icr], (base+"mt_probe_cr"), 50, 0, 200);
      EFF(h_ht_cr[il][icr], (base+"ht_cr"), 50, 0, 400);
      EFF(h_mll_cr[il][icr], (base+"mll_cr"), 30, 0, 300);


    }// end loop over control regions

  }// end loop over lepton types

  #undef PROFILE
  #undef EFFVAR
  #undef EFFVAR2
  #undef EFF
  #undef LABEL

}

/*--------------------------------------------------------------------------------*/
// Process
/*--------------------------------------------------------------------------------*/
void printProgress(const Susy::Event *e, Long64_t counter)
{
  cout << "**** Processing entry " <<counter
       << " run "<<e->run
       << " event "<<e->event << " ****" << endl;
}
//----------------------------------------------------------
Bool_t MeasureFakeRate2::Process(Long64_t entry)
{
  if(m_dbg) cout << "MeasureFakeRate2::Process" << endl;
  GetEntry(entry);
  clearObjects();
  m_chainEntry++;
  increment(n_readin);
  bool lepSf(true), btagSf(false); // computing btag can slow down
  if(m_dbg || m_chainEntry%50000==0) printProgress(nt.evt(), m_chainEntry);
  selectObjects(NtSys_NOM, false, TauID_medium);
  if( !selectEvent() ) return false;
  if(m_signalTaus.size()!=0)      /*increment(n_pass_CRWHSStauv  [m_ET]); else */ return false;
  if( m_baseLeptons.size() >= 2 ) increment(n_pass_atleast2Lep);
  if( m_baseLeptons.size() == 2 ) increment(n_pass_exactly2Lep);
  bool uniqueLepPair(m_baseLeptons.size() == 2);
  if(uniqueLepPair && !susy::passMllMin(m_baseLeptons, 20.0)) return false;
  // Loop over all regions (regardless of data or mc) and only fill relevant data quantitites.
  for(int cr = 0; cr<CR_N; ++cr){
    ControlRegion CR = (ControlRegion) cr;
    m_ch = Ch_all;
    m_probes.clear();
    m_tags.clear();
    bool passCR = false;
    // Data Driven CRs
    if( CR == CR_Real )          passCR = passRealCR(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_SideLow )  passCR = passRealCR(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_SideHigh ) passCR = passRealCR(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_HF )       passCR = passHFCR(m_baseLeptons,m_signalJets2Lep,m_met, CR);
    else if( CR == CR_HF_high )  passCR = passHFCR(m_baseLeptons,m_signalJets2Lep,m_met, CR);
    else if( CR == CR_LFZjet )   passCR = passLFZjetCR(m_baseLeptons,m_signalJets2Lep,m_met);
    else if( CR == CR_LFWjet )   passCR = passLFWjetCR(m_baseLeptons,m_signalJets2Lep,m_met);
    else if( CR == CR_Conv )     passCR = passConvCR(m_baseLeptons,m_signalJets2Lep,m_met);
    else if( CR == CR_CFlip )    passCR = passCFCR(m_baseLeptons,m_signalJets2Lep,m_met);
    // MC Regions
    else if( CR == CR_MCHeavy )  passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCLight )  passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCConv )   passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCQCD )    passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCALL )    passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCReal )   passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    else if( CR == CR_MCNone )   passCR = passMCReg(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    // Remaining signal and control regions
    else                         passCR = passSignalRegion(m_baseLeptons,m_signalJets2Lep,m_met,CR);
    if( passCR ){ // Plot Rates
      for(uint ip=0; ip<m_probes.size(); ++ip){
        if( m_findOptCut ) plotOptimumCuts(m_probes.at(ip), m_signalJets2Lep, m_met, CR);
        else               plotRates(m_probes.at(ip), m_signalJets2Lep, m_met, CR);
      } // end for(ip)
    }// end if Pass Control Region
  }// end loop over Control Regions
  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// Plotting Method
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::plotRates(const Lepton* lep, const JetVector& jets,
				 const Met* met, ControlRegion CR)
{

  LeptonType lt = lep->isEle() ? LT_EL : LT_MU;

  bool pass = m_AltIso ? passAltIso(lep) :
    isSignalLepton(lep, m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
    //isMySignalLepton(lep);


  #define FILL(eff, var)							\
    do{									\
      eff[lt][CR][m_ch]->Fill(pass,m_evtWeight,var);			\
      if(m_ch != Ch_all)						\
	eff[lt][CR][Ch_all]->Fill(pass,m_evtWeight,var);		\
    }while(0)
  #define FILL2(eff, varx, vary)						\
    do{									\
      eff[lt][CR][m_ch]->Fill(pass,m_evtWeight,varx,vary);			\
      if(m_ch != Ch_all)						\
	eff[lt][CR][Ch_all]->Fill(pass,m_evtWeight,varx,vary);		\
    }while(0)

  FILL(h_l_pt, lep->Pt());
  FILL(h_l_pt_coarse, lep->Pt());
  FILL(h_l_eta, fabs(lep->Eta()));
  FILL(h_l_eta_coarse, fabs(lep->Eta()));

  if(nt.evt()->isMC){
    if( isHFLepton(lep) ) FILL(h_l_pt_heavy, lep->Pt());
    else FILL(h_l_pt_others, lep->Pt());
  }

  FILL(h_metrel, m_metRel);
  FILL(h_met, met->Et);
  FILL(h_metrel_coarse, m_metRel);
  FILL(h_met_coarse, met->Et);
  FILL(h_metrel_fine, m_metRel);
  FILL(h_met_fine, met->Et);

  // Number of jets
  FILL(h_njets, jets.size());
  FILL(h_nlightjets, numberOfCLJets(jets));
  FILL(h_nheavyjets, numberOfCBJets(jets));
  if( numberOfCBJets(jets) == 0 )
    FILL(h_nlightjetsNoB, numberOfCLJets(jets));

  // Single bin
  FILL(h_onebin, 0);

  FILL2(h_l_pt_bjet, lep->Pt(), numberOfCBJets(jets));

  FILL2(h_l_pt_eta, lep->Pt(), fabs(lep->Eta()));

  // If the event is MC, save the flavor
  if( nt.evt()->isMC ){
    LeptonSource ls = getLeptonSource(lep);
    if(ls == LS_Real){
      FILL(h_l_type, lep->mcType);
      FILL(h_l_origin, lep->mcOrigin);
    }
    FILL(h_flavor, ls);
    if(ls==LS_HF||ls==LS_LF) FILL(h_flavor, LS_QCD);
    if(ls == LS_HF)   FILL(h_heavy_d0sig, lep->d0Sig(true));
    if(ls == LS_LF)   FILL(h_light_d0sig, lep->d0Sig(true));
    if(ls == LS_Conv) FILL(h_conv_d0sig, lep->d0Sig(true));
  }
  else{
    float d0sig = lep->d0Sig(true);
    FILL(h_heavy_d0sig, d0sig);
    FILL(h_light_d0sig, d0sig);
    FILL(h_conv_d0sig, d0sig);
  }

  // Ht Plots
  float ht = lep->Pt();
  for(uint ij=0; ij<jets.size(); ++ij)
    ht += jets.at(ij)->Pt();
  FILL(h_ht, ht);
  FILL(h_ht_pt, ht/lep->Pt());
  FILL(h_ht_wMet, ht+met->Et);
  FILL(h_ht_pt_wMet, (ht+met->Et)/lep->Pt());

  bool with = isSignalWithEtcone(lep);
  bool without = isSignalWithoutEtcone(lep);
  FILL(h_with_without_Etcone, 0);
  if(with)    FILL(h_with_without_Etcone, 1);
  if(without) FILL(h_with_without_Etcone, 2);


  #undef FILL
  #undef FILL2

}
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::plotCR(const Lepton* tag, const Lepton* probe, const JetVector& jets, const Met* met,
			      ControlRegion CR, CRPLOT CRP)
{



  LeptonType lt = probe->isEle() ? LT_EL : LT_MU;

  //bool pass = m_AltIso ? passAltIso(probe) : isMySignalLepton(probe);
  bool pass = m_AltIso ? passAltIso(probe) : isSignalLepton(probe,
							    m_baseElectrons,
							    m_baseMuons,
							    nt.evt()->nVtx,
							    nt.evt()->isMC);

  #define FILL(eff, var)							\
    do{									\
      eff[lt][CR]->Fill(pass,m_evtWeight,var);			\
    }while(0)

  if(CRP == CRP_mll)       FILL(h_mll_cr, (*tag+*probe).M());
  if(CRP == CRP_met)       FILL(h_met_cr, met->Et);
  if(CRP == CRP_mt_tag)    FILL(h_mt_tag_cr, Mt(tag, met));
  if(CRP == CRP_mt_probe)  FILL(h_mt_probe_cr, Mt(probe, met));

  if(CRP == CRP_ht){
    float ht = tag->Pt() + probe->Pt();
    for(uint i=0; i<jets.size(); ++i) ht += jets.at(i)->Pt();
    FILL(h_ht_cr, ht);
  }


  #undef FILL
}
/*--------------------------------------------------------------------------------*/
void MeasureFakeRate2::plotOptimumCuts(const Lepton* lep, const JetVector& jets, const Met* met,
				       ControlRegion CR)
{

  // Set the Lepton Type
  LeptonType lt = lep->isEle() ? LT_EL : LT_MU;

  // Define preprocessor command to keep clean
  #define FILL(eff, var, pass)			\
    do{						\
      eff[lt][CR]->Fill(pass,m_evtWeight,var);	\
    }while(0)

  // Define preprocessor command to keep clean
  #define FILL2(eff, varx, vary, pass)			\
    do{						\
      eff[lt][CR]->Fill(pass,m_evtWeight,varx,vary);	\
    }while(0)

  // Common var needed:
  float pt = lep->Pt();

  // Need to check everything except PtCone
  bool passOtherIso = isSignalWithoutPtcone(lep);
  float ptcone      = lep->isEle() ? lep->ptcone30 : ((Muon*) lep)->ptcone30ElStyle;
  for(float cut = 0.005; cut<OptIsomax; cut+=0.01)
    //FILL(h_ptcone, cut, (passOtherIso && ptcone/pt < cut));
    FILL(h_ptcone, cut, (ptcone/pt < cut));

  // Need to check everything except EtCone
  passOtherIso = isSignalWithoutEtcone(lep);
  float etcone = lep->isEle() ? ((Electron*) lep)->topoEtcone30Corr : ((Muon*) lep)->etcone30;
  for(float cut = 0.005; cut<OptIsomax; cut+=0.01)
    //FILL(h_etcone, cut, (passOtherIso && etcone/pt < cut));
    FILL(h_etcone, cut, (etcone/pt < cut));

  // Need to check everything except d0Sig
  //passOtherIso = isSignalWithoutd0Sig(lep);
  passOtherIso = isSignalWithoutIP(lep);
  float d0Sig  = fabs(lep->d0Sig(true));
  for(float cut = 0.1; cut<OptSigmax; cut+=0.2)
    FILL(h_d0Sig, cut, ( d0Sig < cut));

  // Check passes all else
  passOtherIso = isSignalWithoutIP(lep);
  float z0sintheta = fabs( lep->z0SinTheta(true) );
  for(float cut = 0.05; cut<OptZ0max; cut+=0.1)
    FILL(h_z0Sig, cut, ( z0sintheta < cut) );

  // Fill the distributions
  FILL(h_dist_relptcone, ptcone/pt, true);
  FILL(h_dist_ptcone, ptcone, true);
  FILL(h_dist_l_pt, pt, true);
  FILL(h_dist_d0Sig, d0Sig, true);
  FILL(h_dist_z0Sig, z0sintheta, true);

  if( m_baseLeptons.size() == 2 )
    FILL2(h_dist_mll_ptcone, Mll(m_baseLeptons[0],m_baseLeptons[1]), ptcone, true);

  #undef FILL2
  #undef FILL

  //
  // TProfiles to check stability
  //

  #define FILLP(prof, varx, vary)		\
    do{							\
      prof[lt][CR]->Fill(varx,vary,m_evtWeight);		\
    }while(0)

  // Stability vars
  float npv = nt.evt()->nVtx;
  float mu  = nt.evt()->avgMu;

  FILLP(p_npv_etcone, npv, getEtcone(lep));
  FILLP(p_npv_ptcone, npv, getPtcone(lep));
  FILLP(p_npv_ptconeElStyle, npv, getPtcone(lep,true));
  FILLP(p_mu_etcone, mu, getEtcone(lep));
  FILLP(p_mu_ptcone, mu, getPtcone(lep));
  FILLP(p_mu_ptconeElStyle, mu, getPtcone(lep,true));


  #undef FILLP

}

/*--------------------------------------------------------------------------------*/
// Control Regions
/*--------------------------------------------------------------------------------*/

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// MC Control region is defined by metrel
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passMCReg(const LeptonVector &leptons,
				 const JetVector &jets,
				 const Met* met,
				 ControlRegion CR)
{

  // Used to measure fake rate in MC in a region
  // that is close to our signal region
  // * 40 < metrel < 100
  // * exactly two leptons
  // * probes are the fake leptons of specific type

  if( !nt.evt()->isMC )     return false;
  if( leptons.size() != 2 ) return false;

  m_ch = getChan(leptons);
  m_metRel = getMetRel(met,leptons,jets);
  //if( !(CR != CR_MCNone && 40 < m_metRel && m_metRel < 100 ) ) return false;
  //if( !(CR != CR_MCNone && 40 < m_metRel  ) ) return false;
  //if( !(CR != CR_MCNone && 20 < met->Et  ) ) return false;
  //if( !(CR != CR_MCNone && 40 < met->Et  ) ) return false;

  for(uint il=0; il<leptons.size(); ++il){

    // Heavy
    if( CR == CR_MCHeavy && isHFLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );

    // Light
    if( CR == CR_MCLight && isLFLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );

    // Conversion
    if( CR == CR_MCConv && isConvLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );

    // None CR
    if( CR == CR_MCNone && isFakeLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );

    // QCD
    if( CR == CR_MCQCD && isQCDLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );

    // All fakes <-- Will include reals
    if( CR == CR_MCALL )
      m_probes.push_back( leptons[il] );

    // Real
    if( CR == CR_MCReal && isRealLepton(leptons[il]) )
      m_probes.push_back( leptons[il] );
  }

  m_evtWeight = getEvtWeight(leptons);

  return true;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// The signal region checks
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passSignalRegion(const LeptonVector &leptons,
					const JetVector &jets,
					const Met* met,
					ControlRegion CR)
{

  if( leptons.size() != 2 ) return false;

  bool passSR = false;
  if( CR == CR_SRmT2a )          passSR = passSRmT2a(leptons,jets,met,false,true);
  else if( CR == CR_SRmT2b )     passSR = passSRmT2b(leptons,jets,met,false,true);
  else if( CR == CR_SRmT2c )     passSR = passSRmT2c(leptons,jets,met,false,true);
  else if( CR == CR_SRWWa )      passSR = passSRWWa(leptons,jets,met,false,true);
  else if( CR == CR_SRWWb )      passSR = passSRWWb(leptons,jets,met,false,true);
  else if( CR == CR_SRWWc )      passSR = passSRWWc(leptons,jets,met,false,true);
  else if( CR == CR_SRZjets )    passSR = passSRZjets(leptons,jets,met,false,true);

  else if( CR == CR_VRSS )         passSR = passVRSS(leptons,jets,met);
  else if( CR == CR_SSInc)         passSR = sameSign(leptons);

  else if( CR == CR_CRWWMet )      passSR = passCRWWMet(leptons,jets,met,false,true);
  else if( CR == CR_CRWWmT2 )      passSR = passCRWWmT2(leptons,jets,met,false,true);
  else if( CR == CR_CRTopMet )     passSR = passCRTopMet(leptons,jets,met,false,true);
  else if( CR == CR_CRTopmT2 )     passSR = passCRTopmT2(leptons,jets,met,false,true);
  else if( CR == CR_CRZVMet )      passSR = passCRZVMet(leptons,jets,met,false,true);
  else if( CR == CR_CRZVmT2_90 )   passSR = passCRZVmT2a(leptons,jets,met,false,true);
  else if( CR == CR_CRZVmT2_120 )  passSR = passCRZVmT2b(leptons,jets,met,false,true);
  else if( CR == CR_CRZVmT2_150 )  passSR = passCRZVmT2c(leptons,jets,met,false,true);
  else if( CR == CR_CRZVmT2_100 )  passSR = passCRZVmT2d(leptons,jets,met,false,true);

  else if( CR == CR_CRTopZjets )   passSR = passCRTopZjets(leptons,jets,met,false,true);
  else if( CR == CR_CRZXZjets )    passSR = passCRZXZjets(leptons,jets,met,false,true);

  else if( CR == CR_PremT2 )       passSR = passCRPremT2(leptons,jets,met);
  else if( CR == CR_SRWHSS )       passSR = passWhSS(leptons, jets, met);

  for(uint i=0; i<leptons.size(); ++i)
    //if( isFakeLepton(leptons[i]) )
    m_probes.push_back( leptons[i] );

  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);
  m_ch = getChan(leptons);
  return passSR;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Real Control Region: Z
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passRealCR(const LeptonVector &leptons,
				  const JetVector &jets,
				  const Met* met,
				  ControlRegion CR)
{

  // Real CRA:
  // * Require Exactly two baseline leptons SF
  // * Require single lepton triggers
  // * Require fall in Z mass window
  // * Require at least one tight lepton for tag

  // Same flavor dilepton
  if( leptons.size() != 2 )  return false;
  if( !sameFlavor(leptons) ) return false;

  // In Z window or side bands
  float mll = Mll(leptons[0],leptons[1]);
  if(CR == CR_Real){
    if( fabs(mll-91.2) > 10 )       return false;
  }
  else if(CR == CR_SideLow){
    if( !(61 < mll && mll < 71) )   return false;
  }
  else if(CR == CR_SideHigh){
    if( !(111 < mll && mll < 121) ) return false;
  }

  // At least one tight lepton
  int nVtx   = nt.evt()->nVtx;
  bool isMC  = nt.evt()->isMC;
  //bool l0sig = m_AltIso ? passAltIso(leptons[0]) : isMySignalLepton(leptons[0]);
  bool l0sig = m_AltIso ? passAltIso(leptons[0]) : isSignalLepton(leptons[0],
								  m_baseElectrons,
								  m_baseMuons,
								  nVtx, isMC);
  //bool l1sig = m_AltIso ? passAltIso(leptons[1]) : isMySignalLepton(leptons[1]);
  bool l1sig = m_AltIso ? passAltIso(leptons[1]) : isSignalLepton(leptons[1],
								  m_baseElectrons,
								  m_baseMuons,
								  nVtx, isMC);
  if( !l0sig && !l1sig ) return false;

  // Pass trigger
  uint l0flag = leptons[0]->trigFlags;
  bool l0trig = leptons[0]->isEle() ? l0flag & TRIG_e24vhi_medium1 : l0flag & TRIG_mu24i_tight;
  uint l1flag = leptons[1]->trigFlags;
  bool l1trig = leptons[1]->isEle() ? l1flag & TRIG_e24vhi_medium1 : l1flag & TRIG_mu24i_tight;

  if( l0sig && l0trig && leptons[0]->Pt() > 25 )
    m_probes.push_back( leptons[1] );

  if( l1sig && l1trig && leptons[1]->Pt() > 25 )
    m_probes.push_back( leptons[0] );

  if( m_probes.size() == 0 ) return false;

  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);

  return true;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Heavy Flavor Control region
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passHFCR(const LeptonVector &leptons,
				const JetVector &jets,
				const Met* met,
				ControlRegion CR)
{

  // This is a heavy flavor tag and probe trying
  // to select b-bbar events by looking at loose
  // muon tag in the event.
  // * Require Muon Stream Only
  // * Exactly one muon overlapping exactly one b-jet
  // * Tag muon Pt > 20 and match to dilepton trigger
  //    -- EF_mu18_tight_e7_medium (elec probe)
  //    -- EF_mu18_tight_mu8_EFFS  (muon probe)
  // * Probe lepton baseline after overlap
  // * Mt(probe,met) < 40 GeV
  // * met < 40 GeV
  // * Veto dimuon events 10 GeV near Z
  // * m(ll) > 40 GeV (probe is muon)

  // Check pre muons and jets
  MuonVector preMuons = getPreMuons(&nt, NtSys_NOM);
  if( preMuons.size() == 0 ) return false;
  if( jets.size() == 0 )     return false;

  // Find tag
  Lepton* tag = NULL;
  int nTags   = 0;
  int nBjets  = 0;
  for(uint ij = 0; ij<jets.size(); ++ij){
    if( !isCentralBJet(jets.at(ij)) ) continue;
    nBjets++;
    for(uint imu=0; imu<preMuons.size(); ++imu){
      if(jets.at(ij)->DeltaR( *preMuons.at(imu) ) < 0.4){
	tag = (Lepton*) preMuons.at(imu);
	nTags++;
      }// end if overlap
    }// end loop over muons
  }// end loop over jets

  if( nTags != 1 )  return false;
  if( nBjets != 1 ) return false;

  // Get the probe
  Lepton* probe = NULL;
  uint nProbes  = 0;
  for(uint il=0; il<leptons.size(); ++il){
    if(leptons.at(il)->isEle()){
      probe = leptons.at(il);
      nProbes++;
    }
    else if( tag != leptons.at(il) ){
      probe = leptons.at(il);
      nProbes++;
    }
  }// end loop over baseline leptons

  if( nProbes != 1 ) return false;

  // Check mll
  plotCR(tag,probe,jets,met,CR_HF,CRP_mll);
  bool isMM = tag->isMu() && probe->isMu();
  if(isMM){
    float mll = Mll(tag,probe);
    float mZ  = 91.2;
    if( fabs(mll-mZ) < 10 ) return false;
    if( mll < 40 )          return false;
  }

  // Check the trigger
  if( tag->Pt() < 20 ) return false;
  uint tagFlag = tag->trigFlags;
  if( !isMM && !((tagFlag & TRIG_mu18_tight) && (tagFlag & TRIG_mu18_tight_e7_medium1)) )
    return false;
  if( isMM && !((tagFlag & TRIG_mu18_tight) && (tagFlag & TRIG_mu18_tight_mu8_EFFS)) )
    return false;

  // Met and Mt Cut
  plotCR(tag,probe,jets,met,CR_HF,CRP_ht);
  plotCR(tag,probe,jets,met,CR_HF,CRP_met);
  if( met->Et > 40 )       return false;

  plotCR(tag,probe,jets,met,CR_HF,CRP_mt_tag);
  plotCR(tag,probe,jets,met,CR_HF,CRP_mt_probe);
  if( CR == CR_HF )
    if( Mt(probe,met) > 40 ) return false;
  if( CR == CR_HF_high )
    if( Mt(probe,met) > 100 ) return false;



  // We have survived so save objects and return
  LeptonVector temp; temp.push_back(tag); temp.push_back(probe);

  m_metRel = getMetRel(met, temp, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(temp, true);

  m_probes.push_back( probe );
  m_tags.push_back( tag );

  return true;


}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Light Flavor Control Region -- Zjet
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passLFZjetCR(const LeptonVector &leptons,
				    const JetVector &jets,
				    const Met* met)
{

  // Light flavor tag and probe using Z+jet events
  // * Same-flavor OS signal leptons 10 GeV of Z peak
  // * Probe lepton of different flavor
  // * M(lll) not within 10 GeV of Z peak
  // * met < 40 GeV
  // * Mt(probe,met) < 50 GeV
  // * Reject events with b-jet
  // * Pass LF Trig function

  // Require exactly 3 leptons
  if( leptons.size() != 3 )       return false;

  // Reject events with b-jet
  if( numberOfCBJets(jets) != 0 ) return false;

  // Use Signal info since it's cleaner
  LeptonVector sigLeps;
  if(m_signalMuons.size() == 2 && m_baseElectrons.size() == 1){
    if(m_signalMuons[0]->q * m_signalMuons[1]->q < 0){
      sigLeps.push_back( (Lepton*) m_signalMuons[0]);
      sigLeps.push_back( (Lepton*) m_signalMuons[1]);
      m_probes.push_back( (Lepton*) m_baseElectrons[0]);
    }
  }
  if(m_baseMuons.size() == 1 && m_signalElectrons.size() == 2){
    if(m_signalElectrons[0]->q * m_signalElectrons[1]->q < 0){
      sigLeps.push_back( (Lepton*) m_signalElectrons[0]);
      sigLeps.push_back( (Lepton*) m_signalElectrons[1]);
      m_probes.push_back( (Lepton*) m_baseMuons[0]);
    }
  }
  if(sigLeps.size() != 2)    return false;

  if( !passLFTrig(sigLeps) ) return false;

  // Mass rejections
  float mZ   = 91.2;
  float mll  = Mll(sigLeps[0],sigLeps[1]);
  float mlll = Mlll(sigLeps[0],sigLeps[1],m_probes[0]);

  if( fabs(mll-mZ) > 10 )  return false;
  if( fabs(mlll-mZ) < 10 ) return false;

  // Mt Cut
  if( Mt(m_probes[0],met) > 50 ) return false;

  // Update and save objects
  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(sigLeps, true) * m_probes[0]->effSF;

  return true;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Light Flavor Control Region -- Wjet
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passLFWjetCR(const LeptonVector &leptons,
				    const JetVector &jets,
				    const Met* met)
{

  // Light Flavor control region for W+jet
  // * Two same-sign same-flavor leptons
  // * Tag is leading probe is subleading
  // * Tag must be signal and match to tight trigger
  // * Veto mll region 70-100 get rid of Zs
  // * Mt(tag,met) > 40
  // * Met/HT > 0.4
  // * Met > 20 GeV
  // * No b jets

    // Select two baseline leptons of SF SS
  // Reject events with bJet
  if( leptons.size() != 2 )     return false;
  if( oppositeFlavor(leptons) ) return false;
  if( oppositeSign(leptons) )   return false;

  // Define tag as leading lepton.
  // This is very nearly true for W+jet
  Lepton* tag = leptons[0];
  Lepton* probe = leptons[1];

  // Require tag to be signal lepton
  //if( !isMySignalLepton(tag) ) return false;
  if( !isSignalLepton(tag,
		      m_baseElectrons,
		      m_baseMuons,
		      nt.evt()->nVtx, nt.evt()->isMC) ) return false;

  // Make sure lepton passes single trigger
  if( tag->isEle() ){
    if( !(tag->Pt() > 25 && tag->matchTrig(TRIG_e24vhi_medium1)) ) return false;
  }
  else if(tag->isMu()){
    if( !(tag->Pt() > 25 && tag->matchTrig(TRIG_mu24i_tight)) ) return false;
  }

  // First cut: Reject Z
  float mll = Mll(tag,probe);
  plotCR(tag,probe,jets,met,CR_LFWjet,CRP_mll);
  //if( fabs(mll - 91.2) < 10 ) return false;
  if( 70 < mll && mll < 100 ) return false;

  // Second Cut: Mt to reduce Z and enhace W+jet
  plotCR(tag,probe,jets,met,CR_LFWjet,CRP_mt_tag);
  plotCR(tag,probe,jets,met,CR_LFWjet,CRP_mt_probe);
  if( Mt(tag,met) < 40 ) return false;

  // Thid cut: MET/HT shown to separate W+jet
  // and all other backgrounds
  float Ht = tag->Pt() + probe->Pt();
  for(uint i=0; i<jets.size(); ++i) Ht += jets.at(i)->Pt();
  plotCR(tag,probe,jets,met,CR_LFWjet,CRP_ht);
  if( met->Et/Ht < 0.4 ) return false;

  // Fourth Cut: Met > 20 to get rid of region with
  // bad agreement..
  plotCR(tag,probe,jets,met,CR_LFWjet,CRP_met);
  if( met-> Et < 20 ) return false;


  // Fifth Cut: Reject events with B jet
  //if( numberOfCBJets(jets) != 0 ) return false;

  // We made it!
  m_metRel = getMetRel(met,leptons,jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons, true);

  m_probes.push_back(probe);

  return true;


}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Conversion Control Region -- Conv
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
bool MeasureFakeRate2::passConvCR(const LeptonVector &leptons,
				  const JetVector &jets,
				  const Met* met)
{

  // Conversion control region to get electron conversion rate
  // * Data only from Muon Stream
  // * Exactly two signal muons of opposite sign
  // * Exactly one baseline electron
  // * M(lll) in enlarged Z window 80-100
  // * Met < 50 GeV
  // * M(mu,mu) > 20
  // * Mt(elec, met) < 40

  ElectronVector preElecs; MuonVector preMuons;
  if( m_signalMuons.size() != 2 )  return false;
  if( m_baseElectrons.size() != 1) return false;
  if( m_baseLeptons.size() != 3)   return false;
  preMuons.push_back( m_signalMuons[0] );
  preMuons.push_back( m_signalMuons[1] );
  preElecs.push_back( m_baseElectrons[0] );


  float mlll = (*preMuons[0] + *preMuons[1] + *preElecs[0]).M();
  if( !(80 < mlll && mlll < 100) )           return false;

  // Opposite sign
  if( preMuons[0]->q * preMuons[1]->q > 0 )  return false;

  // Pt and trigger requirement
  LeptonVector tempL; ElectronVector tempE;
  buildLeptons(tempL, tempE, preMuons);
  if( !passTrigger(tempL) )                  return false;

  // Met < 50
  if( met->Et > 50 )                         return false;

  // M(mu,mu) > 20
  if( (*preMuons[0]+*preMuons[1]).M() < 20 ) return false;

  // Mt(elec, met) < 40
  if( Mt((Lepton*) preElecs[0], met) > 40)   return false;

  m_probes.push_back( (Lepton*) preElecs[0] );

  m_metRel = getMetRel(met, leptons, jets); // Not useful quantity for this region.

  float eff = nt.evt()->isMC ? m_probes[0]->effSF : 1.0;
  m_evtWeight = getEvtWeight(tempL,true,true) * eff;

  return true;


}

/*--------------------------------------------------------------------------------*/
// Pass Charge flip control region
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::passCFCR(const LeptonVector &leptons,
				const JetVector& jets,
				const Met* met)
{

  // Real CRA:
  // * Require Exactly two baseline leptons SF
  // * Require single lepton triggers
  // * Require fall in Z mass window
  // * Require at least one tight lepton for tag
  // * Require to be same-sign

  // Same flavor dilepton
  if( leptons.size() != 2 )  return false;
  if( !sameFlavor(leptons) ) return false;
  if( !sameSign(leptons) )   return false;

  // In Z window or side bands
  float mll = Mll(leptons[0],leptons[1]);
  if( !(70 < mll && mll < 100) ) return false;

  // At least one tight lepton
  int nVtx   = nt.evt()->nVtx;
  bool isMC  = nt.evt()->isMC;
  //bool l0sig = m_AltIso ? passAltIso(leptons[0]) : isMySignalLepton(leptons[0]);
  bool l0sig = m_AltIso ? passAltIso(leptons[0]) : isSignalLepton(leptons[0],
  m_baseElectrons,
  m_baseMuons,
  nVtx, isMC);
  //bool l1sig = m_AltIso ? passAltIso(leptons[1]) : isMySignalLepton(leptons[1]);
  bool l1sig = m_AltIso ? passAltIso(leptons[1]) : isSignalLepton(leptons[1],
  m_baseElectrons,
  m_baseMuons,
  nVtx, isMC);
  if( !l0sig && !l1sig ) return false;

  // Pass trigger
  uint l0flag = leptons[0]->trigFlags;
  bool l0trig = leptons[0]->isEle() ? l0flag & TRIG_e24vhi_medium1 : l0flag & TRIG_mu24i_tight;
  uint l1flag = leptons[1]->trigFlags;
  bool l1trig = leptons[1]->isEle() ? l1flag & TRIG_e24vhi_medium1 : l1flag & TRIG_mu24i_tight;

  if( l0sig && l0trig && leptons[0]->Pt() > 25 )
    m_probes.push_back( leptons[1] );

  if( l1sig && l1trig && leptons[1]->Pt() > 25 )
    m_probes.push_back( leptons[0] );

  if( m_probes.size() == 0 ) return false;

  m_metRel = getMetRel(met, leptons, jets);
  if( nt.evt()->isMC ) m_evtWeight = getEvtWeight(leptons);

  return true;


}

/*--------------------------------------------------------------------------------*/
// Light Flavor Trigger Requirements
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::passLFTrig(const LeptonVector& leps)
{
  // Try LF control region with the same trigger scheme
  // that 3-lepton uses

  if( leps.size() != 2 ) return false;
  uint ie=0;
  uint im=0;
  for(uint i=0; i<leps.size(); ++i){
    if(leps[0]->isEle()) ie++;
    else im++;
  }

  bool passMuons = false;
  uint n1M = 0;
  uint nSym2M = 0;
  uint nAsym2M = 0;
  uint nAsym2M_m18 = 0;
  // Trigger matching for Muon Stream
  for(uint i=0; i<leps.size(); i++){
    if(leps.at(i)->isEle()) continue;
    const Muon* lep = (Muon*) leps[i];
    float pt = lep->Pt();
    // Single muon trigger
    if(pt>25 && lep->matchTrig(TRIG_mu24i_tight)) n1M++;
    // 2m symmetric trigger
    if(pt>14 && lep->matchTrig(TRIG_2mu13)) nSym2M++;
    // 2m asymmetric trigger
    if(lep->matchTrig(TRIG_mu18_tight_mu8_EFFS)){
      nAsym2M++;
      if(pt>18 &&  lep->matchTrig(TRIG_mu18_tight)) nAsym2M_m18++;
    }
  }
  // Logical OR or all Muon Triggers
  if(n1M || nSym2M || nAsym2M || nAsym2M_m18) passMuons = true;

  bool passEgamma = false;
  uint n1E = 0;
  uint nSym2E = 0;
  uint nAsym2E = 0;
  uint nAsym2E_e24 =0;
  // Trigger matching for Electron Stream
  for(uint i=0; i<leps.size(); i++){
    if(!leps.at(i)->isEle()) continue;
    const Electron* lep = (Electron*) leps[i];
    float pt = lep->Pt();
    // Single electron trigger
    if(pt>25 && lep->matchTrig(TRIG_e24vhi_medium1)) n1E++;
    // 2e symmetric trigger
    if(pt>14 && lep->matchTrig(TRIG_2e12Tvh_loose1)) nSym2E++;
    // 2e asymmetric trigger
    if(lep->matchTrig(TRIG_e24vh_medium1_e7_medium1)){
      nAsym2E++;
      if(pt>25 && lep->matchTrig(TRIG_e24vh_medium1)) nAsym2E_e24++;
    }
  }
  // Logical OR or all Electron triggers
  if(n1E || nSym2E || nAsym2E || nAsym2E_e24) passEgamma = true;

  if( !nt.evt()->isMC ){
    DataStream stream = nt.evt()->stream;
    if( im == 2 && stream == Stream_Muons) return passMuons;
    if( ie == 2 && stream == Stream_Egamma) return passEgamma;
    //if( im == 2 )  return passMuons;
    //if( ie == 2 )  return passEgamma;

    return false;
  }

  return (passMuons && im==2) || (passEgamma && ie ==2);
}

/*--------------------------------------------------------------------------------*/
// Alternative Isolations to test
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::passAltIso(const Lepton* lepton)
{

  uint nvtx = nt.evt()->nVtx;
  bool ismc = nt.evt()->isMC;

  // Muons
  if( lepton->isMu() ){
    Muon* mu = (Muon*) lepton;

    if( mu->ptcone30ElStyle/mu->Pt() > 0.12 ) return false;
    if( fabs(mu->d0Sig()) >= MUON_D0SIG_CUT ) return false;
    if( fabs(mu->z0SinTheta()) >= MUON_Z0_SINTHETA_CUT) return false;

    return true;

  }

  // Electrons
  if( lepton->isEle() ){
    Electron* el = (Electron*) lepton;
    bool passNom = isSignalLepton(lepton, m_baseElectrons, m_baseMuons, nvtx, ismc);
    bool passTightd0Sig = fabs(el->d0Sig()) < 3;
    return passNom && passTightd0Sig;

    /*
    Electron* el = (Electron*) lepton;
    if( !el->tightPP ) return false;

    float ptiso = 0.08;
    if( el->ptcone30/el->Pt() >= ptiso ) return false;

    float etiso = 0.09;
    if( elEtTopoConeCorr(el, m_baseElectrons, m_baseMuons, nvtx, ismc)/el->Pt() >= etiso ) return false;

    if( fabs(el->d0Sig()) >= ELECTRON_D0SIG_CUT )            return false;
    if( fabs(el->z0SinTheta()) >= ELECTRON_Z0_SINTHETA_CUT ) return false;

    return true;
    */

  }

  return false;

}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithEtcone(const Lepton* lep)
{

  bool isSignal = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);

  // If Electron, just return
  if( lep->isEle() ) return isSignal;

  // Otherwise muon, add etcone cut:
  return isSignal && ((Muon*) lep)->etcone30/lep->Pt() < 0.18;

}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutEtcone(const Lepton* lep)
{

  // If Electron, just return
  if( lep->isMu() )
    return isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);

  // Otherwise, turn off for electrons
  m_doElEtconeCut = false;
  bool isSignal = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doElEtconeCut = true;

  return isSignal;

}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutPtcone(const Lepton* lep){

  m_doPtconeCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doPtconeCut = true;
  return pass;

}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutd0Sig(const Lepton* lep){

  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;

  float cutVal = lep->isEle() ? ELECTRON_Z0_SINTHETA_CUT : MUON_Z0_SINTHETA_CUT;
  bool passz0  = fabs( lep->z0SinTheta(true) ) < cutVal;

  return pass && passz0;

}
/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutz0Sig(const Lepton* lep){

  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;

  float cutVal = lep->isEle() ? ELECTRON_D0SIG_CUT : MUON_D0SIG_CUT;
  bool passd0  = fabs( lep->d0Sig(true) ) < cutVal;

  return pass && passd0;

}

/*--------------------------------------------------------------------------------*/
bool MeasureFakeRate2::isSignalWithoutIP(const Lepton* lep){

  m_doIPCut = false;
  bool pass = isSignalLepton(lep,m_baseElectrons,m_baseMuons,nt.evt()->nVtx,nt.evt()->isMC);
  m_doIPCut = true;

  return pass;

}

