// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FastTopJets.hh"
#include "Rivet/Projections/FastHiggsJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/PartonicHiggs.hh"
#include "Rivet/Projections/Sphericity.hh"
#include <fastjet/tools/Filter.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/EnergyCorrelator.hh>
#include <fstream>
#include <random>

using namespace fastjet::contrib;


namespace Rivet {
  

  class FatHiggsTagging : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(FatHiggsTagging);
    
    
    /// Book projections and histograms
    void init() {
      
      // get process string 
      _process="UNKNOWN";      
      std::ifstream pfile ("process", std::ifstream::in);      
      getline(pfile, _process);
      
      // Random numbers for simulation of b-tagging efficiencies
      srand (160385);
      
      // Base final state definition
      const FinalState fs(Cuts::pT > 0.1 * GeV && Cuts::abseta < 5.0);
      // charged final states for track jets 
      //agrohsje 
      const ChargedFinalState cfs(Cuts::pT > 0.1 * GeV && Cuts::abseta < 5.0);
            
      // Small-R track jets, 
      //agrohsje include muons ? 
      FastJets trackjets(cfs, FastJets::ANTIKT, 0.2, JetAlg::DECAY_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(trackjets, "trackjets");

      // Large-R jets
      FastJets fatcalojets(fs, FastJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(fatcalojets, "fatcalojets");
      
      // with top tagging 
      FastTopJets fattopjets(fs, FastTopJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(fattopjets, "fattopjets");

      // with higgs and b hadron decay tagging 
      FastHiggsJets fathiggsjets(fs, FastHiggsJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(fathiggsjets, "fathiggsjets");
      
      // identify hadronic top events (no intermediate tau decays) 
      addProjection(PartonicTops(PartonicTops::HADRONIC, false, false), "PartonTops");
            
      //identify bbbar higgs events
      addProjection(PartonicHiggs(PartonicHiggs::BBBAR), "PartonHiggs");
      
      // get bare muons for jet corrections  
      IdentifiedFinalState muons({{PID::MUON, PID::ANTIMUON}}, Cuts::pT > 7 * GeV && Cuts::abseta < 2.5);
      addProjection(muons, "muons");
      
      // get dressed electrons and muons 
      // get photons to dress leptons
      IdentifiedFinalState photons(fs,PID::PHOTON);
      // dressed electrons 
      IdentifiedFinalState emid(fs, {{PID::ELECTRON, PID::POSITRON}});
      PromptFinalState ems(emid, true);
      DressedLeptons dressedelectrons(photons, ems, 0.1, Cuts::pT > 7 * GeV && Cuts::abseta < 2.47, true, true);
      addProjection(dressedelectrons, "dressedelectrons");
      // dressed muons 
      IdentifiedFinalState muid(fs, {{PID::MUON, PID::ANTIMUON}});
      PromptFinalState mus(muid, true);
      DressedLeptons dressedmuons(photons, mus, 0.1, Cuts::pT > 7 * GeV && Cuts::abseta < 2.5, true, true);
      addProjection(dressedmuons, "dressedmuons");

      // neutrinos and DM particles for MET
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      addProjection(neutrinos, "neutrinos");
      
      IdentifiedFinalState darkmatter({{1000022, -1000022}}, Cuts::pT > 0.5 * GeV && Cuts::abseta < 5.0);
      addProjection(darkmatter, "darkmatter");
      
      //Rivet has a list of visible particles. Any particle not in that list is automatically invisible. The missing momentum is calculated as the opposite of the sum of all visible momenta.
      addProjection(MissingMomentum(fs), "MissingET");

      /// Book histogram
      if(_process == "higgs"){
	_h_higgspt        = bookHisto1D("higgspt",            26,  200.,  1500.);
	_h_higgspthbbtag  = bookHisto1D("higgspthbbtag",      26,  200.,  1500.);
      }else if(_process == "ttbar"){
	_h_toppt          = bookHisto1D("toppt",              26,  200.,  1500.);
      }
      _h_fatjetpt         = bookHisto1D("fatjetpt",           26,  200.,  1500.);
      _h_fatjeteta        = bookHisto1D("fatjeteta",          25,  -2.5,    2.5);
      _h_fatjetmass       = bookHisto1D("fatjetmass",         14,    0.,   280.);
      _h_fatjetmasscorr   = bookHisto1D("fatjetmasscorr",     14,    0.,   280.);
      _h_fatjetnsub21     = bookHisto1D("fatjetnsub21",       16,   0.0,    1.6);
      _h_fatjetc2         = bookHisto1D("fatjetc2",           14,   0.0,    0.7);
      _h_fatjetd2         = bookHisto1D("fatjetd2",           25,   0.0,    5.0);
      _h_fatjetntjets     = bookHisto1D("fatjetntjets",        7,  -0.5,    6.5);
      _h_fatjetdphimet    = bookHisto1D("fatjetdphimet",      20,   0.0,    1.0);
      _h_cutflow          = bookHisto1D("cutflow",            10,  -0.5,    9.5);
      _h_weight           = bookHisto1D("weight",            200, -100.,   100.);
      _h_truemet          = bookHisto1D("truemet",            14, 500.0,  1200.);
      _h_met              = bookHisto1D("met",                14, 500.0,  1200.);
      _h_smearedmet       = bookHisto1D("smearedmet",         14, 500.0,  1200.);
      _h_metdiff          = bookHisto1D("metdiff",            20, -100.0, 100.0);
      _h_smearedmetdiff   = bookHisto1D("smearedmetdiff",     20, -100.0, 100.0);
      _h_sphericity       = bookHisto1D("sphericity",         33,  -0.1,    1.0); 
      _h_aplanarity       = bookHisto1D("aplanarity",         36,  -0.1,    0.5); 
      
      _h_higgspt_hbbtagefficiency  = bookScatter2D("higgspt_hbbtagefficiency",  26,  200.,  1500.);

      //_h_deltarfjbh       = bookHisto1D("deltarfjbh",         25,   0.0,    2.5);
      //_h_random             = bookHisto1D("random_numbers",     50,   0.0,    2.0);
      
      // _all=0;
      // _singlebtag=0;
      // _doublebtag=0;
      // _foundhbbbar=0;
      // _allhbbbar=0;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      static random_device rd;
      static mt19937 gen(rd());

      _h_cutflow->fill(0., event.weight());

      // veto events with isolated prompt dressed lepton above certain pT, below some eta
      const vector<DressedLepton> dressedelectrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      const vector<DressedLepton> dressedmuons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      
      // get partonic tops 
      const Particles partontops = apply<ParticleFinder>(event, "PartonTops").particlesByPt();
      
      // get partonic higgs 
      const Particles partonhiggs = apply<ParticleFinder>(event, "PartonHiggs").particlesByPt();
      
      // get muons for jet dressing 
      const Particles& muons = apply<IdentifiedFinalState>(event, "muons").particlesByPt();
      
      // neutrinos and dark matter particles for met 
      const Particles& neutrinos = apply<PromptFinalState>(event, "neutrinos").particlesByPt();
      const Particles& darkmatter = apply<IdentifiedFinalState>(event, "darkmatter").particlesByPt();
      FourMomentum truemet;
      for (const Particle& nu : neutrinos)  truemet += nu.momentum();
      for (const Particle& dm : darkmatter) truemet += dm.momentum();

      const MissingMomentum& mm = applyProjection<MissingMomentum>(event, "MissingET");

      //MET smearing based on https://rivet.hepforge.org/code/dev/SmearingFunctions_8hh_source.html (MET_SMEAR_ATLAS_RUN1) and https://arxiv.org/pdf/1108.5602v2.pdf, Figs 14 and 15
      Vector3 met = mm.vectorEt();
      double set = mm.scalarEt();
      Vector3 smearedmet = met;
      if (met.mod()/GeV < 25*GeV) smearedmet *= 1.05;
      else if (met.mod()/GeV < 40*GeV) smearedmet *= (1.05 - (0.04/15)*(met.mod()/GeV - 25));
      else smearedmet *= 1.01;
      const double resolution = 0.45 * sqrt(set/GeV) * GeV;           //Warum const???
      normal_distribution<> d1(smearedmet.mod(), resolution);
      const double metsmear = max(d1(gen), 0.);
      smearedmet = metsmear * smearedmet.unit();

      if (smearedmet.mod() < 500*GeV) vetoEvent;
      _h_cutflow->fill(1., event.weight());

      // slimmed jets and trimmed fat jets
      const Jets& trackjets  = apply<FastJets>(event, "trackjets").jetsByPt(Cuts::pT > 10*GeV && Cuts::abseta < 2.5);
      
      const Jets& fatcalojets = apply<FastJets>(event, "fatcalojets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);
      
      const Jets& fattopjets = apply<FastTopJets>(event, "fattopjets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);

      const Jets& fathiggsjets = apply<FastHiggsJets>(event, "fathiggsjets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);
      
      // check that ghost tagging doesn't influence number of fat jets in each category 
      if (fatcalojets.size()!= fattopjets.size() || fatcalojets.size()!=fathiggsjets.size()){ 
	std::cout<<"Warning::Ghost tagging affected number of reconstructed jets -> vetoEvent!"<<std::endl;
	vetoEvent;
      }
      // veto if there is no fat jet in the event 
      if (fatcalojets.size()==0) vetoEvent;
      
      //loop over all trackjets and associate to closest fat jet if distance smaller than 1.1  
      //std::cout<<" fat jets in the event " << fatcalojets.size() << " track jets " << trackjets.size() << std::endl;
      vector<unsigned> trackjetfjidx; trackjetfjidx.clear(); 
      for(unsigned itjet = 0 ; itjet < trackjets.size(); itjet++) {
	double drmin(999.);
	unsigned ifjetclose(0);
	for(unsigned ifjet = 0 ; ifjet < fatcalojets.size(); ifjet++) {
	  double tmp=deltaR(trackjets[itjet], fatcalojets[ifjet]);
	  if(tmp < drmin){
	    drmin=tmp;
	    ifjetclose=ifjet;
	  }
	}

	// remove jets with less than 2 tracks 
	if(drmin<1.1 && trackjets[itjet].constituents().size()>1 ) { 
	  trackjetfjidx.push_back(ifjetclose);
	}
	else trackjetfjidx.push_back(999);
      }

      // calculate jet aplanarity/shpericity based on track jets 
      Sphericity sph; sph.calc(trackjets);
      const double sphericity = sph.sphericity(); 
      const double aplanarity = sph.aplanarity();
      
      
      // select good fat jets 
      double beta = 1.0;
      // tau21wta in ATL-PHYS-PUB-2015-035
      NsubjettinessRatio nSub21(2,1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      // energy correlation 
      EnergyCorrelatorC2 c2(beta);
      EnergyCorrelatorD2 d2(beta);
      vector<PseudoJet> trimmedfatjets; trimmedfatjets.clear();
      vector<double> tfjcorrmass; tfjcorrmass.clear();
      vector<int> tfjntracks; tfjntracks.clear();
      for (unsigned ifjet = 0 ; ifjet < fatcalojets.size(); ifjet++) {
	// trimming parameters 
	double rsub = 0.2;
	double fcut = 0.05;
	// from http://fastjet.hepforge.org/svn/contrib/contribs/JetsWithoutJets/tags/1.0.0/example_basic_usage.cc
	fastjet::Filter treetrimmer(rsub, fastjet::SelectorPtFractionMin(fcut));
	PseudoJet pjet = treetrimmer(fatcalojets[ifjet]);
	if ( pjet.perp()<250 )     continue; 
	if ( abs(pjet.eta())>2.0 ) continue;
	//if ( pjet.m()<50 )       continue;
	// loop over track jets to determin fat jet flavor and muon corrections 
	vector <string> trackflavor; trackflavor.clear();
	Particles muoncorr; muoncorr.clear();
	if(trackjets.size()!=trackjetfjidx.size()) 
	  std::cout<<"ERROR::Trackjets and indizes need to have same size!"<<std::endl; 
	for(unsigned itjet = 0 ; itjet < trackjets.size(); itjet++){
	  // check sorting in pT
	  if(itjet<trackjets.size()-1) 
	    if(trackjets[itjet].pT()<trackjets[itjet+1].pT())
	      std::cout<<"ERROR::Tracks not ordered in pT!"<< std::endl;
	  //get jet flavor and muon corrections for fat jet associated track jets 
	  if(trackjetfjidx[itjet]==ifjet) {
	    trackflavor.push_back(getjetflavor(trackjets[itjet]));
	    //get closest muon to correct fat jet mass 
	    double drmin(999.);
	    Particle muonmin;
	    for (Particle muon : muons) {
	      if ( deltaR(muon.momentum(),trackjets[itjet].momentum()) < drmin) {
		drmin = deltaR(muon.momentum(),trackjets[itjet].momentum());
		muonmin = muon;
	      }
	    }
	    //avoid duplicate muon entries 
	    bool usemu=true;
	    for ( Particle muon : muoncorr) if(muon.genParticle()->barcode()==muonmin.genParticle()->barcode()) usemu=false;
	    if (usemu && drmin <0.2) muoncorr.push_back(muonmin);
	  }// end of if track jet associated to fat jet 
	}// end of track jet loop
	// veto non b-tagged jets 
	int btags(0.);
	for (unsigned i=0; i<trackflavor.size(); i++) if ( trackflavor[i] == "bjet" ) btags++;
	if ( btags<2 ) continue;
	// veto trimmed jets with nans in substructure 
	if(std::isnan(nSub21(pjet))||std::isnan(c2(pjet))||std::isnan(d2(pjet))){
	  std::cout<<"WARNING::NAN in jet substructure!"<<std::endl;
	  continue;
	}
	//check b-tagging efficiency/rejection for fat jets 
	// if(_process == "higgs" &&  higgsbbbartagged(fathiggsjets[ifjet], "ALL", Cuts::pT>250*GeV && Cuts::abseta < 2.0, Cuts::pT>5*GeV && Cuts::abseta < 2.5 )){
	//   _all++;
	//   for (unsigned i=0; i<trackflavor.size(); i++) 
	//     std::cout<< "i:" <<  i << " flavor " << trackflavor[i] << std::endl;
	//   if (trackflavor.size()>0)  if(trackflavor[0] == "bjet" ) { _singlebtag++; singlebtag=true; } 
	//   if (trackflavor.size()>1)  if(trackflavor[0] == "bjet" && trackflavor[1] == "bjet") { _doublebtag++; doublebtag=true; }
	// }
	// if(_process == "ttbar" &&  toptagged(fattopjets[ifjet])) { 
	//   _all++;
	//   for (unsigned i=0; i<trackflavor.size(); i++) 
	//     std::cout<< "i:" <<  i << " flavor " << trackflavor[i] << std::endl;
	//   if (trackflavor.size()>0)  if(trackflavor[0] == "bjet" ) { _singlebtag++; singlebtag=true; } 
	//   if (trackflavor.size()>1)  if(trackflavor[0] == "bjet" && trackflavor[1] == "bjet") { _doublebtag++; doublebtag=true; }
	// }
	
	// only plot filled in jet loop, because this should happen only once per higgs, so this is ok to have here 
	if(higgsbbbartagged(fathiggsjets[ifjet], "ALL", Cuts::pT>250*GeV && Cuts::abseta < 2.0, Cuts::pT>5*GeV && Cuts::abseta < 2.5 ))
	  _h_higgspthbbtag->fill(partonhiggs[0].perp()/GeV,  event.weight());  
	
	//Jet smearing assuming 10% mass resolution, see https://indico.in2p3.fr/event/9417/session/2/contribution/5/material/slides/0.pdf and ATL-PHYS-PUB-2015-035
	static const double resolution = 0.1;
	normal_distribution<> d2(1., resolution);
	const double fsmear = max(d2(gen), 0.);

	// push back selected jets and some jet properties 
	trimmedfatjets.push_back(pjet);
	FourMomentum fatjetcorr(pjet.e(), pjet.px(), pjet.py(), pjet.pz());
	for ( Particle muon : muoncorr) fatjetcorr += muon.momentum();
	tfjcorrmass.push_back(fatjetcorr.mass()*fsmear);
	tfjntracks.push_back(trackflavor.size());
      }
      // veto events with less or more than one good fat jet 
      if(trimmedfatjets.size()!=1) vetoEvent;
      _h_cutflow->fill(2., event.weight());
      
      if ((dressedelectrons+dressedmuons).size()>0)  vetoEvent;
      _h_cutflow->fill(3., event.weight());
      
      // fill histograms 
      // partonic top/higgs pT 
      if( _process == "higgs") {
	if ( partonhiggs.size()>0 && partonhiggs[0].children().size()>1 ){
	  if( (partonhiggs[0].perp()>250*GeV && abs(partonhiggs[0].eta()) < 2.0) && 
	      (partonhiggs[0].children()[0].perp()>5*GeV && abs(partonhiggs[0].children()[0].eta()) < 2.5 ) && 
	      (partonhiggs[0].children()[1].perp()>5*GeV && abs(partonhiggs[0].children()[1].eta()) < 2.5 ) ) {
	    _h_higgspt->fill(partonhiggs[0].perp()/GeV,  event.weight());  
	  }	  
	}
      } else if ( _process == "ttbar" ) {
	if (partontops.size() >= 1) if ( partontops[0].perp()>250 * GeV ) _h_toppt->fill(partontops[0].perp()/GeV,  event.weight());
	if (partontops.size() >= 2) if ( partontops[1].perp()>250 * GeV ) _h_toppt->fill(partontops[1].perp()/GeV,  event.weight());
      }
      
      // fat jet related quantities
      _h_fatjetpt->fill(trimmedfatjets[0].perp()/GeV,                                     event.weight());
      _h_fatjeteta->fill(trimmedfatjets[0].eta(),                                         event.weight());
      _h_fatjetmass->fill(trimmedfatjets[0].m()/GeV,                                      event.weight());
      _h_fatjetmasscorr->fill(tfjcorrmass[0]/GeV,                                         event.weight());
      _h_fatjetnsub21->fill(nSub21(trimmedfatjets[0]),                                    event.weight());
      _h_fatjetc2->fill(c2(trimmedfatjets[0]),                                            event.weight());
      _h_fatjetd2->fill(d2(trimmedfatjets[0]),                                            event.weight());
      _h_fatjetntjets->fill(tfjntracks[0],                                                event.weight());
      _h_fatjetdphimet->fill(deltaPhi(trimmedfatjets[0].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());
      
      // event properties 
      _h_weight->fill(event.weight(),                1.);	  
      _h_truemet->fill(truemet.pT()/GeV, event.weight());	  
      _h_met->fill(met.mod()/GeV,         event.weight());	  
      _h_smearedmet->fill(smearedmet.mod()/GeV,event.weight());	  
      _h_metdiff->fill(truemet.pT()/GeV-met.mod()/GeV,event.weight());	  
      _h_smearedmetdiff->fill(met.mod()/GeV-smearedmet.mod()/GeV,event.weight());	  
      _h_aplanarity->fill(aplanarity,    event.weight());	  
      _h_sphericity->fill(sphericity,    event.weight());	  
      
    }//analyze
      
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      double lumi = 3200.; 
      double norm=lumi*crossSection()/picobarn / sumOfWeights();
      
      std::cout<<" lumi*crossSection()/picobarn " << lumi*crossSection()/picobarn << " sum of weiht " << sumOfWeights() << " nrm  "<< norm << std::endl; 
      
      _h_fatjetpt->scaleW(norm*_h_fatjetpt->bin(0).width());
      _h_fatjeteta->scaleW(norm*_h_fatjeteta->bin(0).width());
      _h_fatjetmass->scaleW(norm*_h_fatjetmass->bin(0).width());
      _h_fatjetmasscorr->scaleW(norm*_h_fatjetmasscorr->bin(0).width());
      _h_fatjetnsub21->scaleW(norm*_h_fatjetnsub21->bin(0).width());
      _h_fatjetc2->scaleW(norm*_h_fatjetc2->bin(0).width());
      _h_fatjetd2->scaleW(norm*_h_fatjetd2->bin(0).width());
      _h_fatjetntjets->scaleW(norm*_h_fatjetntjets->bin(0).width());
      _h_fatjetdphimet->scaleW(norm*_h_fatjetdphimet->bin(0).width());
      _h_cutflow->scaleW(norm*_h_cutflow->bin(0).width());
      _h_truemet->scaleW(norm*_h_truemet->bin(0).width());
      _h_met->scaleW(norm*_h_met->bin(0).width());
      _h_smearedmet->scaleW(norm*_h_smearedmet->bin(0).width());
      _h_metdiff->scaleW(norm*_h_metdiff->bin(0).width());
      _h_aplanarity->scaleW(norm*_h_aplanarity->bin(0).width());
      _h_sphericity->scaleW(norm*_h_sphericity->bin(0).width());
      
      if(_process=="higgs"){
	_h_higgspt->scaleW(norm/_h_higgspt->bin(0).width());
	_h_higgspthbbtag->scaleW(norm/_h_higgspthbbtag->bin(0).width());
      	//Rivet::Analysis::efficiency(*_h_higgspthbbtag,  *_h_higgspt,  _h_higgspt_hbbtagefficiency);
      }else if(_process=="ttbar"){
	_h_toppt->scaleW(norm/_h_toppt->bin(0).width());
      }
      
      
      // std::cout<<" Efficiency for double btag="<< _doublebtag/_all << ", rejection=" << _all/_doublebtag << std::endl;
      // std::cout<<" Efficiency for single btag="<< _singlebtag/_all << ", rejection=" << _all/_singlebtag << std::endl;
      
    }
    
  private:
    
    std::string getjetflavor(const Jet& j){
      
      double rndm = rand()/static_cast<double>(RAND_MAX);
      if(!j.bTags().empty()){
	if(rndm<=0.70) return "bjet";
      }else if(!j.cTags().empty()){ 
	if(rndm<=0.18) return "bjet";
	else return "cjet";
      }else{ 
	if(rndm<=0.006) return "bjet";
      }
      return "ljet";
      
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    bool toptagged(const Jet& j, const Cut& ct=Cuts::open()){
      Particles rtn;
      for (const Particle& tp : j.tags()) 
	if (abs(tp.pdgId())==PID::TQUARK && ct->accept(tp)) rtn.push_back(tp);
      return rtn.size()==1; 
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    bool higgsbbbartagged(const Jet& j, string tagged, const Cut& ch=Cuts::open(), const Cut& cb=Cuts::open()){
      Particles rtn; rtn.clear();
      //Particles bpartons; bpartons.clear();
      for (const Particle& tp : j.tags()){ 
	if(Rivet::PID::isHiggs(tp.pdgId()) && ch->accept(tp)) 
	  rtn.push_back(tp);
	if(tagged == "ALL" && Rivet::PID::isHadron(tp.pdgId()) && Rivet::PID::hasBottom(tp.pdgId()) && cb->accept(tp)){
	  rtn.push_back(tp); 
	  //bpartons.push_back(tp);
	  //std::cout<<" tag perp() " << tp.momentum().perp() << " eta() " << tp.momentum().eta() << " phi() " << tp.momentum().phi() << std::endl; 
	}
	//bool found = false;
	// for(Particle tj : _trackjettaggers ) {
	//     if(abs((tj.momentum().perp()-tp.momentum().perp())/tj.momentum().perp())<0.02 && abs((tj.momentum().eta()-tp.momentum().eta())/tj.momentum().eta())<0.02 )
	//       found = true;
	//   }
	//   _allhbbbar++;
	//   if(! found){ 
	//     std::cout<<"hbbbar tagging " << tp.genParticle()->barcode() << " id " << tp.pdgId() << " status " << tp.genParticle()->status() <<" momentum " << tp.momentum().perp() << "  " <<  tp.momentum().eta() <<std::endl;  
	//     for(Particle tj : _trackjettaggers ) 
	//       std::cout<<"track jet tagger " << tj.genParticle()->barcode() << " id " << tj.pdgId() << " status " << tj.genParticle()->status()  <<" momentum " << tj.momentum().perp() << "  " <<  tj.momentum().eta() <<std::endl;  
	//   }else{
	//     _foundhbbbar++;
	//   }
      }
      // if(rtn.size()==3) {
      // 	for( Particle tp: rtn ){
      // 	  if (Rivet::PID::isHadron(tp.pdgId()) && Rivet::PID::hasBottom(tp.pdgId()))
      // 	    _h_deltarfjbh->fill(deltaR(tp.momentum(),j.momentum()),  1.);  
      // 	  std::cout<<"bhadron from higgs " << tp.genParticle()->barcode() << " id " << tp.pdgId() << " status " << tp.genParticle()->status() 
      // 		   <<" momentum " << tp.momentum().perp() << "  " <<  tp.momentum().eta() <<std::endl;
      // 	}
      // }
      //if(bpartons.size() ==2 )
      //std::cout<<" DELTAR " << deltaR(bpartons[0].momentum(), bpartons[1].momentum()) << std::endl;
      if(rtn.size()>3) 
	std::cout<<"Warning::More than 3 particles for higgs tagging due to b hadron ambiguities!"<<std::endl;
      if (tagged == "ALL") return rtn.size()==3; 
      return rtn.size()==1;
    }
    
  private:
    
    std::string  _process;
    // float        _all;
    // float        _singlebtag;
    // float        _doublebtag;
    
    Histo1DPtr   _h_higgspt;
    Histo1DPtr   _h_higgspthbbtag;
    Histo1DPtr   _h_toppt;
    Histo1DPtr   _h_fatjetpt;
    Histo1DPtr   _h_fatjeteta;
    Histo1DPtr   _h_fatjetmass;
    Histo1DPtr   _h_fatjetmasscorr;
    Histo1DPtr   _h_fatjetnsub21;
    Histo1DPtr   _h_fatjetc2;
    Histo1DPtr   _h_fatjetd2;
    Histo1DPtr   _h_fatjetntjets;
    Histo1DPtr   _h_fatjetdphimet;
    Histo1DPtr   _h_cutflow;      
    Histo1DPtr   _h_weight;
    Histo1DPtr   _h_truemet;
    Histo1DPtr   _h_met;
    Histo1DPtr   _h_smearedmet;
    Histo1DPtr   _h_metdiff;
    Histo1DPtr   _h_smearedmetdiff;
    Histo1DPtr   _h_sphericity;
    Histo1DPtr   _h_aplanarity;
    
    // efficiency of hbbbar tag as a function of higgs pT
    Scatter2DPtr _h_higgspt_hbbtagefficiency; 
    
    // Particles    _trackjettaggers; 
    // float        _foundhbbbar;
    // float        _allhbbbar;
    //Histo1DPtr   _h_deltarfjbh;
    //Histo1DPtr     _h_random;
  };
  
  
  DECLARE_RIVET_PLUGIN(FatHiggsTagging);
  
}
// agrohsje taken from ATLAS_2015_I1397637

