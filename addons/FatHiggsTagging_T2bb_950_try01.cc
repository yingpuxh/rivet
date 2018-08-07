// -*- C++ -*-
// // 35.9 fb-1, SUS 16 038
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
#include <math.h>

using namespace fastjet::contrib;


namespace Rivet {
  

  class FatHiggsTagging : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(FatHiggsTagging);
    double rand01() {
      static random_device rd;
      static mt19937 gen(rd());
      return generate_canonical<double, 10>(gen);
    }

    bool isTagging(double eff){
      double rand=rand01();
      cout<<"rand01        "<<rand<<endl;
      return rand<eff;
    } 
    double JET_BTAG_ATLAS_RUN2_MV2C20(const Jet& j) {
    if (j.abseta() > 2.5) return 0;
    if (j.bTagged(Cuts::pT > 5*GeV)) return 0.77;
    if (j.cTagged(Cuts::pT > 5*GeV)) return 1/4.5;
    return 1/140.;
    }
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
            
      // track jets, 
      //agrohsje include muons ? 
      FastJets trackjets(cfs, FastJets::ANTIKT, 0.2, JetAlg::DECAY_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(trackjets, "trackjets");

      // Large-R jets
      //FastJets fatcalojets(fs, FastJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      //addProjection(fatcalojets, "fatcalojets");
     

      //Regular jets
      FastJets calojets(fs, FastJets::ANTIKT, 0.4, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(calojets, "calojets");

 
      // with top tagging 
      FastTopJets regulartopjets(fs, FastTopJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(regulartopjets, "regulartopjets");

      // with higgs and b hadron decay tagging 
      FastHiggsJets regularhiggsjets(fs, FastHiggsJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      addProjection(regularhiggsjets, "regularhiggsjets");
      
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
      DressedLeptons dressedelectrons(photons, ems, 0.1, Cuts::pT > 10 * GeV && Cuts::abseta < 2.5, true, true);
      addProjection(dressedelectrons, "dressedelectrons");
      // dressed muons 
      IdentifiedFinalState muid(fs, {{PID::MUON, PID::ANTIMUON}});
      PromptFinalState mus(muid, true);
      DressedLeptons dressedmuons(photons, mus, 0.1, Cuts::pT > 10 * GeV && Cuts::abseta < 2.5, true, true);
      addProjection(dressedmuons, "dressedmuons");
      //add photons Projection
      //PromptFinalState phos(photons, true);
      //DressedLeptons dressedphotons(photons, phos, 0.1, Cuts::pT > 25 * GeV && Cuts::abseta < 2.5, true, true);
      //addProjection(dressedphotons, "dressedphotons");
      IdentifiedFinalState dressedphotons({{PID::PHOTON}}, Cuts::pT > 25 * GeV && Cuts::abseta < 2.5);
      addProjection(dressedphotons, "dressedphotons");
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
      _h_cutflow          = bookHisto1D("cutflow",           15,  -0.5,    14.5);
      _h_weight           = bookHisto1D("weight",            200, -100.,   100.);
      _h_truemet          = bookHisto1D("truemet",            14, 500.0,  1200.);
      _h_met              = bookHisto1D("met",                14, 500.0,  1200.);
      _h_smearedmet3      = bookHisto1D("smearedmet3",        36, 0.0,  1500.);
      _h_smearedmet4      = bookHisto1D("smearedmet4",        36, 0.0,  1500.);
      _h_metdiff          = bookHisto1D("metdiff",            20, -100.0, 100.0);
      _h_smearedmetdiff   = bookHisto1D("smearedmetdiff",     20, -100.0, 100.0);
      _h_sphericity       = bookHisto1D("sphericity",         33,  -0.1,    1.0); 
      _h_aplanarity       = bookHisto1D("aplanarity",         36,  -0.1,    0.5); 
      _h_AlphaT		  = bookHisto1D("AlphaT",	      36,  0.52,   0.65);
      _h_htmissbefore     = bookHisto1D("htmissbefore",       38,   0.0,   3800.);
      _h_htmissbefore->sumW2();
      _h_htmiss		  = bookHisto1D("htmiss",             20,   700.0,   1200.);
      _h_htmiss->sumW2();
      _h_AzimuthalAngular = bookHisto1D("AzimuthalAngular",   36,   0.5,   M_PI);
      _h_HtmissMETDivision= bookHisto1D("HtmissMETDivision",  15,   0.0,   1.25);
      //_h_truemet0       = bookHisto1D("truemet0",           14, 500.0,  1200.);
      //_h_met0           = bookHisto1D("met0",               14, 500.0,  1200.);
      _h_smearedmet0      = bookHisto1D("smearedmet0",        36, 0.0,  1500.);
      //_h_truemet1       = bookHisto1D("truemet1",           14, 500.0,  1200.);
      //_h_met1           = bookHisto1D("met1",               14, 500.0,  1200.);
      _h_smearedmet1      = bookHisto1D("smearedmet1",        36, 0.0,  1500.);
      //_h_truemet2       = bookHisto1D("truemet2",           14, 500.0,  1200.);
      //_h_met2           = bookHisto1D("met2",               14, 500.0,  1200.);
      _h_smearedmet2      = bookHisto1D("smearedmet2",        36, 0.0,  1500.);
      _h_htbefore         = bookHisto1D("htbefore",           39, 0.0,  3900.);
      _h_htbefore->sumW2();
      _h_ht               = bookHisto1D("ht",    	      20, 700.0,  1200.);
      _h_ht->sumW2();
      _h_httrackjets     = bookHisto1D("httrackjets",         36, 0.0,  2000.);
      _h_dRregularjets0  = bookHisto1D("dRregularjets0",      36,   0.,     1.0);
      _h_dRregularjets1  = bookHisto1D("dRregularjets1",      36,   0.,     1.0);
      _h_dRregularjets2  = bookHisto1D("dRregularjets2",      36,   0.,     1.0);
      _h_dRregularjets3  = bookHisto1D("dRregularjets3",      36,   0.,     1.0);
      _h_dRtrackjets0    = bookHisto1D("dRtrackjets0",        36,   0.,     1.0);
      _h_dRtrackjets1    = bookHisto1D("dRtrackjets1",        36,   0.,     1.0);
      _h_dRtrackjets2    = bookHisto1D("dRtrackjets2",        36,   0.,     1.0);
      _h_dRtrackjets3    = bookHisto1D("dRtrackjets3",        36,   0.,     1.0); 
      _h_regularjetseta0  = bookHisto1D("regularjetseta0",    25,  -3.0,    2.5);
      _h_regularjetseta1  = bookHisto1D("regularjetseta1",    25,  -2.5,    2.5);
      _h_regularjetseta2  = bookHisto1D("regularjetseta2",    25,  -2.5,    2.5);
      _h_regularjetseta3  = bookHisto1D("regularjetseta3",    25,  -2.5,    2.5);
      _h_regularjetspt0   = bookHisto1D("regularjetspt0",     26,   0.,   1500.);
      _h_regularjetspt1   = bookHisto1D("regularjetspt1",     26,   0.,   600.);
      _h_regularjetspt2   = bookHisto1D("regularjetspt2",     26,   0.,   500.);
      _h_regularjetspt3   = bookHisto1D("regularjetspt3",     26,   0.,   400.);
      _h_regularjetsconstituent0  = bookHisto1D("regularjetsconstituent0",    25,  0.,   50.);
      _h_regularjetsconstituent1  = bookHisto1D("regularjetsconstituent1",    25,  0.,   50.);
      _h_regularjetsconstituent2  = bookHisto1D("regularjetsconstituent2",    25,  0.,   50.);
      _h_regularjetsconstituent3  = bookHisto1D("regularjetsconstituent3",    25,  0.,   50.);
      _h_trackjetsconstituent0    = bookHisto1D("trackjetsconstituent0",    25,  0.,    50.);
      _h_trackjetsconstituent1    = bookHisto1D("trackjetsconstituent1",    25,  0.,    50.);
      _h_trackjetsconstituent2    = bookHisto1D("trackjetsconstituent2",    25,  0.,    50.);
      _h_trackjetsconstituent3    = bookHisto1D("trackjetsconstituent3",    25,  0.,    50.);
      _h_trackjetsconstituentsum    = bookHisto1D("trackjetsconstituentsum",    25,  0.,    50.);
      _h_regularjetsconstituentsum  = bookHisto1D("regularjetsconstituentsum",    25,  0.,   50.);
      _h_trackjetseta0  = bookHisto1D("trackjetseta0",    25,  -3.0,    2.5);
      _h_trackjetseta1  = bookHisto1D("trackjetseta1",    25,  -2.5,    2.5);
      _h_trackjetseta2  = bookHisto1D("trackjetseta2",    25,  -2.5,    2.5);
      _h_trackjetseta3  = bookHisto1D("trackjetseta3",    25,  -2.5,    2.5);
      _h_trackjetspt0   = bookHisto1D("trackjetspt0",     26,   0.,   1500.);
      _h_trackjetspt1   = bookHisto1D("trackjetspt1",     26,   0.,   1500.);
      _h_trackjetspt2   = bookHisto1D("trackjetspt2",     26,   0.,   1500.);
      _h_trackjetspt3   = bookHisto1D("trackjetspt3",     26,   0.,   1500.);
      _h_trackjetssize    = bookHisto1D("trackjetssize",       20,   -5.,   15.);
      _h_regularjetssize  = bookHisto1D("regularjetssize",     1,   2.,   3.);
      _h_bjetssizetry  = bookHisto1D("bjetssizetry",     5,   0.,   5.);
      _h_bjetssize        = bookHisto1D("bjetssize",           1,   1.,   2.);
      _h_bjetssizebefore        = bookHisto1D("bjetssizebefore",           6,   0.,   6.);
      _h_regularjetssizebefore  = bookHisto1D("regularjetssizebefore",     14,   0.,   14.);
      //_h_regularjet0dphimet      = bookHisto1D("regularjet0dphimet",          20,   0.0,    1.0);
      //_h_regularjet1dphimet      = bookHisto1D("regularjet1dphimet",          20,   0.0,    1.0);
      //_h_regularjet2dphimet      = bookHisto1D("regularjet2dphimet",          20,   0.0,    1.0);
      //_h_regularjet3dphimet      = bookHisto1D("regularjet3dphimet",          20,   0.0,    1.0);
      _h_regularjetsht             = bookHisto1D("regularjetsht",               26,    0.,  1500.);
      _h_higgspt_hbbtagefficiency  = bookScatter2D("higgspt_hbbtagefficiency",  26,  200.,  1500.);
      _h_htmht                     = bookHisto2D("htmht",   20,   700.0,   1200.,  20,  700.0, 1200.);
      _h_htmhtnb0                     = bookHisto2D("htmhtnb0",   20, 700.0,  2700.,  20,  700.0, 1200.);
      _h_htmhtnb1                     = bookHisto2D("htmhtnb1",   20, 700.0,  2700.,  20,  700.0, 1200.);
      _h_htmhtnb2                     = bookHisto2D("htmhtnb2",   20, 700.0,  2700.,  20,  700.0, 1200.);
      _h_htmhtnb3                     = bookHisto2D("htmhtnb3",   20, 700.0,  2700.,  20,  700.0, 1200.);
      _h_nbnjet                    = bookHisto2D("nbnjet",  1,   1.,   2.,  1,   2.,   3.);


      //_h_deltarfjbh              = bookHisto1D("deltarfjbh",                  25,   0.0,    2.5);
      //_h_random                  = bookHisto1D("random_numbers",              50,   0.0,    2.0);
      
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
      // _h_smearedmet0->fill(smearedmet.mod()/GeV,event.weight());

      // veto events with isolated prompt dressed lepton above certain pT
      const vector<DressedLepton> dressedelectrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      const vector<DressedLepton> dressedmuons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      if (dressedelectrons.size()>0) vetoEvent;
      if (dressedmuons.size()>0) vetoEvent;
      _h_cutflow->fill(1., event.weight());
    //  const vector<DressedLepton> dressedphotons = apply<DressedLeptons>(event, "dressedphotons").dressedLeptons();     
      const Jets& trackjets  = apply<FastJets>(event, "trackjets").jetsByPt(Cuts::pT > 10*GeV && Cuts::abseta < 5.0);

      const Jets& calojets  = apply<FastJets>(event, "calojets").jetsByPt(Cuts::pT > 40*GeV && Cuts::abseta < 5.0);//&& Cuts::abseta < 3.0);
      //matching tracjets and keep the number of calojets which will be useful in determing jet flavour. The matching track jets are ordered in pT.   
      vector<unsigned> trackjetfjidx1; trackjetfjidx1.clear();
      Jets jss;
      //matching calojets and keep the matching calojets in pT. 
      Jets js;
      for(unsigned itjet = 0 ; itjet < calojets.size(); itjet++) {
         double drmin(999.);
         unsigned ifjetclose(0);
         unsigned itjetclose(0);
         for(unsigned ifjet = 0 ; ifjet < trackjets.size(); ifjet++) {
          double tmp=deltaR(trackjets[ifjet], calojets[itjet]);
          if(tmp < drmin){
            drmin=tmp;
            itjetclose=itjet;
            ifjetclose=ifjet;
          }
         }
         if (drmin!=999){
           cout<<"                              Before track jet matching (before drmin<0.4)"<<endl;
           cout<<"calojets.pT        "<<calojets[itjet].pT()<<endl;
           cout<<"trackjets.pT        "<<trackjets[ifjetclose].pT()<<endl;
           cout<<"calojets.phi      "<<calojets[itjet].phi()<<endl;
           cout<<"trackjets.phi      "<<trackjets[ifjetclose].phi()<<endl;
           cout<<"calojets.eta        "<<calojets[itjet].eta()<<endl;
           cout<<"trackjets.eta        "<<trackjets[ifjetclose].eta()<<endl;
         }
          if(drmin<0.4){
           jss.push_back(trackjets[ifjetclose]);
           trackjetfjidx1.push_back(itjetclose);
           js.push_back(calojets[itjetclose]);
           cout<<"                       After track calo jets matching"<<endl;
           cout<<"               drmin  "<<drmin<<endl;
           cout<<"calojets.pT        "<<calojets[itjet].pT()<<endl;
           cout<<"trackjets.pT        "<<trackjets[ifjetclose].pT()<<endl;
           cout<<"calojets.phi   "<<calojets[itjet].phi()<<endl;
           cout<<"trackjets.phi      "<<trackjets[ifjetclose].phi()<<endl;
           cout<<"calojets.eta        "<<calojets[itjet].eta()<<endl;
           cout<<"trackjets.eta        "<<trackjets[ifjetclose].eta()<<endl;
           cout<<"                       "<<endl;        //js.size() contains all calojets with at leaset trackjets and drmin<0.4
          }
      }
      cout<<"calojets.size		"<<js.size()<<endl;
      cout<<"trackjets.size		"<<jss.size()<<endl;


      
       //Jets js;
       //foreach (const Jet& j, calojets){
       //js.push_back(j);
       //}
      _h_cutflow->fill(2., event.weight());
       //get photons selections
      const Particles photons = apply<IdentifiedFinalState>(event, "dressedphotons").particlesByPt(); 
      for(unsigned itjet = 0 ; itjet < photons.size(); itjet++){
       double drmin=100000.00;
       for(unsigned ifjet = 0 ; ifjet < js.size(); ifjet++){
         double tmp = deltaR(photons[itjet],js[ifjet]);
         if (tmp<drmin)
         {drmin=tmp;
         }
       }
       if (drmin>0.4)
       {vetoEvent;
       break;
       }
      }
      _h_cutflow->fill(3., event.weight());
      //draw ht, htmiss and njets before any selections
      double  htbefore=0.0;
      FourMomentum sum_jetbefore;
      for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
       if (js[itjet].perp()/GeV>40 && js[itjet].abseta()<3.0){
         htbefore+=js[itjet].pT();
         Jet j = js[itjet];
         sum_jetbefore += js[itjet].momentum();
       }
      }
      _h_htbefore->fill(htbefore,                                                             event.weight());
      _h_htmissbefore->fill(sum_jetbefore.pT(),                                               event.weight());
      int countjets=0;
      for(unsigned itjet = 0 ; itjet < js.size(); itjet++){
      if (js[itjet].perp()/GeV>40 && js[itjet].abseta()<3.0)
       countjets=countjets+1;
      }
      _h_regularjetssizebefore->fill(countjets,                                               event.weight());


      //draw number of bjets before any selections
//      vector <string> trackflavors; trackflavors.clear();
      vector <string> trackflavors; trackflavors.clear();
      int btag=0;
      for (unsigned ifjet = 0 ; ifjet < js.size(); ifjet++) {
       if(jss.size()!=trackjetfjidx1.size())
       std::cout<<"ERROR::Trackjets and indizes need to have same size!"<<std::endl;
       if (js[ifjet].perp()/GeV>40 && js[ifjet].abseta()<3.0){
       int isbtagg(0.0);
       double eff(0.0);
       for(unsigned itjet = 0 ; itjet < jss.size(); itjet++){ 
         if(trackjetfjidx1[itjet]==ifjet) {
         trackflavors.push_back(getjetflavor(jss[itjet]));
         break;
         }   
       }
       for (unsigned i=0; i<trackflavors.size(); i++)
          if ( trackflavors[i] == "bjet" ){
            eff = 0.77;
            isbtagg=isTagging(eff);  
            if (isbtagg!=0)  btag=btag+1;
          }
       trackflavors.clear();
       }
      }
      _h_bjetssizebefore->fill(btag,                                                             event.weight());
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
      //fill the MET beofre any selections
      _h_smearedmet0->fill(smearedmet.mod()/GeV,event.weight());      
 
      //const Jets& fatcalojets = apply<FastJets>(event, "fatcalojets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);
      
      const Jets& regulartopjets = apply<FastTopJets>(event, "regulartopjets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);

      const Jets& regularhiggsjets = apply<FastHiggsJets>(event, "regularhiggsjets").jetsByPt(Cuts::pT > 250*GeV && Cuts::abseta < 2.0);
      
      // check that ghost tagging doesn't influence number of fat jets in each category 
      if (regulartopjets.size()!=regularhiggsjets.size()){ 
	std::cout<<"Warning::Ghost tagging affected number of reconstructed jets -> vetoEvent!"<<std::endl;
	vetoEvent;
      }
      // veto if there is no (fat jet or) regular jet in the event 
      //if (fatcalojets.size()==0) vetoEvent;
      //if (js.size()==0) vetoEvent;
      // forward jets acceptance criteria
//  SUS_15_005: forward jet veto
      //for(unsigned itjet = 0 ; itjet < js.size(); itjet++){
      // if(js[itjet].perp()/GeV>40 && js[itjet].abseta()>2.4)
      // vetoEvent;
     // }
     //_h_cutflow->fill(4., event.weight());
     //end forward jet veto
      Sphericity sph; sph.calc(jss);
      const double sphericity = sph.sphericity(); 
      const double aplanarity = sph.aplanarity();
      //jet acceptance criteria
      for (unsigned itjet=0;itjet<js.size();itjet++){
        if(itjet<js.size()-1)
        if(js[itjet].pT()<js[itjet+1].pT()){
         std::cout<<"ERROR::js not ordered in pT!"<< std::endl;
         //cout<<"js[itjet].pT       "<<js[itjet].pT()<<endl;
         //cout<<"js[itjet+1].pT       "<<js[itjet+1].pT()<<endl;
         //cout<<"js[itjet]       "<<itjet<<endl;
      } }
      if (js.size()<2) vetoEvent;
      _h_cutflow->fill(4., event.weight());
// missing CHF selection
      _h_cutflow->fill(5., event.weight());
      if (js[0].perp()/GeV<100) vetoEvent;
     _h_cutflow->fill(6., event.weight());
     //std::cout<<" regular jets in the event " << calojets.size() << " track jets " << trackjets.size()<< std::endl;
    //calculate ht and htmiss               
    double  ht=0.0;
    FourMomentum sum_jet;
    for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
     ht+=js[itjet].pT();
     Jet j = js[itjet];
     sum_jet += js[itjet].momentum();
    }
    double httrackjets=0.0;
    for (unsigned itjet = 0 ; itjet < jss.size(); itjet++){
     httrackjets+=jss[itjet].pT();
    }
    //cout<<"ht: "<<ht<<endl;
    //cout<<"sum_jet: "<<sum_jet.pT()<<endl;
    //   Jets js;
    // foreach (const Jet& j, calojets) {
    //   js.push_back(j); 

    //}
    // criterias for Ht, ht miss and met
    if (ht < 200) vetoEvent;
    _h_cutflow->fill(7., event.weight());

    if (sum_jet.pT() < 200) vetoEvent;
    _h_cutflow->fill(8., event.weight());
    for(unsigned itjet = 0 ; itjet < js.size(); itjet++){
       if(js[itjet].perp()/GeV>40 && js[itjet].abseta()>2.4)
              vetoEvent;
    }
    _h_cutflow->fill(9., event.weight());
    if (sum_jet.pT()/smearedmet.mod() > 1.25) vetoEvent;
    _h_cutflow->fill(10., event.weight());
    //if (sum_jet.pT()/smearedmet.mod() > 1.25) vetoEvent;
    //_h_cutflow->fill(3., event.weight());
    //_h_smearedmet3->fill(smearedmet.mod()/GeV,event.weight());   
    // select good  calojets 
    double beta = 1.0;
    // tau21wta in ATL-PHYS-PUB-2015-035
    NsubjettinessRatio nSub21(2,1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    // energy correlation 
    EnergyCorrelatorC2 c2(beta);
    EnergyCorrelatorD2 d2(beta);
    //vector<double> tfjcorrmass; tfjcorrmass.clear();
    vector<int> tfjntracks; tfjntracks.clear();
    vector <string> trackflavor; trackflavor.clear(); 
    int btags=0;
    for (unsigned ifjet = 0 ; ifjet < js.size(); ifjet++) {
     // loop over track jets to determin fat jet flavor and muon corrections 
      //vector <string> trackflavor; trackflavor.clear();
      Particles muoncorr; muoncorr.clear();
      if(jss.size()!=trackjetfjidx1.size())
      std::cout<<"ERROR::Trackjets and indizes need to have same size!"<<std::endl;
      int isbtag(0.0);
      double eff(0.0);
      for(unsigned itjet = 0 ; itjet < jss.size(); itjet++){
        if(trackjetfjidx1[itjet]==ifjet) {
         //eff = JET_BTAG_ATLAS_RUN2_MV2C20(jss[itjet]);
         //isbtag=isTagged(eff);
        trackflavor.push_back(getjetflavor(jss[itjet]));
        break;
        }
      }// end of track jet loop
      for (unsigned i=0; i<trackflavor.size(); i++){
          if ( trackflavor[i] == "bjet" ){
           eff = 0.77;
           isbtag=isTagging(eff);
           if (isbtag!=0)  btags=btags+1;
          }
      }
      //Jet smearing assuming 10% mass resolution, see https://indico.in2p3.fr/event/9417/session/2/contribution/5/material/slides/0.pdf and ATL-PHYS-PUB-2015-035
      static const double resolution = 0.1;
      normal_distribution<> d2(1., resolution);
      //const double fsmear = max(d2(gen), 0.);
      //push back some jet properties and selected jets
      FourMomentum regularjetcorr(js[ifjet].perp(), js[ifjet].px(), js[ifjet].py(), js[ifjet].pz());
      //for ( Particle muon : muoncorr) regularjetcorr += muon.momentum();
      //tfjcorrmass.push_back(regularjetcorr.mass()*fsmear);
      tfjntracks.push_back(trackflavor.size());
      trackflavor.clear();
    }
    //Jets js;
    //foreach (const Jet& j, calojets) {
    //js.push_back(j);
    //}
    //Jets jss;
    //foreach (const Jet& j, trackjets) {
    //jss.push_back(j);
    //}
    //Signal Region Criteria 
    double drmax=0.0;
    //criteria for htmiss and calojets -> azimuthal angle
    //Alpha_T only suits for multijets
    //Alpha_T for dijets
    double AlphaT = 0.0;
    if (js.size()==2)
    {
      if (js[0].pT()<js[1].pT())
      {cout<<"calojets not ordered in pT"<<endl;
      double mass= sqrt((js[0].pT()+js[1].pT())*(js[0].pT()+js[1].pT())-(js[0].px()+js[1].px())*(js[0].px()+js[1].px())-(js[0].py()+js[1].py())*(js[0].py()+js[1].py()));
      AlphaT= js[0].pT()/mass;
      }
      else
      {double mass= sqrt((js[0].pT()+js[1].pT())*(js[0].pT()+js[1].pT())-(js[0].px()+js[1].px())*(js[0].px()+js[1].px())-(js[0].py()+js[1].py())*(js[0].py()+js[1].py()));
      AlphaT= js[1].pT()/mass;
      }
    }
    else if (js.size()>2)
    {
    for (unsigned itjet=0;itjet<js.size()-1;itjet++)
      {if (js[itjet].pT()<js[itjet+1].pT())
      cout<<"js pT not in ordered."<<endl;
      }
    double ht1=js[0].pT();
    double ht2=js[1].pT();
    double sumpx=0.0; 
    double sumpy=0.0;
    for (unsigned itjet=0;itjet<js.size();itjet++)
      {sumpx=sumpx+js[itjet].px();
      sumpy=sumpy+js[itjet].py();
    }
    for (unsigned itjet=2;itjet<js.size();itjet++)
      {if (ht1<ht2)
      ht1=ht1+js[itjet].pT();
      else if (ht>ht2)
      ht2=ht2+js[itjet].pT();
      else ht1=ht1+js[itjet].pT(); 
    }
    double deltaht=fabs(ht1-ht2);
    AlphaT= 0.5*(ht-deltaht)/sqrt(ht*ht-sumpx*sumpx-sumpy*sumpy);
    }
    if (ht)
    {if (ht >200 && ht< 250)
      {if (AlphaT<0.65) vetoEvent;}
    if (ht>250&& ht< 300)
      {if (AlphaT<0.60) vetoEvent;}
    if (ht>300 && ht < 350)
      {if (AlphaT<0.55) vetoEvent;}
    if (ht>350 &&ht <400)
      {if (AlphaT<0.53) vetoEvent;}
    if (ht >400 && ht <900)
      {if (AlphaT<0.52) vetoEvent;}
    //if (ht>600 && ht<800)
    //  {if (AlphaT< 0.52) vetoEvent;}
    if(ht > 900) 
      {if (AlphaT < 0) vetoEvent;}
    }
    _h_cutflow->fill(11., event.weight());
    for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
      double tmp=deltaPhi(sum_jet-js[itjet].momentum(), js[itjet]);
      if(tmp > drmax){
      drmax=tmp;
      }
    }
    double angle = M_PI-0.5;
    if (drmax > angle) vetoEvent;
    _h_cutflow->fill(12., event.weight());

   vector <string> trackflavorst; trackflavorst.clear();
   int btagt=0;
    for (unsigned ifjet = 0 ; ifjet < js.size(); ifjet++) {
       int isbtagg(0.0);
       double eff(0.0);
       for(unsigned itjet = 0 ; itjet < jss.size(); itjet++){
         if(trackjetfjidx1[itjet]==ifjet) {
         trackflavorst.push_back(getjetflavor(jss[itjet]));
         break;
         }
       }
       for (unsigned i=0; i<trackflavorst.size(); i++){
         if ( trackflavorst[i] == "bjet" ){
           eff = 0.77;
           isbtagg=isTagging(eff);
           if ( isbtagg!=0 ) btagt=btagt+1; 
         }
       }
       trackflavorst.clear(); 
    }
    _h_bjetssizetry->fill(btagt,                                                             event.weight());



    //if (sum_jet.pT()/smearedmet.mod() > 1.25) vetoEvent;
    //_h_cutflow->fill(4., event.weight());
   //checking whether all the calojets can be targetted by trackjets,if not, discrepency++
   int discrepency=0;
   for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
     double drmin(999.);
     for(unsigned ifjet = 0 ; ifjet < jss.size(); ifjet++) {
        double tmp=deltaR(jss[ifjet], js[itjet]);
        if(tmp < drmin)
         drmin=tmp;
     }
     if(drmin>0.4){
       discrepency=discrepency+1;
       cout<<"          Trackjets fail to decorate this calojet."<<endl;
       cout<<"calojets number   "<<itjet<<endl;
       cout<<"calojets.pT "<<js[itjet].pT()<<endl;
       cout<<"calojets.eta "<<js[itjet].eta()<<endl;
       cout<<"drmin between trackjets and calojets              "<<drmin<<endl;
    }}
    if (js.size()!=jss.size()){
     cout<<"		calojets.size()	"<<js.size()<<endl;
     cout<<"		trackjets.size()	"<<jss.size()<<endl;
     for (unsigned itjet = 0 ; itjet < js.size(); itjet++){
      cout<<"calojets.pT	"<<js[itjet].pT()<<endl;
      cout<<"calojets.phi        "<<js[itjet].phi()<<endl;
      cout<<"calojets.abseta        "<<js[itjet].abseta()<<endl;}
      for (unsigned itjet = 0 ; itjet < jss.size(); itjet++){
      cout<<"trackjets.pT        "<<jss[itjet].pT()<<endl;
      cout<<"trackjets.phi        "<<jss[itjet].phi()<<endl;
      cout<<"trackjets.abseta        "<<jss[itjet].abseta()<<endl;}
    } 
    std::cout<<" 		calo jets in the event " << js.size() << " 		track jets " << jss.size()<< std::endl; 
    cout<<"ht: "<<ht<<endl;
    cout<<"htmiss: "<<sum_jet.pT()<<endl;
    cout<<"trackjets.size: "<<jss.size()<<endl;
    cout<<"calojets.size: "<<js.size()<<endl;
    cout<<"ANGLE    "<<angle<<endl;
    cout<<"drmax    "<<drmax<<endl;
    cout<<"Alpha_T	"<<AlphaT<<endl;
    cout<<"discrepency	"<<discrepency<<endl;
    cout<<"b_jets	"<<btags<<endl;


    for(unsigned itjet = 0 ; itjet < jss.size(); itjet++) {
     _h_trackjetsconstituentsum->fill(jss[itjet].constituents().size(),          event.weight());
     if (itjet==0)
      _h_trackjetsconstituent0->fill(jss[itjet].constituents().size(),		event.weight());
     if (itjet==1)
      _h_trackjetsconstituent1->fill(jss[itjet].constituents().size(),          event.weight());
     if (itjet==2)
      _h_trackjetsconstituent2->fill(jss[itjet].constituents().size(),          event.weight());
     if (itjet==3)
      _h_trackjetsconstituent3->fill(jss[itjet].constituents().size(),          event.weight());
    }
    for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
     _h_regularjetsconstituentsum->fill(js[itjet].constituents().size(),          event.weight());
     if (itjet==0)
      _h_regularjetsconstituent0->fill(js[itjet].constituents().size(),          event.weight());
     if (itjet==1)
      _h_regularjetsconstituent1->fill(js[itjet].constituents().size(),          event.weight());
     if (itjet==2)
      _h_regularjetsconstituent2->fill(js[itjet].constituents().size(),          event.weight());
     if (itjet==3)
      _h_regularjetsconstituent3->fill(js[itjet].constituents().size(),          event.weight());
    }
    //cout<<"jet 0 pT"<<js[0].perp()<<endl;
    //cout<<"jet 0 eta"<<js[0].eta()<<endl; 
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
    for(unsigned itjet = 0 ; itjet < jss.size(); itjet++) {
       double drminb=10000.0;
       for(unsigned ifjet = 0 ; ifjet < calojets.size(); ifjet++) {
        double  tmp=deltaR(jss[itjet], calojets[ifjet]);
        if (tmp <drminb)
        drminb=tmp;
       } 
       if (itjet==0)
       _h_dRtrackjets0->fill(drminb,	event.weight());
       if (itjet==1)
       _h_dRtrackjets1->fill(drminb,      event.weight());
       if (itjet==2)
       _h_dRtrackjets2->fill(drminb,      event.weight());
       if (itjet==3)
       _h_dRtrackjets3->fill(drminb,      event.weight());
    }
    for(unsigned itjet = 0 ; itjet < calojets.size(); itjet++) {
       double drminb=100000.0;
       for(unsigned ifjet = 0 ; ifjet < jss.size(); ifjet++) {
        double tmp=deltaR(jss[ifjet], calojets[itjet]);
        if (tmp <drminb)
        drminb=tmp;
       }
       if (itjet==0)
       _h_dRregularjets0->fill(drminb,      event.weight());
       if (itjet==1)
       _h_dRregularjets1->fill(drminb,      event.weight());
       if (itjet==2)
       _h_dRregularjets2->fill(drminb,      event.weight());
       if (itjet==3)
       _h_dRregularjets3->fill(drminb,      event.weight());
    }
      
      // calojet[0] related quantities
     // _h_fatjetpt->fill(js[0].perp()/GeV,                                     event.weight());
      //_h_fatjeteta->fill(js[0].eta(),                                         event.weight());
      //_h_fatjetmass->fill(calojets[0].mass()/GeV,                             event.weight());
      //_h_fatjetmasscorr->fill(tfjcorrmass[0]/GeV,                             event.weight());
      //_h_fatjetnsub21->fill(nSub21(calojets[0]),                              event.weight());
     // _h_fatjetc2->fill(c2(calojets[0]),                                      event.weight());
      //_h_fatjetd2->fill(d2(calojets[0]),                                      event.weight());
      //_h_fatjetntjets->fill(tfjntracks[0],                                    event.weight());
      //_h_fatjetdphimet->fill(deltaPhi(calojets[0].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());
     // for(unsigned itjet = 0 ; itjet < calojets.size(); itjet++) {
     // _h_regularjetsht->fill(calojets[itjet].perp()/GeV,                      event.weight());
     //}
    _h_ht->fill(ht,    						         	event.weight()); 
    _h_httrackjets->fill(httrackjets,						event.weight());
    _h_htmiss->fill(sum_jet.pT(),						event.weight());
    _h_htmht->fill(ht,sum_jet.pT(),        					event.weight());
    _h_nbnjet->fill(btags,js.size(),						event.weight()); 
    _h_HtmissMETDivision->fill(sum_jet.pT()/smearedmet.mod(),		        event.weight());	
    _h_AlphaT->fill(AlphaT,							event.weight());
    _h_AzimuthalAngular->fill(M_PI-drmax,					event.weight());
    _h_regularjetssize->fill(js.size(),                                         event.weight());
    _h_trackjetssize->fill(jss.size(), 	                                        event.weight());
    _h_bjetssize->fill(btags,                                                   event.weight());
    
    if (btags==0)
    _h_htmhtnb0->fill(ht,sum_jet.pT(),                                             event.weight());
    else if (btags==1)
    _h_htmhtnb1->fill(ht,sum_jet.pT(),                                             event.weight());
    else if (btags==2)    
    _h_htmhtnb2->fill(ht,sum_jet.pT(),                                             event.weight());
    else
    _h_htmhtnb3->fill(ht,sum_jet.pT(),                                             event.weight());

    if (js.size()==1)
    {_h_regularjetseta0->fill(js[0].eta(),                                      event.weight());
    _h_regularjetspt0->fill(js[0].perp()/GeV,                                   event.weight());
    //_h_regularjet0dphimet->fill(deltaPhi(calojets[itjet].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());}
    }
    else if (js.size()==2)
    {
    _h_regularjetseta0->fill(js[0].eta(),                                       event.weight());
    _h_regularjetspt0->fill(js[0].perp()/GeV,                                   event.weight());
    _h_regularjetseta1->fill(js[1].eta(),                                       event.weight());
    _h_regularjetspt1->fill(js[1].perp()/GeV,                                   event.weight());
    //_h_regularjet1dphimet->fill(deltaPhi(calojets[itjet].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());}
    }
    else if (js.size()==3)
    {
    _h_regularjetseta0->fill(js[0].eta(),                                       event.weight());
    _h_regularjetspt0->fill(js[0].perp()/GeV,                                   event.weight());
    _h_regularjetseta1->fill(js[1].eta(),                                       event.weight());
    _h_regularjetspt1->fill(js[1].perp()/GeV,                                   event.weight());
    _h_regularjetseta2->fill(js[2].eta(),                                       event.weight());
    _h_regularjetspt2->fill(js[2].perp()/GeV,                                   event.weight());
    //_h_regularjet2dphimet->fill(deltaPhi(calojets[itjet].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());}
    }
    //if (js.size()==4)
   else
   {
   _h_regularjetseta0->fill(js[0].eta(),                                        event.weight());
   _h_regularjetspt0->fill(js[0].perp()/GeV,                                    event.weight());
   _h_regularjetseta1->fill(js[1].eta(),                                        event.weight());
   _h_regularjetspt1->fill(js[1].perp()/GeV,                                    event.weight());
   _h_regularjetseta2->fill(js[2].eta(),                                        event.weight());
   _h_regularjetspt2->fill(js[2].perp()/GeV,                                    event.weight());
   _h_regularjetseta3->fill(js[3].eta(),                                        event.weight());
   _h_regularjetspt3->fill(js[3].perp()/GeV,                                    event.weight());
   //_h_regularjet3dphimet->fill(deltaPhi(calojets[itjet].phi_02pi(),met.phi(ZERO_2PI))/PI, event.weight());}
   }
   //if (itjet>3)
   //break;

   if (jss.size()==1)
   {_h_trackjetseta0->fill(jss[0].eta(),                                      event.weight());
    _h_trackjetspt0->fill(jss[0].perp()/GeV,                                   event.weight());
   }
   else if (jss.size()==2)
   {
    _h_trackjetseta0->fill(jss[0].eta(),                                       event.weight());
    _h_trackjetspt0->fill(jss[0].perp()/GeV,                                   event.weight());
    _h_trackjetseta1->fill(jss[1].eta(),                                       event.weight());
    _h_trackjetspt1->fill(jss[1].perp()/GeV,                                   event.weight());
   }
   else if (jss.size()==3)
   {
    _h_trackjetseta0->fill(jss[0].eta(),                                       event.weight());
    _h_trackjetspt0->fill(jss[0].perp()/GeV,                                   event.weight());
    _h_trackjetseta1->fill(jss[1].eta(),                                       event.weight());
    _h_trackjetspt1->fill(jss[1].perp()/GeV,                                   event.weight());
    _h_trackjetseta2->fill(jss[2].eta(),                                       event.weight());
    _h_trackjetspt2->fill(jss[2].perp()/GeV,                                   event.weight());
   }
   else 
   {
   _h_trackjetseta0->fill(jss[0].eta(),                                        event.weight());
   _h_trackjetspt0->fill(jss[0].perp()/GeV,                                    event.weight());
   _h_trackjetseta1->fill(jss[1].eta(),                                        event.weight());
   _h_trackjetspt1->fill(jss[1].perp()/GeV,                                    event.weight());
   _h_trackjetseta2->fill(jss[2].eta(),                                        event.weight());
   _h_trackjetspt2->fill(jss[2].perp()/GeV,                                    event.weight());
   _h_trackjetseta3->fill(jss[3].eta(),                                        event.weight());
   _h_trackjetspt3->fill(jss[3].perp()/GeV,                                    event.weight());
   }
   // event properties 
   _h_weight->fill(event.weight(),               					    1.);	  
   _h_truemet->fill(truemet.pT()/GeV, 						event.weight());	  
   _h_met->fill(met.mod()/GeV,         						event.weight());	  
   _h_smearedmet4->fill(smearedmet.mod()/GeV,		         		event.weight());	  
   _h_metdiff->fill(truemet.pT()/GeV-met.mod()/GeV,				event.weight());	  
   _h_smearedmetdiff->fill(met.mod()/GeV-smearedmet.mod()/GeV,			event.weight());	  
   _h_aplanarity->fill(aplanarity,    						event.weight());	  
   _h_sphericity->fill(sphericity,   					        event.weight());	  
      
   }//analyze
      
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      double lumi = 35900.; 
      double norm=lumi*crossSection()/picobarn / sumOfWeights();
      
      std::cout<<" lumi*crossSection()/picobarn " << lumi*crossSection()/picobarn << " sum of weiht " << sumOfWeights() << " nrm  "<< norm << std::endl; 
      
      // _h_fatjetpt->scaleW(norm*_h_fatjetpt->bin(0).width());
      // _h_fatjeteta->scaleW(norm*_h_fatjeteta->bin(0).width());
      // _h_fatjetmass->scaleW(norm*_h_fatjetmass->bin(0).width());
      // _h_fatjetmasscorr->scaleW(norm*_h_fatjetmasscorr->bin(0).width());
      // _h_fatjetnsub21->scaleW(norm*_h_fatjetnsub21->bin(0).width());
      // _h_fatjetc2->scaleW(norm*_h_fatjetc2->bin(0).width());
      // _h_fatjetd2->scaleW(norm*_h_fatjetd2->bin(0).width());
      // _h_fatjetntjets->scaleW(norm*_h_fatjetntjets->bin(0).width());
      // _h_fatjetdphimet->scaleW(norm*_h_fatjetdphimet->bin(0).width());
      // _h_cutflow->scaleW(norm*_h_cutflow->bin(0).width());
      // _h_truemet->scaleW(norm*_h_truemet->bin(0).width());
      // _h_met->scaleW(norm*_h_met->bin(0).width());
      // _h_smearedmet->scaleW(norm*_h_smearedmet->bin(0).width());
      // _h_metdiff->scaleW(norm*_h_metdiff->bin(0).width());
      // _h_aplanarity->scaleW(norm*_h_aplanarity->bin(0).width());
      // _h_sphericity->scaleW(norm*_h_sphericity->bin(0).width());
      
      _h_fatjetpt->scaleW(norm);
      _h_fatjeteta->scaleW(norm);
      _h_fatjetmass->scaleW(norm);
      _h_fatjetmasscorr->scaleW(norm);
      _h_fatjetnsub21->scaleW(norm);
      _h_fatjetc2->scaleW(norm);
      _h_fatjetd2->scaleW(norm);
      _h_fatjetntjets->scaleW(norm);
      _h_fatjetdphimet->scaleW(norm);
      _h_cutflow->scaleW(norm);
      _h_truemet->scaleW(norm);
      _h_met->scaleW(norm);
      _h_smearedmet3->scaleW(norm);
      _h_smearedmet4->scaleW(norm);
      _h_metdiff->scaleW(norm);
      _h_htbefore->scaleW(norm);
      _h_htmissbefore->scaleW(norm);
      _h_ht->scaleW(norm);
      _h_httrackjets->scaleW(norm);
      _h_aplanarity->scaleW(norm);
      _h_sphericity->scaleW(norm);
      _h_HtmissMETDivision->scaleW(norm);
      _h_htmiss->scaleW(norm);
      _h_AlphaT->scaleW(norm);
      _h_AzimuthalAngular->scaleW(norm);
      // _h_truemet2->scaleW(norm);
      //_h_met2->scaleW(norm);
      _h_smearedmet2->scaleW(norm);
      //_h_truemet1->scaleW(norm);
      // _h_met1->scaleW(norm);
      _h_smearedmet1->scaleW(norm);
      // _h_truemet0->scaleW(norm);
      //_h_met0->scaleW(norm);
      _h_trackjetssize->scaleW(norm);
      _h_regularjetssizebefore->scaleW(norm);
      _h_regularjetssize->scaleW(norm);
      _h_bjetssizetry->scaleW(norm);
      _h_bjetssize->scaleW(norm);
      _h_bjetssizebefore->scaleW(norm);
      _h_htmht->scaleW(norm);
      _h_htmhtnb0->scaleW(norm);
      _h_htmhtnb1->scaleW(norm);
      _h_htmhtnb2->scaleW(norm);
      _h_htmhtnb3->scaleW(norm);
      _h_nbnjet->scaleW(norm);
      _h_smearedmet0->scaleW(norm);
      _h_dRregularjets0->scaleW(norm);
      _h_dRregularjets1->scaleW(norm);
      _h_dRregularjets2->scaleW(norm);
      _h_dRregularjets3->scaleW(norm);
      _h_dRtrackjets0->scaleW(norm);
      _h_dRtrackjets1->scaleW(norm);
      _h_dRtrackjets2->scaleW(norm);
      _h_dRtrackjets3->scaleW(norm);
      _h_regularjetseta0->scaleW(norm);
      _h_regularjetseta1->scaleW(norm);
      _h_regularjetseta2->scaleW(norm);
      _h_regularjetseta3->scaleW(norm);
      _h_regularjetspt0->scaleW(norm);
      _h_regularjetspt1->scaleW(norm);
      _h_regularjetspt2->scaleW(norm);
      _h_regularjetspt3->scaleW(norm);
      _h_trackjetseta0->scaleW(norm);
      _h_trackjetseta1->scaleW(norm);
      _h_trackjetseta2->scaleW(norm);
      _h_trackjetseta3->scaleW(norm);
      _h_trackjetspt0->scaleW(norm);
      _h_trackjetspt1->scaleW(norm);
      _h_trackjetspt2->scaleW(norm);
      _h_trackjetspt3->scaleW(norm);
      //_h_regularjet0dphimet->scaleW(norm);
      //_h_regularjet1dphimet->scaleW(norm);
      //_h_regularjet2dphimet->scaleW(norm);
      //_h_regularjet3dphimet->scaleW(norm);
      _h_regularjetsht->scaleW(norm);
      _h_regularjetsconstituent0->scaleW(norm);
      _h_regularjetsconstituent1->scaleW(norm);
      _h_regularjetsconstituent2->scaleW(norm);
      _h_regularjetsconstituent3->scaleW(norm);
      _h_trackjetsconstituent0->scaleW(norm);
      _h_trackjetsconstituent1->scaleW(norm);
      _h_trackjetsconstituent2->scaleW(norm);
      _h_trackjetsconstituent3->scaleW(norm);
      _h_trackjetsconstituentsum->scaleW(norm);
      _h_regularjetsconstituentsum->scaleW(norm);

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
    //  Particles bpartons; bpartons.clear();
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
    Histo1DPtr   _h_smearedmet3;
    Histo1DPtr   _h_smearedmet4;
    Histo1DPtr   _h_metdiff;
    Histo1DPtr   _h_smearedmetdiff;
    Histo1DPtr   _h_sphericity;
    Histo1DPtr   _h_aplanarity;
    Histo1DPtr   _h_ht;
    Histo1DPtr   _h_httrackjets;
    Histo1DPtr   _h_htmiss;
    Histo1DPtr   _h_htbefore;
    Histo1DPtr   _h_htmissbefore;
    Histo1DPtr   _h_trackjetssize;
    Histo1DPtr   _h_regularjetssizebefore;
    Histo1DPtr   _h_regularjetssize;
    Histo1DPtr   _h_bjetssizetry;
    Histo1DPtr   _h_bjetssize;
    Histo1DPtr   _h_bjetssizebefore;
    Histo2DPtr   _h_htmht;
    Histo2DPtr   _h_htmhtnb0;
    Histo2DPtr   _h_htmhtnb1;
    Histo2DPtr   _h_htmhtnb2;
    Histo2DPtr   _h_htmhtnb3;
    Histo2DPtr   _h_nbnjet;
    Histo1DPtr   _h_AlphaT;
    Histo1DPtr   _h_AzimuthalAngular;
    Histo1DPtr   _h_HtmissMETDivision;
    //Histo1DPtr   _h_truemet2;
    //Histo1DPtr   _h_met2;
    Histo1DPtr   _h_smearedmet2;
    //Histo1DPtr   _h_truemet1;
    //Histo1DPtr   _h_met1;
    Histo1DPtr   _h_smearedmet1;
    //Histo1DPtr   _h_truemet0;
    //Histo1DPtr   _h_met0;
    Histo1DPtr   _h_smearedmet0;
    Histo1DPtr   _h_dRtrackjets0;
    Histo1DPtr   _h_dRtrackjets1;
    Histo1DPtr   _h_dRtrackjets2;
    Histo1DPtr   _h_dRtrackjets3;
    Histo1DPtr   _h_dRregularjets0;
    Histo1DPtr   _h_dRregularjets1;
    Histo1DPtr   _h_dRregularjets2;
    Histo1DPtr   _h_dRregularjets3;
    Histo1DPtr   _h_regularjetseta0;
    Histo1DPtr   _h_regularjetseta1;
    Histo1DPtr   _h_regularjetseta2;
    Histo1DPtr   _h_regularjetseta3;
    Histo1DPtr   _h_regularjetspt0;
    Histo1DPtr   _h_regularjetspt1;
    Histo1DPtr   _h_regularjetspt2;
    Histo1DPtr   _h_regularjetspt3;
    Histo1DPtr   _h_trackjetseta0;
    Histo1DPtr   _h_trackjetseta1;
    Histo1DPtr   _h_trackjetseta2;
    Histo1DPtr   _h_trackjetseta3;
    Histo1DPtr   _h_trackjetspt0;
    Histo1DPtr   _h_trackjetspt1;
    Histo1DPtr   _h_trackjetspt2;
    Histo1DPtr   _h_trackjetspt3;
    Histo1DPtr   _h_trackjetsconstituent0;
    Histo1DPtr   _h_trackjetsconstituent1;
    Histo1DPtr   _h_trackjetsconstituent2;
    Histo1DPtr   _h_trackjetsconstituent3;
    Histo1DPtr   _h_regularjetsconstituent0;
    Histo1DPtr   _h_regularjetsconstituent1;
    Histo1DPtr   _h_regularjetsconstituent2;
    Histo1DPtr   _h_regularjetsconstituent3;
    Histo1DPtr   _h_trackjetsconstituentsum;
    Histo1DPtr   _h_regularjetsconstituentsum;

    //Histo1DPtr   _h_regularjet0dphimet;
    //Histo1DPtr   _h_regularjet1dphimet;
    //Histo1DPtr   _h_regularjet2dphimet;
    //Histo1DPtr   _h_regularjet3dphimet;
    Histo1DPtr   _h_regularjetsht;    
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

