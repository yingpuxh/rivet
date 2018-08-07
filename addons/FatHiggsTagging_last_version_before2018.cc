// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/AnalysisHandler.hh"
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
      vector <double> bin={0.0, 0.2, 0.3, 0.6, 0.8, 1.0};
      _h_variablebinwideth = bookHisto1D("variablebinwidth", bin); 
      _h_fatjetpt         = bookHisto1D("fatjetpt",           26,  200.,  1500.);
      _h_fatjeteta        = bookHisto1D("fatjeteta",          25,  -2.5,    2.5);
      _h_fatjetmass       = bookHisto1D("fatjetmass",         14,    0.,   280.);
      _h_fatjetmasscorr   = bookHisto1D("fatjetmasscorr",     14,    0.,   280.);
      _h_fatjetnsub21     = bookHisto1D("fatjetnsub21",       16,   0.0,    1.6);
      _h_fatjetc2         = bookHisto1D("fatjetc2",           14,   0.0,    0.7);
      _h_fatjetd2         = bookHisto1D("fatjetd2",           25,   0.0,    5.0);
      _h_fatjetntjets     = bookHisto1D("fatjetntjets",        7,  -0.5,    6.5);
      _h_fatjetdphimet    = bookHisto1D("fatjetdphimet",      20,   0.0,    1.0);
      _h_cutflow          = bookHisto1D("cutflow",           316,  -0.5,    315.5);
      _h_cutflownjet4          = bookHisto1D("cutflownjet4",           307,  -0.5,    306.5);
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
      _h_htmiss		  = bookHisto1D("htmiss",             40,   0.0,   4000.);
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
      _h_ht               = bookHisto1D("ht",    	      40, 0.0,  4000.);
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
      _h_regularjetssize  = bookHisto1D("regularjetssize",     4,   1.,   5.);
      _h_bjetssizetry  = bookHisto1D("bjetssizetry",     5,   0.,   5.);
      _h_bjetssize        = bookHisto1D("bjetssize",           3,   0.,   3.);
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
      //_h_ch1bkg          = bookHisto1D("ch1bkg",           4,  -0.5,    3.5);
      vector <double> sbin={130., 225., 275., 325., 375., 425., 475., 525., 575., 625., 675., 725., 775., 9999.};
      _h_ch1sig= bookHisto1D("ch1sig", sbin);      
      _h_ch1try= bookHisto1D("ch1try",sbin);
      _h_ch1data = bookHisto1D("data_obs_ch1", sbin);
      _h_ch1bkg = bookHisto1D("ch1bkg", sbin);
      _h_ch1sigbefore          = bookHisto1D("ch1sigbefore",           4,  -0.5,    3.5);
      vector <double> sssbin={130., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 9999.};
      _h_ch2bkg          = bookHisto1D("ch2bkg",           sssbin);
      _h_ch2sig          = bookHisto1D("ch2sig",           sssbin);
      _h_ch2data = bookHisto1D("data_obs_ch2", sssbin);
      vector <double> ssbin={130., 175., 225., 275., 325., 375., 425., 475., 9999.};
      _h_ch3bkg          = bookHisto1D("ch3bkg",           ssbin);
      _h_ch3sig          = bookHisto1D("ch3sig",           ssbin);
      _h_ch3data = bookHisto1D("data_obs_ch3", ssbin);
      vector <double> ch4bin={130., 225., 275., 325., 375., 425., 475., 525., 575., 9999.};
      _h_ch4bkg          = bookHisto1D("ch4bkg",           ch4bin);
      _h_ch4sig          = bookHisto1D("ch4sig",           ch4bin);
      _h_ch4data = bookHisto1D("data_obs_ch4", ch4bin);
      vector <double> ch5bin={130., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 9999.}; 
      _h_ch5bkg          = bookHisto1D("ch5bkg",           ch5bin);
      _h_ch5sig          = bookHisto1D("ch5sig",           ch5bin);
      _h_ch5data = bookHisto1D("data_obs_ch5", ch5bin);
      _h_ch6data = bookHisto1D("data_obs_ch6", sssbin);
      _h_ch6bkg          = bookHisto1D("ch6bkg",           sssbin);
      _h_ch6sig          = bookHisto1D("ch6sig",           sssbin);
      vector <double> ch7bin={130., 9999.}; 
      _h_ch7bkg          = bookHisto1D("ch7bkg",           ch7bin);
      _h_ch7sig          = bookHisto1D("ch7sig",           ch7bin);
      _h_ch7data = bookHisto1D("data_obs_ch7", ch7bin);
      _h_ch8data = bookHisto1D("data_obs_ch8", ch7bin);
      _h_ch8bkg          = bookHisto1D("ch8bkg",           ch7bin);
      _h_ch8sig          = bookHisto1D("ch8sig",           ch7bin);
      vector <double> ch9bin={130., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 9999.}; 
      _h_ch9bkg          = bookHisto1D("ch9bkg",           ch9bin);
      _h_ch9sig          = bookHisto1D("ch9sig",           ch9bin);
      _h_ch9data = bookHisto1D("data_obs_ch9", ch9bin);
      vector <double> ch10bin={130., 200., 250., 300., 350., 400., 450., 500., 550.,9999.};
      _h_ch10bkg          = bookHisto1D("ch10bkg",           ch10bin);
      _h_ch10sig          = bookHisto1D("ch10sig",           ch10bin);
      _h_ch10data = bookHisto1D("data_obs_ch10", ch10bin);
      vector <double> ch11bin={130., 225., 275., 325., 375., 425., 475., 525., 575., 625., 675., 725., 9999.};
      _h_ch11bkg          = bookHisto1D("ch11bkg",           ch11bin);
      _h_ch11sig          = bookHisto1D("ch11sig",           ch11bin);
      _h_ch11data = bookHisto1D("data_obs_ch11", ch11bin);
      _h_ch12data = bookHisto1D("data_obs_ch12", sssbin);
      _h_ch12bkg          = bookHisto1D("ch12bkg",           sssbin);
      _h_ch12sig          = bookHisto1D("ch12sig",           sssbin);
      vector <double> ch13bin={130., 175., 225., 275., 325., 375., 425., 9999.};
      _h_ch13bkg          = bookHisto1D("ch13bkg",           ch13bin);
      _h_ch13sig          = bookHisto1D("ch13sig",           ch13bin);
      _h_ch13data = bookHisto1D("data_obs_ch13", ch13bin);
      vector <double> ch14bin={130., 225., 275., 325., 375., 425., 475., 9999.}; 
      _h_ch14data = bookHisto1D("data_obs_ch14", ch14bin);
      _h_ch14bkg          = bookHisto1D("ch14bkg",           ch14bin);
      _h_ch14sig          = bookHisto1D("ch14sig",           ch14bin);
      vector <double> ch15bin={130., 250., 300., 350., 400., 450., 500., 575., 625., 9999.}; 
      _h_ch15bkg          = bookHisto1D("ch15bkg",           ch15bin);
      _h_ch15sig          = bookHisto1D("ch15sig",           ch15bin);
      _h_ch15data = bookHisto1D("data_obs_ch15", ch15bin);
      vector <double> ch16bin={130., 325, 375., 425., 475., 600., 700., 800., 9999.};  
      _h_ch16bkg          = bookHisto1D("ch16bkg",           ch16bin);
      _h_ch16sig          = bookHisto1D("ch16sig",           ch16bin);
      _h_ch16data = bookHisto1D("data_obs_ch16", ch16bin);
      _h_ch17data = bookHisto1D("data_obs_ch17", ch9bin);
      _h_ch17bkg          = bookHisto1D("ch17bkg",           ch9bin);
      _h_ch17sig          = bookHisto1D("ch17sig",           ch9bin);
      _h_ch18data = bookHisto1D("data_obs_ch18", ch9bin);
      _h_ch18bkg          = bookHisto1D("ch18bkg",           ch9bin);
      _h_ch18sig          = bookHisto1D("ch18sig",           ch9bin);
      vector <double> ch19bin={130., 225., 275., 325., 375., 425., 475., 525., 600., 700., 800., 9999.}; 
      _h_ch19bkg          = bookHisto1D("ch19bkg",           ch19bin);
      _h_ch19sig          = bookHisto1D("ch19sig",           ch19bin);
      _h_ch19data = bookHisto1D("data_obs_ch19", ch19bin);
      _h_ch20data = bookHisto1D("data_obs_ch20", ch7bin);
      _h_ch20bkg          = bookHisto1D("ch20bkg",           ch7bin);
      _h_ch20sig          = bookHisto1D("ch20sig",           ch7bin);
      vector <double> ch21bin={130., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 9999.}; 
      _h_ch21bkg          = bookHisto1D("ch21bkg",           ch21bin);
      _h_ch21sig          = bookHisto1D("ch21sig",           ch21bin);
      _h_ch21data = bookHisto1D("data_obs_ch21", ch21bin);
      vector <double> ch22bin={130., 175., 225., 275., 325., 375., 425., 475., 525., 575., 650., 800., 9999.};
      _h_ch22bkg          = bookHisto1D("ch22bkg",           ch22bin);
      _h_ch22sig          = bookHisto1D("ch22sig",           ch22bin);
      _h_ch22data = bookHisto1D("data_obs_ch22", ch22bin);

      //_h_deltarfjbh              = bookHisto1D("deltarfjbh",                  25,   0.0,    2.5);
      //_h_random                  = bookHisto1D("random_numbers",              50,   0.0,    2.0);
      
      // _all=0;
      // _singlebtag=0;
      // _doublebtag=0;
      // _foundhbbbar=0;
      // _allhbbbar=0;
      _h_variablebinwideth->fill(0., 2.0);
      _h_variablebinwideth->fill(0.2,4.0);
      //_h_ch1sig->fill(575., 3.0);
      //_h_ch1sig->fill(725., 2.0);
      _h_ch1data->fill(130., 2);
      _h_ch1bkg->fill(130., 1.900000);
      _h_ch1data->fill(225., 4);
      _h_ch1bkg->fill(225., 3.380000);
      _h_ch1data->fill(275., 3);
      _h_ch1bkg->fill(275., 3.960000);
      _h_ch1data->fill(325., 8);
      _h_ch1bkg->fill(325., 4.040000);
      _h_ch1data->fill(375., 3);
      _h_ch1bkg->fill(375., 3.500000);
      _h_ch1data->fill(425., 3);
      _h_ch1bkg->fill(425., 3.060000);
      _h_ch1data->fill(475., 2);
      _h_ch1bkg->fill(475., 3.270000);
      _h_ch1data->fill(525., 5);
      _h_ch1bkg->fill(525., 3.840000);
      _h_ch1data->fill(575., 12);
      _h_ch1bkg->fill(575., 8.060000);
      _h_ch1data->fill(625., 11);
      _h_ch1bkg->fill(625., 6.930000);
      _h_ch1data->fill(675., 2);
      _h_ch1bkg->fill(675., 4.200000);
      _h_ch1data->fill(725., 3);
      _h_ch1bkg->fill(725., 1.810000);
      _h_ch1data->fill(775., 0);
      _h_ch1bkg->fill(775., 0.525000);
      _h_ch2data->fill(130., 6);
      _h_ch2bkg->fill(130., 2.680000);
      _h_ch2data->fill(250., 6);
      _h_ch2bkg->fill(250., 6.290000);
      _h_ch2data->fill(300., 8);
      _h_ch2bkg->fill(300., 6.640000);
      _h_ch2data->fill(350., 4);
      _h_ch2bkg->fill(350., 6.470000);
      _h_ch2data->fill(400., 4);
      _h_ch2bkg->fill(400., 5.800000);
      _h_ch2data->fill(450., 5);
      _h_ch2bkg->fill(450., 4.660000);
      _h_ch2data->fill(500., 2);
      _h_ch2bkg->fill(500., 3.330000);
      _h_ch2data->fill(550., 3);
      _h_ch2bkg->fill(550., 3.380000);
      _h_ch2data->fill(600., 4);
      _h_ch2bkg->fill(600., 2.990000);
      _h_ch2data->fill(650., 3);
      _h_ch2bkg->fill(650., 2.410000);
      _h_ch2data->fill(700., 2);
      _h_ch2bkg->fill(700., 1.720000);
      _h_ch2data->fill(750., 1);
      _h_ch2bkg->fill(750., 2.010000);
      _h_ch2data->fill(800., 9);
      _h_ch2bkg->fill(800., 9.600000);
      _h_ch3data->fill(130., 4);
      _h_ch3bkg->fill(130., 1.550000);
      _h_ch3data->fill(175., 6);
      _h_ch3bkg->fill(175., 3.420000);
      _h_ch3data->fill(225., 3);
      _h_ch3bkg->fill(225., 3.070000);
      _h_ch3data->fill(275., 6);
      _h_ch3bkg->fill(275., 3.470000);
      _h_ch3data->fill(325., 3);
      _h_ch3bkg->fill(325., 5.540000);
      _h_ch3data->fill(375., 13);
      _h_ch3bkg->fill(375., 8.400000);
      _h_ch3data->fill(425., 4);
      _h_ch3bkg->fill(425., 5.090000);
      _h_ch3data->fill(475., 0);
      _h_ch3bkg->fill(475., 0.608000);
      _h_ch4data->fill(130., 0);
      _h_ch4bkg->fill(130., 0.434000);
      _h_ch4data->fill(225., 0);
      _h_ch4bkg->fill(225., 0.900000);
      _h_ch4data->fill(275., 1);
      _h_ch4bkg->fill(275., 1.280000);
      _h_ch4data->fill(325., 1);
      _h_ch4bkg->fill(325., 0.712000);
      _h_ch4data->fill(375., 1);
      _h_ch4bkg->fill(375., 0.968000);
      _h_ch4data->fill(425., 0);
      _h_ch4bkg->fill(425., 0.982000);
      _h_ch4data->fill(475., 2);
      _h_ch4bkg->fill(475., 3.560000);
      _h_ch4data->fill(525., 0);
      _h_ch4bkg->fill(525., 2.010000);
      _h_ch4data->fill(575., 0);
      _h_ch4bkg->fill(575., 0.985000);
      _h_ch5data->fill(130., 0);
      _h_ch5bkg->fill(130., 0.367000);
      _h_ch5data->fill(250., 0);
      _h_ch5bkg->fill(250., 0.387000);
      _h_ch5data->fill(300., 0);
      _h_ch5bkg->fill(300., 0.442000);
      _h_ch5data->fill(350., 0);
      _h_ch5bkg->fill(350., 0.187000);
      _h_ch5data->fill(400., 0);
      _h_ch5bkg->fill(400., 0.322000);
      _h_ch5data->fill(450., 1);
      _h_ch5bkg->fill(450., 0.424000);
      _h_ch5data->fill(500., 2);
      _h_ch5bkg->fill(500., 0.294000);
      _h_ch5data->fill(550., 1);
      _h_ch5bkg->fill(550., 0.688000);
      _h_ch5data->fill(600., 0);
      _h_ch5bkg->fill(600., 0.760000);
      _h_ch5data->fill(650., 0);
      _h_ch5bkg->fill(650., 0.420000);
      _h_ch5data->fill(700., 0);
      _h_ch5bkg->fill(700., 0.534000);
      _h_ch5data->fill(750., 0);
      _h_ch5bkg->fill(750., 0.108000);
      _h_ch6data->fill(130., 0);
      _h_ch6bkg->fill(130., 0.219000);
      _h_ch6data->fill(250., 0);
      _h_ch6bkg->fill(250., 0.517000);
      _h_ch6data->fill(300., 0);
      _h_ch6bkg->fill(300., 0.497000);
      _h_ch6data->fill(350., 0);
      _h_ch6bkg->fill(350., 0.479000);
      _h_ch6data->fill(400., 0);
      _h_ch6bkg->fill(400., 0.299000);
      _h_ch6data->fill(450., 1);
      _h_ch6bkg->fill(450., 0.352000);
      _h_ch6data->fill(500., 1);
      _h_ch6bkg->fill(500., 0.252000);
      _h_ch6data->fill(550., 0);
      _h_ch6bkg->fill(550., 0.214000);
      _h_ch6data->fill(600., 0);
      _h_ch6bkg->fill(600., 0.201000);
      _h_ch6data->fill(650., 0);
      _h_ch6bkg->fill(650., 0.273000);
      _h_ch6data->fill(700., 0);
      _h_ch6bkg->fill(700., 0.239000);
      _h_ch6data->fill(750., 0);
      _h_ch6bkg->fill(750., 0.137000);
      _h_ch6data->fill(800., 0);
      _h_ch6bkg->fill(800., 0.711000);
      _h_ch7data->fill(130., 3);
      _h_ch7bkg->fill(130., 1.310000);
      _h_ch8data->fill(130., 0);
      _h_ch8bkg->fill(130., 0.269000);
      _h_ch9data->fill(130., 1);
      _h_ch9bkg->fill(130., 0.078600);
      _h_ch9data->fill(200., 3);
      _h_ch9bkg->fill(200., 4.460000);
      _h_ch9data->fill(250., 7);
      _h_ch9bkg->fill(250., 9.470000);
      _h_ch9data->fill(300., 11);
      _h_ch9bkg->fill(300., 9.740000);
      _h_ch9data->fill(350., 13);
      _h_ch9bkg->fill(350., 8.700000);
      _h_ch9data->fill(400., 10);
      _h_ch9bkg->fill(400., 7.620000);
      _h_ch9data->fill(450., 1);
      _h_ch9bkg->fill(450., 6.430000);
      _h_ch9data->fill(500., 6);
      _h_ch9bkg->fill(500., 5.720000);
      _h_ch9data->fill(550., 5);
      _h_ch9bkg->fill(550., 4.760000);
      _h_ch9data->fill(600., 8);
      _h_ch9bkg->fill(600., 3.870000);
      _h_ch9data->fill(650., 4);
      _h_ch9bkg->fill(650., 3.450000);
      _h_ch9data->fill(700., 3);
      _h_ch9bkg->fill(700., 2.790000);
      _h_ch9data->fill(750., 2);
      _h_ch9bkg->fill(750., 2.670000);
      _h_ch9data->fill(800., 5);
      _h_ch9bkg->fill(800., 8.270000);
      _h_ch10data->fill(130., 5);
      _h_ch10bkg->fill(130., 2.280000);
      _h_ch10data->fill(200., 10);
      _h_ch10bkg->fill(200., 5.800000);
      _h_ch10data->fill(250., 4);
      _h_ch10bkg->fill(250., 6.330000);
      _h_ch10data->fill(300., 3);
      _h_ch10bkg->fill(300., 4.640000);
      _h_ch10data->fill(350., 3);
      _h_ch10bkg->fill(350., 4.230000);
      _h_ch10data->fill(400., 6);
      _h_ch10bkg->fill(400., 3.720000);
      _h_ch10data->fill(450., 1);
      _h_ch10bkg->fill(450., 3.830000);
      _h_ch10data->fill(500., 0);
      _h_ch10bkg->fill(500., 2.170000);
      _h_ch10data->fill(550., 0);
      _h_ch10bkg->fill(550., 0.742000);
      _h_ch11data->fill(130., 1);
      _h_ch11bkg->fill(130., 0.401000);
      _h_ch11data->fill(225., 0);
      _h_ch11bkg->fill(225., 2.050000);
      _h_ch11data->fill(275., 1);
      _h_ch11bkg->fill(275., 2.730000);
      _h_ch11data->fill(325., 3);
      _h_ch11bkg->fill(325., 3.100000);
      _h_ch11data->fill(375., 3);
      _h_ch11bkg->fill(375., 2.680000);
      _h_ch11data->fill(425., 3);
      _h_ch11bkg->fill(425., 1.890000);
      _h_ch11data->fill(475., 2);
      _h_ch11bkg->fill(475., 2.130000);
      _h_ch11data->fill(525., 3);
      _h_ch11bkg->fill(525., 1.870000);
      _h_ch11data->fill(575., 1);
      _h_ch11bkg->fill(575., 1.870000);
      _h_ch11data->fill(625., 0);
      _h_ch11bkg->fill(625., 1.150000);
      _h_ch11data->fill(675., 1);
      _h_ch11bkg->fill(675., 0.600000);
      _h_ch11data->fill(725., 0);
      _h_ch11bkg->fill(725., 0.214000);
      _h_ch12data->fill(130., 0);
      _h_ch12bkg->fill(130., 0.859000);
      _h_ch12data->fill(250., 1);
      _h_ch12bkg->fill(250., 1.280000);
      _h_ch12data->fill(300., 0);
      _h_ch12bkg->fill(300., 1.750000);
      _h_ch12data->fill(350., 2);
      _h_ch12bkg->fill(350., 1.160000);
      _h_ch12data->fill(400., 1);
      _h_ch12bkg->fill(400., 1.190000);
      _h_ch12data->fill(450., 1);
      _h_ch12bkg->fill(450., 1.090000);
      _h_ch12data->fill(500., 0);
      _h_ch12bkg->fill(500., 1.050000);
      _h_ch12data->fill(550., 2);
      _h_ch12bkg->fill(550., 0.534000);
      _h_ch12data->fill(600., 0);
      _h_ch12bkg->fill(600., 0.368000);
      _h_ch12data->fill(650., 1);
      _h_ch12bkg->fill(650., 0.416000);
      _h_ch12data->fill(700., 1);
      _h_ch12bkg->fill(700., 0.454000);
      _h_ch12data->fill(750., 1);
      _h_ch12bkg->fill(750., 0.340000);
      _h_ch12data->fill(800., 0);
      _h_ch12bkg->fill(800., 1.150000);
      _h_ch13data->fill(130., 1);
      _h_ch13bkg->fill(130., 1.280000);
      _h_ch13data->fill(175., 2);
      _h_ch13bkg->fill(175., 3.280000);
      _h_ch13data->fill(225., 7);
      _h_ch13bkg->fill(225., 4.200000);
      _h_ch13data->fill(275., 2);
      _h_ch13bkg->fill(275., 2.700000);
      _h_ch13data->fill(325., 0);
      _h_ch13bkg->fill(325., 2.790000);
      _h_ch13data->fill(375., 3);
      _h_ch13bkg->fill(375., 1.590000);
      _h_ch13data->fill(425., 1);
      _h_ch13bkg->fill(425., 0.190000);
      _h_ch14data->fill(130., 1);
      _h_ch14bkg->fill(130., 0.698000);
      _h_ch14data->fill(225., 1);
      _h_ch14bkg->fill(225., 0.666000);
      _h_ch14data->fill(275., 0);
      _h_ch14bkg->fill(275., 1.280000);
      _h_ch14data->fill(325., 1);
      _h_ch14bkg->fill(325., 0.575000);
      _h_ch14data->fill(375., 1);
      _h_ch14bkg->fill(375., 0.765000);
      _h_ch14data->fill(425., 0);
      _h_ch14bkg->fill(425., 0.637000);
      _h_ch14data->fill(475., 1);
      _h_ch14bkg->fill(475., 0.460000);
      _h_ch15data->fill(130., 0);
      _h_ch15bkg->fill(130., 0.147000);
      _h_ch15data->fill(250., 0);
      _h_ch15bkg->fill(250., 0.157000);
      _h_ch15data->fill(300., 0);
      _h_ch15bkg->fill(300., 0.090200);
      _h_ch15data->fill(350., 0);
      _h_ch15bkg->fill(350., 0.165000);
      _h_ch15data->fill(400., 0);
      _h_ch15bkg->fill(400., 0.104000);
      _h_ch15data->fill(450., 0);
      _h_ch15bkg->fill(450., 0.117000);
      _h_ch15data->fill(500., 0);
      _h_ch15bkg->fill(500., 0.108000);
      _h_ch15data->fill(575., 0);
      _h_ch15bkg->fill(575., 0.213000);
      _h_ch15data->fill(625., 1);
      _h_ch15bkg->fill(625., 0.102000);
      _h_ch16data->fill(130., 0);
      _h_ch16bkg->fill(130., 0.166000);
      _h_ch16data->fill(325., 0);
      _h_ch16bkg->fill(325., 0.294000);
      _h_ch16data->fill(375., 0);
      _h_ch16bkg->fill(375., 0.065700);
      _h_ch16data->fill(425., 0);
      _h_ch16bkg->fill(425., 0.174000);
      _h_ch16data->fill(475., 0);
      _h_ch16bkg->fill(475., 0.130000);
      _h_ch16data->fill(600., 0);
      _h_ch16bkg->fill(600., 0.134000);
      _h_ch16data->fill(700., 0);
      _h_ch16bkg->fill(700., 0.171000);
      _h_ch16data->fill(800., 1);
      _h_ch16bkg->fill(800., 0.180000);
      _h_ch17data->fill(130., 1);
      _h_ch17bkg->fill(130., 1.260000);
      _h_ch17data->fill(200., 12);
      _h_ch17bkg->fill(200., 6.720000);
      _h_ch17data->fill(250., 6);
      _h_ch17bkg->fill(250., 9.990000);
      _h_ch17data->fill(300., 9);
      _h_ch17bkg->fill(300., 9.440000);
      _h_ch17data->fill(350., 8);
      _h_ch17bkg->fill(350., 7.950000);
      _h_ch17data->fill(400., 3);
      _h_ch17bkg->fill(400., 7.610000);
      _h_ch17data->fill(450., 2);
      _h_ch17bkg->fill(450., 5.820000);
      _h_ch17data->fill(500., 4);
      _h_ch17bkg->fill(500., 4.410000);
      _h_ch17data->fill(550., 7);
      _h_ch17bkg->fill(550., 3.410000);
      _h_ch17data->fill(600., 3);
      _h_ch17bkg->fill(600., 3.140000);
      _h_ch17data->fill(650., 6);
      _h_ch17bkg->fill(650., 2.420000);
      _h_ch17data->fill(700., 2);
      _h_ch17bkg->fill(700., 2.010000);
      _h_ch17data->fill(750., 2);
      _h_ch17bkg->fill(750., 1.730000);
      _h_ch17data->fill(800., 3);
      _h_ch17bkg->fill(800., 5.380000);
      _h_ch18data->fill(130., 0);
      _h_ch18bkg->fill(130., 0.336000);
      _h_ch18data->fill(200., 1);
      _h_ch18bkg->fill(200., 1.670000);
      _h_ch18data->fill(250., 1);
      _h_ch18bkg->fill(250., 1.870000);
      _h_ch18data->fill(300., 2);
      _h_ch18bkg->fill(300., 1.810000);
      _h_ch18data->fill(350., 1);
      _h_ch18bkg->fill(350., 1.550000);
      _h_ch18data->fill(400., 0);
      _h_ch18bkg->fill(400., 1.520000);
      _h_ch18data->fill(450., 2);
      _h_ch18bkg->fill(450., 1.010000);
      _h_ch18data->fill(500., 1);
      _h_ch18bkg->fill(500., 0.746000);
      _h_ch18data->fill(550., 0);
      _h_ch18bkg->fill(550., 0.713000);
      _h_ch18data->fill(600., 0);
      _h_ch18bkg->fill(600., 0.451000);
      _h_ch18data->fill(650., 0);
      _h_ch18bkg->fill(650., 0.526000);
      _h_ch18data->fill(700., 1);
      _h_ch18bkg->fill(700., 0.206000);
      _h_ch18data->fill(750., 0);
      _h_ch18bkg->fill(750., 0.352000);
      _h_ch18data->fill(800., 1);
      _h_ch18bkg->fill(800., 1.700000);
      _h_ch19data->fill(130., 0);
      _h_ch19bkg->fill(130., 0.271000);
      _h_ch19data->fill(225., 0);
      _h_ch19bkg->fill(225., 0.379000);
      _h_ch19data->fill(275., 0);
      _h_ch19bkg->fill(275., 0.405000);
      _h_ch19data->fill(325., 0);
      _h_ch19bkg->fill(325., 0.723000);
      _h_ch19data->fill(375., 1);
      _h_ch19bkg->fill(375., 0.455000);
      _h_ch19data->fill(425., 0);
      _h_ch19bkg->fill(425., 0.133000);
      _h_ch19data->fill(475., 1);
      _h_ch19bkg->fill(475., 0.239000);
      _h_ch19data->fill(525., 0);
      _h_ch19bkg->fill(525., 0.241000);
      _h_ch19data->fill(600., 0);
      _h_ch19bkg->fill(600., 0.151000);
      _h_ch19data->fill(700., 0);
      _h_ch19bkg->fill(700., 0.148000);
      _h_ch19data->fill(800., 0);
      _h_ch19bkg->fill(800., 0.307000);
      _h_ch20data->fill(130., 0);
      _h_ch20bkg->fill(130., 0.086100);
      _h_ch21data->fill(130., 1);
      _h_ch21bkg->fill(130., 0.261000);
      _h_ch21data->fill(150., 3);
      _h_ch21bkg->fill(150., 2.350000);
      _h_ch21data->fill(200., 3);
      _h_ch21bkg->fill(200., 4.050000);
      _h_ch21data->fill(250., 3);
      _h_ch21bkg->fill(250., 4.320000);
      _h_ch21data->fill(300., 2);
      _h_ch21bkg->fill(300., 3.150000);
      _h_ch21data->fill(350., 2);
      _h_ch21bkg->fill(350., 2.180000);
      _h_ch21data->fill(400., 1);
      _h_ch21bkg->fill(400., 2.150000);
      _h_ch21data->fill(450., 0);
      _h_ch21bkg->fill(450., 1.250000);
      _h_ch21data->fill(500., 3);
      _h_ch21bkg->fill(500., 1.110000);
      _h_ch21data->fill(550., 1);
      _h_ch21bkg->fill(550., 0.863000);
      _h_ch21data->fill(600., 0);
      _h_ch21bkg->fill(600., 0.700000);
      _h_ch21data->fill(650., 2);
      _h_ch21bkg->fill(650., 0.379000);
      _h_ch21data->fill(700., 0);
      _h_ch21bkg->fill(700., 0.555000);
      _h_ch21data->fill(750., 0);
      _h_ch21bkg->fill(750., 0.290000);
      _h_ch21data->fill(800., 0);
      _h_ch21bkg->fill(800., 0.708000);
      _h_ch22data->fill(130., 2);
      _h_ch22bkg->fill(130., 0.718000);
      _h_ch22data->fill(175., 0);
      _h_ch22bkg->fill(175., 0.958000);
      _h_ch22data->fill(225., 2);
      _h_ch22bkg->fill(225., 1.260000);
      _h_ch22data->fill(275., 5);
      _h_ch22bkg->fill(275., 1.380000);
      _h_ch22data->fill(325., 3);
      _h_ch22bkg->fill(325., 0.821000);
      _h_ch22data->fill(375., 2);
      _h_ch22bkg->fill(375., 0.755000);
      _h_ch22data->fill(425., 1);
      _h_ch22bkg->fill(425., 0.314000);
      _h_ch22data->fill(475., 1);
      _h_ch22bkg->fill(475., 0.061900);
      _h_ch22data->fill(525., 0);
      _h_ch22bkg->fill(525., 0.147000);
      _h_ch22data->fill(575., 0);
      _h_ch22bkg->fill(575., 0.161000);
      _h_ch22data->fill(650., 0);
      _h_ch22bkg->fill(650., 0.436000);
      _h_ch22data->fill(800., 0);
      _h_ch22bkg->fill(800., 0.146000);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      static random_device rd;
      static mt19937 gen(rd());

      //_h_cutflow->fill(0., event.weight());
      // _h_smearedmet0->fill(smearedmet.mod()/GeV,event.weight());

      // veto events with isolated prompt dressed lepton above certain pT
      const vector<DressedLepton> dressedelectrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      const vector<DressedLepton> dressedmuons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      if (dressedelectrons.size()>0) vetoEvent;
      if (dressedmuons.size()>0) vetoEvent;
      //_h_cutflow->fill(1., event.weight());
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
      //_h_cutflow->fill(2., event.weight());
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
      //_h_cutflow->fill(3., event.weight());
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
      for(unsigned itjet = 0 ; itjet < js.size(); itjet++){
       if(js[itjet].perp()/GeV>40 && js[itjet].abseta()>3.0)
       vetoEvent;
      }
      //_h_cutflow->fill(4., event.weight());
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
     // _h_cutflow->fill(5., event.weight());
      if (js[0].perp()/GeV<100) vetoEvent;
      //_h_cutflow->fill(6., event.weight());
      if (js[0].abseta()>2.5) vetoEvent;
     //_h_cutflow->fill(7., event.weight());
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
   // _h_cutflow->fill(8., event.weight());

    if (sum_jet.pT() < 130) vetoEvent;
   // _h_cutflow->fill(9., event.weight());
    if (sum_jet.pT()/smearedmet.mod() > 1.25) vetoEvent;
   // _h_cutflow->fill(10., event.weight());
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
    if (ht >400 && ht <600)
      {if (AlphaT<0.52) vetoEvent;}
    if (ht>600 && ht<800)
      {if (AlphaT< 0.52) vetoEvent;}
    if(ht > 800) 
      {if (AlphaT < 0) vetoEvent;}
    }
    //_h_cutflow->fill(11., event.weight());
    for(unsigned itjet = 0 ; itjet < js.size(); itjet++) {
      double tmp=deltaPhi(sum_jet-js[itjet].momentum(), js[itjet]);
      if(tmp > drmax){
      drmax=tmp;
      }
    }
    double angle = M_PI-0.5;
    if (drmax > angle) vetoEvent;
    //_h_cutflow->fill(12., event.weight());

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
  if(js.size()==1){
    if(btags==0){
      if(ht >200 && ht< 250)
      _h_cutflow->fill(0., event.weight()); 
      else if(ht >250 && ht< 300)
      _h_cutflow->fill(1., event.weight());
      else if(ht >300 && ht< 350)
      _h_cutflow->fill(2., event.weight());
      else if(ht >350 && ht< 400)
      _h_cutflow->fill(3., event.weight());
      else if(ht >400 && ht< 500)
      _h_cutflow->fill(4., event.weight());
      else if(ht >500 && ht< 600)
      _h_cutflow->fill(5., event.weight());
      else _h_cutflow->fill(6., event.weight());}
    else if(btags==1){
      if(ht >200 && ht< 250)
      _h_cutflow->fill(7., event.weight());
      else if(ht >250 && ht< 300)
      _h_cutflow->fill(8., event.weight());
      else if(ht >300 && ht< 350)
      _h_cutflow->fill(9., event.weight());
      else if(ht >350 && ht< 400)
      _h_cutflow->fill(10., event.weight());
      else if(ht >400 && ht< 500)
      _h_cutflow->fill(11., event.weight());
      else _h_cutflow->fill(12., event.weight());}       
  }
  if(js.size()==2){
    if(btags==0){
      if(ht >200 && ht< 250){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(13., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(14., event.weight());
        else _h_cutflow->fill(15., event.weight());}  
      else if(ht >250 && ht< 300){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(16., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(17., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(18., event.weight());
        else _h_cutflow->fill(19., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(20., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(21., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(22., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(23., event.weight());
        else _h_cutflow->fill(24., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(25., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(26., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(27., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(28., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(29., event.weight());
        else _h_cutflow->fill(30., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(31., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(32., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(33., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(34., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(35., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(36., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(37., event.weight());
        else _h_cutflow->fill(38., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(39., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(40., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(41., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(42., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(43., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(44., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(45., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflow->fill(46., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflow->fill(47., event.weight());
        else _h_cutflow->fill(48., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<225){
        _h_ch1sig->fill(130., event.weight());
        _h_ch1try->fill(130.);
        _h_cutflow->fill(49., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_cutflow->fill(50., event.weight());
        _h_ch1try->fill(225.);
        _h_ch1sig->fill(225., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_cutflow->fill(51., event.weight());
        _h_ch1try->fill(275.);
        _h_ch1sig->fill(275., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_cutflow->fill(52., event.weight());
        _h_ch1try->fill(325.);
        _h_ch1sig->fill(325., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch1sig->fill(375., event.weight());
        _h_ch1try->fill(375.);
        _h_cutflow->fill(53., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch1sig->fill(425., event.weight());
        _h_ch1try->fill(425.);
        _h_cutflow->fill(54., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<525){
        _h_ch1sig->fill(475., event.weight());
        _h_ch1try->fill(475.);
        _h_cutflow->fill(55., event.weight());}
        else if(sum_jet.pT()>525 && sum_jet.pT()<575){
        _h_ch1sig->fill(525., event.weight());
        _h_ch1try->fill(525.);
        _h_cutflow->fill(56., event.weight());}
        else if(sum_jet.pT()>575 && sum_jet.pT()<625){
        _h_ch1sigbefore->fill(0., event.weight());
        _h_ch1try->fill(575.);
        _h_ch1sig->fill(575., event.weight());
        _h_cutflow->fill(57., event.weight());}
        else if(sum_jet.pT()>625 && sum_jet.pT()<675){
        _h_ch1sig->fill(625., event.weight());
        _h_ch1try->fill(625.);
        _h_cutflow->fill(58., event.weight());}
        else if(sum_jet.pT()>675 && sum_jet.pT()<725){
        _h_cutflow->fill(59., event.weight());
        _h_ch1try->fill(675.);
        _h_ch1sig->fill(675., event.weight());}
        else if(sum_jet.pT()>725 && sum_jet.pT()<775){
        _h_cutflow->fill(60., event.weight());
        _h_ch1sigbefore->fill(3., event.weight());
        _h_ch1try->fill(725.);
        _h_ch1sig->fill(725., event.weight());
        }
        else{_h_ch1sig->fill(775., event.weight());
          _h_ch1try->fill(775.);
          _h_cutflow->fill(61., event.weight());}}
      else {
        if(sum_jet.pT()<250){
        _h_ch2sig->fill(130., event.weight());
        _h_cutflow->fill(62., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch2sig->fill(250., event.weight());
        _h_cutflow->fill(63., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch2sig->fill(300., event.weight());
        _h_cutflow->fill(64., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch2sig->fill(350., event.weight());
        _h_cutflow->fill(65., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch2sig->fill(400., event.weight());
        _h_cutflow->fill(66., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch2sig->fill(450., event.weight());
        _h_cutflow->fill(67., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch2sig->fill(500., event.weight());
        _h_cutflow->fill(68., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch2sig->fill(550., event.weight());
        _h_cutflow->fill(69., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_ch2sig->fill(600., event.weight());
        _h_cutflow->fill(70., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_ch2sig->fill(650., event.weight());
        _h_cutflow->fill(71., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        //_h_ch2sig->fill(0., event.weight());
        _h_ch2sig->fill(700., event.weight());
        _h_cutflow->fill(72., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        //_h_ch2sig->fill(1., event.weight());
        _h_ch2sig->fill(750., event.weight());
        _h_cutflow->fill(73., event.weight());}
        else{ //_h_ch2sig->fill(2., event.weight());
        _h_ch2sig->fill(800., event.weight());
        _h_cutflow->fill(74., event.weight());}}
      }

    if(btags==1){
      if(ht >200 && ht< 250){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(75., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(76., event.weight());
        else _h_cutflow->fill(77., event.weight());}  
      else if(ht >250 && ht< 300){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(78., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(79., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(80., event.weight());
        else _h_cutflow->fill(81., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(82., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(83., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(84., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(85., event.weight());
        else _h_cutflow->fill(86., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(87., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(88., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(89., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(90., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(91., event.weight());
        else _h_cutflow->fill(92., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175){
        _h_ch3sig->fill(130., event.weight());
        _h_cutflow->fill(93., event.weight());}
        else if(sum_jet.pT()>175 && sum_jet.pT()<225){
        _h_ch3sig->fill(175., event.weight());
        _h_cutflow->fill(94., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch3sig->fill(225., event.weight());
        _h_cutflow->fill(95., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch3sig->fill(275., event.weight());
        _h_cutflow->fill(96., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch3sig->fill(325., event.weight());
        _h_cutflow->fill(97., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch3sig->fill(375., event.weight());
        _h_cutflow->fill(98., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_cutflow->fill(99., event.weight());
        _h_ch3sig->fill(425., event.weight());}
        else{_h_ch3sig->fill(475., event.weight());
        _h_cutflow->fill(100., event.weight());}}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<225){
        _h_ch4sig->fill(130., event.weight());
        _h_cutflow->fill(101., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch4sig->fill(225., event.weight());
        _h_cutflow->fill(102., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch4sig->fill(275., event.weight());
        _h_cutflow->fill(103., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch4sig->fill(325., event.weight());
        _h_cutflow->fill(104., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch4sig->fill(375., event.weight());
        _h_cutflow->fill(105., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch4sig->fill(425., event.weight());
        _h_cutflow->fill(106., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<525){
        _h_ch4sig->fill(475., event.weight());
        _h_cutflow->fill(107., event.weight());}
        else if(sum_jet.pT()>525 && sum_jet.pT()<575){
        _h_ch4sig->fill(525., event.weight());
        _h_cutflow->fill(108., event.weight());}
        else{_h_ch4sig->fill(575., event.weight());
         _h_cutflow->fill(109., event.weight());}}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<250){
        _h_ch5sig->fill(130., event.weight());
        _h_cutflow->fill(110., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch5sig->fill(250., event.weight());
        _h_cutflow->fill(111., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch5sig->fill(300., event.weight());
        _h_cutflow->fill(112., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch5sig->fill(350., event.weight());
        _h_cutflow->fill(113., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch5sig->fill(400., event.weight());
        _h_cutflow->fill(114., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch5sig->fill(450., event.weight());
        _h_cutflow->fill(115., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch5sig->fill(500., event.weight());
        _h_cutflow->fill(116., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch5sig->fill(550., event.weight());
        _h_cutflow->fill(117., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_ch5sig->fill(600., event.weight());
        _h_cutflow->fill(118., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_cutflow->fill(119., event.weight());
        _h_ch5sig->fill(650., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_cutflow->fill(120., event.weight());
        _h_ch5sig->fill(700., event.weight());}
        else{ _h_cutflow->fill(121., event.weight());
        _h_ch5sig->fill(750., event.weight());}}
      else {
        if(sum_jet.pT()<250){
        _h_ch6sig->fill(130., event.weight());
        _h_cutflow->fill(122., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch6sig->fill(250., event.weight());
        _h_cutflow->fill(123., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch6sig->fill(300., event.weight());
        _h_cutflow->fill(124., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch6sig->fill(350., event.weight());
        _h_cutflow->fill(125., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch6sig->fill(400., event.weight());
        _h_cutflow->fill(126., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch6sig->fill(450., event.weight());
        _h_cutflow->fill(127., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch6sig->fill(500., event.weight());
        _h_cutflow->fill(128., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_cutflow->fill(129., event.weight());
        _h_ch6sig->fill(550., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_cutflow->fill(130., event.weight());
        _h_ch6sig->fill(600., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_cutflow->fill(131., event.weight());
        _h_ch6sig->fill(650., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_cutflow->fill(132., event.weight());
        _h_ch6sig->fill(700., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_cutflow->fill(133., event.weight());
        _h_ch6sig->fill(750., event.weight());}
        else{ _h_cutflow->fill(134., event.weight());
        _h_ch6sig->fill(800., event.weight());}}
      }

    if(btags==2){
      if(ht >200 && ht< 250){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(135., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(136., event.weight());
        else _h_cutflow->fill(137., event.weight());}  
      else if(ht >250 && ht< 300){
        if(sum_jet.pT()<225)
        _h_cutflow->fill(138., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(139., event.weight());
        else _h_cutflow->fill(140., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<250)
        _h_cutflow->fill(141., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflow->fill(142., event.weight());
        else _h_cutflow->fill(143., event.weight());}   
      else if(ht >350 && ht< 400){
        _h_cutflow->fill(144., event.weight());}
      else if(ht >400 && ht< 500){
        _h_ch7sig->fill(130., event.weight());
        _h_cutflow->fill(145., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<425)
        _h_cutflow->fill(146., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<550)
        _h_cutflow->fill(147., event.weight());
        else _h_cutflow->fill(148., event.weight());}
      else {_h_ch8sig->fill(130., event.weight());
        _h_cutflow->fill(149., event.weight());}
      }

  }

  if(js.size()==3){
    if(btags==0){
      if(ht >200 && ht< 250){
        _h_cutflow->fill(150., event.weight());}  
      else if(ht >250 && ht< 300){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(151., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(152., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(153., event.weight());
        else _h_cutflow->fill(154., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(155., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(156., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(157., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(158., event.weight());
        else _h_cutflow->fill(159., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(160., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(161., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(162., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(163., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(164., event.weight());
        else _h_cutflow->fill(165., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(166., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(167., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(168., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(169., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(170., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(171., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(172., event.weight());
        else _h_cutflow->fill(173., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(174., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(175., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(176., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(177., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(178., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(179., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(180., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflow->fill(181., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflow->fill(182., event.weight());
        else _h_cutflow->fill(183., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<225)
        _h_cutflow->fill(184., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(185., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(186., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(187., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(188., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(189., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflow->fill(190., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflow->fill(191., event.weight());
        else if(sum_jet.pT()>575 && sum_jet.pT()<625)
        _h_cutflow->fill(192., event.weight());
        else if(sum_jet.pT()>625 && sum_jet.pT()<675)
        _h_cutflow->fill(193., event.weight());
        else if(sum_jet.pT()>675 && sum_jet.pT()<725)
        _h_cutflow->fill(194., event.weight());
        else if(sum_jet.pT()>725 && sum_jet.pT()<775)
        _h_cutflow->fill(195., event.weight());
        else _h_cutflow->fill(196., event.weight());}
      else {
        if(sum_jet.pT()<200){
        _h_ch9sig->fill(130., event.weight());
        _h_cutflow->fill(197., event.weight());}
        else if(sum_jet.pT()>200 && sum_jet.pT()<250){
        _h_ch9sig->fill(200., event.weight());
        _h_cutflow->fill(198., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch9sig->fill(250., event.weight());
        _h_cutflow->fill(199., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch9sig->fill(300., event.weight());
        _h_cutflow->fill(200., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch9sig->fill(350., event.weight());
        _h_cutflow->fill(201., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch9sig->fill(400., event.weight());
        _h_cutflow->fill(202., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch9sig->fill(450., event.weight());
        _h_cutflow->fill(203., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch9sig->fill(500., event.weight());
        _h_cutflow->fill(204., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch9sig->fill(550., event.weight());
        _h_cutflow->fill(205., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_ch9sig->fill(600., event.weight());
        _h_cutflow->fill(206., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_cutflow->fill(207., event.weight());
        _h_ch9sig->fill(650., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_ch9sig->fill(700., event.weight());
        _h_cutflow->fill(208., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_ch9sig->fill(750., event.weight());
        _h_cutflow->fill(209., event.weight());}
        else{ _h_cutflow->fill(210., event.weight());
        _h_ch9sig->fill(800., event.weight());}}
      }

    if(btags==1){
      if(ht >200 && ht< 250){
        _h_cutflow->fill(211., event.weight());}  
      else if(ht >250 && ht< 300){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(212., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(213., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(214., event.weight());
        else _h_cutflow->fill(215., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(216., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(217., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(218., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(219., event.weight());
        else _h_cutflow->fill(220., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(221., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(222., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(223., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(224., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(225., event.weight());
        else _h_cutflow->fill(226., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175)
        _h_cutflow->fill(227., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflow->fill(228., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflow->fill(229., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflow->fill(230., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflow->fill(231., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflow->fill(232., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflow->fill(233., event.weight());
        else _h_cutflow->fill(234., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<200){
        _h_ch10sig->fill(130., event.weight());
        _h_cutflow->fill(235., event.weight());}
        else if(sum_jet.pT()>200 && sum_jet.pT()<250){
        _h_ch10sig->fill(200., event.weight());
        _h_cutflow->fill(236., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch10sig->fill(250., event.weight());
        _h_cutflow->fill(237., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch10sig->fill(300., event.weight());
        _h_cutflow->fill(238., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch10sig->fill(350., event.weight());
        _h_cutflow->fill(239., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch10sig->fill(400., event.weight());
        _h_cutflow->fill(240., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch10sig->fill(450., event.weight());
        _h_cutflow->fill(241., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_cutflow->fill(242., event.weight());
        _h_ch10sig->fill(500., event.weight());}
        else{_h_ch10sig->fill(550., event.weight());
         _h_cutflow->fill(243., event.weight());}}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<225){
        _h_ch11sig->fill(130., event.weight());
        _h_cutflow->fill(244., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch11sig->fill(225., event.weight());
        _h_cutflow->fill(245., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch11sig->fill(275., event.weight());
        _h_cutflow->fill(246., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch11sig->fill(325., event.weight());
        _h_cutflow->fill(247., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch11sig->fill(375., event.weight());
        _h_cutflow->fill(248., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch11sig->fill(425., event.weight());
        _h_cutflow->fill(249., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<525){
        _h_ch11sig->fill(475., event.weight());
        _h_cutflow->fill(250., event.weight());}
        else if(sum_jet.pT()>525 && sum_jet.pT()<575){
        _h_ch11sig->fill(525., event.weight());
        _h_cutflow->fill(251., event.weight());}
        else if(sum_jet.pT()>575 && sum_jet.pT()<625){
        _h_cutflow->fill(252., event.weight());
        _h_ch11sig->fill(575., event.weight());}
        else if(sum_jet.pT()>625 && sum_jet.pT()<675){
        _h_cutflow->fill(253., event.weight());
        _h_ch11sig->fill(625., event.weight());}
        else if(sum_jet.pT()>675 && sum_jet.pT()<725){
        _h_cutflow->fill(254., event.weight());
        _h_ch11sig->fill(675., event.weight());}
        else{ _h_cutflow->fill(255., event.weight());
        _h_ch11sig->fill(725., event.weight());}}
      else {
        if(sum_jet.pT()<250){
        _h_ch12sig->fill(130., event.weight());
        _h_cutflow->fill(256., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch12sig->fill(250., event.weight());
        _h_cutflow->fill(257., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch12sig->fill(300., event.weight());
        _h_cutflow->fill(258., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch12sig->fill(350., event.weight());
        _h_cutflow->fill(259., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch12sig->fill(400., event.weight());
        _h_cutflow->fill(260., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch12sig->fill(450., event.weight());
        _h_cutflow->fill(261., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch12sig->fill(500., event.weight());
        _h_cutflow->fill(262., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_cutflow->fill(263., event.weight());
        _h_ch12sig->fill(550., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_cutflow->fill(264., event.weight());
        _h_ch12sig->fill(600., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_cutflow->fill(265., event.weight());
        _h_ch12sig->fill(650., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_cutflow->fill(266., event.weight());
        _h_ch12sig->fill(700., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_cutflow->fill(267., event.weight());
        _h_ch12sig->fill(750., event.weight());}
        else{ _h_cutflow->fill(268., event.weight());
        _h_ch12sig->fill(800., event.weight());}}
      }

    if(btags==2){
      if(ht >250 && ht< 300){
        if(sum_jet.pT()<200)
        _h_cutflow->fill(269., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflow->fill(270., event.weight());
        else _h_cutflow->fill(271., event.weight());}
      else if(ht >300 && ht< 350){
        if(sum_jet.pT()<150)
        _h_cutflow->fill(272., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflow->fill(273., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflow->fill(274., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflow->fill(275., event.weight());
        else _h_cutflow->fill(276., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<150)
        _h_cutflow->fill(277., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflow->fill(278., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflow->fill(279., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflow->fill(280., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflow->fill(281., event.weight());
        else _h_cutflow->fill(282., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175){
        _h_ch13sig->fill(130., event.weight());
        _h_cutflow->fill(283., event.weight());}
        else if(sum_jet.pT()>175 && sum_jet.pT()<225){
        _h_ch13sig->fill(175., event.weight());
        _h_cutflow->fill(284., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch13sig->fill(225., event.weight());
        _h_cutflow->fill(285., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch13sig->fill(275., event.weight());
        _h_cutflow->fill(286., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch13sig->fill(325., event.weight());
        _h_cutflow->fill(287., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch13sig->fill(375., event.weight());
        _h_cutflow->fill(288., event.weight());}
        else{ _h_cutflow->fill(289., event.weight());
        _h_ch13sig->fill(425., event.weight());}}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<225){
        _h_ch14sig->fill(130., event.weight());
        _h_cutflow->fill(290., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch14sig->fill(225., event.weight());
        _h_cutflow->fill(291., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch14sig->fill(275., event.weight());
        _h_cutflow->fill(292., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch14sig->fill(325., event.weight());
        _h_cutflow->fill(293., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch14sig->fill(375., event.weight());
        _h_cutflow->fill(294., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch14sig->fill(425., event.weight());
        _h_cutflow->fill(295., event.weight());}
        else{ _h_cutflow->fill(296., event.weight());
        _h_ch14sig->fill(475., event.weight());}}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<250){
        _h_ch15sig->fill(130., event.weight());
        _h_cutflow->fill(297., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch15sig->fill(250., event.weight());
        _h_cutflow->fill(298., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch15sig->fill(300., event.weight());
        _h_cutflow->fill(299., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch15sig->fill(350., event.weight());
        _h_cutflow->fill(300., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch15sig->fill(400., event.weight());
        _h_cutflow->fill(301., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch15sig->fill(450., event.weight());
        _h_cutflow->fill(302., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<575){
        _h_ch15sig->fill(500., event.weight());
        _h_cutflow->fill(303., event.weight());}
        else if(sum_jet.pT()>575 && sum_jet.pT()<625){
        _h_ch15sig->fill(575., event.weight());
        _h_cutflow->fill(304., event.weight());}
        else{ _h_cutflow->fill(305., event.weight());
        _h_ch15sig->fill(625., event.weight());}}
      else {
        if(sum_jet.pT()<325){
        _h_ch16sig->fill(130., event.weight());
        _h_cutflow->fill(306., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch16sig->fill(325., event.weight());
        _h_cutflow->fill(307., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch16sig->fill(375., event.weight());
        _h_cutflow->fill(308., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch16sig->fill(425., event.weight());
        _h_cutflow->fill(309., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<600){
        _h_ch16sig->fill(475., event.weight());
        _h_cutflow->fill(310., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<700){
        _h_ch16sig->fill(600., event.weight());
        _h_cutflow->fill(311., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<800){
        _h_ch16sig->fill(700., event.weight());
        _h_cutflow->fill(312., event.weight());}
        else{ _h_cutflow->fill(313., event.weight());
        _h_ch16sig->fill(800., event.weight());}}
      }

    if(btags==3){
      if(ht >250 && ht< 300){
        _h_cutflow->fill(314., event.weight());}  
      else if(ht >400){
        _h_cutflow->fill(315., event.weight());}
      }

  }

  if(js.size()==4){
    if(btags==0){
      if(ht >300 && ht< 350){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(0., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(1., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(2., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(3., event.weight());
        else _h_cutflownjet4->fill(4., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(5., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(6., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(7., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(8., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(9., event.weight());
        else _h_cutflownjet4->fill(10., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(11., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(12., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(13., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(14., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(15., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(16., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(17., event.weight());
        else _h_cutflownjet4->fill(18., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(19., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(20., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(21., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(22., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(23., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(24., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(25., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(26., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflownjet4->fill(27., event.weight());
        else _h_cutflownjet4->fill(28., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(29., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(30., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(31., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(32., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(33., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(34., event.weight());
        else if(sum_jet.pT()>450 && sum_jet.pT()<500)
        _h_cutflownjet4->fill(35., event.weight());
        else if(sum_jet.pT()>500 && sum_jet.pT()<550)
        _h_cutflownjet4->fill(36., event.weight());
        else if(sum_jet.pT()>550 && sum_jet.pT()<600)
        _h_cutflownjet4->fill(37., event.weight());
        else if(sum_jet.pT()>600 && sum_jet.pT()<650)
        _h_cutflownjet4->fill(38., event.weight());
        else if(sum_jet.pT()>650 && sum_jet.pT()<700)
        _h_cutflownjet4->fill(39., event.weight());
        else if(sum_jet.pT()>700 && sum_jet.pT()<750)
        _h_cutflownjet4->fill(40., event.weight());
        else _h_cutflownjet4->fill(41., event.weight());}
      else {
        if(sum_jet.pT()<200){
        _h_ch17sig->fill(130., event.weight());
        _h_cutflownjet4->fill(42., event.weight());}
        else if(sum_jet.pT()>200 && sum_jet.pT()<250){
        _h_ch17sig->fill(200., event.weight()); 
        _h_cutflownjet4->fill(43., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch17sig->fill(250., event.weight());
        _h_cutflownjet4->fill(44., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch17sig->fill(300., event.weight());
        _h_cutflownjet4->fill(45., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch17sig->fill(350., event.weight());
        _h_cutflownjet4->fill(46., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch17sig->fill(400., event.weight());
        _h_cutflownjet4->fill(47., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch17sig->fill(450., event.weight());
        _h_cutflownjet4->fill(48., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch17sig->fill(500., event.weight());
        _h_cutflownjet4->fill(49., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch17sig->fill(550., event.weight());
        _h_cutflownjet4->fill(50., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_ch17sig->fill(600., event.weight());
        _h_cutflownjet4->fill(51., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_ch17sig->fill(650., event.weight());
        _h_cutflownjet4->fill(52., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_ch17sig->fill(700., event.weight());
        _h_cutflownjet4->fill(53., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_ch17sig->fill(750., event.weight());
        _h_cutflownjet4->fill(54., event.weight());}
        else{ _h_cutflownjet4->fill(55., event.weight());
        _h_ch17sig->fill(800., event.weight());}}
      }

    if(btags==1){
      if(ht >300 && ht< 350){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(56., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(57., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(58., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(59., event.weight());
        else _h_cutflownjet4->fill(60., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(61., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(62., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(63., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(64., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(65., event.weight());
        else _h_cutflownjet4->fill(66., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(67., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(68., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(69., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(70., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(71., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(72., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(73., event.weight());
        else _h_cutflownjet4->fill(74., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(75., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(76., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(77., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(78., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(79., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(80., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(81., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(82., event.weight());
        else _h_cutflownjet4->fill(83., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(84., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(85., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(86., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(87., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(88., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(89., event.weight());
        else if(sum_jet.pT()>450 && sum_jet.pT()<500)
        _h_cutflownjet4->fill(90., event.weight());
        else if(sum_jet.pT()>500 && sum_jet.pT()<550)
        _h_cutflownjet4->fill(91., event.weight());
        else if(sum_jet.pT()>550 && sum_jet.pT()<600)
        _h_cutflownjet4->fill(92., event.weight());
        else if(sum_jet.pT()>600 && sum_jet.pT()<650)
        _h_cutflownjet4->fill(93., event.weight());
        else if(sum_jet.pT()>650 && sum_jet.pT()<700)
        _h_cutflownjet4->fill(94., event.weight());
        else _h_cutflownjet4->fill(95., event.weight());}
      else {
        if(sum_jet.pT()<200){
        _h_ch18sig->fill(130., event.weight());
        _h_cutflownjet4->fill(96., event.weight());}
        else if(sum_jet.pT()>200 && sum_jet.pT()<250){
        _h_ch18sig->fill(200., event.weight());
        _h_cutflownjet4->fill(97., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch18sig->fill(250., event.weight());
        _h_cutflownjet4->fill(98., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch18sig->fill(300., event.weight());
        _h_cutflownjet4->fill(99., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch18sig->fill(350., event.weight());
        _h_cutflownjet4->fill(100., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch18sig->fill(400., event.weight());
        _h_cutflownjet4->fill(101., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch18sig->fill(450., event.weight());
        _h_cutflownjet4->fill(102., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch18sig->fill(500., event.weight());
        _h_cutflownjet4->fill(103., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch18sig->fill(550., event.weight());
        _h_cutflownjet4->fill(104., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_cutflownjet4->fill(105., event.weight());
        _h_ch18sig->fill(600., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_ch18sig->fill(650., event.weight());
        _h_cutflownjet4->fill(106., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_cutflownjet4->fill(107., event.weight());
        _h_ch18sig->fill(700., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_cutflownjet4->fill(108., event.weight());
        _h_ch18sig->fill(750., event.weight());}
        else{ _h_cutflownjet4->fill(109., event.weight());
        _h_ch18sig->fill(800., event.weight());}}
      }

    if(btags==2){
      if(ht >300 && ht< 350){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(110., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(111., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(112., event.weight());
        else _h_cutflownjet4->fill(113., event.weight());}   
      else if(ht >350 && ht< 400){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(114., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(115., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(116., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(117., event.weight());
        else _h_cutflownjet4->fill(118., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(119., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(120., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(121., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(122., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(123., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(124., event.weight());
        else _h_cutflownjet4->fill(125., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(126., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(127., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(128., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(129., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(130., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(131., event.weight());
        else _h_cutflownjet4->fill(132., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<275)
        _h_cutflownjet4->fill(133., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(134., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(135., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(136., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(137., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(138., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<600)
        _h_cutflownjet4->fill(139., event.weight());
        else if(sum_jet.pT()>600 && sum_jet.pT()<650)
        _h_cutflownjet4->fill(140., event.weight());
        else _h_cutflownjet4->fill(141., event.weight());}
      else {
        if(sum_jet.pT()<225){
        _h_ch19sig->fill(130., event.weight());
        _h_cutflownjet4->fill(142., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch19sig->fill(225., event.weight());
        _h_cutflownjet4->fill(143., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch19sig->fill(275., event.weight());
        _h_cutflownjet4->fill(144., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch19sig->fill(325., event.weight());
        _h_cutflownjet4->fill(145., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch19sig->fill(375., event.weight());
        _h_cutflownjet4->fill(146., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch19sig->fill(425., event.weight());
        _h_cutflownjet4->fill(147., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<525){
        _h_ch19sig->fill(475., event.weight());
        _h_cutflownjet4->fill(148., event.weight());}
        else if(sum_jet.pT()>525 && sum_jet.pT()<600){
        _h_ch19sig->fill(525., event.weight()); 
        _h_cutflownjet4->fill(149., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<700){
        _h_cutflownjet4->fill(150., event.weight());
        _h_ch19sig->fill(600., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<800){
        _h_ch19sig->fill(700., event.weight());
        _h_cutflownjet4->fill(151., event.weight());}
        else{_h_ch19sig->fill(800., event.weight());
         _h_cutflownjet4->fill(152., event.weight());}}
      }

    if(btags==3){
      if(ht >300 && ht< 350){
        _h_cutflownjet4->fill(153., event.weight());}  
      else if(ht >350 && ht< 400){
        _h_cutflownjet4->fill(154., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(155., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(156., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(157., event.weight());        
        else _h_cutflownjet4->fill(158., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<375)
        _h_cutflownjet4->fill(159., event.weight());
        else _h_cutflownjet4->fill(160., event.weight());}
      else if(ht >600 && ht< 800){
        _h_cutflownjet4->fill(161., event.weight());}
      else {_h_ch20sig->fill(130., event.weight());
        _h_cutflownjet4->fill(162., event.weight());}
      }
  }

  if(js.size()==5){
    if(btags==0){
      if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(163., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(164., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(165., event.weight());
        else _h_cutflownjet4->fill(166., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(167., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(168., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(169., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(170., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(171., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(172., event.weight());
        else _h_cutflownjet4->fill(173., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(174., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(175., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(176., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(177., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(178., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(179., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(180., event.weight());
        else if(sum_jet.pT()>450 && sum_jet.pT()<500)
        _h_cutflownjet4->fill(181., event.weight());
        else _h_cutflownjet4->fill(182., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<225)
        _h_cutflownjet4->fill(183., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(184., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(185., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(186., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(187., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(188., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(189., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflownjet4->fill(190., event.weight());
        else if(sum_jet.pT()>575 && sum_jet.pT()<625)
        _h_cutflownjet4->fill(191., event.weight());
        else if(sum_jet.pT()>625 && sum_jet.pT()<675)
        _h_cutflownjet4->fill(192., event.weight());
        else if(sum_jet.pT()>675 && sum_jet.pT()<725)
        _h_cutflownjet4->fill(193., event.weight());
        else _h_cutflownjet4->fill(194., event.weight());}
      else {
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(195., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(196., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(197., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(198., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(199., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(200., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(201., event.weight());
        else if(sum_jet.pT()>450 && sum_jet.pT()<500)
        _h_cutflownjet4->fill(202., event.weight());
        else if(sum_jet.pT()>500 && sum_jet.pT()<550)
        _h_cutflownjet4->fill(203., event.weight());
        else if(sum_jet.pT()>550 && sum_jet.pT()<600)
        _h_cutflownjet4->fill(204., event.weight());
        else if(sum_jet.pT()>600 && sum_jet.pT()<650)
        _h_cutflownjet4->fill(205., event.weight());
        else if(sum_jet.pT()>650 && sum_jet.pT()<700)
        _h_cutflownjet4->fill(206., event.weight());
        else if(sum_jet.pT()>700 && sum_jet.pT()<750)
        _h_cutflownjet4->fill(207., event.weight());
        else if(sum_jet.pT()>750 && sum_jet.pT()<800)
        _h_cutflownjet4->fill(208., event.weight());
        else _h_cutflownjet4->fill(209., event.weight());}
      }

    if(btags==1){
      if(ht >350 && ht< 400){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(210., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(211., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(212., event.weight());
        else _h_cutflownjet4->fill(213., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<150)
        _h_cutflownjet4->fill(214., event.weight());
        else if(sum_jet.pT()>150 && sum_jet.pT()<200)
        _h_cutflownjet4->fill(215., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(216., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(217., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(218., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(219., event.weight());
        else _h_cutflownjet4->fill(220., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(221., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(222., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(223., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(224., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(225., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(226., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(227., event.weight());
        else _h_cutflownjet4->fill(228., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(229., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(230., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(231., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(232., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(233., event.weight());
        else if(sum_jet.pT()>400 && sum_jet.pT()<450)
        _h_cutflownjet4->fill(234., event.weight());
        else if(sum_jet.pT()>450 && sum_jet.pT()<500)
        _h_cutflownjet4->fill(235., event.weight());
        else if(sum_jet.pT()>500 && sum_jet.pT()<550)
        _h_cutflownjet4->fill(236., event.weight());
        else if(sum_jet.pT()>550 && sum_jet.pT()<600)
        _h_cutflownjet4->fill(237., event.weight());
        else if(sum_jet.pT()>600 && sum_jet.pT()<650)
        _h_cutflownjet4->fill(238., event.weight());
        else _h_cutflownjet4->fill(239., event.weight());}
      else {
        if(sum_jet.pT()<150){
        _h_ch21sig->fill(130., event.weight());
        _h_cutflownjet4->fill(240., event.weight());}
        else if(sum_jet.pT()>150 && sum_jet.pT()<200){
        _h_ch21sig->fill(150., event.weight());
        _h_cutflownjet4->fill(241., event.weight());}
        else if(sum_jet.pT()>200 && sum_jet.pT()<250){
        _h_ch21sig->fill(200., event.weight());
        _h_cutflownjet4->fill(242., event.weight());}
        else if(sum_jet.pT()>250 && sum_jet.pT()<300){
        _h_ch21sig->fill(250., event.weight());
        _h_cutflownjet4->fill(243., event.weight());}
        else if(sum_jet.pT()>300 && sum_jet.pT()<350){
        _h_ch21sig->fill(300., event.weight());
        _h_cutflownjet4->fill(244., event.weight());}
        else if(sum_jet.pT()>350 && sum_jet.pT()<400){
        _h_ch21sig->fill(350., event.weight());
        _h_cutflownjet4->fill(245., event.weight());}
        else if(sum_jet.pT()>400 && sum_jet.pT()<450){
        _h_ch21sig->fill(400., event.weight());
        _h_cutflownjet4->fill(246., event.weight());}
        else if(sum_jet.pT()>450 && sum_jet.pT()<500){
        _h_ch21sig->fill(450., event.weight());
        _h_cutflownjet4->fill(247., event.weight());}
        else if(sum_jet.pT()>500 && sum_jet.pT()<550){
        _h_ch21sig->fill(500., event.weight());
        _h_cutflownjet4->fill(248., event.weight());}
        else if(sum_jet.pT()>550 && sum_jet.pT()<600){
        _h_ch21sig->fill(550., event.weight());
        _h_cutflownjet4->fill(249., event.weight());}
        else if(sum_jet.pT()>600 && sum_jet.pT()<650){
        _h_ch21sig->fill(600., event.weight());
        _h_cutflownjet4->fill(250., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<700){
        _h_ch21sig->fill(650., event.weight());
        _h_cutflownjet4->fill(251., event.weight());}
        else if(sum_jet.pT()>700 && sum_jet.pT()<750){
        _h_ch21sig->fill(700., event.weight());
        _h_cutflownjet4->fill(252., event.weight());}
        else if(sum_jet.pT()>750 && sum_jet.pT()<800){
        _h_ch21sig->fill(750., event.weight());
        _h_cutflownjet4->fill(253., event.weight());}
        else{ _h_cutflownjet4->fill(254., event.weight());
        _h_ch21sig->fill(800., event.weight());}}
      }

    if(btags==2){
      if(ht >350 && ht< 400){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(255., event.weight());
        else _h_cutflownjet4->fill(256., event.weight());}
      else if(ht >400 && ht< 500){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(257., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(258., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(259., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(260., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(261., event.weight());
        else _h_cutflownjet4->fill(262., event.weight());}
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(263., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<250)
        _h_cutflownjet4->fill(264., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(265., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(266., event.weight());
        else if(sum_jet.pT()>350 && sum_jet.pT()<400)
        _h_cutflownjet4->fill(267., event.weight());
        else _h_cutflownjet4->fill(268., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<225)
        _h_cutflownjet4->fill(269., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(270., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(271., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(272., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(273., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<475)
        _h_cutflownjet4->fill(274., event.weight());
        else if(sum_jet.pT()>475 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(275., event.weight());
        else _h_cutflownjet4->fill(276., event.weight());}
      else {
        if(sum_jet.pT()<175){
        _h_ch22sig->fill(130., event.weight());
        _h_cutflownjet4->fill(277., event.weight());}
        else if(sum_jet.pT()>175 && sum_jet.pT()<225){
        _h_ch22sig->fill(175., event.weight());
        _h_cutflownjet4->fill(278., event.weight());}
        else if(sum_jet.pT()>225 && sum_jet.pT()<275){
        _h_ch22sig->fill(225., event.weight());
        _h_cutflownjet4->fill(279., event.weight());}
        else if(sum_jet.pT()>275 && sum_jet.pT()<325){
        _h_ch22sig->fill(275., event.weight());
        _h_cutflownjet4->fill(280., event.weight());}
        else if(sum_jet.pT()>325 && sum_jet.pT()<375){
        _h_ch22sig->fill(325., event.weight());
        _h_cutflownjet4->fill(281., event.weight());}
        else if(sum_jet.pT()>375 && sum_jet.pT()<425){
        _h_ch22sig->fill(375., event.weight());
        _h_cutflownjet4->fill(282., event.weight());}
        else if(sum_jet.pT()>425 && sum_jet.pT()<475){
        _h_ch22sig->fill(425., event.weight());
        _h_cutflownjet4->fill(283., event.weight());}
        else if(sum_jet.pT()>475 && sum_jet.pT()<525){
        _h_ch22sig->fill(475., event.weight());
        _h_cutflownjet4->fill(284., event.weight());}
        else if(sum_jet.pT()>525 && sum_jet.pT()<575){
        _h_ch22sig->fill(525., event.weight());
        _h_cutflownjet4->fill(285., event.weight());}
        else if(sum_jet.pT()>575 && sum_jet.pT()<650){
        _h_ch22sig->fill(575., event.weight());
        _h_cutflownjet4->fill(286., event.weight());}
        else if(sum_jet.pT()>650 && sum_jet.pT()<800){
        _h_ch22sig->fill(650., event.weight());
        _h_cutflownjet4->fill(287., event.weight());}
        else{ _h_cutflownjet4->fill(288., event.weight());
        _h_ch22sig->fill(800., event.weight());}}
      }

    if(btags==3){
      if(ht >400 && ht< 500){
        _h_cutflownjet4->fill(289., event.weight());}  
      else if(ht >500 && ht< 600){
        if(sum_jet.pT()<175)
        _h_cutflownjet4->fill(290., event.weight());
        else if(sum_jet.pT()>175 && sum_jet.pT()<225)
        _h_cutflownjet4->fill(291., event.weight());
        else if(sum_jet.pT()>225 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(292., event.weight());        
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(293., event.weight());
        else _h_cutflownjet4->fill(294., event.weight());}
      else if(ht >600 && ht< 800){
        if(sum_jet.pT()<250)
        _h_cutflownjet4->fill(295., event.weight());
        else if(sum_jet.pT()>250 && sum_jet.pT()<300)
        _h_cutflownjet4->fill(296., event.weight());
        else if(sum_jet.pT()>300 && sum_jet.pT()<350)
        _h_cutflownjet4->fill(297., event.weight());
        else _h_cutflownjet4->fill(298., event.weight());}
      else {
        if(sum_jet.pT()<200)
        _h_cutflownjet4->fill(299., event.weight());
        else if(sum_jet.pT()>200 && sum_jet.pT()<275)
        _h_cutflownjet4->fill(300., event.weight());
        else if(sum_jet.pT()>275 && sum_jet.pT()<325)
        _h_cutflownjet4->fill(301., event.weight());
        else if(sum_jet.pT()>325 && sum_jet.pT()<375)
        _h_cutflownjet4->fill(302., event.weight());
        else if(sum_jet.pT()>375 && sum_jet.pT()<425)
        _h_cutflownjet4->fill(303., event.weight());
        else if(sum_jet.pT()>425 && sum_jet.pT()<525)
        _h_cutflownjet4->fill(304., event.weight());
        else if(sum_jet.pT()>525 && sum_jet.pT()<575)
        _h_cutflownjet4->fill(305., event.weight());
        else _h_cutflownjet4->fill(306., event.weight());}
      }
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
      
      double lumi = 2300.; 
      double norm=lumi*crossSection()/picobarn / sumOfWeights();
      double normtry=lumi*crossSection()/picobarn;
      
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
      _h_cutflownjet4->scaleW(norm);
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
      _h_ch1sigbefore->scaleW(norm);
      _h_ch1try->scaleW(normtry);
      _h_ch1sig->scaleW(norm);
      _h_ch2sig->scaleW(norm);
      _h_ch3sig->scaleW(norm);
      _h_ch4sig->scaleW(norm);
      _h_ch5sig->scaleW(norm);
      _h_ch6sig->scaleW(norm);
      _h_ch7sig->scaleW(norm);
      _h_ch8sig->scaleW(norm);
      _h_ch9sig->scaleW(norm);
      _h_ch10sig->scaleW(norm);
      _h_ch11sig->scaleW(norm);
      _h_ch12sig->scaleW(norm);
      _h_ch13sig->scaleW(norm);
      _h_ch14sig->scaleW(norm);
      _h_ch15sig->scaleW(norm);
      _h_ch16sig->scaleW(norm);
      _h_ch17sig->scaleW(norm);
      _h_ch18sig->scaleW(norm);
      _h_ch19sig->scaleW(norm);
      _h_ch20sig->scaleW(norm);
      _h_ch21sig->scaleW(norm);
      _h_ch22sig->scaleW(norm);



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
    Histo1DPtr   _h_variablebinwideth; 
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
    Histo1DPtr   _h_cutflownjet4;
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
    Histo1DPtr   _h_ch1sigbefore;
    Histo1DPtr   _h_ch1sig;
    Histo1DPtr   _h_ch1try;
    Histo1DPtr   _h_ch1data;
    Histo1DPtr   _h_ch2data;
    Histo1DPtr   _h_ch3data;
    Histo1DPtr   _h_ch4data;
    Histo1DPtr   _h_ch5data;
    Histo1DPtr   _h_ch6data;
    Histo1DPtr   _h_ch7data;
    Histo1DPtr   _h_ch8data;
    Histo1DPtr   _h_ch9data;
    Histo1DPtr   _h_ch10data;
    Histo1DPtr   _h_ch11data;
    Histo1DPtr   _h_ch12data;
    Histo1DPtr   _h_ch13data;
    Histo1DPtr   _h_ch14data;
    Histo1DPtr   _h_ch15data;
    Histo1DPtr   _h_ch16data;
    Histo1DPtr   _h_ch17data;
    Histo1DPtr   _h_ch18data;
    Histo1DPtr   _h_ch19data;
    Histo1DPtr   _h_ch20data;
    Histo1DPtr   _h_ch21data;
    Histo1DPtr   _h_ch22data;
    Histo1DPtr   _h_ch1bkg;
    Histo1DPtr   _h_ch2sig;
    Histo1DPtr   _h_ch2bkg;
    Histo1DPtr   _h_ch3sig;
    Histo1DPtr   _h_ch3bkg;
    Histo1DPtr   _h_ch4sig;
    Histo1DPtr   _h_ch4bkg;
    Histo1DPtr   _h_ch5sig;
    Histo1DPtr   _h_ch5bkg;
    Histo1DPtr   _h_ch6sig;
    Histo1DPtr   _h_ch6bkg;
    Histo1DPtr   _h_ch7sig;
    Histo1DPtr   _h_ch7bkg;
    Histo1DPtr   _h_ch8sig;
    Histo1DPtr   _h_ch8bkg;
    Histo1DPtr   _h_ch9sig;
    Histo1DPtr   _h_ch9bkg;
    Histo1DPtr   _h_ch10sig;
    Histo1DPtr   _h_ch10bkg;
    Histo1DPtr   _h_ch11sig;
    Histo1DPtr   _h_ch11bkg;
    Histo1DPtr   _h_ch12sig;
    Histo1DPtr   _h_ch12bkg;
    Histo1DPtr   _h_ch13sig;
    Histo1DPtr   _h_ch13bkg;
    Histo1DPtr   _h_ch14sig;
    Histo1DPtr   _h_ch14bkg;
    Histo1DPtr   _h_ch15sig;
    Histo1DPtr   _h_ch15bkg;
    Histo1DPtr   _h_ch16sig;
    Histo1DPtr   _h_ch16bkg;
    Histo1DPtr   _h_ch17sig;
    Histo1DPtr   _h_ch17bkg;
    Histo1DPtr   _h_ch18sig;
    Histo1DPtr   _h_ch18bkg;
    Histo1DPtr   _h_ch19sig;
    Histo1DPtr   _h_ch19bkg;
    Histo1DPtr   _h_ch20sig;
    Histo1DPtr   _h_ch20bkg;
    Histo1DPtr   _h_ch21sig;
    Histo1DPtr   _h_ch21bkg;
    Histo1DPtr   _h_ch22sig;
    Histo1DPtr   _h_ch22bkg;

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

