/**					\
 * DMAnalysisTreeMaker			\
 * \
 * Produces analysis trees from edm-ntuples adding extra variables for resolved and unresolved tops\
 * For custom systematics scenarios\
 * \
 * \\Author A. Orso M. Iorio\
 * \
 * \
 *\\version  $Id:\
 * \
 * \
*/ 

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/VectorUtil.h>
#include "./MT2Utility.h"
#include "./mt2w_bisect.h"
#include "./mt2bl_bisect.h"
#include "./Mt2Com_bisect.h"
//#include "./EquationSolver.h"
#include "./DMTopVariables.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <TLorentzVector.h>
#include <TMVA/Reader.h>
#include <string>
#include <iostream>
//#include "TopTagger/Resolved/interface/KinematicFitter.hh"

//using namespace reco;
using namespace edm;
using namespace std;

namespace LHAPDF
{
void initPDFSet(int nset, const std::string &filename, int member = 0);
int numberPDF(int nset);
void usePDFMember(int nset, int member);
double xfx(int nset, double x, double Q, int fl);
double getXmin(int nset, int member);
double getXmax(int nset, int member);
double getQ2min(int nset, int member);
double getQ2max(int nset, int member);
void extrapolate(bool extrapolate = true);
}

class  DMAnalysisTreeMaker : public edm::EDAnalyzer 
{
public:
  explicit DMAnalysisTreeMaker( const edm::ParameterSet & );   

private:

  virtual void analyze(const edm::Event &, const edm::EventSetup & );
  virtual void beginRun(const edm::Run  &, const edm::EventSetup & );
  virtual void endJob();
  vector<string> additionalVariables(string);
  string makeName(string label,string pref,string var);
  string makeBranchNameCat(string label, string cat, string pref, string var);
  string makeBranchName(string label, string pref, string var);
  void initializePdf(string centralpdf,string variationpdf);
  void initTreeWeightHistory(bool useLHEW);
  void getEventPdf();
  int eventFlavour(bool getFlavour, int nb, int nc,int nudsg);
  bool flavourFilter(string ch, int nb, int nc,int nudsg); 

  void initCategoriesSize(string label);
  void setCategorySize(string label, string category, size_t size);
  void fillCategory(string label, string category, int pos_nocat, int pos_cat);

  double getWPtWeight(double ptW);
  double getZPtWeight(double ptZ);

  double getWEWKPtWeight(double ptW);
  double getZEWKPtWeight(double ptZ);
  double getAPtWeight(double ptA);
  double getTopPtWeight(double ptT,double ptTbar, bool extrap = false);
  bool getEventTriggers();
  bool getMETFilters();
  void getEventLHEWeights();

  bool isEWKID(int id);
  //
  // Set up MVA reader
  //
  // spectator variables, not used for MVA evaluation
  int isSig, b_mis, w_mis, wb_mis;
  float mtop;
  // MVA input variables
  //  float bdt_qgid1, bdt_qgid2;
  //  float bdt_dphij1b, bdt_dphij2b, bdt_drj1b, bdt_drj2b;
  //  float bdt_bjcsv, bdt_jet1csv, bdt_jet2csv;
  //  float bdt_prob;
  //  TMVA::Reader res_topmvaReader;

  double jetUncertainty(double pt, double eta, string syst);
  double jetUncertainty8(double pt, double eta, string syst);
  double smear(double pt, double genpt, double eta, string syst);
  double getEffectiveArea(string particle, double eta);
  double resolSF(double eta, string syst);
  double getScaleFactor(double pt, double eta, double partonFlavour, string syst);
  double pileUpSF(string syst);
  double nInitEvents;
  float nTightJets;

  bool isInVector(std::vector<std::string> v, std::string s);
  bool isMCWeightName(std::string s);
  std::vector<edm::ParameterSet > physObjects;
  std::vector<edm::InputTag > variablesFloat, variablesInt, singleFloat,  singleInt;
  
  std::vector<edm::InputTag > variablesDouble, singleDouble;
  
  edm::EDGetTokenT< LHEEventProduct > t_lhes_;
  edm::EDGetTokenT< GenEventInfoProduct > t_genprod_;

  edm::EDGetTokenT< std::vector<string> > t_triggerNames_;
  edm::EDGetTokenT< std::vector<float> > t_triggerBits_;
  edm::EDGetTokenT< std::vector<int> > t_triggerPrescales_;
  edm::EDGetTokenT<std::vector<string> > t_triggerNamesR_;
  
  edm::EDGetTokenT< unsigned int > t_lumiBlock_;
  edm::EDGetTokenT< unsigned int > t_runNumber_;
  edm::EDGetTokenT< ULong64_t > t_eventNumber_;
  
  edm::EDGetTokenT< std::vector<string> > t_metNames_;
  edm::EDGetTokenT< std::vector<float> > t_metBits_;

  edm::EDGetTokenT< std::vector<float> > t_pvZ_,t_pvChi2_,t_pvRho_;
  edm::EDGetTokenT< std::vector<int> > t_pvNdof_;

  edm::EDGetTokenT< double > t_Rho_;
  edm::EDGetTokenT<int> t_ntrpu_;

  edm::EDGetTokenT< std::vector<float> > jetAK8vSubjetIndex0;
  edm::EDGetTokenT< std::vector<float> > jetAK8vSubjetIndex1;

  edm::EDGetTokenT< std::vector<float> > genPartID ;
  edm::EDGetTokenT< std::vector<float> > genPartStatus;
  edm::EDGetTokenT< std::vector<float> > genPartMom0ID; 
  edm::EDGetTokenT< std::vector<float> > genPartPt;
  edm::EDGetTokenT< std::vector<float> > genPartPhi;
  edm::EDGetTokenT< std::vector<float> > genPartEta;
  edm::EDGetTokenT< std::vector<float> > genPartE;

  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;
  
  //edm::EDGetTokenT<reco::GenParticleCollection>(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles")));
  edm::EDGetTokenT<reco::GenParticleCollection> t_genParticleCollection_;

  TH1D * nInitEventsHisto;
  TTree * treesBase;
  map<string, TTree * > trees;
  std::vector<string> names;
  std::vector<string> systematics;
  map< string , float[100] > vfloats_values;
  map< string , int[100] > vints_values;
  map< string , vector<string> > obj_to_floats,obj_to_ints, obj_to_doubles;
  map< string , string > obs_to_obj;
  map< string , string > obj_to_pref;
  map< string , std::vector<string> > obj_cats;

  map< string , double[100] > vdoubles_values;
  map< string , double[100] > vdouble_values;
  map< string , double > double_values;


  map< string , float > float_values;
  map< string , int > int_values;
  //  map< string , ULong64_t > int_values;
  map< string , int > sizes;

  map< string , bool > got_label; 
  map< string , int > max_instances; 
  map< int, int > subj_jet_map;

  map<string, edm::Handle<std::vector<float> > > h_floats;
  map<string, edm::Handle<std::vector<int> > > h_ints;
  map<string, edm::Handle<float> > h_float;
  map<string, edm::Handle<int> >h_int;
  
  map<string, edm::Handle<std::vector<double> > > h_doubles;
  map<string, edm::Handle<double> > h_double;
  
  map<string, edm::EDGetTokenT< std::vector<float> >  > t_floats;
  map<string, edm::EDGetTokenT< std::vector<int> > > t_ints;
  map<string, edm::EDGetTokenT<float>  > t_float;
  map<string, edm::EDGetTokenT<int> >t_int;
  map<string, edm::EDGetTokenT<std::vector<double> > > t_doubles;
  map<string, edm::EDGetTokenT<double> > t_double;
  

  string gen_label,  mu_label, ele_label, jets_label, boosted_tops_label, boosted_tops_subjets_label, met_label, photon_label;//metNoHF_label

  bool getPartonW, getPartonTop, doWReweighting, doTopReweighting;
  bool getParticleWZ, getWZFlavour;
  
  //Do resolved top measurement:
  bool doResolvedTopHad,doResolvedTopSemiLep;
  int max_leading_jets_for_top;
  int max_bjets_for_top;
  int max_genparticles;
  int n0;

  bool recalculateEA;

  //MC info:
  edm::ParameterSet channelInfo;
  std::string channel;
  double crossSection, originalEvents;
  bool useLHEWeights, useLHE, useTriggers,cutOnTriggers, useMETFilters, addPV;
  bool addLHAPDFWeights;
  string centralPdfSet,variationPdfSet;
  std::vector<string> SingleElTriggers, SingleMuTriggers, PhotonTriggers, hadronicTriggers,metFilters, HadronPFHT900Triggers, HadronPFHT800Triggers, HadronPFJet450Triggers;
  int maxPdf, maxWeights;
  edm::Handle<LHEEventProduct > lhes;
  edm::Handle<GenEventInfoProduct> genprod;

  //Trigger info
  edm::Handle<std::vector<float> > triggerBits;
  edm::Handle<std::vector<string> > triggerNames;
  edm::Handle<std::vector<int> > triggerPrescales;
  edm::Handle<std::vector<string> > triggerNamesR ;
  
  edm::Handle<unsigned int> lumiBlock;
  edm::Handle<unsigned int> runNumber;
  edm::Handle<ULong64_t> eventNumber;

  edm::Handle<std::vector<float> > metBits;
  edm::Handle<std::vector<string> > metNames;
  edm::InputTag metNames_;
  
  edm::Handle<std::vector<float> > pvZ,pvChi2,pvRho;
  edm::Handle<std::vector<int> > pvNdof;

  edm::Handle<std::vector<float> > ak8jetvSubjetIndex0;
  edm::Handle<std::vector<float> > ak8jetvSubjetIndex1;
  
  float nPV;
  edm::Handle<int> ntrpu;

  //JEC info
  bool changeJECs;
  bool isData, applyRes;
  bool isV2;
  edm::Handle<double> rho;
  double Rho;
  edm::Handle<reco::GenParticleCollection> genParticles;   
  
  //edm::Handle<double> Rho;
  std::vector<double> jetScanCuts;
  std::vector<JetCorrectorParameters> jecPars;
  JetCorrectorParameters *jecParsL1, *jecParsL1RC, *jecParsL2, *jecParsL3, *jecParsL2L3Residuals;
  JetCorrectionUncertainty *jecUnc;
  FactorizedJetCorrector *jecCorr;

  std::vector<JetCorrectorParameters> jecPars8;
  JetCorrectorParameters *jecParsL18, *jecParsL1RC8, *jecParsL28, *jecParsL38, *jecParsL2L3Residuals8;
  JetCorrectionUncertainty *jecUnc8;
  FactorizedJetCorrector *jecCorr8;

  bool isFirstEvent;
  //Do preselection
  bool doPreselection;
  
  class BTagWeight
  {
  private:
    int minTags;
    int maxTags;
  public:
    struct JetInfo
    {
      JetInfo(float mceff, float datasf) : eff(mceff), sf(datasf) {}
      float eff;
      float sf;
    };
    BTagWeight():
      minTags(0), maxTags(0)
    {
      ;
    }
    BTagWeight(int jmin, int jmax) :
      minTags(jmin) , maxTags(jmax) {}
    bool filter(int t);
    float weight(vector<JetInfo> jets, int tags);
    float weightWithVeto(vector<JetInfo> jetsTags, int tags, vector<JetInfo> jetsVetoes, int vetoes);
  };
  vector<BTagWeight::JetInfo> jsfscsvt, 
    jsfscsvt_b_tag_up, 
    jsfscsvt_b_tag_down, 
    jsfscsvt_mistag_up, 
    jsfscsvt_mistag_down;

  vector<BTagWeight::JetInfo> jsfscsvm, 
    jsfscsvm_b_tag_up, 
    jsfscsvm_b_tag_down, 
    jsfscsvm_mistag_up, 
    jsfscsvm_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvl, 
    jsfscsvl_b_tag_up, 
    jsfscsvl_b_tag_down, 
    jsfscsvl_mistag_up, 
    jsfscsvl_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvt_subj, 
    jsfscsvt_subj_b_tag_up, 
    jsfscsvt_subj_b_tag_down, 
    jsfscsvt_subj_mistag_up, 
    jsfscsvt_subj_mistag_down;

  vector<BTagWeight::JetInfo> jsfscsvm_subj, 
    jsfscsvm_subj_b_tag_up, 
    jsfscsvm_subj_b_tag_down, 
    jsfscsvm_subj_mistag_up, 
    jsfscsvm_subj_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvl_subj, 
    jsfscsvl_subj_b_tag_up, 
    jsfscsvl_subj_b_tag_down, 
    jsfscsvl_subj_mistag_up, 
    jsfscsvl_subj_mistag_down;
  
  BTagWeight b_csvt_0_tags= BTagWeight(0,1000),
    b_csvt_1_tag= BTagWeight(1,1000),
    b_csvt_1_2_tags= BTagWeight(2,1000),
    b_csvt_2_tags= BTagWeight(3,1000);

  double b_weight_csvt_0_tags,
    b_weight_csvt_1_tag,
    b_weight_csvt_1_2_tags,
    b_weight_csvt_2_tags;
  double b_weight_csvt_0_tags_mistag_up,
    b_weight_csvt_1_tag_mistag_up,
    b_weight_csvt_1_2_tags_mistag_up,
    b_weight_csvt_2_tags_mistag_up;
  double b_weight_csvt_0_tags_mistag_down,
    b_weight_csvt_1_tag_mistag_down,
    b_weight_csvt_1_2_tags_mistag_down,
    b_weight_csvt_2_tags_mistag_down;
  double b_weight_csvt_0_tags_b_tag_down,
    b_weight_csvt_1_tag_b_tag_down,
    b_weight_csvt_1_2_tags_b_tag_down,
    b_weight_csvt_2_tags_b_tag_down;
  double b_weight_csvt_0_tags_b_tag_up,
    b_weight_csvt_1_tag_b_tag_up,
    b_weight_csvt_1_2_tags_b_tag_up,
    b_weight_csvt_2_tags_b_tag_up;
  
  BTagWeight b_csvm_0_tags= BTagWeight(0,1000),
    b_csvm_1_tag= BTagWeight(1,1000),
    b_csvm_1_2_tags= BTagWeight(2,1000),
    b_csvm_2_tags= BTagWeight(3,1000);
  
  double b_weight_csvm_0_tags,
    b_weight_csvm_1_tag,
    b_weight_csvm_1_2_tags,
    b_weight_csvm_2_tags;
  double b_weight_csvm_0_tags_mistag_up,
    b_weight_csvm_1_tag_mistag_up,
    b_weight_csvm_1_2_tags_mistag_up,
    b_weight_csvm_2_tags_mistag_up;
  double b_weight_csvm_0_tags_mistag_down,
    b_weight_csvm_1_tag_mistag_down,
    b_weight_csvm_1_2_tags_mistag_down,
    b_weight_csvm_2_tags_mistag_down;
  double b_weight_csvm_0_tags_b_tag_down,
    b_weight_csvm_1_tag_b_tag_down,
    b_weight_csvm_1_2_tags_b_tag_down,
    b_weight_csvm_2_tags_b_tag_down;
  double b_weight_csvm_0_tags_b_tag_up,
    b_weight_csvm_1_tag_b_tag_up,
    b_weight_csvm_1_2_tags_b_tag_up,
    b_weight_csvm_2_tags_b_tag_up;
  
  BTagWeight b_subj_csvm_0_tags= BTagWeight(0,1000),
    b_subj_csvm_1_tag= BTagWeight(1,1000),
    b_subj_csvm_1_2_tags= BTagWeight(0,1),
    b_subj_csvm_2_tags= BTagWeight(1,1);

  double b_weight_subj_csvm_0_tags,
    b_weight_subj_csvm_1_tag,
    b_weight_subj_csvm_1_2_tags,
    b_weight_subj_csvm_2_tags;
  double b_weight_subj_csvm_0_tags_mistag_up,
    b_weight_subj_csvm_1_tag_mistag_up,
    b_weight_subj_csvm_1_2_tags_mistag_up,
    b_weight_subj_csvm_2_tags_mistag_up;
  double b_weight_subj_csvm_0_tags_mistag_down,
    b_weight_subj_csvm_1_tag_mistag_down,
    b_weight_subj_csvm_1_2_tags_mistag_down,
    b_weight_subj_csvm_2_tags_mistag_down;
  double b_weight_subj_csvm_0_tags_b_tag_down,
    b_weight_subj_csvm_1_tag_b_tag_down,
    b_weight_subj_csvm_1_2_tags_b_tag_down,
    b_weight_subj_csvm_2_tags_b_tag_down;
  double b_weight_subj_csvm_0_tags_b_tag_up,
    b_weight_subj_csvm_1_tag_b_tag_up,
    b_weight_subj_csvm_1_2_tags_b_tag_up,
    b_weight_subj_csvm_2_tags_b_tag_up;
  
  BTagWeight b_csvl_0_tags= BTagWeight(0,1000),
    b_csvl_1_tag= BTagWeight(1,1000),
    b_csvl_1_2_tags= BTagWeight(0,1),
    b_csvl_2_tags= BTagWeight(1,1);
  
  double b_weight_csvl_0_tags_mistag_up,
    b_weight_csvl_1_tag_mistag_up,
    b_weight_csvl_1_2_tags_mistag_up,
    b_weight_csvl_2_tags_mistag_up;
  double b_weight_csvl_0_tags_mistag_down,
    b_weight_csvl_1_tag_mistag_down,
    b_weight_csvl_1_2_tags_mistag_down,
    b_weight_csvl_2_tags_mistag_down;
  double b_weight_csvl_0_tags_b_tag_down,
    b_weight_csvl_1_tag_b_tag_down,
    b_weight_csvl_1_2_tags_b_tag_down,
    b_weight_csvl_2_tags_b_tag_down;
  double b_weight_csvl_0_tags_b_tag_up,
    b_weight_csvl_1_tag_b_tag_up,
    b_weight_csvl_1_2_tags_b_tag_up,
    b_weight_csvl_2_tags_b_tag_up;
  double b_weight_csvl_0_tags,
    b_weight_csvl_1_tag,
    b_weight_csvl_1_2_tags,
    b_weight_csvl_2_tags;
  
  BTagWeight b_subj_csvl_0_tags= BTagWeight(0,1000),
    b_subj_csvl_1_tag= BTagWeight(1,1000),
    b_subj_csvl_1_2_tags= BTagWeight(2,1000),
    b_subj_csvl_2_tags= BTagWeight(3,1000);
  
  double b_weight_subj_csvl_0_tags_mistag_up,
    b_weight_subj_csvl_1_tag_mistag_up,
    b_weight_subj_csvl_1_2_tags_mistag_up,
    b_weight_subj_csvl_2_tags_mistag_up;
  double b_weight_subj_csvl_0_tags_mistag_down,
    b_weight_subj_csvl_1_tag_mistag_down,
    b_weight_subj_csvl_1_2_tags_mistag_down,
    b_weight_subj_csvl_2_tags_mistag_down;
  double b_weight_subj_csvl_0_tags_b_tag_down,
    b_weight_subj_csvl_1_tag_b_tag_down,
    b_weight_subj_csvl_1_2_tags_b_tag_down,
    b_weight_subj_csvl_2_tags_b_tag_down;
  double b_weight_subj_csvl_0_tags_b_tag_up,
    b_weight_subj_csvl_1_tag_b_tag_up,
    b_weight_subj_csvl_1_2_tags_b_tag_up,
    b_weight_subj_csvl_2_tags_b_tag_up;
  double b_weight_subj_csvl_0_tags,
    b_weight_subj_csvl_1_tag,
    b_weight_subj_csvl_1_2_tags,
    b_weight_subj_csvl_2_tags;
  
  double MCTagEfficiency(string algo, int flavor, double pt, double eta);
  double MCTagEfficiencySubjet(string algo, int flavor, double pt, double eta); 
  double TagScaleFactor(string algo, int flavor, string syst,double pt);
  double TagScaleFactorSubjet(string algo, int flavor, string syst,double pt);
  
  //
  bool doBTagSF;
  bool doPU;
  
  string season;
  //  season = dataPUFile_;
  
  string distr;
  
  
};


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig){
  
  mu_label = iConfig.getParameter<std::string >("muLabel");
  ele_label = iConfig.getParameter<std::string >("eleLabel");
  jets_label = iConfig.getParameter<std::string >("jetsLabel");
  photon_label = iConfig.getParameter<std::string >("photonLabel");
  boosted_tops_label = iConfig.getParameter<std::string >("boostedTopsLabel");
  boosted_tops_subjets_label = iConfig.getParameter<std::string >("boostedTopsSubjetsLabel");
  met_label = iConfig.getParameter<std::string >("metLabel");
  gen_label = iConfig.getParameter<std::string >("genLabel");
  physObjects = iConfig.template getParameter<std::vector<edm::ParameterSet> >("physicsObjects");
  
  channelInfo = iConfig.getParameter<edm::ParameterSet >("channelInfo"); // The physics of the channel, e.g. the cross section, #original events, etc.
  channel = channelInfo.getParameter<std::string>("channel");
  crossSection = channelInfo.getParameter<double>("crossSection");
  originalEvents = channelInfo.getParameter<double>("originalEvents");

  doPreselection = iConfig.getUntrackedParameter<bool>("doPreselection",true);
  doPU = iConfig.getUntrackedParameter<bool>("doPU",true);

  useLHEWeights = channelInfo.getUntrackedParameter<bool>("useLHEWeights",false);
  useLHE = channelInfo.getUntrackedParameter<bool>("useLHE",false);
  addLHAPDFWeights = channelInfo.getUntrackedParameter<bool>("addLHAPDFWeights",false);

  getPartonW = channelInfo.getUntrackedParameter<bool>("getPartonW",false);
  getParticleWZ = channelInfo.getUntrackedParameter<bool>("getParticleWZ",false);
  getPartonTop = channelInfo.getUntrackedParameter<bool>("getPartonTop",false);
  doWReweighting = channelInfo.getUntrackedParameter<bool>("doWReweighting",false);

  getWZFlavour = channelInfo.getUntrackedParameter<bool>("getWZFlavour",false);

  doResolvedTopSemiLep = iConfig.getUntrackedParameter<bool>("doResolvedTopSemiLep",false);
  doResolvedTopHad = iConfig.getUntrackedParameter<bool>("doResolvedTopHad",false);

  edm::InputTag genprod_ = iConfig.getParameter<edm::InputTag>( "genprod" );
  t_genprod_ = consumes<GenEventInfoProduct>( genprod_ );
  
  useTriggers = iConfig.getUntrackedParameter<bool>("useTriggers",true);
  cutOnTriggers = iConfig.getUntrackedParameter<bool>("cutOnTriggers",true);

  edm::InputTag ak8jetvSubjetIndex0_ = iConfig.getParameter<edm::InputTag>("ak8jetvSubjetIndex0");
  jetAK8vSubjetIndex0  = consumes< std::vector<float> >( ak8jetvSubjetIndex0_);
  edm::InputTag ak8jetvSubjetIndex1_ = iConfig.getParameter<edm::InputTag>("ak8jetvSubjetIndex1");
  jetAK8vSubjetIndex1 = consumes< std::vector<float> >( ak8jetvSubjetIndex1_);

  edm::InputTag lumiBlock_ = iConfig.getParameter<edm::InputTag>("lumiBlock");
  t_lumiBlock_ = consumes< unsigned int >( lumiBlock_ );
  edm::InputTag runNumber_ = iConfig.getParameter<edm::InputTag>("runNumber");
  t_runNumber_ = consumes< unsigned int >( runNumber_ );
  edm::InputTag eventNumber_ = iConfig.getParameter<edm::InputTag>("eventNumber");
  t_eventNumber_ = consumes< ULong64_t >( eventNumber_ );
  
  if(useTriggers){
     edm::InputTag triggerBits_ = iConfig.getParameter<edm::InputTag>("triggerBits");
    t_triggerBits_ = consumes< std::vector<float> >( triggerBits_ );
    edm::InputTag triggerPrescales_ = iConfig.getParameter<edm::InputTag>("triggerPrescales");
    t_triggerPrescales_ = consumes< std::vector<int> >( triggerPrescales_ );

    //testing to read from Run
    edm::InputTag triggerNamesR_ = iConfig.getParameter<edm::InputTag>("triggerNames");
    t_triggerNamesR_ = mayConsume< std::vector<string>, edm::InRun>(edm::InputTag("TriggerUserData","triggerNameTree"));

    SingleElTriggers= channelInfo.getParameter<std::vector<string> >("SingleElTriggers");
    SingleMuTriggers= channelInfo.getParameter<std::vector<string> >("SingleMuTriggers");
    PhotonTriggers= channelInfo.getParameter<std::vector<string> >("PhotonTriggers");
    hadronicTriggers= channelInfo.getParameter<std::vector<string> >("hadronicTriggers");
    HadronPFHT900Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFHT900Triggers");
    HadronPFHT800Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFHT800Triggers");
    HadronPFJet450Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFJet450Triggers");
  }
  useMETFilters = iConfig.getUntrackedParameter<bool>("useMETFilters",true);

  if(useMETFilters){
    metFilters = channelInfo.getParameter<std::vector<string> >("metFilters");
    edm::InputTag metBits_ = iConfig.getParameter<edm::InputTag>("metBits");
    t_metBits_ = consumes< std::vector<float> >( metBits_ );
    metNames_ = iConfig.getParameter<edm::InputTag>("metNames");
    t_metNames_ = consumes< std::vector<string>, edm::InRun >( metNames_ );
  }
  
  addPV = iConfig.getUntrackedParameter<bool>("addPV",true);
  changeJECs = iConfig.getUntrackedParameter<bool>("changeJECs",false);
  recalculateEA = iConfig.getUntrackedParameter<bool>("recalculateEA",true);

  isData = iConfig.getUntrackedParameter<bool>("isData",false);
  applyRes = iConfig.getUntrackedParameter<bool>("applyRes",false);
  isV2= iConfig.getUntrackedParameter<bool>("isV2",false);
  
  t_Rho_ = consumes<double>( edm::InputTag( "fixedGridRhoFastjetAll" ) ) ;
  t_genParticleCollection_ = consumes<reco::GenParticleCollection>(edm::InputTag("filteredPrunedGenParticles"));

  if(addPV || changeJECs){
    edm::InputTag pvZ_ = iConfig.getParameter<edm::InputTag >("vertexZ");
    t_pvZ_ = consumes< std::vector<float> >( pvZ_ );
    edm::InputTag pvChi2_ = iConfig.getParameter<edm::InputTag >("vertexChi2");
    t_pvChi2_ = consumes< std::vector<float> >( pvChi2_ );
    edm::InputTag pvRho_ = iConfig.getParameter<edm::InputTag >("vertexRho");
    t_pvRho_ = consumes< std::vector<float> >( pvRho_ );
    edm::InputTag pvNdof_ = iConfig.getParameter<edm::InputTag >("vertexNdof");
    t_pvNdof_ = consumes< std::vector< int > >( pvNdof_ );
  }

  if (doPU)
    t_ntrpu_ = consumes< int >( edm::InputTag( "eventUserData","puNtrueInt" ) );

  
  maxWeights = 9;
  if(useLHEWeights){
    maxWeights = channelInfo.getUntrackedParameter<int>("maxWeights",9);//Usually we do have 9 weights for the scales, might vary depending on the lhe
  }
  if(addLHAPDFWeights){
    centralPdfSet = channelInfo.getUntrackedParameter<string>("pdfSet","NNPDF");
    variationPdfSet = channelInfo.getUntrackedParameter<string>("pdfSet","NNPDF");
    initializePdf(centralPdfSet,variationPdfSet);

  }
  if(doResolvedTopHad){
    max_leading_jets_for_top  = iConfig.getUntrackedParameter<int>("maxLeadingJetsForTop",8);//Take the 8 leading jets for the top permutations
  }
  if(doResolvedTopSemiLep){
    max_bjets_for_top  = iConfig.getUntrackedParameter<int>("maxBJetsForTop",2);//Take the 8 leading jets for the top permutations
  }
  if(!isData){//G
    max_genparticles  = iConfig.getUntrackedParameter<int>("maxGenPart",30);
   }
  systematics = iConfig.getParameter<std::vector<std::string> >("systematics");

  jetScanCuts = iConfig.getParameter<std::vector<double> >("jetScanCuts");


  
  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();

  bool addNominal=false;
  for (size_t s = 0; s<systematics.size();++s){
    if(systematics.at(s).find("noSyst")!=std::string::npos) {
      addNominal=true;
      break;
    }
  }
  if(systematics.size()==0){
    addNominal=true;
    systematics.push_back("noSyst");
  }//In case there's no syst specified, do the nominal scenario
  //addNominal=true;
  Service<TFileService> fs;
  TFileDirectory DMTrees;// = fs->mkdir( "systematics_trees" );

  if(addNominal){
    DMTrees = fs->mkdir( "systematics_trees" );
  }

  trees["noSyst"] =  new TTree((channel+"__noSyst").c_str(),(channel+"__noSyst").c_str());

  nInitEventsHisto = new TH1D("initialEvents","initalEvents",10,0,10);

  edm::InputTag lhes_ = iConfig.getParameter<edm::InputTag>( "lhes" );
  t_lhes_ = consumes< LHEEventProduct >( lhes_ );
  
  for (;itPsets!=physObjects.end();++itPsets){ 
    int maxI = itPsets->getUntrackedParameter< int >("maxInstances",10);
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI");
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("singleD"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    string namelabel = itPsets->getParameter< string >("label");
    string nameprefix = itPsets->getParameter< string >("prefix");
    bool saveBaseVariables = itPsets->getUntrackedParameter<bool>("saveBaseVariables",true);
    bool saveNoCat = itPsets->getUntrackedParameter<bool>("saveNoCat",true);

    std::vector<std::string > categories = itPsets->getParameter<std::vector<std::string> >("categories");
    std::vector<std::string > toSave= itPsets->getParameter<std::vector<std::string> >("toSave");
    
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();

    for(size_t sc = 0; sc< categories.size() ;++sc){
      string category = categories.at(sc);
      obj_cats[namelabel].push_back(category);
    }
    stringstream max_instance_str;
    max_instance_str<<maxI;
    max_instances[namelabel]=maxI;
    string nameobs = namelabel;
    string prefix = nameprefix;

    for (;itF != variablesFloat.end();++itF){
      
      string name=itF->instance()+"_"+itF->label();
      string nameinstance=itF->instance();
      string nameshort=itF->instance();
      
      string nametobranch = makeBranchName(namelabel,prefix,nameinstance);
      
      name = nametobranch;
      nameshort = nametobranch;
    
      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itF->instance()))) trees["noSyst"]->Branch(nameshort.c_str(), &vfloats_values[name],(nameshort+"["+max_instance_str.str()+"]/F").c_str());
      names.push_back(name);
      obj_to_floats[namelabel].push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;

      t_floats[ name ] = consumes< std::vector<float> >( *itF );
      
      for(size_t sc = 0; sc< categories.size() ;++sc){
	string category = categories.at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameinstance);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
      }
    }
  
    for (;itI != variablesInt.end();++itI){
      string name=itI->instance()+"_"+itI->label();
      string nameshort=itF->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;

      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itI->instance())) ) trees["noSyst"]->Branch(nameshort.c_str(), &vints_values[name],(nameshort+"["+max_instance_str.str()+"]/I").c_str());
      for(size_t sc = 0; sc< categories.size() ;++sc){
	string category = categories.at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameshort);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/I").c_str());
      }

      names.push_back(name);
      obj_to_ints[namelabel].push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;

      t_ints[ name ] = consumes< std::vector<int> >( *itI );

    }
    
    if (variablesFloat.size()>0){
      string nameshortv = namelabel;
      vector<string> extravars = additionalVariables(nameshortv);
      for(size_t addv = 0; addv < extravars.size();++addv){
	string name = nameshortv+"_"+extravars.at(addv);
	if (saveNoCat && (saveBaseVariables || isInVector(toSave, extravars.at(addv)) || isInVector(toSave, "allExtra") ) )trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+max_instance_str.str()+"]/F").c_str());
	for(size_t sc = 0; sc< categories.size() ;++sc){
	  string category = categories.at(sc);
	  string nametobranchcat = nameshortv+category+"_"+extravars.at(addv);
	  string namecat = nametobranchcat;
	  cout << "extra var "<< extravars.at(addv)<<endl;
	  cout << " namecat "<< namecat<< endl;
	  if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
	}

	obj_to_floats[namelabel].push_back(name);
	obs_to_obj[name] = nameobs;
	obj_to_pref[nameobs] = prefix;
      }
    }
    names.push_back(nameobs);
    cout << "size part: nameobs is  "<< nameobs<<endl;
    if(saveNoCat) trees["noSyst"]->Branch((nameobs+"_size").c_str(), &sizes[nameobs]);
    for(size_t sc = 0; sc< categories.size() ;++sc){
      string category = categories.at(sc);
      trees["noSyst"]->Branch((nameobs+category+"_size").c_str(), &sizes[nameobs+category]);
    }
    
    //Initialize single pset objects
     for (;itsF != singleFloat.end();++itsF){
      string name=itsF->instance()+itsF->label();
      string nameshort=itsF->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_float[ name ] = consumes< float >( *itsF );
      if(saveBaseVariables|| isInVector(toSave,itsF->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &float_values[name]);
    }
 
    for (;itsD != singleDouble.end();++itsD){
      string name=itsD->instance()+itsD->label();
      string nameshort=itsD->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_double[ name ] = consumes< double >( *itsD );
      if(saveBaseVariables|| isInVector(toSave,itsD->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &double_values[name]);
    }
    for (;itsI != singleInt.end();++itsI){
      string name=itsI->instance()+itsI->label();
      string nameshort=itsI->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_int[ name ] = consumes< int >( *itsI );
      if(saveBaseVariables|| isInVector(toSave,itsI->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &int_values[name]);
    }
  }
  if(doResolvedTopSemiLep){
    string nameshortv= "resolvedTopSemiLep";
    vector<string> extravarstop = additionalVariables(nameshortv);
    double max_instances_top = max_bjets_for_top*max((int)(max_instances[ele_label]),(int)(max_instances[mu_label]) );
    max_instances[nameshortv]=max_instances_top;
    stringstream mtop;
    mtop << max_instances_top;
    cout << " max instances top is "<< max_instances_top << " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }
  if(doResolvedTopHad){
    string nameshortv= "resolvedTopHad";
    vector<string> extravarstop = additionalVariables(nameshortv);
    double max_instances_top = TMath::Binomial(min((int)(max_instances[jets_label]),max_leading_jets_for_top),4);
    max_instances[nameshortv]=max_instances_top;
    stringstream mtop;
    mtop << max_instances_top;
    cout << " max instances top is "<< max_instances_top<< " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }
  
  if(!isData){//G
    string nameshortv= "genPart";
    vector<string> extravarstop = additionalVariables(nameshortv);
    //double max_instances_gen = max_genparticles*max((int)(max_instances[gen_label]),(int)(max_instances[gen_label]) );
    double max_instances_gen = 30*max((int)0,(int)1);
    max_instances[nameshortv]=max_instances_gen;
    stringstream mtop;
    mtop << max_instances_gen;
    //cout << " max instances top is "<< max_instances_top<< " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }

  string nameshortv= "Event";
  vector<string> extravars = additionalVariables(nameshortv);
  for(size_t addv = 0; addv < extravars.size();++addv){
    string name = nameshortv+"_"+extravars.at(addv);

    if (name.find("EventNumber")!=std::string::npos){
      std::cout<<"=====================sto riempendo il branch event number"<<std::endl;
      trees["noSyst"]->Branch(name.c_str(), &double_values[name],(name+"/D").c_str());
    }
    else trees["noSyst"]->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
  }
  
  //Prepare the trees cloning all branches and setting the correct names/titles:
  if(!addNominal){
    DMTrees = fs->mkdir( "systematics_trees" );
  }
  
  trees["EventHistory"] =  new TTree("EventHistory","EventHistory");
  trees["WeightHistory"] =  new TTree("WeightHistory","WeightHistory");
  trees["EventHistory"]->Branch("initialEvents",&nInitEvents);

  for(size_t s=0;s< systematics.size();++s){
    std::string syst  = systematics.at(s);
    if(syst=="noSyst")continue;
    trees[syst]= trees["noSyst"]->CloneTree();
    trees[syst]->SetName((channel+"__"+syst).c_str());
    trees[syst]->SetTitle((channel+"__"+syst).c_str());
  }

  initTreeWeightHistory(useLHEWeights);

  string L1Name = "Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt"; //
  string L1RCName = "Spring16_25nsV6_MC_L1RC_AK4PFchs.txt"; 
  string L2Name = "Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt";
  string L3Name = "Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt";
  string L2L3ResName = "Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt";
  if(isData){
    L1Name   = "Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    L1RCName = "Spring16_25nsV6_DATA_L1RC_AK4PFchs.txt";  
    L2Name   = "Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt";
    L3Name   = "Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    L2L3ResName = "Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt";
  }

  if(isV2){
      L1Name = "Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt";
      L2Name = "Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt";
      L3Name = "Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt";
      L2L3ResName = "Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt"; 
  }

  jecParsL1  = new JetCorrectorParameters(L1Name);
  jecParsL2  = new JetCorrectorParameters(L2Name);
  jecParsL3  = new JetCorrectorParameters(L3Name);
  jecParsL2L3Residuals  = new JetCorrectorParameters(L2L3ResName);
  jecPars.push_back(*jecParsL1);
  jecPars.push_back(*jecParsL2);
  jecPars.push_back(*jecParsL3);
  if(isData)jecPars.push_back(*jecParsL2L3Residuals);

  jecCorr = new FactorizedJetCorrector(jecPars);
  jecUnc  = new JetCorrectionUncertainty(*(new JetCorrectorParameters("Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt", "Total")));

  //JEC on AK8

  string L1Name8 = "Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt"; 
  string L1RCName8 = "Spring16_25nsV6_MC_L1RC_AK8PFchs.txt"; 
  string L2Name8 = "Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt";
  string L3Name8 = "Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt";
  string L2L3ResName8 = "Spring16_25nsV6_MC_L2L3Residual_AK8PFchs.txt";
  if(isData){
    L1Name8   = "Spring16_25nsV6_DATA_L1FastJet_AK8PFchs.txt";
    L1RCName8 = "Spring16_25nsV6_DATA_L1RC_AK8PFchs.txt";  
    L2Name8   = "Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt";
    L3Name8   = "Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt";
    L2L3ResName8 = "Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt";
  }

  if(isV2){
      L1Name8 = "Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt";
      L2Name8 = "Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt";
      L3Name8 = "Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt";
      L2L3ResName8 = "Spring16_25nsV6_MC_L2L3Residual_AK8PFchs.txt"; 
  }

  jecParsL18  = new JetCorrectorParameters(L1Name8);
  jecParsL28  = new JetCorrectorParameters(L2Name8);
  jecParsL38  = new JetCorrectorParameters(L3Name8);
  jecParsL2L3Residuals8  = new JetCorrectorParameters(L2L3ResName8);
  jecPars8.push_back(*jecParsL18);
  jecPars8.push_back(*jecParsL28);
  jecPars8.push_back(*jecParsL38);
  if(isData)jecPars8.push_back(*jecParsL2L3Residuals8);

  jecCorr8 = new FactorizedJetCorrector(jecPars8);
  jecUnc8  = new JetCorrectionUncertainty(*(new JetCorrectorParameters("Spring16_25nsV6_DATA_UncertaintySources_AK8PFchs.txt", "Total")));

  
  isFirstEvent = true;
  doBTagSF= true;
  if(isData)doPU= false;
  
  season = "Summer11";
  
  distr = "pileUpDistr" + season + ".root";

}


void DMAnalysisTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    iRun.getByLabel(edm::InputTag("TriggerUserData","triggerNameTree"), triggerNamesR);
    iRun.getByLabel(metNames_, metNames);
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //cout << "trigger test tname "<< tname <<endl; 
    }
}


void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();

  nInitEventsHisto->Fill(0.1);
  nInitEvents+=1;
  //cout << "opening analyze" << endl;
  // event info
  iEvent.getByToken(t_lumiBlock_,lumiBlock );
  iEvent.getByToken(t_runNumber_,runNumber );
  iEvent.getByToken(t_eventNumber_,eventNumber );

  //  cout << "starting gen" << endl;
  if(useLHE){
    iEvent.getByToken(t_lhes_, lhes);
  }
  if(addLHAPDFWeights){
    iEvent.getByToken(t_genprod_, genprod);
  }
  
  vector<TLorentzVector> genlep,gentop,genantitop,genw,gennu,genz,gena;
  vector<TLorentzVector> pare,parebar,parmu,parmubar,parnumu,parnumubar,partau,partaubar,parnue,parnuebar,parnutau,parnutaubar,parz,parw;

  //Parton-level info                                                                                                                               
  if(isData==true){
    getPartonW=false;
    getParticleWZ=false;
    getPartonTop=false;
    doWReweighting=false;
    doTopReweighting=false;
    useLHE=false;
    useLHEWeights=false;
  }
  
  if(getPartonW || getPartonTop || doWReweighting || doTopReweighting){
    if(!useLHE)return;
    genlep.clear();
    gentop.clear();
    genantitop.clear();
    genw.clear();
    genz.clear();
    gena.clear();
    gennu.clear();
    pare.clear();parebar.clear();parmu.clear();parmubar.clear();parnumu.clear();parnumubar.clear();parnue.clear();parnuebar.clear();parnutau.clear();parnutaubar.clear();parz.clear();parw.clear();

    float_values["Event_Z_EW_Weight"]= 1.0;
    float_values["Event_W_EW_Weight"]= 1.0;
    float_values["Event_Z_QCD_Weight"]= 1.0;
    float_values["Event_W_QCD_Weight"]= 1.0;
    float_values["Event_Z_Weight"]= 1.0;
    float_values["Event_W_Weight"]= 1.0;
    float_values["Event_T_Weight"]= 1.0;
    float_values["Event_T_Ext_Weight"]= 1.0;
    size_t nup=lhes->hepeup().NUP;

    for( size_t i=0;i<nup;++i){
      //      cout << " particle number " << i << endl;
      int id = lhes->hepeup().IDUP[i];
      float px = lhes->hepeup().PUP[i][0];
      float py = lhes->hepeup().PUP[i][1];
      float pz = lhes->hepeup().PUP[i][2];
      float energy = lhes->hepeup().PUP[i][3];
    
      TLorentzVector vec;
      math::XYZTLorentzVector part = math::XYZTLorentzVector(px, py, pz, energy);
      float pt = part.pt();
      float phi = part.phi();
      float eta = part.eta();
      
      if(pt>0){
	vec.SetPtEtaPhiE(pt, eta, phi, energy);
	
	if(abs (id) == 11 || abs (id) == 13 || abs(id) == 15){
	  genlep.push_back(vec);
	}
	if(abs (id) == 12 || abs (id) == 14 || abs(id) == 16){
	  gennu.push_back(vec);
	}
	if(id == 6 ){
	  gentop.push_back(vec);
	}
	if(id == -6 ){
	  genantitop.push_back(vec);
	}
	if(abs (id) == 24 ){
	  genw.push_back(vec);
	}
	if(abs (id) == 23 ){
	  genz.push_back(vec);
	}
	if(abs (id) == 22 ){
	  gena.push_back(vec);
	}
      }      
    }
    
    if(getPartonTop && gentop.size()==1){
      float_values["Event_T_Pt"]= gentop.at(0).Pt();
      float_values["Event_T_Eta"]= gentop.at(0).Eta();
      float_values["Event_T_Phi"]= gentop.at(0).Phi();
      float_values["Event_T_E"]= gentop.at(0).Energy();
      float_values["Event_T_Mass"]= gentop.at(0).M();
    }
    if(getPartonTop && genantitop.size()==1){
      float_values["Event_Tbar_Pt"]= genantitop.at(0).Pt();
      float_values["Event_Tbar_Eta"]= genantitop.at(0).Eta();
      float_values["Event_Tbar_Phi"]= genantitop.at(0).Phi();
      float_values["Event_Tbar_E"]= genantitop.at(0).Energy();
      float_values["Event_Tbar_Mass"]= genantitop.at(0).M();
      
    }
    if((getPartonW || doWReweighting )) {
      if(genw.size()==1){
	float_values["Event_W_Pt"]= genw.at(0).Pt();
	float_values["Event_W_Eta"]= genw.at(0).Eta();
	float_values["Event_W_Phi"]= genw.at(0).Phi();
	float_values["Event_W_E"]= genw.at(0).Energy();
	float_values["Event_W_Mass"]= genw.at(0).M();	
	
	double ptW = genw.at(0).Pt();
	double wweight = getWPtWeight(ptW);			
	float_values["Event_W_QCD_Weight"]= wweight;
      }
      else (float_values["Event_W_QCD_Weight"]=1.0);
    }
    if((getPartonW || doWReweighting )){ 
      if(genz.size()==1){
	float_values["Event_Z_Pt"]= genz.at(0).Pt();
	float_values["Event_Z_Eta"]= genz.at(0).Eta();
	float_values["Event_Z_Phi"]= genz.at(0).Phi();
	float_values["Event_Z_E"]= genz.at(0).Energy();
	float_values["Event_Z_Mass"]= genz.at(0).M();	
	
	double ptW = genz.at(0).Pt();
	double wweight = getZPtWeight(ptW);			
	cout << wweight << endl;
	float_values["Event_Z_QCD_Weight"]= wweight;
      }
      else (float_values["Event_Z_QCD_Weight"]=1.0);
    }
    //cout << "getPartonW is \t" << getPartonW << "and doWReweighting is\t" << doWReweighting << endl;
    if((getPartonW || doWReweighting ) ) {
      if(gena.size()==1){       
      float_values["Event_a_Pt"]= gena.at(0).Pt();
      float_values["Event_a_Eta"]= gena.at(0).Eta();
      float_values["Event_a_Phi"]= gena.at(0).Phi();
      float_values["Event_a_E"]= gena.at(0).Energy();
      float_values["Event_a_Mass"]= gena.at(0).M();	
      
      double ptW = gena.at(0).Pt();
      double wweight = getAPtWeight(ptW);			
      //cout << wweight << endl;
      float_values["Event_a_Weight"]= wweight;
      }
      else (float_values["Event_a_Weight"]=1.0);
    }
    if( (getPartonTop || doTopReweighting)) {
      if (gentop.size()==1 && genantitop.size()==1 && getPartonTop){
	double ptT = gentop.at(0).Pt();
	double ptTbar = genantitop.at(0).Pt();
	double tweight = getTopPtWeight(ptT,ptTbar);			
	double tweightext = getTopPtWeight(ptT,ptTbar,true);			
	float_values["Event_T_Weight"]= tweight;
	float_values["Event_T_Ext_Weight"]= tweightext;
	//cout << tweight << endl;
      }
      else {(float_values["Event_T_Weight"]=1.0);
	(float_values["Event_T_Ext_Weight"]=1.0);}
    }
    if(useLHEWeights){
      getEventLHEWeights();
    }
    //  }
  //if(1>0){
    //G
    //cout << "starting gen"<<endl;
    //if(!isData) {//G
    if(0>1){ 
    iEvent.getByToken(t_genParticleCollection_ , genParticles);
    
      //const bool print = false;
      int momId=-999;
      int momStatus=-999;
  
      std::vector< int > MomId;
      std::vector< int > MomStatus;
      for(size_t i=0; i<genParticles->size(); ++i) {
	  const reco::GenParticle& p = (*genParticles)[i];
	  momId = p.numberOfMothers() ? p.mother()->pdgId() : 0;                                                                                                                            
	  momStatus = p.numberOfMothers() ? p.mother()->status() : 0;    
	  MomId.push_back(momId);
	  MomStatus.push_back(momStatus);
      }

      cout << "starting gen 2" << endl;
      for(size_t i=0;/*i<genParticles->size()*/(size_t)max_instances[gen_label]/*max((int)genParticles->size(),(int)30)*/;++i){
	string nameshortv= "genPart";
	string pref = obj_to_pref[gen_label];

	const reco::GenParticle& p = (*genParticles)[i];
	cout << "genPart" << endl;
	/*vfloats_values[makeName(gen_label,pref,"Pt")][i]=(float)((p.p4()).Pt());
	vfloats_values[makeName(gen_label,pref,"Eta")][i]=(float)((p.p4()).Eta());
	vfloats_values[makeName(gen_label,pref,"Phi")][i]=(float)((p.p4()).Phi());
	vfloats_values[makeName(gen_label,pref,"E")][i]=(float)((p.p4()).E());
	vfloats_values[makeName(gen_label,pref,"Status")][i]=(int)(p.status());
	vfloats_values[makeName(gen_label,pref,"Id")][i]=(int)(p.pdgId());
	*/
	cout << "genPart mom" << endl;
	int momId=-999;
	//int momStatus =-999;
	/*momId = p.numberOfMothers() ? p.mother()->pdgId() : 0;
	momStatus = p.numberOfMothers() ? p.mother()->status() : 0;

	vfloats_values[makeName(gen_label,pref,"Mom0Id")][i]=(float)(momId);
	vfloats_values[makeName(gen_label,pref,"Mom0Status")][i]=(float)(momStatus);
	*/
	/*for(size_t j = 0, n=p.numberOfDaughters(); j<n; ++j) {
	  vfloats_values[makeName(gen_label,pref,"dauId1")][i]=(float)(p.daughter(j)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus1")][i]=(float)(p.daughter(j)->status());
	  vfloats_values[makeName(gen_label,pref,"dauId2")][i]=(float)(p.daughter(j+1)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus2")][i]=(float)(p.daughter(j+1)->status());
	  vfloats_values[makeName(gen_label,pref,"dauId3")][i]=(float)(p.daughter(j+1)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus3")][i]=(float)(p.daughter(j+1)->status());
	  }*/
	cout << "ewk" << endl;
	if(p.pdgId()==momId && isEWKID(p.pdgId()) && getParticleWZ){
	  //cout << "inside gen loop" << endl;
	  TLorentzVector vec;
	  vec.SetPtEtaPhiE((p.p4()).Pt(),(p.p4()).Eta(),(p.p4()).Phi(), (p.p4()).E());
	  if(p.pdgId()==11  && p.status()==1){
	    pare.push_back(vec);	  }
	  if(p.pdgId()==-11  && p.status()==1){
	    parebar.push_back(vec);	  }
	  if(p.pdgId()==13 && p.status()==1){
	    parmu.push_back(vec);	  }
	  if(p.pdgId()==-13 && p.status()==1){
	    parmubar.push_back(vec);	  }

	  if(p.pdgId()==15  && p.status()==2){
	    partau.push_back(vec);  }
	  if(p.pdgId()==-15  && p.status()==2){
	    partaubar.push_back(vec);  }
	  
	  if(p.pdgId()==12  && p.status()==1){
	    parnue.push_back(vec);	  }
	  if(p.pdgId()==14  && p.status()==1){
	    parnumu.push_back(vec);	  }
	  if(p.pdgId()==16  && p.status()==1){
	    parnutau.push_back(vec);	  }

	  if(p.pdgId()==-12  && p.status()==1){
	    parnuebar.push_back(vec);	  }
	  if(p.pdgId()==-14  && p.status()==1){
	    parnumubar.push_back(vec);	  }
	  if(p.pdgId()==-16  && p.status()==1){
	    parnutaubar.push_back(vec);	  }
	  
	  if(abs(p.pdgId())==23 && p.status()==62){
	    double wweight = getZEWKPtWeight((p.p4()).Pt());
	    float_values["Event_Z_EW_Weight"]= wweight;
	    //    cout << "w weight:\t" << wweight << endl;
	  }
	  if(abs(p.pdgId())==24 && p.status()==62){
	    double wweight = getWEWKPtWeight((p.p4()).Pt());
	    cout << "third if ok"<< endl;
	    float_values["Event_W_EW_Weight"]= wweight;
	     cout << "z weight:\t" << wweight << endl;
	     }   
	  // cout << "closing ewk loop"<< endl;
	}
	//cout << "closing genloop" << endl;
      }
      //cout << "starting Z" << endl;
      //Z
      if(parmu.size()>0&& parmubar.size()>0){parz.push_back(parmu.at(0)+parmubar.at(0)) ;}
      if(pare.size()>0&& parebar.size()>0){parz.push_back(pare.at(0)+parebar.at(0)) ;}
      if(partau.size()>0&& partaubar.size()>0){parz.push_back(partau.at(0)+partaubar.at(0)) ;}
      if(parnumu.size()>0&& parnumubar.size()>0){parz.push_back(parnumu.at(0)+parnumubar.at(0)) ;}
      if(parnue.size()>0&& parnuebar.size()>0){parz.push_back(parnue.at(0)+parnuebar.at(0)) ;}
      if(parnutau.size()>0&& parnutaubar.size()>0){parz.push_back(parnutau.at(0)+parnutaubar.at(0)) ;}
      if(   float_values["Event_Z_EW_Weight"] ==1 &&parz.size()>0 )    float_values["Event_Z_EW_Weight"]= getZEWKPtWeight(parz.at(0).Pt());
      if(   float_values["Event_Z_QCD_Weight"] ==1 &&parz.size()>0 )    float_values["Event_Z_QCD_Weight"]= getWPtWeight(parz.at(0).Pt());

      //W
      //      cout << "starting W" << endl;
      if(parmu.size()>0&& parnumubar.size()>0){parw.push_back(parmu.at(0)+parnumubar.at(0)) ;}

      if(pare.size()>0&& parnuebar.size()>0){parw.push_back(pare.at(0)+parnuebar.at(0)) ;}

      if(partau.size()>0&& parnutaubar.size()>0){parw.push_back(partau.at(0)+parnutaubar.at(0)) ;}
      
      if(parnumu.size()>0&& parmubar.size()>0){parw.push_back(parnumu.at(0)+parmubar.at(0)) ;}

      if(parnue.size()>0&& parebar.size()>0){parw.push_back(parnue.at(0)+parebar.at(0)) ;}


      if(parnutau.size()>0&& partaubar.size()>0){parw.push_back(parnutau.at(0)+partaubar.at(0)) ;}

      if(   float_values["Event_W_EW_Weight"] ==1 &&parw.size()>0 )    {
	float_values["Event_W_EW_Weight"]= getWEWKPtWeight(parw.at(0).Pt());
      }
      if(   float_values["Event_W_QCD_Weight"] ==1 &&parw.size()>0 )    float_values["Event_W_QCD_Weight"]= getWPtWeight(parw.at(0).Pt());
    }
    float_values["Event_W_Weight"]= float_values["Event_W_EW_Weight"]*float_values["Event_W_QCD_Weight"];
    float_values["Event_Z_Weight"]= float_values["Event_Z_EW_Weight"]*float_values["Event_Z_QCD_Weight"];   
    //cout << "ending mega if" << endl;
  }
  trees["WeightHistory"]->Fill();

  //Part 3: filling the additional variables

  //  cout << "starting boosted and trigger" << endl;
  //boosted part
  iEvent.getByToken(jetAK8vSubjetIndex0, ak8jetvSubjetIndex0);
  iEvent.getByToken(jetAK8vSubjetIndex1, ak8jetvSubjetIndex1);
 
  //Part 0: trigger preselection
 if(useTriggers){
    //Trigger names are retrieved from the run tree
    iEvent.getByToken(t_triggerBits_,triggerBits );
    iEvent.getByToken(t_triggerPrescales_,triggerPrescales );
    bool triggerOr = getEventTriggers();
    
    if(isFirstEvent){
      for(size_t bt = 0; bt < triggerNamesR->size();++bt){
	std::string tname = triggerNamesR->at(bt);
	cout << "trigger test tname "<< tname << " passes "<< triggerBits->at(bt)<< endl;
      }
    }
    
    if(cutOnTriggers && !triggerOr) return;
  }

  if(useMETFilters){
   iEvent.getByToken(t_metBits_,metBits );
    //    iEvent.getByToken(t_metNames_,metNames );
    if(isFirstEvent){
      for(size_t bt = 0; bt < metNames->size();++bt){
	std::string tname = metNames->at(bt);
	//cout << "test tname "<< tname << " passes "<< metBits->at(bt)<< endl;
      }
    }
    getMETFilters();
  }
  
  if(changeJECs || recalculateEA){
    iEvent.getByToken(t_Rho_ ,rho);
    Rho = *rho; 
  }
  float_values["Event_Rho"] = (double)Rho;
  if(isFirstEvent){
    isFirstEvent = false;
  }
    
  if( addPV || changeJECs){
    iEvent.getByToken(t_pvZ_,pvZ);
    iEvent.getByToken(t_pvChi2_,pvChi2);
    iEvent.getByToken(t_pvNdof_,pvNdof);
    iEvent.getByToken(t_pvRho_,pvRho);
    nPV = pvZ->size();
  }
  
  //  cout << "starting vectors" << endl;
  //Part 1 taking the obs values from the edm file
  for (;itPsets!=physObjects.end();++itPsets){ 
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI"); 
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();
    string namelabel = itPsets->getParameter< string >("label");
    string nameprefix = itPsets->getParameter< string >("prefix");
    size_t maxInstance=(size_t)max_instances[namelabel];


    variablesDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesD"); 
    singleDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("singleD"); 

    std::vector<edm::InputTag >::const_iterator itD = variablesDouble.begin();
    std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
    
  
    //Vectors of floats
    for (;itF != variablesFloat.end();++itF){
      string varname=itF->instance();
      
      string name = makeBranchName(namelabel,nameprefix,varname);
      
      //string namelabel;
      float tmp =1.0;
      iEvent.getByToken(t_floats[name] ,h_floats[name]);
      //iEvent.getByLabel(*(itF),h_floats[name]);
      //      cout << "name "<< name <<endl;
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_floats[name]->size()){tmp = h_floats[name]->at(fi);}
	else { tmp = -9999.; }
	//	cout << " setting name "<< name<< " at instance "<< fi <<" to value "<< tmp <<endl;
	vfloats_values[name][fi]=tmp;
	for (size_t sc=0;sc< obj_cats[namelabel].size();++sc){
	  string category = obj_cats[namelabel].at(sc);
	  string namecat = makeBranchNameCat(namelabel,category,nameprefix,varname);
	  vfloats_values[namecat][fi]=-9999.;
	}
      }
      sizes[namelabel]=h_floats[name]->size();
    }

    for (;itD != variablesDouble.end();++itD){
      string varname=itD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      //string namelabel;
      float tmp =1.0;
      iEvent.getByToken(t_doubles[name] ,h_doubles[name]);
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_doubles[name]->size()){tmp = h_doubles[name]->at(fi);}
	else { tmp = -9999.; }
	vdoubles_values[name][fi]=tmp;
      }
      sizes[namelabel]=h_doubles[name]->size();
    }
    
    //Vectors of ints
    for (;itI != variablesInt.end();++itI){
      string varname=itI->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      int tmp = 1;
      iEvent.getByToken(t_ints[name] ,h_ints[name]);
      //iEvent.getByLabel(*(itI),h_ints[name]);
      for (size_t fi = 0;fi < maxInstance;++fi){
	if(fi <h_ints[name]->size()){tmp = h_ints[name]->at(fi);}
	else { tmp = -9999.; }
	vints_values[name][fi]=tmp;
      }
    }  
    
    //    std::cout << " checkpoint ints"<<endl;
    //Single floats/ints
    for (;itsF != singleFloat.end();++itsF){
      string varname=itsF->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_float[name],h_float[name]);
      //iEvent.getByLabel(*(itsF),h_float[name]);
      float_values[name]=*h_float[name];
    }

    for (;itsD != singleDouble.end();++itsD){
      string varname=itsD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_double[name] ,h_double[name]);
      //iEvent.getByLabel(*(itsD),h_double[name]);
      double_values[name]=*h_double[name];
    }
    for (;itsI != singleInt.end();++itsI){
      string varname=itsI->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_int[name],h_int[name]);
      //iEvent.getByLabel(*(itsI),h_int[name]);
      int_values[name]=*h_int[name];
    }
    //    std::cout << " checkpoint singles"<<endl;
  }

  //  std::cout << " checkpoint part 1"<<endl;


  //Part 2: selection and analysis-level changes
  //This might change for each particular systematics, 
  //e.g. for each jet energy scale variation, for MET or even lepton energy scale variations


  vector<TLorentzVector> photons;
  vector<TLorentzVector> electrons;
  vector<TLorentzVector> muons;
  vector<TLorentzVector> leptons;
  vector<TLorentzVector> loosemuons;
  vector<TLorentzVector> jets;
  vector<TLorentzVector> jetsnob;
  vector<TLorentzVector> bjets;
  vector<TLorentzVector> subjets;
  vector<TLorentzVector> topjets;
  vector<TLorentzVector> type1topjets;
  vector<TLorentzVector> type2topjets;
  vector<TLorentzVector> resolvedtops;

  vector<float> leptonsCharge;

  vector<int> flavors;
  //  cout << "starting loop systematics" << endl;

  for (size_t s = 0; s< systematics.size();++s){
    
    int nb=0,nc=0,nudsg=0;

    int ncsvl_tags=0,ncsvt_tags=0,ncsvm_tags=0;//, njets_tottag=0;
    int ncsvl_subj_tags=0,ncsvm_subj_tags=0;
    getEventTriggers();

    photons.clear();
    leptons.clear();
    electrons.clear();
    muons.clear();
    loosemuons.clear();
    jets.clear();
    jetsnob.clear();
    bjets.clear();
    type2topjets.clear();
    type1topjets.clear();
    topjets.clear();
    resolvedtops.clear();
    string syst = systematics.at(s);
    nTightJets=0;

    jsfscsvt.clear();
    jsfscsvt_b_tag_up.clear(); 
    jsfscsvt_b_tag_down.clear(); 
    jsfscsvt_mistag_up.clear(); 
    jsfscsvt_mistag_down.clear();

    jsfscsvm.clear(); 
    jsfscsvm_b_tag_up.clear(); 
    jsfscsvm_b_tag_down.clear(); 
    jsfscsvm_mistag_up.clear(); 
    jsfscsvm_mistag_down.clear();
    
    jsfscsvl.clear(); 
    jsfscsvl_b_tag_up.clear(); 
    jsfscsvl_b_tag_down.clear(); 
    jsfscsvl_mistag_up.clear();
    jsfscsvl_mistag_down.clear();

    jsfscsvm_subj.clear(); 
    jsfscsvm_subj_b_tag_up.clear(); 
    jsfscsvm_subj_b_tag_down.clear(); 
    jsfscsvm_subj_mistag_up.clear(); 
    jsfscsvm_subj_mistag_down.clear();
    
    jsfscsvl_subj.clear(); 
    jsfscsvl_subj_b_tag_up.clear(); 
    jsfscsvl_subj_b_tag_down.clear(); 
    jsfscsvl_subj_mistag_up.clear();
    jsfscsvl_subj_mistag_down.clear();


    //---------------- Soureek Adding PU Info ------------------------------
    //if(doPU_){
    //  iEvent.getByToken(t_ntrpu_,ntrpu);
    //  nTruePU=*ntrpu;
    //  getPUSF();
    //}
    
    int mapBJets[20], mapEle[20], mapMu[20];
    for(int i = 0; i<20;++i){
      mapBJets[i]=-1; mapEle[i]=-1; mapMu[i]=-1;
    } 
    int lepidx=0;
    int bjetidx=0;
    //cout << "starting photons" << endl;
    //Photons
    for(int ph = 0;ph < max_instances[photon_label] ;++ph){
      string pref = obj_to_pref[photon_label];
      
      
      float pt = vfloats_values[makeName(photon_label,pref,"Pt")][ph];
      float eta = vfloats_values[makeName(photon_label,pref,"Eta")][ph];      
      
      float sieie = vfloats_values[makeName(photon_label,pref,"SigmaIEtaIEta")][ph];      
      float hoe = vfloats_values[makeName(photon_label,pref,"HoverE")][ph];      
      
      float abseta = fabs(eta);

      float pho_isoC  = vfloats_values[makeName(photon_label,pref,"ChargedHadronIso")][ph];      
      float pho_isoP  = vfloats_values[makeName(photon_label,pref,"NeutralHadronIso")][ph];      
      float pho_isoN     =  vfloats_values[makeName(photon_label,pref,"PhotonIso")][ph];      

      float pho_isoCea  = vfloats_values[makeName(photon_label,pref,"ChargedHadronIsoEAcorrected")][ph];      
      float pho_isoPea  = vfloats_values[makeName(photon_label,pref,"PhotonIsoEAcorrected")][ph];      
      float pho_isoNea     =  vfloats_values[makeName(photon_label,pref,"NeutralHadronIsoEAcorrected")][ph];      

      if(recalculateEA){
	pho_isoCea     = std::max( double(0.0) ,(pho_isoC - Rho*getEffectiveArea("ch_hadrons",abseta)));
	pho_isoPea     = std::max( double(0.0) ,(pho_isoP - Rho*getEffectiveArea("photons",abseta)));
	pho_isoNea     = std::max( double(0.0) ,(pho_isoN - Rho*getEffectiveArea("neu_hadrons",abseta)));
      }
      
      bool isBarrel = (abseta<1.479);
      bool isEndcap = (abseta>1.479 && abseta < 2.5);

      vfloats_values[photon_label+"_isLooseSpring15"][ph]=0.0;
      vfloats_values[photon_label+"_isMediumSpring15"][ph]=0.0;
      vfloats_values[photon_label+"_isTightSpring15"][ph]=0.0;
 
      if(isBarrel){

	if( sieie < 0.0103 &&   hoe < 0.05 &&   pho_isoCea < 2.44 &&   pho_isoNea < (2.57+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 1.31 &&   pho_isoNea < (0.60+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 0.91 &&   pho_isoNea < (0.33+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
      if(isEndcap){
	if( sieie < 0.0277 &&   hoe < 0.05 &&   pho_isoCea < 1.84 &&   pho_isoNea < (4.00+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 1.25 &&   pho_isoNea < (1.65+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 0.65 &&   pho_isoNea < (0.93+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
      //WORKING POINTS TO BE CHECKED
    
    }
    //G                                                                                                                                              
    //cout << "starting muons"<<endl;
    //    cout << "starting muons" << endl;
    //Muons
    for(int mu = 0;mu < max_instances[mu_label] ;++mu){
      string pref = obj_to_pref[mu_label];
      float isTight = vfloats_values[makeName(mu_label,pref,"IsTightMuon")][mu];
      float isLoose = vfloats_values[makeName(mu_label,pref,"IsLooseMuon")][mu];
      float isMedium = vfloats_values[makeName(mu_label,pref,"IsMediumMuon")][mu];
      float isSoft = vfloats_values[makeName(mu_label,pref,"IsSoftMuon")][mu];

      float pt = vfloats_values[makeName(mu_label,pref,"Pt")][mu];
      float eta = vfloats_values[makeName(mu_label,pref,"Eta")][mu];
      float phi = vfloats_values[makeName(mu_label,pref,"Phi")][mu];
      float energy = vfloats_values[makeName(mu_label,pref,"E")][mu];
      float iso = vfloats_values[makeName(mu_label,pref,"Iso04")][mu];
      
      float muCharge = vfloats_values[makeName(mu_label,pref,"Charge")][mu];

      if(isMedium>0 && pt> 30 && abs(eta) < 2.1 && iso <0.25){ 
	++float_values["Event_nMediumMuons"];
	TLorentzVector muon;
	muon.SetPtEtaPhiE(pt, eta, phi, energy);
	muons.push_back(muon);
	leptons.push_back(muon);
	flavors.push_back(13);
	
	leptonsCharge.push_back(muCharge);
	
	mapMu[lepidx]=mu; 
	if(isInVector(obj_cats[mu_label],"Medium")){
	  fillCategory(mu_label,"Medium",mu,float_values["Event_nMediumMuons"]-1);
	}
	++lepidx;
      }
      
      if(isInVector(obj_cats[mu_label],"Medium")){
	sizes[mu_label+"Medium"]=(int)float_values["Event_nMediumMuons"];
      }
      
      if(isLoose>0 && pt> 10 && abs(eta) < 2.4 && iso<0.25){
	if(isInVector(obj_cats[mu_label],"Loose")){
	  ++float_values["Event_nLooseMuons"];
	  if(isInVector(obj_cats[mu_label],"Loose")){
	    fillCategory(mu_label,"Loose",mu,float_values["Event_nLooseMuons"]-1);
	  }
	}
      }
      if(isInVector(obj_cats[mu_label],"Loose")){
	sizes[mu_label+"Loose"]=(int)float_values["Event_nLooseMuons"];
      }


      if(isTight>0 && pt> 10 && abs(eta) < 2.4 && iso<0.25){
        if(isInVector(obj_cats[mu_label],"Tight")){
          ++float_values["Event_nTightMuons"];
          if(isInVector(obj_cats[mu_label],"Tight")){
            fillCategory(mu_label,"Tight",mu,float_values["Event_nTightMuons"]-1);
          }
        }
      }
      if(isInVector(obj_cats[mu_label],"Tight")){
        sizes[mu_label+"Tight"]=(int)float_values["Event_nTightMuons"];
      }
      
      if(isSoft>0 && pt> 30 && abs(eta) < 2.4){
	++float_values["Event_nSoftMuons"]; 
      }
      if(isLoose>0 && pt > 15){
	TLorentzVector muon;
	muon.SetPtEtaPhiE(pt, eta, phi, energy);
	loosemuons.push_back(muon);
      }
    }
    //G                                                                                                                                              
    //cout << "starting electrons"<<endl;
    //    cout << "starting electrons" << endl;
    //Electrons:
    for(int el = 0;el < max_instances[ele_label] ;++el){
      string pref = obj_to_pref[ele_label];
      float pt = vfloats_values[makeName(ele_label,pref,"Pt")][el];
      float isTight = vfloats_values[makeName(ele_label,pref,"isTight")][el];
      float isLoose = vfloats_values[makeName(ele_label,pref,"isLoose")][el];
      float isMedium = vfloats_values[makeName(ele_label,pref,"isMedium")][el];
      float isVeto = vfloats_values[makeName(ele_label,pref,"isVeto")][el];

      isTight = vfloats_values[makeName(ele_label,pref,"vidTight")][el];
      isLoose = vfloats_values[makeName(ele_label,pref,"vidLoose")][el];
      isMedium = vfloats_values[makeName(ele_label,pref,"vidMedium")][el];
      isVeto = vfloats_values[makeName(ele_label,pref,"vidVeto")][el];
      float eta = vfloats_values[makeName(ele_label,pref,"Eta")][el];
      float scEta = vfloats_values[makeName(ele_label,pref,"scEta")][el];
      float phi = vfloats_values[makeName(ele_label,pref,"Phi")][el];
      float energy = vfloats_values[makeName(ele_label,pref,"E")][el];      
      float iso = vfloats_values[makeName(ele_label,pref,"Iso03")][el];

      float elCharge = vfloats_values[makeName(ele_label,pref,"Charge")][el];

      bool passesDRmu = true;
      bool passesTightCuts = false;
      if(fabs(scEta)<=1.479){
	passesTightCuts = isTight >0.0 && iso < 0.0588 ;
      } //is barrel electron
      if (fabs(scEta)>1.479){
	passesTightCuts = isTight >0.0 && iso < 0.0571 ;
      }

      if(pt> 30 && fabs(eta) < 2.1 && passesTightCuts){
	TLorentzVector ele;
	ele.SetPtEtaPhiE(pt, eta, phi, energy);	
	double minDR=999;
	double deltaR = 999;
	for (size_t m = 0; m < (size_t)loosemuons.size(); ++m){
	  deltaR = ele.DeltaR(loosemuons[m]);
	  minDR = min(minDR, deltaR);
	}
	if(!loosemuons.size()) minDR=999;
	if(minDR>0.1){ 
	  electrons.push_back(ele); 
	  flavors.push_back(11);
	  leptons.push_back(ele);
	  
	  leptonsCharge.push_back(elCharge);

	  ++float_values["Event_nTightElectrons"];
	  mapEle[lepidx]=el;
	  ++lepidx;
	  if(isInVector(obj_cats[ele_label],"Tight")){
	    fillCategory(ele_label,"Tight",el,float_values["Event_nTightElectrons"]-1);
	  }
	}
	else {passesDRmu = false;}
      }
      if(isInVector(obj_cats[ele_label],"Tight")){
	sizes[ele_label+"Tight"]=(int)float_values["Event_nTightElectrons"];
      }

      if(isLoose>0 && pt> 30 && fabs(eta) < 2.5){
	++float_values["Event_nLooseElectrons"];

      }

      if(isMedium>0 && pt> 30 && fabs(eta) < 2.5){
	++float_values["Event_nMediumElectrons"]; 
      }
      
      if(isVeto>0 && pt> 10 && fabs(eta) < 2.5 ){
	if((fabs(scEta)<=1.479 && (iso<0.175)) 
	   || ((fabs(scEta)>1.479) && (iso<0.159))){
	  ++float_values["Event_nVetoElectrons"]; 
	  if(isInVector(obj_cats[ele_label],"Veto")){
	    fillCategory(ele_label,"Veto",el,float_values["Event_nVetoElectrons"]-1);
	  }
	}
      }
      if(isInVector(obj_cats[ele_label],"Veto")){
	sizes[ele_label+"Veto"]=(int)float_values["Event_nVetoElectrons"];
      }
      
      vfloats_values[ele_label+"_PassesDRmu"][el]=(float)passesDRmu;
    } 
    int firstidx=-1, secondidx=-1;
    double maxpt=0.0;

    for(size_t l =0; l< leptons.size();++l){
      //double maxpt= 0.0;
      double lpt= leptons.at(l).Pt();
      if(lpt>maxpt){maxpt = lpt;firstidx=l;}
      
    }

    maxpt=0.0;
    for(size_t l =0; l< leptons.size();++l){
      //double maxpt= 0.0;
      double lpt= leptons.at(l).Pt();
      if(lpt>maxpt&&firstidx!=(int)l){maxpt = lpt;secondidx=l;}
    }
    if(firstidx>-1){
      float_values["Event_Lepton1_Pt"]=leptons.at(firstidx).Pt(); 
      float_values["Event_Lepton1_Phi"]=leptons.at(firstidx).Phi(); 
      float_values["Event_Lepton1_Eta"]=leptons.at(firstidx).Eta(); 
      float_values["Event_Lepton1_E"]=leptons.at(firstidx).E(); 
      float_values["Event_Lepton1_Flavour"]=flavors.at(firstidx);

      float_values["Event_Lepton1_Charge"]=leptonsCharge.at(firstidx);

    }
    if(secondidx>-1){
      float_values["Event_Lepton2_Pt"]=leptons.at(secondidx).Pt(); 
      float_values["Event_Lepton2_Phi"]=leptons.at(secondidx).Phi(); 
      float_values["Event_Lepton2_Eta"]=leptons.at(secondidx).Eta(); 
      float_values["Event_Lepton2_E"]=leptons.at(secondidx).E(); 
      float_values["Event_Lepton2_Flavour"]=flavors.at(secondidx);

      float_values["Event_Lepton2_Charge"]=leptonsCharge.at(secondidx);

    }
    //G                                                                                                                                              
    //cout << "starting jets"<<endl;
    //    cout << "starting jets" << endl;
    //Jets:
    double Ht=0;
    double corrMetPx =0;
    double corrMetPy =0;
    double DUnclusteredMETPx=0.0;
    double DUnclusteredMETPy=0.0;

    string prefm = obj_to_pref[met_label];
    float metZeroCorrPt = vfloats_values[makeName(met_label,prefm,"UncorrPt")][0];
    float metZeroCorrPhi = vfloats_values[makeName(met_label,prefm,"UncorrPhi")][0];
    float metZeroCorrY = metZeroCorrPt*sin(metZeroCorrPhi);
    float metZeroCorrX = metZeroCorrPt*cos(metZeroCorrPhi);

    for(int j = 0;j < max_instances[jets_label] ;++j){
      string pref = obj_to_pref[jets_label];
      float pt = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float ptzero = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float genpt = vfloats_values[makeName(jets_label,pref,"GenJetPt")][j];
      float eta = vfloats_values[makeName(jets_label,pref,"Eta")][j];
      float phi = vfloats_values[makeName(jets_label,pref,"Phi")][j];
      float energy = vfloats_values[makeName(jets_label,pref,"E")][j];
     
      float ptCorr = -9999;
      float energyCorr = -9999;
      float smearfact = -9999;

      float jecscale = vfloats_values[makeName(jets_label,pref,"jecFactor0")][j];
      float area = vfloats_values[makeName(jets_label,pref,"jetArea")][j];
      
      float juncpt=0.;
      float junce=0.;
      if(pt>0){
	  TLorentzVector jetUncorr;

	  jetUncorr.SetPtEtaPhiE(pt,eta,phi,energy);
	  
	  jetUncorr= jetUncorr*jecscale;
	  DUnclusteredMETPx+=jetUncorr.Pt()*cos(phi);
	  DUnclusteredMETPy+=jetUncorr.Pt()*sin(phi);
	  
	  juncpt=jetUncorr.Perp();
	  junce=jetUncorr.E();

	  if(changeJECs){
	   
	    jecCorr->setJetPhi(jetUncorr.Phi());
	    jecCorr->setJetEta(jetUncorr.Eta());
	    jecCorr->setJetE(jetUncorr.E());
	    jecCorr->setJetPt(jetUncorr.Perp());
	    jecCorr->setJetA(area);
	    jecCorr->setRho(Rho);
	    jecCorr->setNPV(nPV);
	    double recorr =  jecCorr->getCorrection();
	    jetUncorr = jetUncorr *recorr;
	    pt = jetUncorr.Pt();
	    eta = jetUncorr.Eta();
	    energy = jetUncorr.Energy();
	    phi = jetUncorr.Phi();
	  }

	smearfact = smear(pt, genpt, eta, syst);
	ptCorr = pt * smearfact;
	energyCorr = energy * smearfact;
	float unc = jetUncertainty(ptCorr,eta,syst);
	ptCorr = ptCorr * (1 + unc);
	energyCorr = energyCorr * (1 + unc);
	//cout << "pt uncor: " << pt << "pt corr: " << ptCorr << endl;
	corrMetPx -=(cos(phi)*(ptCorr-ptzero));
        corrMetPy -=(sin(phi)*(ptCorr-ptzero));
	
      }
          
      float csv = vfloats_values[makeName(jets_label,pref,"CSVv2")][j];
      float partonFlavour = vfloats_values[makeName(jets_label,pref,"PartonFlavour")][j];
      int flavor = int(partonFlavour);
  
      if(getWZFlavour){
	if(abs(flavor)==5){++nb;}
	else{ 
	  if(abs(flavor)==4){++nc;}
	  else {++nudsg;}
	}
      }
 
      vfloats_values[jets_label+"_NoCorrPt"][j]=juncpt;
      vfloats_values[jets_label+"_NoCorrE"][j]=junce;

      vfloats_values[jets_label+"_CorrPt"][j]=ptCorr;
      vfloats_values[jets_label+"_CorrE"][j]=energyCorr;
      vfloats_values[jets_label+"_CorrEta"][j]=eta;
      vfloats_values[jets_label+"_CorrPhi"][j]=phi;

      bool isCSVT = csv  > 0.9535;
      bool isCSVM = csv  > 0.8484;
      bool isCSVL = csv  > 0.5426;
      vfloats_values[jets_label+"_IsCSVT"][j]=isCSVT;
      vfloats_values[jets_label+"_IsCSVM"][j]=isCSVM;
      vfloats_values[jets_label+"_IsCSVL"][j]=isCSVL;
      
      float bsf = getScaleFactor(ptCorr,eta,partonFlavour,"noSyst");
      float bsfup = getScaleFactor(ptCorr,eta,partonFlavour,"up");
      float bsfdown = getScaleFactor(ptCorr,eta,partonFlavour,"down");
      
      vfloats_values[jets_label+"_BSF"][j]=bsf;
      vfloats_values[jets_label+"_BSFUp"][j]=bsfup;
      vfloats_values[jets_label+"_BSFDown"][j]=bsfdown;
      
      bool passesID = true;
      
      if(!(jecscale*energy > 0))passesID = false;
      else{
        float neuMulti = vfloats_values[makeName(jets_label,pref,"neutralMultiplicity")][j];
        float chMulti = vfloats_values[makeName(jets_label,pref,"chargedMultiplicity")][j];
        float chHadEnFrac = vfloats_values[makeName(jets_label,pref,"chargedHadronEnergyFrac")][j];
        float chEmEnFrac = vfloats_values[makeName(jets_label,pref,"chargedEmEnergyFrac")][j];
        float neuEmEnFrac = vfloats_values[makeName(jets_label,pref,"neutralEmEnergyFrac")][j];
        float neuHadEnFrac = vfloats_values[makeName(jets_label,pref,"neutralHadronEnergyFrac")][j];
	float numConst = chMulti + neuMulti;

        if(fabs(eta)<=2.7){
          passesID =  (neuHadEnFrac<0.90 && neuEmEnFrac<0.90 && numConst>1) && ( (abs(eta)<=2.4 && chHadEnFrac>0 && chMulti>0 && chEmEnFrac<0.99) || abs(eta)>2.4);
        }
	else if(fabs(eta)>2.7 && fabs(eta)<=3.0){
	  passesID = neuMulti <2 && neuEmEnFrac < 0.90;
	}
        else if(fabs(eta)>3){
          passesID = neuEmEnFrac<0.90 && neuMulti>10 ;
        }
      }
      
      vfloats_values[jets_label+"_PassesID"][j]=(float)passesID;
      
      //Remove overlap with tight electrons/muons
      double minDR=9999;
      double minDRThrEl=0.3;
      double minDRThrMu=0.4;
      bool passesDR=true;
 
      for (size_t e = 0; e < (size_t)electrons.size(); ++e){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(electrons.at(e).Pt(),electrons.at(e).Eta(),electrons.at(e).Phi(),electrons.at(e).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrEl)passesDR = false;
      }
      for (size_t m = 0; m < (size_t)muons.size(); ++m){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(muons.at(m).Pt(),muons.at(m).Eta(),muons.at(m).Phi(),muons.at(m).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	//minDR = min(minDR,deltaR(muons.at(m) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrMu)passesDR = false;
      }
      
      vfloats_values[jets_label+"_MinDR"][j]=minDR;
      vfloats_values[jets_label+"_PassesDR"][j]=(float)passesDR;
      
      vfloats_values[jets_label+"_IsTight"][j]=0.0;
      vfloats_values[jets_label+"_IsLoose"][j]=0.0;
	
      if(passesID && passesDR) vfloats_values[jets_label+"_IsLoose"][j]=1.0;

      if(passesID && passesDR && pt>50 && abs(eta)<2.4){
	Ht+=pt;
      }

      float_values["Event_Ht"] = (float)Ht;
      
      for (size_t ji = 0; ji < (size_t)jetScanCuts.size(); ++ji){
	stringstream j_n;
	double jetval = jetScanCuts.at(ji);
	j_n << "Cut" <<jetval;
	bool passesCut = ( ptCorr > jetval && fabs(eta) < 4.);

	if(!passesID || !passesCut || !passesDR) continue;
	if(ji==0){
	  vfloats_values[jets_label+"_IsTight"][j]=1.0;
	  TLorentzVector jet;
	  jet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	  jets.push_back(jet);

	  if(!isCSVM)     jetsnob.push_back(jet);

	}

	if(passesCut &&  passesID && passesDR){	
	  float_values["Event_nJets"+j_n.str()]+=1;if(ji==0){
	    nTightJets+=1;
	    if(isInVector(obj_cats[jets_label],"Tight")){
	      fillCategory(jets_label,"Tight",j,float_values["Event_nJets"+j_n.str()]-1);
	    }
	  }
	}

	if(passesCut &&  passesID && passesDR){
	  double csvteff = MCTagEfficiency("csvt",flavor, ptCorr, eta);
	  double sfcsvt = TagScaleFactor("csvt", flavor, "noSyst", ptCorr);
	  
	  double csvleff = MCTagEfficiency("csvl",flavor,ptCorr, eta);
	  double sfcsvl = TagScaleFactor("csvl", flavor, "noSyst", ptCorr);

	  double csvmeff = MCTagEfficiency("csvm",flavor,ptCorr, eta);
	  double sfcsvm = TagScaleFactor("csvm", flavor, "noSyst", ptCorr);
	  //cout << "csvmeff: "  << csvmeff << " sfcsvm: " << sfcsvm << endl; 

	  double sfcsvt_mistag_up = TagScaleFactor("csvt", flavor, "mistag_up", ptCorr);
	  double sfcsvl_mistag_up = TagScaleFactor("csvl", flavor, "mistag_up", ptCorr);
	  double sfcsvm_mistag_up = TagScaleFactor("csvm", flavor, "mistag_up", ptCorr);

	  double sfcsvt_mistag_down = TagScaleFactor("csvt", flavor, "mistag_down", ptCorr);
	  double sfcsvl_mistag_down = TagScaleFactor("csvl", flavor, "mistag_down", ptCorr);
	  double sfcsvm_mistag_down = TagScaleFactor("csvm", flavor, "mistag_down", ptCorr);

	  double sfcsvt_b_tag_down = TagScaleFactor("csvt", flavor, "b_tag_down", ptCorr);
	  double sfcsvl_b_tag_down = TagScaleFactor("csvl", flavor, "b_tag_down", ptCorr);
	  double sfcsvm_b_tag_down = TagScaleFactor("csvm", flavor, "b_tag_down", ptCorr);
	  
	  double sfcsvt_b_tag_up = TagScaleFactor("csvt", flavor, "b_tag_up", ptCorr);
	  double sfcsvl_b_tag_up = TagScaleFactor("csvl", flavor, "b_tag_up", ptCorr);
	  double sfcsvm_b_tag_up = TagScaleFactor("csvm", flavor, "b_tag_up", ptCorr);


	  jsfscsvt.push_back(BTagWeight::JetInfo(csvteff, sfcsvt));
	  jsfscsvl.push_back(BTagWeight::JetInfo(csvleff, sfcsvl));
	  jsfscsvm.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm));

	  jsfscsvt_mistag_up.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_mistag_up));
	  jsfscsvl_mistag_up.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_mistag_up));
	  jsfscsvm_mistag_up.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_mistag_up));

	  jsfscsvt_b_tag_up.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_b_tag_up));
	  jsfscsvl_b_tag_up.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_b_tag_up));
	  jsfscsvm_b_tag_up.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_b_tag_up));

	  jsfscsvt_mistag_down.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_mistag_down));
	  jsfscsvl_mistag_down.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_mistag_down));
	  jsfscsvm_mistag_down.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_mistag_down));

	  jsfscsvt_b_tag_down.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_b_tag_down));
	  jsfscsvl_b_tag_down.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_b_tag_down));
	  jsfscsvm_b_tag_down.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_b_tag_down));

	  //cout<< "jet :# "<< j << "sf is " <<sfcsvt << " eff is "<<csvteff<<endl;
	}
	
	if(isCSVT && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) {
	  float_values["Event_nCSVTJets"+j_n.str()]+=1.0;
	  ncsvt_tags +=1;
	}

	if(isCSVL && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) { 
	  ncsvl_tags +=1;
	}
	if(isCSVM && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) { 
	  float_values["Event_nCSVMJets"+j_n.str()]+=1.0;
	  if(ji==0){
	    ncsvm_tags +=1;
	    //	    cout << " jet "<< j<< " isBJET "<<endl;
	    TLorentzVector bjet;
	    bjet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	    bjets.push_back(bjet);
	    mapBJets[bjetidx]=j;
	    ++bjetidx;
	  }
	}
	
	if(isCSVL && passesCut &&  passesID && passesDR && abs(eta) < 2.4) float_values["Event_nCSVLJets"+j_n.str()]+=1;
	
      }
    }
 
    if(isInVector(obj_cats[jets_label],"Tight")){
      sizes[jets_label+"Tight"]=(int)nTightJets;
    }
    
    float_values["Event_eventFlavour"]=eventFlavour(getWZFlavour, nb, nc, nudsg);
    if(syst.find("unclusteredMet")!= std::string::npos ){
      
      DUnclusteredMETPx=metZeroCorrX+DUnclusteredMETPx;
      DUnclusteredMETPy=metZeroCorrY+DUnclusteredMETPy;
      
      double signmet = 1.0; 
      if(syst.find("down")!=std::string::npos) signmet=-1.0;
      corrMetPx -=signmet*DUnclusteredMETPx*0.1;
      corrMetPy -=signmet*DUnclusteredMETPy*0.1;
    }
    //G                                                                                                                                              
    //cout << "starting met"<<endl;
    //    cout << "starting met" << endl;
    //Met and mt
    string pref = obj_to_pref[met_label];
    float metpt = vfloats_values[makeName(met_label,pref,"Pt")][0];
    float metphi = vfloats_values[makeName(met_label,pref,"Phi")][0];

    
    float metPy = metpt*sin(metphi);
    float metPx = metpt*cos(metphi);
    metPx+=corrMetPx; metPy+=corrMetPy; // add JEC/JER contribution
    
    //Correcting the pt
    float metptCorr = sqrt(metPx*metPx + metPy*metPy);
    vfloats_values[met_label+"_CorrPt"][0]=metptCorr;

    //Correcting the phi
    float metphiCorr = metphi;
    if(metPx<0){
      if(metPy>0)metphiCorr = atan(metPy/metPx)+3.141592;
      if(metPy<0)metphiCorr = atan(metPy/metPx)-3.141592;
    }
    else  metphiCorr = (atan(metPy/metPx));

    vfloats_values[met_label+"_CorrPhi"][0]=metphiCorr;

    //Preselection part

    if(doPreselection){
      bool passes = true;
      bool metCondition = (metptCorr >100.0 && Ht > 400.);

      float lep1phi = float_values["Event_Lepton1_Phi"];
      float lep1pt = float_values["Event_Lepton1_Pt"];

      float lep2phi = float_values["Event_Lepton2_Phi"];
      float lep2pt = float_values["Event_Lepton2_Pt"];

      float lep1px = 0.0;
      float lep1py = 0.0;
      float lep2px = 0.0;
      float lep2py = 0.0;
      float metPxLep=metPx;
      float metPyLep=metPy;
      
      if(lep1pt>0.0){
	lep1px = lep1pt*cos(lep1phi);
	lep1py = lep1pt*sin(lep1phi);

	if(lep2pt>0.0){
	  lep2px = lep2pt*cos(lep2phi);
	  lep2py = lep2pt*sin(lep2phi);
	}	
      } 
    
      metPxLep+=lep1px;
      metPyLep+=lep1py;

      metPxLep+=lep2px;
      metPyLep+=lep2py;
      
      double metLep=sqrt(metPxLep*metPxLep+metPyLep*metPyLep);
      
      metPxLep+=lep1px;

      metCondition = metCondition || (metLep>100.0);
	
      passes = passes && metCondition;
      passes = passes && nTightJets>=3;
      if (!passes ) {
	//Reset event weights/#objects
	string nameshortv= "Event";
	vector<string> extravars = additionalVariables(nameshortv);
	for(size_t addv = 0; addv < extravars.size();++addv){
	  string name = nameshortv+"_"+extravars.at(addv);
	  if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;

	}
	continue;
      }
    }
    
    TLorentzVector lepton;
    
    if( ( (electrons.size()==1 && muons.size()==0 ) || (muons.size()==1 && electrons.size()==0) ) && bjets.size()>0 ){
      if(electrons.size()==1) lepton = electrons[0];
      else if(muons.size()==1) lepton = muons[0];
      
      TVector2 met( metptCorr*cos(metphiCorr), metptCorr*sin(metphiCorr));
      float phi_lmet = fabs(deltaPhi(lepton.Phi(), metphiCorr) );
      float mt = sqrt(2* lepton.Pt() * metptCorr * ( 1- cos(phi_lmet)));
      float_values["Event_mt"] = (float)mt;
      Mt2Com_bisect *Mt2cal = new Mt2Com_bisect();
      double Mt2w = Mt2cal->calculateMT2w(jetsnob,bjets,lepton, met,"MT2w");
      float_values["Event_Mt2w"] = (float)Mt2w;    
    }
    
    //G                                                                                                                                              
    //cout << "starting ak8"<<endl;
    //    cout << "starting AK8 jets" << endl;
    for(int s = 0;s < min(max_instances[boosted_tops_subjets_label],sizes[boosted_tops_subjets_label]) ;++s){
      string pref = obj_to_pref[boosted_tops_subjets_label];
      float pt  = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Pt")][s];
      float eta = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Eta")][s];
      float phi = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Phi")][s];
      float e   = vfloats_values[makeName(boosted_tops_subjets_label,pref,"E")][s];

      float partonFlavourSubjet = vfloats_values[makeName(boosted_tops_subjets_label,pref,"PartonFlavour")][s];
      int flavorSubjet = int(partonFlavourSubjet);

      float bsfsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"noSyst");
      float bsfupsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"up");
      float bsfdownsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"down");
      
      vfloats_values[boosted_tops_subjets_label+"_BSF"][s]=bsfsubj;
      vfloats_values[boosted_tops_subjets_label+"_BSFUp"][s]=bsfupsubj;
      vfloats_values[boosted_tops_subjets_label+"_BSFDown"][s]=bsfdownsubj;
      
      TLorentzVector subjet;
      subjet.SetPtEtaPhiE(pt, eta, phi, e);       
      double minDR=999;
      float subjcsv = vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][s];
     
      bool isCSVM = (subjcsv>0.8484);

      vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][s] = (float)subjcsv;
      
      if(subjcsv>0.5426 && fabs(eta) < 2.4) {
	ncsvl_subj_tags +=1;
      }
      
      if(subjcsv>0.8484 && fabs(eta) < 2.4) {
	ncsvm_subj_tags +=1;
      }
      //cout << "medium csv subj tag is: " << ncsvm_subj_tags << endl;
      double csvleff_subj = MCTagEfficiencySubjet("csvl",flavorSubjet,pt, eta);
      double sfcsvl_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "noSyst", pt);
      
      double csvmeff_subj = MCTagEfficiencySubjet("csvm",flavorSubjet,pt, eta);
      double sfcsvm_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "noSyst", pt);

      double sfcsvl_mistag_up_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "mistag_up", pt);
      double sfcsvm_mistag_up_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "mistag_up", pt);

      double sfcsvl_mistag_down_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "mistag_down", pt);
      double sfcsvm_mistag_down_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "mistag_down", pt);

      double sfcsvl_b_tag_down_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "b_tag_down", pt);
      double sfcsvm_b_tag_down_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "b_tag_down", pt);
      
      double sfcsvl_b_tag_up_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "b_tag_up", pt);
      double sfcsvm_b_tag_up_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "b_tag_up", pt);     

      jsfscsvl_subj.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_subj));
      jsfscsvm_subj.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_subj));
      jsfscsvl_subj_mistag_up.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_mistag_up_subj));
      jsfscsvm_subj_mistag_up.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_mistag_up_subj));
      
      jsfscsvl_subj_b_tag_up.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_b_tag_up_subj));
      jsfscsvm_subj_b_tag_up.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_b_tag_up_subj));
      
      jsfscsvl_subj_mistag_down.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_mistag_down_subj));
      jsfscsvm_subj_mistag_down.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_mistag_down_subj));
      
      jsfscsvl_subj_b_tag_down.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_b_tag_down_subj));
      jsfscsvm_subj_b_tag_down.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_b_tag_down_subj));
           
      
      for(int t = 0;t < min(max_instances[boosted_tops_label],sizes[boosted_tops_label]) ;++t){
	  
	  float ptj = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
	  if (ptj<0.0)continue;
	  float etaj = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
	  float phij = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
	  float ej = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
	  TLorentzVector topjet;
	  topjet.SetPtEtaPhiE(ptj, etaj, phij, ej);       
	  
	  float DR = subjet.DeltaR(topjet); 
	    if(DR < minDR){
	    minDR = DR;
	    subj_jet_map[s]=t;
	    }
      }
	size_t tm = subj_jet_map[s];
	if(isCSVM)vfloats_values[boosted_tops_label+"_nCSVM"][tm]+=1;
	vfloats_values[boosted_tops_label+"_nJ"][tm]+=1;
    }
    
    for(int t = 0;t < max_instances[boosted_tops_label] ;++t){
      string pref = obj_to_pref[boosted_tops_label];
      float prunedMass   = vfloats_values[makeName(boosted_tops_label,pref,"prunedMass")][t];
      float softDropMass = vfloats_values[makeName(boosted_tops_label,pref,"softDropMass")][t];

      float topPt        = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
      float topEta = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
      //float topPhi = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
      float topE = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
      float tau1         = vfloats_values[makeName(boosted_tops_label,pref,"tau1")][t];
      float tau2         = vfloats_values[makeName(boosted_tops_label,pref,"tau2")][t];
      float tau3         = vfloats_values[makeName(boosted_tops_label,pref,"tau3")][t];

      float genpt8 = vfloats_values[makeName(boosted_tops_label,pref,"GenJetPt")][t];
      //float genpt8 = vfloats_values[jets_label+"_GenJetPt"][i];
      //float genpt = vfloats_values[makeName(jets_label,pref,"GenJetPt")][j];

      float ptCorr8 = -9999;
      float energyCorr8 = -9999;
      float smearfact8 = -9999;
      //cout << genpt8 << endl;    
      if(topPt>0){
	smearfact8 = smear(topPt, genpt8, topEta, syst);
	
	//cout << genpt8 << endl;
	ptCorr8 = topPt * smearfact8;
	energyCorr8 = topE * smearfact8;
	float unc8 = jetUncertainty8(ptCorr8,topEta,syst);
	//cout << unc8 << endl;
	ptCorr8 = ptCorr8 * (1 + unc8);
	energyCorr8 = energyCorr8 * (1 + unc8);
	//cout <<"pt uncor: " << topPt << "pt corr: " << ptCorr8 << endl;
	//cout << topPt << " " << ptCorr8 << endl;
      }
      //cout << topPt << " " << ptCorr8 << endl;
      
      vfloats_values[boosted_tops_label+"_CorrPt"][t]=ptCorr8;
      vfloats_values[boosted_tops_label+"_CorrE"][t]=energyCorr8;
      
      float tau3OVERtau2 = (tau2!=0. ? tau3/tau2 : 9999.);
      float tau2OVERtau1 = (tau1!=0. ? tau2/tau1 : 9999.);

      vfloats_values[makeName(boosted_tops_label,pref,"tau3OVERtau2")][t]=(float)tau3OVERtau2;
      vfloats_values[makeName(boosted_tops_label,pref,"tau2OVERtau1")][t]=(float)tau2OVERtau1;
      
      math::PtEtaPhiELorentzVector p4bestTop;
      math::PtEtaPhiELorentzVector p4bestB;

      /*math::PtEtaPhiELorentzVector j12;
      math::PtEtaPhiELorentzVector j13;
      math::PtEtaPhiELorentzVector j23;

      math::PtEtaPhiELorentzVector j1;
      math::PtEtaPhiELorentzVector j2;
      math::PtEtaPhiELorentzVector j3;*/
      
      int indexv0 = vfloats_values[makeName(boosted_tops_label,pref,"vSubjetIndex0")][t];
      int indexv1 = vfloats_values[makeName(boosted_tops_label,pref,"vSubjetIndex1")][t];

      int nCSVsubj = 0;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0] > 0.8484) ++nCSVsubj;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] > 0.8484) ++nCSVsubj;

      int nCSVsubj_tm = 0;

      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0] > 0.5426 && vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0]<0.8484) ++nCSVsubj_tm;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] > 0.5426 && vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] < 0.8484) ++nCSVsubj_tm;

      int nCSVsubj_t = 0;
      if(vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0]<0.5426) ++nCSVsubj_t;
      if(vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1]< 0.5426) ++nCSVsubj_t;
      
      vfloats_values[makeName(boosted_tops_label,pref,"nCSVsubj")][t]=(float)nCSVsubj;
      vfloats_values[makeName(boosted_tops_label,pref,"nCSVsubj_tm")][t]=(float)nCSVsubj_tm;
      /*
      j1 = math::PtEtaPhiELorentzVector(vfloats_values[makeName(boosted_tops_subjets_label,pref,"Pt")][index0],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Eta")][index0],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Phi")][index0],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"E")][index0]);
      j2 = math::PtEtaPhiELorentzVector(vfloats_values[makeName(boosted_tops_subjets_label,pref,"Pt")][index1],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Eta")][index1],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Phi")][index1],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"E")][index1]);
      j3 = math::PtEtaPhiELorentzVector(vfloats_values[makeName(boosted_tops_subjets_label,pref,"Pt")][index2],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Eta")][index2],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"Phi")][index2],
					 vfloats_values[makeName(boosted_tops_subjets_label,pref,"E")][index2]);

      j12 = j1 + j2;
      j13 = j1 + j3;
      j23 = j2 + j3;
      */
      // double m12 = j12.M();
      // double m13 = j13.M();
      // double m23 = j23.M();

      // double topmass = 172.3;
      // double wmass = 80.4;

      // to removed unused masses and ti add mass of ak8 jet and ptjet

      //bool isTop = ( ( (trimmedMass  <= 250 && trimmedMass >= 140) or (filteredMass <= 250 && filteredMass >= 140) or (prunedMass   <= 250 && prunedMass >= 140) )
      bool isTop = ( ( softDropMass <= 220 && softDropMass >=105 )
		     and (tau3OVERtau2 <= 0.81 )
		     and ( topPt > 500)
		     );

       bool isW = ( ( prunedMass <= 95 && prunedMass >=65 )
		    and (tau2OVERtau1 <= 0.6) //0.5
		    and ( topPt > 200)
		    );
    
      math::PtEtaPhiELorentzVector p4ak8;
   
       if (isW) {
	p4ak8 = math::PtEtaPhiELorentzVector(vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"E")][t]);

	float bestTopMass = 0.;
	float dRmin = 2.6;
	
	for(int i = 0; i < sizes[jets_label]; ++i){
	 
	  math::PtEtaPhiELorentzVector p4ak4 = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][i], 
									    vfloats_values[jets_label+"_CorrEta"][i], 
									    vfloats_values[jets_label+"_CorrPhi"][i], 
									    vfloats_values[jets_label+"_CorrE"][i] );
	  double dR = ROOT::Math::VectorUtil::DeltaR(p4ak8,p4ak4);

	  if (dR <= 0.8 or dR > 2.5) continue;

	  math::PtEtaPhiELorentzVector p4top = (p4ak8+p4ak4);
	  float topMass = p4top.mass();
	  if ( dR < dRmin ) {
	    dRmin = dR;
	    bestTopMass = topMass;
	    p4bestB= p4ak4;
	  }
	}

	if (bestTopMass > 250 or bestTopMass < 140) isW=false; 
      }

      if(isW){
	p4bestTop = p4bestB + p4ak8;
      }
      if(isTop) p4bestTop = p4ak8;
      

      vfloats_values[makeName(boosted_tops_label,pref,"isType2")][t]=(float)isW;
      vfloats_values[makeName(boosted_tops_label,pref,"isType1")][t]=(float)isTop;

      if(isW || isTop){
	TLorentzVector topjet;
	float ptj  = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
	float etaj = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
	float phij = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
	float ej   = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
	
	vfloats_values[makeName(boosted_tops_label,pref,"TopPt")][t]   = p4bestTop.pt();
	vfloats_values[makeName(boosted_tops_label,pref,"TopEta")][t]  = p4bestTop.eta();
	vfloats_values[makeName(boosted_tops_label,pref,"TopPhi")][t]  = p4bestTop.phi();
	vfloats_values[makeName(boosted_tops_label,pref,"TopE")][t]    = p4bestTop.energy();
	vfloats_values[makeName(boosted_tops_label,pref,"TopMass")][t] = p4bestTop.mass();
	if (isW)
	  vfloats_values[makeName(boosted_tops_label,pref,"TopWMass")][t] = vfloats_values[makeName(boosted_tops_label,pref,"prunedMass")][t];


	topjet.SetPtEtaPhiE(ptj, etaj, phij, ej);       
	if(isW){
	  float_values["Event_nType2TopJets"]+=1;
	  type2topjets.push_back(topjet);
	}
	if(isTop){
	  float_values["Event_nType1TopJets"]+=1;
	  type1topjets.push_back(topjet);
	}
      }
    }
    
    //    int nTightLeptons = electrons.size()+muons.size();
    size_t cat = 0;

    size_t ni = 9;
    cat+= 100000*(min(ni,muons.size()));
    cat+= 10000*(min(ni,electrons.size()));
    cat+= 1000*(min(ni,jets.size()));
    cat+= 100*(min(ni,bjets.size()));
    cat+= 10*(min(ni,type2topjets.size()));
    cat+= 1*(min(ni,type1topjets.size()));

    float_values["Event_category"]=(float)cat;
    string namelabel= "resolvedTopSemiLep";
    //Resolved tops Semileptonic:
    if(doResolvedTopSemiLep){
      sizes[namelabel]=0;
      if(((electrons.size()==1 && muons.size()==0 ) || (muons.size()==1 && electrons.size()==0)) &&  bjets.size()>0 &&   jets.size()>0){
	//	if(jets.size()==2 && bjets.size()==1) cout << " check this one "<<endl;
	TopUtilities topUtils;
	size_t t = 0;
	//	cout << " size b "<< bjets.size()<< " size l  "<< leptons.size() << " size 0 "<< sizes[namelabel]<<endl ;
	for(size_t b =0; b<bjets.size();++b){
	  for(size_t l =0; l<leptons.size();++l){
	    //	  double metPx= 1.0, metPy =1.0;
	    math::PtEtaPhiELorentzVector topSemiLep = topUtils.top4Momentum(leptons.at(l), bjets.at(b),metPx , metPy);
	    if(t > (size_t)max_instances[namelabel])continue;
	    vfloats_values[namelabel+"_Pt"][t]=topSemiLep.pt();
	    vfloats_values[namelabel+"_Eta"][t]=topSemiLep.eta();
	    vfloats_values[namelabel+"_Phi"][t]=topSemiLep.phi();
	    vfloats_values[namelabel+"_E"][t]=topSemiLep.energy();
	    vfloats_values[namelabel+"_Mass"][t]=topSemiLep.mass();
	    vfloats_values[namelabel+"_MT"][t]= topUtils.topMtw(leptons.at(l),bjets.at(b),metPx,metPy);
	    vfloats_values[namelabel+"_LBMPhi"][t]=deltaPhi((leptons.at(l)+bjets.at(b)).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_LMPhi"][t]= deltaPhi(leptons.at(l).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_BMPhi"][t]= deltaPhi(bjets.at(b).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_TMPhi"][t]= deltaPhi(topSemiLep.Phi(), metphiCorr);
	    vfloats_values[namelabel+"_LBPhi"][t]= deltaPhi(leptons.at(l).Phi(), bjets.at(b).Phi());

	    vfloats_values[namelabel+"_IndexB"][t]= mapBJets[b];
	    float lidx = -9999;
	    float flav = -9999;
	    if(mapMu[l]!=-1){lidx=mapMu[l]; flav = 11;}
	    if(mapEle[l]!=-1){lidx=mapEle[l]; flav = 13;}
	    if(mapMu[l]!=-1 && mapEle[l]!=-1) cout<<" what is going on? " <<endl;
	    
	    vfloats_values[namelabel+"_IndexL"][t]= lidx;
	    vfloats_values[namelabel+"_LeptonFlavour"][t]= flav;
	    
	    ++t;
	    ++sizes[namelabel];
	  }
	}
      }
      if((bjets.size()==0 || leptons.size()==0)){
	vector<string> extravars = additionalVariables(namelabel);
	for(size_t t =0;t<(size_t)max_instances[namelabel];++t){
	  for(size_t addv = 0; addv < extravars.size();++addv){
	    vfloats_values[namelabel+"_"+extravars.at(addv)][t]=-9999;
	    //vfloats_values[namelabel+"_"+extravars.at(addv)][t]=0.0;
	  }
	}
      }
      //      cout << "after all loop size is "<<  sizes[namelabel]<<endl;
    }
    //Resolved tops Fullhadronic:
   
    if(doResolvedTopHad || doResolvedTopSemiLep) {
      //      if(getwobjets) cout << " test 1 "<<endl;
      string namelabel= "resolvedTopHad";
      sizes[namelabel]=0;
      
      // ==================== Implementing kinematic fitter ==========================

      if(jets.size()>2 && bjets.size()>0   ){
      
	int maxJetLoop = min((int)(max_instances[jets_label]),max_leading_jets_for_top);
	maxJetLoop = min(maxJetLoop,sizes[jets_label]);
	string pref = obj_to_pref[jets_label];
	size_t t = 0;

	TLorentzVector jet1, jet2, jet3;
	
	
	for(int i = 0; i < maxJetLoop ;++i){
	  for(int j = i+1; j < maxJetLoop ;++j){
	    for(int k = j+1; k < maxJetLoop ;++k){
		      
	      if(vfloats_values[jets_label+"_CorrPt"][i]> 0. && vfloats_values[jets_label+"_CorrPt"][j]> 0. && vfloats_values[jets_label+"_CorrPt"][k]> 0.){
	      
	     
	      //
	      // Kinematic fitter
	      //
		//	      KinematicFitter fitter;
	      
	      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      // Stuff that goes into event loop
	      //
	      //*********************************
	      //
	      // Example tri-jet combination
	      // The MVA training is done such that the "b-jet" is the jet in the triplet
	      // that has the highest CSVv2+IVF value. I've called this one "jet3" in the example.
	      //

	      float  csv_i= (float)vfloats_values[jets_label+"_CSVv2"][i];
	      float  csv_j= (float)vfloats_values[jets_label+"_CSVv2"][j];
	      float  csv_k= (float)vfloats_values[jets_label+"_CSVv2"][k];
	      
	      vector<float> csv = {csv_i, csv_j, csv_k};
	      vector<int> idx = {i, j, k};
	   
	      ++t;
      
	      size_t tt = 0;
	      //	cout << " test 3 "<<endl;
	      string pref = obj_to_pref[jets_label];
	      
	      bool isIBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][i] && (fabs(vfloats_values[jets_label+"_CorrEta"][i]) < 2.4);
	      bool isITight= (bool)vfloats_values[jets_label+"_IsTight"][i];
	      
	      bool isJBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][j] && (fabs(vfloats_values[jets_label+"_CorrEta"][j]) < 2.4);
	      bool isJTight=(bool)vfloats_values[jets_label+"_IsTight"][j];
		
		//for(int k = j+1; k < maxJetLoop ;++k){
		int nBJets =0;
	
		bool isKBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][k] && (fabs(vfloats_values[jets_label+"_CorrEta"][k]) < 2.4);
		bool isKTight=(bool)vfloats_values[jets_label+"_IsTight"][k];
		
		if(isIBJet)nBJets++;
		if(isJBJet)nBJets++;
		if(isKBJet)nBJets++;
	
		if(nBJets !=1  || !(isITight && isJTight && isKTight) ) continue;
		math::PtEtaPhiELorentzVector p4i = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][i], vfloats_values[jets_label+"_CorrEta"][i], vfloats_values[jets_label+"_CorrPhi"][i], vfloats_values[jets_label+"_CorrE"][i] );
		math::PtEtaPhiELorentzVector p4j = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][j], vfloats_values[jets_label+"_CorrEta"][j], vfloats_values[jets_label+"_CorrPhi"][j], vfloats_values[jets_label+"_CorrE"][j] );
		math::PtEtaPhiELorentzVector p4k = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][k], vfloats_values[jets_label+"_CorrEta"][k], vfloats_values[jets_label+"_CorrPhi"][k], vfloats_values[jets_label+"_CorrE"][k] );
		math::PtEtaPhiELorentzVector topHad = (p4i+p4j)+p4k;
		vfloats_values[namelabel+"_Pt"][tt]=topHad.pt(); 
		vfloats_values[namelabel+"_Eta"][tt]=topHad.eta();
		vfloats_values[namelabel+"_Phi"][tt]=topHad.phi();
		vfloats_values[namelabel+"_E"][tt]=topHad.e();
		vfloats_values[namelabel+"_Mass"][tt]=topHad.mass();
		if(isKBJet){
		  vfloats_values[namelabel+"_IndexB"][tt]= k; vfloats_values[namelabel+"_IndexJ1"][tt]= i;  vfloats_values[namelabel+"_IndexJ2"][tt]= j; 
		}
		if(isJBJet){
		  vfloats_values[namelabel+"_IndexB"][tt]= j; vfloats_values[namelabel+"_IndexJ1"][tt]= i;  vfloats_values[namelabel+"_IndexJ2"][tt]= k; }
		if(isIBJet){
		  //		cout << " isibjet"<<endl;
		  vfloats_values[namelabel+"_IndexB"][tt]= i; vfloats_values[namelabel+"_IndexJ1"][tt]= j;  vfloats_values[namelabel+"_IndexJ2"][tt]= k; }

		math::PtEtaPhiELorentzVector p4w, p4b ;
		float massdrop_=0., deltaRjets=0.;
		if(isIBJet){ p4w = p4j+p4k; p4b= p4i;
		  massdrop_ = max( p4j.mass(),  p4k.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4j,  p4k);}
		if(isJBJet){p4w = p4i+p4k; p4b= p4j;
		  massdrop_ = max( p4i.mass(),  p4k.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4i,  p4k);}
		if(isKBJet){p4w = p4j+p4i; p4b= p4k;
		  massdrop_ = max( p4i.mass(),  p4j.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4i,  p4j);}
		vfloats_values[namelabel+"_Pt"][tt]=topHad.pt(); 
		vfloats_values[namelabel+"_Eta"][tt]=topHad.eta();
		vfloats_values[namelabel+"_Phi"][tt]=topHad.phi();
		vfloats_values[namelabel+"_E"][tt]=topHad.e();
		vfloats_values[namelabel+"_Mass"][tt]=topHad.mass();
		vfloats_values[namelabel+"_WMass"][tt]=p4w.mass();
		vfloats_values[namelabel+"_massDrop"][tt]=massdrop_*deltaRjets;
		//cout<<"Mass drop: "<<massdrop_*deltaRjets<<endl;
		
		if(topHad.mass()<0. || p4w.mass()<0){
		  float genpti = vfloats_values[jets_label+"_GenJetPt"][i];
		  float genptj = vfloats_values[jets_label+"_GenJetPt"][j];
		  float genptk = vfloats_values[jets_label+"_GenJetPt"][k];
		  
		  float pti = vfloats_values[jets_label+"_Pt"][i];
		  float ptj = vfloats_values[jets_label+"_Pt"][j];
		  float ptk = vfloats_values[jets_label+"_Pt"][k];
		  
		  std::cout<<"Top Mass: "<<topHad.mass()<<std::endl;
		  std::cout<<"W Mass: "<<p4w.mass()<<std::endl;
		  std::cout<<"Pt1: "<<p4i.pt()<<" gen "<< genpti<<" pt0 "<<pti <<" eta "<< p4i.eta()<< " phi "<< p4i.phi()<< " e "<< p4i.e()<<std::endl;
		  std::cout<<"Pt2: "<<p4j.pt()<< " gen "<< genptj<< " pt0 "<<ptj <<" eta "<< p4j.eta()<< " phi "<< p4j.phi()<< " e "<< p4j.e()<<std::endl;
		  std::cout<<"Pt3> "<<p4k.pt()<< " gen "<< genptk<<" pt0 "<<ptk<< " eta "<< p4k.eta()<< " phi "<< p4k.phi()<< " e "<< p4k.e()<<std::endl;
		}
		if(t > (size_t)max_instances[namelabel])continue;
		vfloats_values[namelabel+"_WMPhi"][tt]=deltaPhi(p4w.Phi(), metphiCorr);
		vfloats_values[namelabel+"_TMPhi"][tt]= deltaPhi(topHad.Phi(), metphiCorr);
		vfloats_values[namelabel+"_BMPhi"][tt]= deltaPhi(p4b.Phi(), metphiCorr);
		vfloats_values[namelabel+"_WBPhi"][tt]= deltaPhi(p4w.Phi(), p4b.Phi());

		++tt;
		++sizes[namelabel];
		
		//  } //end of if statement on the number of bjets
	      }//end if statement on jets
	    }// end 3rd loop on jets
	  }// end 2st loop on jets
	}// end 1st loop on jets
	
      }// end if statement (at least 3 jets)
      // =============================================================================
      
      //	    }
      //	  }	
      //	}
      //}//end if on number of jets and bjets
      
      //      if(getwobjets)cout << " nresolvedtophad "<< sizes["resolvedTopHad"]<<endl;
      //      cout << " namelabel? "<< namelabel<< endl;
      if(sizes[namelabel]==0){
	vector<string> extravars = additionalVariables(namelabel);
	for(size_t t =0;t<(size_t)max_instances[namelabel];++t){
	  //	  if(t> 45)continue;
	  for(size_t addv = 0; addv < extravars.size();++addv){
	    //	    cout << "resetting variable "<< extravars.at(addv)<< " name "<< namelabel+"_"+extravars.at(addv)<<" t is "<< t << endl;
	    vfloats_values[namelabel+"_"+extravars.at(addv)][t]=-9999;
	  }
	}
      }
    }
    //G                                                                                                                                              
    //cout << "starting btag"<<endl;
    //    cout << "starting btagging" << endl;
    //BTagging part
    if(doBTagSF){
    //CSVT
      //0 tags
      b_weight_csvt_0_tags = b_csvt_0_tags.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_0_tags_mistag_up = b_csvt_0_tags.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_0_tags_mistag_down = b_csvt_0_tags.weight(jsfscsvt_mistag_down, ncsvt_tags);  
      b_weight_csvt_0_tags_b_tag_up = b_csvt_0_tags.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_0_tags_b_tag_down = b_csvt_0_tags.weight(jsfscsvt_b_tag_down, ncsvt_tags);

      //1
      b_weight_csvt_1_tag = b_csvt_1_tag.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_1_tag_mistag_up = b_csvt_1_tag.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_1_tag_mistag_down = b_csvt_1_tag.weight(jsfscsvt_mistag_down, ncsvt_tags);  
      b_weight_csvt_1_tag_b_tag_up = b_csvt_1_tag.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_1_tag_b_tag_down = b_csvt_1_tag.weight(jsfscsvt_b_tag_down, ncsvt_tags);
      
      //2 tags
      b_weight_csvt_2_tags = b_csvt_2_tags.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_2_tags_mistag_up = b_csvt_2_tags.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_2_tags_mistag_down = b_csvt_2_tags.weight(jsfscsvt_mistag_down, ncsvt_tags);  
      b_weight_csvt_2_tags_b_tag_up = b_csvt_2_tags.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_2_tags_b_tag_down = b_csvt_2_tags.weight(jsfscsvt_b_tag_down, ncsvt_tags);  
      
      //1-2 tags
      b_weight_csvt_1_2_tags = b_csvt_1_2_tags.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_1_2_tags_b_tag_up = b_csvt_1_2_tags.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_1_2_tags_b_tag_down = b_csvt_1_2_tags.weight(jsfscsvt_b_tag_down, ncsvt_tags);  
      b_weight_csvt_1_2_tags_mistag_up = b_csvt_1_2_tags.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_1_2_tags_mistag_down = b_csvt_1_2_tags.weight(jsfscsvt_mistag_down, ncsvt_tags);  

    //CSVM
      //cout << "subj tags is: " << ncsvm_subj_tags << endl;
      //0 tags
      b_weight_csvm_0_tags = b_csvm_0_tags.weight(jsfscsvm, ncsvm_tags);  
      b_weight_csvm_0_tags_mistag_up = b_csvm_0_tags.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_0_tags_mistag_down = b_csvm_0_tags.weight(jsfscsvm_mistag_down, ncsvm_tags);  
      b_weight_csvm_0_tags_b_tag_up = b_csvm_0_tags.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_0_tags_b_tag_down = b_csvm_0_tags.weight(jsfscsvm_b_tag_down, ncsvm_tags);  
      
      b_weight_subj_csvm_0_tags = b_csvm_0_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_mistag_up = b_csvm_0_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_mistag_down = b_csvm_0_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_b_tag_up = b_csvm_0_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_b_tag_down = b_csvm_0_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
   
      //1 tag
      b_weight_csvm_1_tag = b_csvm_1_tag.weight(jsfscsvm, ncsvm_tags);  
      b_weight_csvm_1_tag_mistag_up = b_csvm_1_tag.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_1_tag_mistag_down = b_csvm_1_tag.weight(jsfscsvm_mistag_down, ncsvm_tags);  
      b_weight_csvm_1_tag_b_tag_up = b_csvm_1_tag.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_1_tag_b_tag_down = b_csvm_1_tag.weight(jsfscsvm_b_tag_down, ncsvm_tags);  

      b_weight_subj_csvm_1_tag = b_csvm_1_tag.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_mistag_up = b_csvm_1_tag.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_mistag_down = b_csvm_1_tag.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_b_tag_up = b_csvm_1_tag.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_b_tag_down = b_csvm_1_tag.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
      //2 tags
      b_weight_csvm_2_tags = b_csvm_2_tags.weight(jsfscsvm, ncsvm_tags);  
      b_weight_csvm_2_tags_mistag_up = b_csvm_2_tags.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_2_tags_mistag_down = b_csvm_2_tags.weight(jsfscsvm_mistag_down, ncsvm_tags);  
      b_weight_csvm_2_tags_b_tag_up = b_csvm_2_tags.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_2_tags_b_tag_down = b_csvm_2_tags.weight(jsfscsvm_b_tag_down, ncsvm_tags);  

      b_weight_subj_csvm_2_tags = b_csvm_2_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_mistag_up = b_csvm_2_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_mistag_down = b_csvm_2_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_b_tag_up = b_csvm_2_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_b_tag_down = b_csvm_2_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
      //1-2 tags
      b_weight_csvm_1_2_tags = b_csvm_1_2_tags.weight(jsfscsvm, ncsvm_tags);  
      //giorgia
      //cout << "ncsvm_tag: " << ncsvm_tags <<  " b_weight_csvm_0_tags: " << b_weight_csvm_0_tags << " b_weight_csvm_1_2_tags: " << b_weight_csvm_1_2_tags << endl;
      b_weight_csvm_1_2_tags_b_tag_up = b_csvm_1_2_tags.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_1_2_tags_b_tag_down = b_csvm_1_2_tags.weight(jsfscsvm_b_tag_down, ncsvm_tags);  
      b_weight_csvm_1_2_tags_mistag_up = b_csvm_1_2_tags.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_1_2_tags_mistag_down = b_csvm_1_2_tags.weight(jsfscsvm_mistag_down, ncsvm_tags);  

      b_weight_subj_csvm_1_2_tags = b_csvm_1_2_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_2_tags_mistag_up = b_csvm_1_2_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_2_tags_mistag_down = b_csvm_1_2_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_2_tags_b_tag_up = b_csvm_1_2_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_2_tags_b_tag_down = b_csvm_1_2_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
    //CSVL
      //0 tags
      b_weight_csvl_0_tags = b_csvl_0_tags.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_0_tags_mistag_up = b_csvl_0_tags.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_0_tags_mistag_down = b_csvl_0_tags.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      b_weight_csvl_0_tags_b_tag_up = b_csvl_0_tags.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_0_tags_b_tag_down = b_csvl_0_tags.weight(jsfscsvl_b_tag_down, ncsvl_tags);  

     b_weight_subj_csvl_0_tags = b_csvl_0_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_mistag_up = b_csvl_0_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_mistag_down = b_csvl_0_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_b_tag_up = b_csvl_0_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_b_tag_down = b_csvl_0_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      //1 tag
      b_weight_csvl_1_tag = b_csvl_1_tag.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_1_tag_mistag_up = b_csvl_1_tag.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_1_tag_mistag_down = b_csvl_1_tag.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      b_weight_csvl_1_tag_b_tag_up = b_csvl_1_tag.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_1_tag_b_tag_down = b_csvl_1_tag.weight(jsfscsvl_b_tag_down, ncsvl_tags);  

      b_weight_subj_csvl_1_tag = b_csvl_1_tag.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_mistag_up = b_csvl_1_tag.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_mistag_down = b_csvl_1_tag.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_b_tag_up = b_csvl_1_tag.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_b_tag_down = b_csvl_1_tag.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      //2 tags
      b_weight_csvl_2_tags = b_csvl_2_tags.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_2_tags_mistag_up = b_csvl_2_tags.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_2_tags_mistag_down = b_csvl_2_tags.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      b_weight_csvl_2_tags_b_tag_up = b_csvl_2_tags.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_2_tags_b_tag_down = b_csvl_2_tags.weight(jsfscsvl_b_tag_down, ncsvl_tags);  

      b_weight_subj_csvl_2_tags = b_csvl_2_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_mistag_up = b_csvl_2_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_mistag_down = b_csvl_2_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_b_tag_up = b_csvl_2_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_b_tag_down = b_csvl_2_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      //1-2 tags
      b_weight_csvl_1_2_tags = b_csvl_1_2_tags.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_1_2_tags_b_tag_up = b_csvl_1_2_tags.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_1_2_tags_b_tag_down = b_csvl_1_2_tags.weight(jsfscsvl_b_tag_down, ncsvl_tags);  
      b_weight_csvl_1_2_tags_mistag_up = b_csvl_1_2_tags.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_1_2_tags_mistag_down = b_csvl_1_2_tags.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      
      b_weight_subj_csvl_1_2_tags = b_csvl_1_2_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_2_tags_mistag_up = b_csvl_1_2_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_2_tags_mistag_down = b_csvl_1_2_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_2_tags_b_tag_up = b_csvl_1_2_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_2_tags_b_tag_down = b_csvl_1_2_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      float_values["Event_bWeight0CSVL_subj"]=b_weight_subj_csvl_0_tags;
      float_values["Event_bWeight1CSVL_subj"]=b_weight_subj_csvl_1_tag;
      float_values["Event_bWeight2CSVL_subj"]=b_weight_subj_csvl_2_tags;
      float_values["Event_bWeight1_2CSVL_subj"]=b_weight_subj_csvl_1_2_tags;
      
      float_values["Event_bWeight0CSVM_subj"]=b_weight_subj_csvm_0_tags;
      float_values["Event_bWeight1CSVM_subj"]=b_weight_subj_csvm_1_tag;
      float_values["Event_bWeight2CSVM_subj"]=b_weight_subj_csvm_2_tags;
      float_values["Event_bWeight1_2CSVM_subj"]=b_weight_subj_csvm_1_2_tags;
      //cout << "b_weight_subj_csvm_1_2_tags: " << b_weight_subj_csvm_1_2_tags << endl; 
      //Mistag
      float_values["Event_bWeightMisTagUp0CSVL_subj"]=b_weight_subj_csvl_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVL_subj"]=b_weight_subj_csvl_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVL_subj"]=b_weight_subj_csvl_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1_2CSVL_subj"]=b_weight_subj_csvl_1_2_tags_mistag_up;
      
      float_values["Event_bWeightMisTagUp0CSVM_subj"]=b_weight_subj_csvm_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVM_subj"]=b_weight_subj_csvm_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVM_subj"]=b_weight_subj_csvm_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1_2CSVM_subj"]=b_weight_subj_csvm_1_2_tags_mistag_up;
    
      float_values["Event_bWeightMisTagDown0CSVL_subj"]=b_weight_subj_csvl_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVL_subj"]=b_weight_subj_csvl_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVL_subj"]=b_weight_subj_csvl_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1_2CSVL_subj"]=b_weight_subj_csvl_1_2_tags_mistag_down;
      
      float_values["Event_bWeightMisTagDown0CSVM_subj"]=b_weight_subj_csvm_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVM_subj"]=b_weight_subj_csvm_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVM_subj"]=b_weight_subj_csvm_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1_2CSVM_subj"]=b_weight_subj_csvm_1_2_tags_mistag_down;
      
      //Btag
      float_values["Event_bWeightBTagUp0CSVL_subj"]=b_weight_subj_csvl_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVL_subj"]=b_weight_subj_csvl_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVL_subj"]=b_weight_subj_csvl_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1_2CSVL_subj"]=b_weight_subj_csvl_1_2_tags_b_tag_up;
      
      float_values["Event_bWeightBTagUp0CSVM_subj"]=b_weight_subj_csvm_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVM_subj"]=b_weight_subj_csvm_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVM_subj"]=b_weight_subj_csvm_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1_2CSVM_subj"]=b_weight_subj_csvm_1_2_tags_b_tag_up;

      float_values["Event_bWeightBTagDown0CSVL_subj"]=b_weight_subj_csvl_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVL_subj"]=b_weight_subj_csvl_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVL_subj"]=b_weight_subj_csvl_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1_2CSVL_subj"]=b_weight_subj_csvl_1_2_tags_b_tag_down;
      
      float_values["Event_bWeightBTagDown0CSVM_subj"]=b_weight_subj_csvm_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVM_subj"]=b_weight_subj_csvm_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVM_subj"]=b_weight_subj_csvm_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1_2CSVM_subj"]=b_weight_subj_csvm_1_2_tags_b_tag_down;
      
      /////AK4
      
      float_values["Event_bWeight0CSVL"]=b_weight_csvl_0_tags;
      float_values["Event_bWeight1CSVL"]=b_weight_csvl_1_tag;
      float_values["Event_bWeight2CSVL"]=b_weight_csvl_2_tags;
      float_values["Event_bWeight1_2CSVL"]=b_weight_csvl_1_2_tags;
      
      float_values["Event_bWeight0CSVM"]=b_weight_csvm_0_tags;
      float_values["Event_bWeight1CSVM"]=b_weight_csvm_1_tag;
      float_values["Event_bWeight2CSVM"]=b_weight_csvm_2_tags;
      float_values["Event_bWeight1_2CSVM"]=b_weight_csvm_1_2_tags;
      //cout << "b_weight_csvm_1_2_tags: " << b_weight_csvm_1_2_tags << endl;

      float_values["Event_bWeight0CSVT"]=b_weight_csvt_0_tags;
      float_values["Event_bWeight1CSVT"]=b_weight_csvt_1_tag;
      float_values["Event_bWeight2CSVT"]=b_weight_csvt_2_tags;
      float_values["Event_bWeight1_2CSVT"]=b_weight_csvt_1_2_tags;
      
   
      //Mistag
      float_values["Event_bWeightMisTagUp0CSVL"]=b_weight_csvl_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVL"]=b_weight_csvl_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVL"]=b_weight_csvl_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1_2CSVL"]=b_weight_csvl_1_2_tags_mistag_up;
      
      float_values["Event_bWeightMisTagUp0CSVM"]=b_weight_csvm_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVM"]=b_weight_csvm_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVM"]=b_weight_csvm_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1_2CSVM"]=b_weight_csvm_1_2_tags_mistag_up;
    
      float_values["Event_bWeightMisTagUp0CSVT"]=b_weight_csvt_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVT"]=b_weight_csvt_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVT"]=b_weight_csvt_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1_2CSVT"]=b_weight_csvt_1_2_tags_mistag_up;
      
      
      float_values["Event_bWeightMisTagDown0CSVL"]=b_weight_csvl_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVL"]=b_weight_csvl_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVL"]=b_weight_csvl_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1_2CSVL"]=b_weight_csvl_1_2_tags_mistag_down;
      
      float_values["Event_bWeightMisTagDown0CSVM"]=b_weight_csvm_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVM"]=b_weight_csvm_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVM"]=b_weight_csvm_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1_2CSVM"]=b_weight_csvm_1_2_tags_mistag_down;
      
      float_values["Event_bWeightMisTagDown0CSVT"]=b_weight_csvt_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVT"]=b_weight_csvt_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVT"]=b_weight_csvt_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1_2CSVT"]=b_weight_csvt_1_2_tags_mistag_down;

      //Btag
      float_values["Event_bWeightBTagUp0CSVL"]=b_weight_csvl_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVL"]=b_weight_csvl_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVL"]=b_weight_csvl_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1_2CSVL"]=b_weight_csvl_1_2_tags_b_tag_up;
      
      float_values["Event_bWeightBTagUp0CSVM"]=b_weight_csvm_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVM"]=b_weight_csvm_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVM"]=b_weight_csvm_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1_2CSVM"]=b_weight_csvm_1_2_tags_b_tag_up;
      
      float_values["Event_bWeightBTagUp0CSVT"]=b_weight_csvt_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVT"]=b_weight_csvt_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVT"]=b_weight_csvt_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1_2CSVT"]=b_weight_csvt_1_2_tags_b_tag_up;
      
      float_values["Event_bWeightBTagDown0CSVL"]=b_weight_csvl_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVL"]=b_weight_csvl_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVL"]=b_weight_csvl_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1_2CSVL"]=b_weight_csvl_1_2_tags_b_tag_down;
      
      float_values["Event_bWeightBTagDown0CSVM"]=b_weight_csvm_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVM"]=b_weight_csvm_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVM"]=b_weight_csvm_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1_2CSVM"]=b_weight_csvm_1_2_tags_b_tag_down;
      
      float_values["Event_bWeightBTagDown0CSVT"]=b_weight_csvt_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVT"]=b_weight_csvt_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVT"]=b_weight_csvt_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1_2CSVT"]=b_weight_csvt_1_2_tags_b_tag_down;
    }
    
    
    float LHEWeightSign=1.0;
    if(useLHE){
      //LHE and luminosity weights:
      float weightsign = lhes->hepeup().XWGTUP;
      float_values["Event_LHEWeight"]=weightsign;
      LHEWeightSign = weightsign/fabs(weightsign);
      float_values["Event_LHEWeightSign"]=LHEWeightSign;
     }
    float weightLumi = crossSection/originalEvents;
    float_values["Event_weight"]=weightLumi*LHEWeightSign;
    
    //Part 3: filling the additional variables

    if(useLHEWeights){
      getEventLHEWeights();
    }
    if(addLHAPDFWeights){
      getEventPdf();
    }
    
    if(doPU){
      iEvent.getByToken(t_ntrpu_,ntrpu);
      int nTruePV=*ntrpu;
      float_values["Event_nTruePV"]=(float)(nTruePV);
    }

    if(addPV){
      float nGoodPV = 0.0;
      for (size_t v = 0; v < pvZ->size();++v){
	bool isGoodPV = (
			 fabs(pvZ->at(v)) < 24.0 &&
			 pvNdof->at(v) > 4.0 &&
			 pvRho->at(v) <2.0
			 );
	if (isGoodPV)nGoodPV+=1.0;
      }	
      float_values["Event_nGoodPV"]=(float)(nGoodPV);
     float_values["Event_nPV"]=(float)(nPV);
    }

    //technical event informationx
    double_values["Event_EventNumber"]=*eventNumber;
    float_values["Event_LumiBlock"]=*lumiBlock;
    float_values["Event_RunNumber"]=*runNumber;
 
    trees[syst]->Fill();
    
    //Reset event weights/#objects
    string nameshortv= "Event";
    vector<string> extravars = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravars.size();++addv){
      string name = nameshortv+"_"+extravars.at(addv);

      if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;
    }
  }
  for(int t = 0;t < max_instances[boosted_tops_label] ;++t){
    vfloats_values[boosted_tops_label+"_nCSVM"][t]=0;
    vfloats_values[boosted_tops_label+"_nJ"][t]=0;
  }
   
}

bool DMAnalysisTreeMaker::flavourFilter(string ch, int nb, int nc, int nl)
{

  if (ch == "WJets_wbb" || ch == "ZJets_wbb") return (nb > 0 );
  if (ch == "WJets_wcc" || ch == "ZJets_wcc") return (nb == 0 && nc > 0);
  if (ch == "WJets_wlight" || ch == "ZJets_wlight") return (nb == 0 && nc == 0);
  return true;
}



int DMAnalysisTreeMaker::eventFlavour(bool getFlavour, int nb, int nc, int nl)
{
  if (!getFlavour) return 0;
  else
    {
      if ( flavourFilter("WJets_wlight", nb, nc, nl) ) return 1;
      if ( flavourFilter("WJets_wcc", nb, nc, nl) ) return 2;
      if ( flavourFilter("WJets_wbb", nb, nc, nl) ) return 3;
    }
  return 0;
}


string DMAnalysisTreeMaker::makeBranchName(string label, string pref, string var){ //substitutes the "pref" word with "label+'_'"
  string outVar = var;
  size_t prefPos=outVar.find(pref);
  size_t prefLength = pref.length();
  outVar.replace(prefPos,prefLength,label+"_");
  return outVar;
}

string DMAnalysisTreeMaker::makeBranchNameCat(string label, string cat, string pref, string var){
  return makeBranchName (label+cat,pref,var);
}

string DMAnalysisTreeMaker::makeName(string label,string pref,string var){
  return label+"_"+var;
}

void DMAnalysisTreeMaker::initCategoriesSize(string label){
  for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
    string category = obj_cats[label].at(sc);
    sizes[label+category]=0;
  }
}

void DMAnalysisTreeMaker::setCategorySize(string label, string category, size_t size){
    sizes[label+category]=size;
}
void DMAnalysisTreeMaker::fillCategory(string label, string category, int pos_nocat, int pos_cat){
  //  for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
  //    string category = obj_cats.at(sc);
  for (size_t obj =0; obj< obj_to_floats[label].size(); ++obj){
    string var = obj_to_floats[label].at(obj);
    string varCat = makeBranchNameCat(label,category,label+"_",var);
    //    cout << " var "<< var << " varcat "<< varCat<<endl;
    vfloats_values[varCat][pos_cat]= vfloats_values[var][pos_nocat];
  }
}


vector<string> DMAnalysisTreeMaker::additionalVariables(string object){
  vector<string> addvar;
  bool ismuon=object.find("muon")!=std::string::npos;
  bool isphoton=object.find("photon")!=std::string::npos;
  bool iselectron=object.find("electron")!=std::string::npos;
  bool ismet=object.find("met")!=std::string::npos;
  bool isjet=object.find("jet")!=std::string::npos && object.find("AK4")!=std::string::npos;
  bool isak8=object.find("jet")!=std::string::npos && object.find("AK8")!=std::string::npos && object.find("sub")==std::string::npos;
  bool isak8subjet=object.find("jet")!=std::string::npos && object.find("AK8")!=std::string::npos && object.find("sub")!=std::string::npos;
  bool isevent=object.find("Event")!=std::string::npos;
  bool isResolvedTopHad=object.find("resolvedTopHad")!=std::string::npos;
  bool isResolvedTopSemiLep=object.find("resolvedTopSemiLep")!=std::string::npos;
  //  bool isgen = object.find("genPart")!=std::string::npos;

  if(ismuon || iselectron){
    addvar.push_back("SFTrigger");
  }
  
  if(isphoton){
    addvar.push_back("isLooseSpring15");
    addvar.push_back("isMediumSpring15");
    addvar.push_back("isTightSpring15");
  }
  if(iselectron){
    addvar.push_back("PassesDRmu");
  }
  if(ismet){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrPhi");
  }
  if(isjet){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrE");
    addvar.push_back("NoCorrPt");
    addvar.push_back("NoCorrE");
    addvar.push_back("MinDR");
    addvar.push_back("IsCSVT");
    addvar.push_back("IsCSVM");
    addvar.push_back("IsCSVL");
    addvar.push_back("PassesID");
    addvar.push_back("PassesDR");
    addvar.push_back("CorrMass");
    addvar.push_back("IsTight");
    addvar.push_back("IsLoose");
  }
  if(isak8){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrE");
    addvar.push_back("isType1");
    addvar.push_back("isType2");
    addvar.push_back("TopPt");
    addvar.push_back("TopEta");
    addvar.push_back("TopPhi");
    addvar.push_back("TopE");
    addvar.push_back("TopMass");
    addvar.push_back("TopWMass");
    addvar.push_back("nJ");
    addvar.push_back("nCSVM");
    addvar.push_back("CSVv2");
    addvar.push_back("nCSVsubj");
    addvar.push_back("nCSVsubj_tm");
    addvar.push_back("tau3OVERtau2");
    addvar.push_back("tau2OVERtau1");

  }  
  if(1>0){
    addvar.push_back("Pt");
    addvar.push_back("Eta");
    addvar.push_back("Phi");
    addvar.push_back("E");
    addvar.push_back("Id");
    addvar.push_back("Status");
    addvar.push_back("MomId");
    addvar.push_back("Mom0Status");
    /*addvar.push_back("dauId1");
    addvar.push_back("dauStatus1");
    addvar.push_back("dauId2");
    addvar.push_back("dauStatus2");
    addvar.push_back("dauId3");
    addvar.push_back("dauStatus3");*/
  }
  if(isak8subjet){
    ;//    addvar.push_back("CorrPt");
  }

  if(isResolvedTopHad ){
    addvar.push_back("Pt");    addvar.push_back("Eta");    addvar.push_back("Phi");    addvar.push_back("E"); addvar.push_back("Mass");  addvar.push_back("massDrop");
    addvar.push_back("WMass"); addvar.push_back("BMPhi");  addvar.push_back("WMPhi");  addvar.push_back("TMPhi");  addvar.push_back("WBPhi");
    addvar.push_back("IndexB");    addvar.push_back("IndexJ1");    addvar.push_back("IndexJ2");  addvar.push_back("IndexB_MVA");    addvar.push_back("IndexJ1_MVA");    addvar.push_back("IndexJ2_MVA");  addvar.push_back("MVA");  addvar.push_back("WMassPreFit"); addvar.push_back("MassPreFit"); addvar.push_back("PtPreFit"); addvar.push_back("EtaPreFit"); addvar.push_back("PhiPreFit"); addvar.push_back("BMPhiPreFit");  addvar.push_back("WMPhiPreFit");  addvar.push_back("TMPhiPreFit");  addvar.push_back("WBPhiPreFit"); addvar.push_back("WMassPostFit"); addvar.push_back("MassPostFit"); addvar.push_back("PtPostFit"); addvar.push_back("EtaPostFit"); addvar.push_back("PhiPostFit"); addvar.push_back("BMPhiPostFit");  addvar.push_back("WMPhiPostFit");  addvar.push_back("TMPhiPostFit");  addvar.push_back("WBPhiPostFit");  addvar.push_back("FitProb");  addvar.push_back("DPhiJet1b"); addvar.push_back("DPhiJet2b"); addvar.push_back("DRJet1b"); addvar.push_back("DRJet2b"); 

  }

  if(isResolvedTopSemiLep ){
    addvar.push_back("Pt");    addvar.push_back("Eta");    addvar.push_back("Phi");    addvar.push_back("E"); addvar.push_back("Mass");   
    addvar.push_back("MT");    addvar.push_back("LBMPhi");    addvar.push_back("LMPhi");    addvar.push_back("LBPhi");     addvar.push_back("BMPhi");  addvar.push_back("TMPhi"); 
    addvar.push_back("IndexL");    addvar.push_back("LeptonFlavour");    addvar.push_back("IndexB");
  }
  
  if(isevent){
    addvar.push_back("weight");
    addvar.push_back("nTightMuons");
    addvar.push_back("nSoftMuons");
    addvar.push_back("nLooseMuons");
    addvar.push_back("nMediumMuons");    
    addvar.push_back("nTightElectrons");
    addvar.push_back("nMediumElectrons");
    addvar.push_back("nLooseElectrons");
    addvar.push_back("nVetoElectrons");
    addvar.push_back("nElectronsSF");
    addvar.push_back("mt");
    addvar.push_back("Mt2w");
    addvar.push_back("category");
    addvar.push_back("nMuonsSF");
    addvar.push_back("nCSVTJets");
    addvar.push_back("nCSVMJets");
    addvar.push_back("nCSVLJets");
    addvar.push_back("nTightJets");
    addvar.push_back("nLooseJets");
    addvar.push_back("nType1TopJets");
    addvar.push_back("nType2TopJets");
    addvar.push_back("Ht");
    addvar.push_back("nGoodPV");
    addvar.push_back("nPV");
    addvar.push_back("nTruePV");
    addvar.push_back("Rho");

    addvar.push_back("bWeight0CSVT");
    addvar.push_back("bWeight1CSVT");
    addvar.push_back("bWeight2CSVT");
    addvar.push_back("bWeight1_2CSVT");

    addvar.push_back("bWeight0CSVM");
    addvar.push_back("bWeight1CSVM");
    addvar.push_back("bWeight2CSVM");
    addvar.push_back("bWeight1_2CSVM");

    addvar.push_back("bWeight0CSVL");
    addvar.push_back("bWeight1CSVL");
    addvar.push_back("bWeight2CSVL");
    addvar.push_back("bWeight1_2CSVL");


    addvar.push_back("bWeightMisTagDown0CSVT");
    addvar.push_back("bWeightMisTagDown1CSVT");
    addvar.push_back("bWeightMisTagDown2CSVT");
    addvar.push_back("bWeightMisTagDown1_2CSVT");

    addvar.push_back("bWeightMisTagDown0CSVM");
    addvar.push_back("bWeightMisTagDown1CSVM");
    addvar.push_back("bWeightMisTagDown2CSVM");
    addvar.push_back("bWeightMisTagDown1_2CSVM");

    addvar.push_back("bWeightMisTagDown0CSVL");
    addvar.push_back("bWeightMisTagDown1CSVL");
    addvar.push_back("bWeightMisTagDown2CSVL");
    addvar.push_back("bWeightMisTagDown1_2CSVL");

    addvar.push_back("bWeightMisTagUp0CSVT");
    addvar.push_back("bWeightMisTagUp1CSVT");
    addvar.push_back("bWeightMisTagUp2CSVT");
    addvar.push_back("bWeightMisTagUp1_2CSVT");

    addvar.push_back("bWeightMisTagUp0CSVM");
    addvar.push_back("bWeightMisTagUp1CSVM");
    addvar.push_back("bWeightMisTagUp2CSVM");
    addvar.push_back("bWeightMisTagUp1_2CSVM");

    addvar.push_back("bWeightMisTagUp0CSVL");
    addvar.push_back("bWeightMisTagUp1CSVL");
    addvar.push_back("bWeightMisTagUp2CSVL");
    addvar.push_back("bWeightMisTagUp1_2CSVL");

    addvar.push_back("bWeightBTagUp0CSVT");
    addvar.push_back("bWeightBTagUp1CSVT");
    addvar.push_back("bWeightBTagUp2CSVT");
    addvar.push_back("bWeightBTagUp1_2CSVT");

    addvar.push_back("bWeightBTagUp0CSVM");
    addvar.push_back("bWeightBTagUp1CSVM");
    addvar.push_back("bWeightBTagUp2CSVM");
    addvar.push_back("bWeightBTagUp1_2CSVM");

    addvar.push_back("bWeightBTagUp0CSVL");
    addvar.push_back("bWeightBTagUp1CSVL");
    addvar.push_back("bWeightBTagUp2CSVL");
    addvar.push_back("bWeightBTagUp1_2CSVL");

    addvar.push_back("bWeightBTagDown0CSVT");
    addvar.push_back("bWeightBTagDown1CSVT");
    addvar.push_back("bWeightBTagDown2CSVT");
    addvar.push_back("bWeightBTagDown1_2CSVT");

    addvar.push_back("bWeightBTagDown0CSVM");
    addvar.push_back("bWeightBTagDown1CSVM");
    addvar.push_back("bWeightBTagDown2CSVM");
    addvar.push_back("bWeightBTagDown1_2CSVM");

    addvar.push_back("bWeightBTagDown0CSVL");
    addvar.push_back("bWeightBTagDown1CSVL");
    addvar.push_back("bWeightBTagDown2CSVL");
    addvar.push_back("bWeightBTagDown1_2CSVL");

    //subjets

    addvar.push_back("bWeight0CSVM_subj");
    addvar.push_back("bWeight1CSVM_subj");
    addvar.push_back("bWeight2CSVM_subj");
    addvar.push_back("bWeight1_2CSVM_subj");

    addvar.push_back("bWeight0CSVL_subj");
    addvar.push_back("bWeight1CSVL_subj");
    addvar.push_back("bWeight2CSVL_subj");
    addvar.push_back("bWeight1_2CSVL_subj");

    addvar.push_back("bWeightMisTagDown0CSVM_subj");
    addvar.push_back("bWeightMisTagDown1CSVM_subj");
    addvar.push_back("bWeightMisTagDown2CSVM_subj");
    addvar.push_back("bWeightMisTagDown1_2CSVM_subj");

    addvar.push_back("bWeightMisTagDown0CSVL_subj");
    addvar.push_back("bWeightMisTagDown1CSVL_subj");
    addvar.push_back("bWeightMisTagDown2CSVL_subj");
    addvar.push_back("bWeightMisTagDown1_2CSVL_subj");

    addvar.push_back("bWeightMisTagUp0CSVM_subj");
    addvar.push_back("bWeightMisTagUp1CSVM_subj");
    addvar.push_back("bWeightMisTagUp2CSVM_subj");
    addvar.push_back("bWeightMisTagUp1_2CSVM_subj");

    addvar.push_back("bWeightMisTagUp0CSVL_subj");
    addvar.push_back("bWeightMisTagUp1CSVL_subj");
    addvar.push_back("bWeightMisTagUp2CSVL_subj");
    addvar.push_back("bWeightMisTagUp1_2CSVL_subj");

    addvar.push_back("bWeightBTagUp0CSVM_subj");
    addvar.push_back("bWeightBTagUp1CSVM_subj");
    addvar.push_back("bWeightBTagUp2CSVM_subj");
    addvar.push_back("bWeightBTagUp1_2CSVM_subj");

    addvar.push_back("bWeightBTagUp0CSVL_subj");
    addvar.push_back("bWeightBTagUp1CSVL_subj");
    addvar.push_back("bWeightBTagUp2CSVL_subj");
    addvar.push_back("bWeightBTagUp1_2CSVL_subj");

    addvar.push_back("bWeightBTagDown0CSVM_subj");
    addvar.push_back("bWeightBTagDown1CSVM_subj");
    addvar.push_back("bWeightBTagDown2CSVM_subj");
    addvar.push_back("bWeightBTagDown1_2CSVM_subj");

    addvar.push_back("bWeightBTagDown0CSVL_subj");
    addvar.push_back("bWeightBTagDown1CSVL_subj");
    addvar.push_back("bWeightBTagDown2CSVL_subj");
    addvar.push_back("bWeightBTagDown1_2CSVL_subj");

    addvar.push_back("T_Pt");
    addvar.push_back("T_Eta");
    addvar.push_back("T_Phi");
    addvar.push_back("T_E");

    addvar.push_back("Tbar_Pt");
    addvar.push_back("Tbar_Eta");
    addvar.push_back("Tbar_Phi");
    addvar.push_back("Tbar_E");

    addvar.push_back("W_Pt");
    addvar.push_back("W_Eta");
    addvar.push_back("W_Phi");
    addvar.push_back("W_E");

    addvar.push_back("Z_Pt");
    addvar.push_back("Z_Eta");
    addvar.push_back("Z_Phi");
    addvar.push_back("Z_E");

    addvar.push_back("Z_QCD_Weight");
    addvar.push_back("W_QCD_Weight");
    
    addvar.push_back("Z_Weight");
    addvar.push_back("W_Weight");

    addvar.push_back("Z_EW_Weight");
    addvar.push_back("W_EW_Weight");

    addvar.push_back("T_Weight");
    addvar.push_back("T_Ext_Weight");
    addvar.push_back("eventFlavour");

    addvar.push_back("Lepton1_Pt");
    addvar.push_back("Lepton1_Eta");
    addvar.push_back("Lepton1_Phi");
    addvar.push_back("Lepton1_E");
    addvar.push_back("Lepton1_Charge");
    addvar.push_back("Lepton1_Flavour");

    addvar.push_back("Lepton2_Pt");
    addvar.push_back("Lepton2_Eta");
    addvar.push_back("Lepton2_Phi");
    addvar.push_back("Lepton2_E");
    addvar.push_back("Lepton2_Charge");
    addvar.push_back("Lepton2_Flavour");
    
    addvar.push_back("LHEWeightSign");
    //addvar.push_back("LHEWeightAVG");
    addvar.push_back("LHEWeight");
    addvar.push_back("EventNumber");
    addvar.push_back("LumiBlock");
    addvar.push_back("RunNumber");
    
    for (size_t j = 0; j < (size_t)jetScanCuts.size(); ++j){
      stringstream j_n;
      double jetval = jetScanCuts.at(j);
      j_n << "Cut" <<jetval;
      
      addvar.push_back("nJets"+j_n.str());
      addvar.push_back("nCSVTJets"+j_n.str());
      addvar.push_back("nCSVMJets"+j_n.str());
      addvar.push_back("nCSVLJets"+j_n.str());
      addvar.push_back("bWeight1CSVTWeight"+j_n.str());
      addvar.push_back("bWeight2CSVTWeight"+j_n.str());
      addvar.push_back("bWeight1CSVMWeight"+j_n.str());
      addvar.push_back("bWeight2CSVMWeight"+j_n.str());
      addvar.push_back("bWeight1CSVLWeight"+j_n.str());
      addvar.push_back("bWeight2CSVLWeight"+j_n.str());
      addvar.push_back("bWeight0CSVLWeight"+j_n.str());
      addvar.push_back("bWeight0CSVLWeight"+j_n.str());
    }

    if(useLHEWeights){
      for (size_t w = 0; w <= (size_t)maxWeights; ++w)  {
	stringstream w_n;
	w_n << w;
	addvar.push_back("LHEWeight"+w_n.str());
      }
    }
    if(addLHAPDFWeights){
      for (size_t p = 1; p <= (size_t)maxPdf; ++p)  {
	//cout << " pdf # " << pdf_weights[p - 1] << " test "<< endl; 
	stringstream w_n;
	w_n << p;
	addvar.push_back("PDFWeight" + w_n.str());
      }
    }
    if(useMETFilters){
      for (size_t lt = 0; lt < metFilters.size(); ++lt)  {
	string trig = metFilters.at(lt);
	addvar.push_back("passes"+trig);
      }
      addvar.push_back("passesMETFilters");
    }
    if(useTriggers){
      for (size_t lt = 0; lt < SingleElTriggers.size(); ++lt)  {
	string trig = SingleElTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < SingleMuTriggers.size(); ++lt)  {
	string trig = SingleMuTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < PhotonTriggers.size(); ++lt)  {
	string trig = PhotonTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t ht = 0; ht < hadronicTriggers.size(); ++ht)  {
	string trig = hadronicTriggers.at(ht);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFHT900Triggers.size(); ++lt)  {
        string trig = HadronPFHT900Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFHT800Triggers.size(); ++lt)  {
        string trig = HadronPFHT800Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFJet450Triggers.size(); ++lt)  {
        string trig = HadronPFJet450Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      
      addvar.push_back("passesSingleElTriggers");
      addvar.push_back("passesSingleMuTriggers");
      addvar.push_back("passesPhotonTriggers");
      addvar.push_back("passesHadronicTriggers");
      addvar.push_back("passesHadronPFHT900Triggers");
      addvar.push_back("passesHadronPFHT800Triggers");
      addvar.push_back("passesHadronPFJet450Triggers");
    }
  }
  
  return addvar;
}

void DMAnalysisTreeMaker::initializePdf(string central, string varied){

    if(central == "CT") {  LHAPDF::initPDFSet(1, "cteq66.LHgrid"); }
    if(central == "CT10") {  LHAPDF::initPDFSet(1, "CT10.LHgrid"); }
    if(central == "CT10f4") {  LHAPDF::initPDFSet(1, "CT10f4.LHgrid"); }
    if(central == "NNPDF") { LHAPDF::initPDFSet(1, "NNPDF21_100.LHgrid");  }
    if(central == "MSTW") { LHAPDF::initPDFSet(1, "MSTW2008nlo68cl.LHgrid");  }

    if(varied == "CT") {  LHAPDF::initPDFSet(2, "cteq66.LHgrid"); maxPdf = 44; }
    if(varied == "CT10") {  LHAPDF::initPDFSet(2, "CT10.LHgrid"); maxPdf = 52; }
    if(varied == "CT10f4") {  LHAPDF::initPDFSet(2, "CT10f4.LHgrid"); maxPdf = 52; }
    if(varied == "NNPDF") { LHAPDF::initPDFSet(2, "NNPDF21_100.LHgrid");  maxPdf = 100; }
    if(varied == "MSTW") { LHAPDF::initPDFSet(2, "MSTW2008nlo68cl.LHgrid"); maxPdf = 40; }

}

double DMAnalysisTreeMaker::getWEWKPtWeight(double ptW){ // TO DO?, in aggrement with ttDM PAS, Table 6
  //EWK

  if(ptW<150.)return 0.980859;
  if(ptW>=150. && ptW <200.)return 0.962119;
  if(ptW>=200. && ptW <250.)return 0.944429;
  if(ptW>=250. && ptW <300.)return 0.927686;
  if(ptW>=300. && ptW <350.)return 0.911802;
  if(ptW>=350. && ptW <400.)return 0.8967;
  if(ptW>=400. && ptW <500.)return 0.875368;
  if(ptW>=500. && ptW <600.)return 0.849097;
  if(ptW>=600. && ptW <1000)return 0.792159;
  return 1.0;
}

double DMAnalysisTreeMaker::getZEWKPtWeight(double ptW){ //TO DO?, in aggrement with ttDM PAS, Table 6
  if(ptW<150.)return 0.984525;
  if(ptW>=150. && ptW <200.)return 0.969079;
  if(ptW>=200. && ptW <250.)return 0.954627;
  if(ptW>=250. && ptW <300.)return 0.941059;
  if(ptW>=300. && ptW <350.)return 0.928284;
  if(ptW>=350. && ptW <400.)return 0.91622;
  if(ptW>=400. && ptW <500.)return 0.899312;
  if(ptW>=500. && ptW <600.)return 0.878693;
  if(ptW>=600. && ptW <1000)return 0.834718;
  return 1.0;

}

double DMAnalysisTreeMaker::getWPtWeight(double ptW){ //TO DO?, in aggrement with ttDM PAS, Table 5
  //QCD
  
  if(ptW<150.)return 1.89123;
  if(ptW>=150. && ptW <200.)return 1.70414;
  if(ptW>=200. && ptW <250.)return 1.60726;
  if(ptW>=250. && ptW <300.)return 1.57206;
  if(ptW>=300. && ptW <350.)return 1.51689;
  if(ptW>=350. && ptW <400.)return 1.4109;
  if(ptW>=400. && ptW <500.)return 1.30758;
  if(ptW>=500. && ptW <600.)return 1.32046;
  if(ptW>=600. && ptW <1000)return 1.26853;
  return 1.0;
}

double DMAnalysisTreeMaker::getAPtWeight(double ptW){ //TO DO, ma che e'?
  if(ptW<150.)return 1.24087;
  if(ptW>=150. && ptW <200.)return 1.55807;
  if(ptW>=200. && ptW <250.)return 1.51043;
  if(ptW>=250. && ptW <300.)return 1.47333;
  if(ptW>=300. && ptW <350.)return 1.43497;
  if(ptW>=350. && ptW <400.)return 1.37846;
  if(ptW>=400. && ptW <500.)return 1.29202;
  if(ptW>=500. && ptW <600.)return 1.31414;
  if(ptW>=600.)return 1.20454;
  return 1.0;
}

double DMAnalysisTreeMaker::getZPtWeight(double ptW){//TO DO?, in aggrement with ttDM PAS, Table 5
  if(ptW<150.)return 1.685005;
  if(ptW>=150. && ptW <200.)return 1.552560;
  if(ptW>=200. && ptW <250.)return 1.522595;
  if(ptW>=250. && ptW <300.)return 1.520624;
  if(ptW>=300. && ptW <350.)return 1.432282;
  if(ptW>=350. && ptW <400.)return 1.457417;
  if(ptW>=400. && ptW <500.)return 1.368499;
  if(ptW>=500. && ptW <600.)return 1.358024;
  if(ptW>=600.)return 1.164847;;
  return 1.0;
}

double DMAnalysisTreeMaker::getTopPtWeight(double ptT, double ptTbar, bool extrap){ //TO DO
  if((ptT>0.0 && ptTbar>0.0) ){
    if (extrap || (ptT<=400.0 && ptTbar <=400.0)){
      double a = 0.0615;
      double b = -0.0005;
      double sfT = exp(a+b*ptT);
      double sfTbar = exp(a+b*ptTbar);
    return sqrt(sfT*sfTbar); 
    }
  }
  return 1.0;
}

bool DMAnalysisTreeMaker::getMETFilters(){
  bool METFilterAND=true;
  for(size_t mf =0; mf< metFilters.size();++mf){
    string fname = metFilters.at(mf);
    for(size_t bt = 0; bt < metNames->size();++bt){
      std::string tname = metNames->at(bt);
      //      cout << "test tname "<<endl;
      if(tname.find(fname)!=std::string::npos){
	METFilterAND = METFilterAND && (metBits->at(bt)>0);
	float_values["Event_passes"+fname]=metBits->at(bt);
      }
    }
  }
  float_values["Event_passesMETFilters"]=(float)METFilterAND;
  return METFilterAND;
}

bool DMAnalysisTreeMaker::getEventTriggers(){
  bool eleOR=false, muOR=false, hadronOR=false, phOR=false, H900OR=false, H800OR=false, J450OR=false;
  for(size_t lt =0; lt< SingleElTriggers.size();++lt){
    string lname = SingleElTriggers.at(lt);
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(lname)!=std::string::npos){
	eleOR = muOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t bt = 0; bt < triggerNamesR->size();++bt){
    std::string tname = triggerNamesR->at(bt);
    //    std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;
  }
  
  for(size_t lt =0; lt< SingleMuTriggers.size();++lt){
    string lname = SingleMuTriggers.at(lt);
    //    std::cout << lname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;
      if(tname.find(lname)!=std::string::npos){
	//	cout << " matches "<<endl;
	muOR = muOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }
  for(size_t pt =0; pt< PhotonTriggers.size();++pt){
    string pname = PhotonTriggers.at(pt);
    //std::cout << pname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(pname)!=std::string::npos){
	phOR = phOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+pname]=triggerBits->at(bt);
	float_values["Event_prescale"+pname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t ht =0; ht< hadronicTriggers.size();++ht){
    string hname = hadronicTriggers.at(ht);
    //std::cout << hname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(hname)!=std::string::npos){
	//bool before = hadronOR;
	hadronOR = hadronOR || (triggerBits->at(bt)>0);
	//bool after = hadronOR;
	//if(before != after){
	//std::cout<< "hadr name "<< hname << std::endl;
	//std::cout<< "trig name "<< tname << std::endl;
	//std::cout << "hadronOR before " << before << std::endl;
	//std::cout << "hadronOR after " << after << std::endl;
	//}
	//hadronOR = hadronOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+hname]=triggerBits->at(bt);
	float_values["Event_prescale"+hname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t lt =0; lt< HadronPFHT900Triggers.size();++lt){
    string lname = HadronPFHT900Triggers.at(lt);
    //    std::cout << lname << std::endl;                                                                                                                                                
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                         
      if(tname.find(lname)!=std::string::npos){
        //      cout << " matches "<<endl;                                                                                                                                                
        H900OR = H900OR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
        float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }
  for(size_t lt =0; lt< HadronPFHT800Triggers.size();++lt){
    string lname = HadronPFHT800Triggers.at(lt);
    //    std::cout << lname << std::endl;                                                                                                                                         
                                                                                                                                                                                    
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                   
                                                                                                                                                                                       
      if(tname.find(lname)!=std::string::npos){
	//      cout << " matches "<<endl;                                                                                                                                        
                                                                                                                                                                                   
	H800OR = H800OR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t lt =0; lt< HadronPFJet450Triggers.size();++lt){
    string lname = HadronPFJet450Triggers.at(lt);
    //    std::cout << lname << std::endl;                                                                                                                                            
                                                                                                                                                                                     
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                     
                                                                                                                                                                                      
      if(tname.find(lname)!=std::string::npos){
	//      cout << " matches "<<endl;                                                                                                                                                                                                                                                                                                                                     
	J450OR = J450OR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  float_values["Event_passesSingleElTriggers"]=(float)eleOR;
  float_values["Event_passesSingleMuTriggers"]=(float)muOR;
  float_values["Event_passesPhotonTriggers"]=(float)phOR;
  float_values["Event_passesHadronicTriggers"]=(float)hadronOR;
  float_values["Event_passesHadronicPFHT900Triggers"]=(float)H900OR;
  float_values["Event_passesHadronicPFHT800Triggers"]=(float)H800OR;
  float_values["Event_passesHadronicPFJet450Triggers"]=(float)J450OR;
  return (eleOR || muOR || hadronOR || phOR || H900OR || H800OR || J450OR);
}


void DMAnalysisTreeMaker::getEventPdf(){

  //  std::cout << " getting pdf "<<endl;

  double scalePDF = genprod->pdf()->scalePDF;
  double x1 =  genprod->pdf()->x.first;
  double x2 =  genprod->pdf()->x.second;
  int id1 =  genprod->pdf()->id.first;
  int id2 =  genprod->pdf()->id.second;

  //  std::cout << " maxpdf "<< maxPdf << " accessing x1 " << x1<< id1<<std::endl;


  LHAPDF::usePDFMember(1, 0);
  double xpdf1 = LHAPDF::xfx(1, x1, scalePDF, id1);
  double xpdf2 = LHAPDF::xfx(1, x2, scalePDF, id2);
  double w0 = xpdf1 * xpdf2;
  int maxPDFCount = maxPdf;

  std::cout << "entering pdf loop" <<std::endl;
  for (int p = 1; p <= maxPdf; ++p)
    {
      
      if ( p > maxPDFCount ) continue;
      LHAPDF::usePDFMember(2, p);
      double xpdf1_new = LHAPDF::xfx(2, x1, scalePDF, id1);
      double xpdf2_new = LHAPDF::xfx(2, x2, scalePDF, id2);
      double pweight = xpdf1_new * xpdf2_new / w0;
      stringstream w_n;
      w_n << p;
      float_values["PDFWeight"+w_n.str()]= pweight;
    }
  
}


void DMAnalysisTreeMaker::getEventLHEWeights(){
  //  std::cout << " in weight "<<endl;
  size_t wgtsize=  lhes->weights().size();
  //  std::cout << "weight size "<< wgtsize<<endl;
  for (size_t i = 0; i <  wgtsize; ++i)  {
    if (i<= (size_t)maxWeights){ 
      stringstream w_n;
      w_n << i;

      float ww = (float)lhes->weights().at(i).wgt;
      
      //      cout << "ww # " << i<< "is "<<ww <<endl;
      //      cout << "id  is "<< std::string(lhes->weights().at(i).id.data()) <<endl;
      //      cout <<" floatval before "<< float_values["Event_LHEWeight"+w_n.str()]<<endl;

      float_values["Event_LHEWeight"+w_n.str()]= ww;
      //if(i>=11)float_values["Event_LHEWeightAVG"]+= ww;

      //      cout <<" floatval after "<< float_values["Event_LHEWeight"+w_n.str()]<<endl;

    }
    //    float_values["Event_LHEWeightAVG"]+= ww;

    //    else cout << "WARNING! there are " << wgtsize << " weights, and you accomodated for only "<< maxWeights << " weights, check your configuration file/your lhe!!!"<<endl;
  }
  
}

void DMAnalysisTreeMaker::initTreeWeightHistory(bool useLHEW){
  cout << " preBranch "<<endl;

  trees["WeightHistory"]->Branch("Event_Z_EW_Weight",&float_values["Event_Z_EW_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_EW_Weight",&float_values["Event_W_EW_Weight"]);
  trees["WeightHistory"]->Branch("Event_Z_QCD_Weight",&float_values["Event_Z_QCD_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_QCD_Weight",&float_values["Event_W_QCD_Weight"]);
  trees["WeightHistory"]->Branch("Event_Z_Weight",&float_values["Event_Z_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_Weight",&float_values["Event_W_Weight"]);

  trees["WeightHistory"]->Branch("Event_T_Weight",&float_values["Event_T_Weight"]);
  trees["WeightHistory"]->Branch("Event_T_Ext_Weight",&float_values["Event_T_Ext_Weight"]);
  trees["WeightHistory"]->Branch("Event_T_Pt",&float_values["Event_T_Pt"]);
  trees["WeightHistory"]->Branch("Event_Tbar_Pt",&float_values["Event_Tbar_Pt"]);

  trees["WeightHistory"]->Branch("Event_W_Pt",&float_values["Event_W_Pt"]);
  trees["WeightHistory"]->Branch("Event_Z_Pt",&float_values["Event_Z_Pt"]);
  
  cout << " preBranch Weight "<<endl;

  //  size_t wgtsize=  lhes->weights().size();
  //  std::cout << "weight size "<< wgtsize<<endl;
  if(useLHEW){
    for (size_t w = 1; w <=  (size_t)maxWeights; ++w)  {
      stringstream w_n;
      w_n << w;
      string name = "Event_LHEWeight"+w_n.str();
      cout << " pre single w # "<< w <<endl;
      trees["WeightHistory"]->Branch(name.c_str(),&float_values[name],(name+"/F").c_str());
      //trees["noSyst"]->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
    }
  }
}


double DMAnalysisTreeMaker::smear(double pt, double genpt, double eta, string syst){
  double resolScale = resolSF(fabs(eta), syst);
  double smear =1.0;
  if(genpt>0) smear = std::max((double)(0.0), (double)(pt + (pt - genpt) * resolScale) / pt);
  return  smear;
}

double DMAnalysisTreeMaker::resolSF(double eta, string syst)//DONE
{//from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  double fac = 0.;
  if (syst == "jer__up")fac = 1.;
  if (syst == "jer__down")fac = -1.;
  if (eta <= 0.5)                       return 0.122 + (0.026 * fac);
  else if ( eta > 0.5 && eta <= 0.8 )   return 0.167 + (0.048 * fac);
  else if ( eta > 0.8 && eta <= 1.1 )   return 0.168 + (0.046 * fac);
  else if ( eta > 1.1 && eta <= 1.3 )   return 0.029 + (0.066 * fac);
  else if ( eta > 1.3 && eta <= 1.7 )   return 0.115 + (0.003 * fac);
  else if ( eta > 1.7 && eta <= 1.9 )   return 0.041 + (0.0062 * fac);
  else if ( eta > 1.9 && eta <= 2.1 )   return 0.167 + (0.086 * fac);
  else if ( eta > 2.1 && eta <= 2.3 )   return 0.094 + (0.093 *fac);
  else if ( eta > 2.3 && eta <= 2.5 )   return 0.168 + (0.120 *fac);
  else if ( eta > 2.5 && eta <= 2.8 )   return 0.266 + (0.132 *fac);
  else if ( eta > 2.8 && eta <= 3.0 )   return 0.595 + (0.175 *fac);
  else if ( eta > 3.0 && eta <= 3.2 )   return 0.998 + (0.066 *fac);
  else if ( eta > 3.2 && eta <= 5.0 )   return 0.226 + (0.145 *fac);
  return 0.1;
 }

double DMAnalysisTreeMaker::getEffectiveArea(string particle, double eta){ //TO DO
  double aeta = fabs(eta);
  if(particle=="photon"){
    if(aeta<1.0)return 0.0725;
    if(aeta<1.479 && aeta >1.0)return 0.0604;
    if(aeta<2.0 && aeta >1.479)return 0.0320;
    if(aeta<2.2 && aeta >2.0)return 0.0512;
    if(aeta<2.3 && aeta >2.2)return 0.0766;
    if(aeta<2.4 && aeta >2.3)return 0.0949;
    if(aeta>2.4)return 0.1160;
  }
  if(particle=="ch_hadron"){
    if(aeta<1.0)return 0.0157;
    if(aeta<1.479 && aeta >1.0)return 0.0143;
    if(aeta<2.0 && aeta >1.479)return 0.0115;
    if(aeta<2.2 && aeta >2.0)return 0.0094;
    if(aeta<2.3 && aeta >2.2)return 0.0095;
    if(aeta<2.4 && aeta >2.3)return 0.0068;
    if(aeta>2.4)return 0.0053;
  }
  if(particle=="neu_hadron"){
    if(aeta<1.0)return 0.0143;
    if(aeta<1.479 && aeta >1.0)return 0.0210;
    if(aeta<2.0 && aeta >1.479)return 0.0147;
    if(aeta<2.2 && aeta >2.0)return 0.0082;
    if(aeta<2.3 && aeta >2.2)return 0.0124;
    if(aeta<2.4 && aeta >2.3)return 0.0186;
    if(aeta>2.4)return 0.0320;
  }
  return 0.0;


};
double DMAnalysisTreeMaker::jetUncertainty(double ptCorr, double eta, string syst)
{
  if(ptCorr<0)return ptCorr;
  if(syst == "jes__up" || syst == "jes__down"){
    double fac = 1.;
    if (syst == "jes__down")fac = -1.;
    jecUnc->setJetEta(eta);
    jecUnc->setJetPt(ptCorr);
    double JetCorrection = jecUnc->getUncertainty(true);
    return JetCorrection*fac;
  }
  return 0.0;
}

double DMAnalysisTreeMaker::jetUncertainty8(double ptCorr, double eta, string syst)
{
  if(ptCorr<0)return ptCorr;
  if(syst == "jes__up" || syst == "jes__down"){
    double fac = 1.;
    if (syst == "jes__down")fac = -1.;
    jecUnc8->setJetEta(eta);
    jecUnc8->setJetPt(ptCorr);
    double JetCorrection = jecUnc8->getUncertainty(true);
    return JetCorrection*fac;
  }
  return 0.0;
}

double DMAnalysisTreeMaker::getScaleFactor(double ptCorr,double etaCorr,double partonFlavour, string syst){
  return 1.0;
}


bool DMAnalysisTreeMaker::isMCWeightName(string s){
  
  if(s=="Z_Weight")return true;
  if(s=="W_Weight")return true;

  if(s=="Z_QCD_Weight")return true;
  if(s=="W_QCD_Weight")return true;

  if(s=="Z_EW_Weight")return true;
  if(s=="W_EW_Weight")return true;
  
  if(s=="T_Weight")return true;
  if(s=="T_Ext_Weight")return true;
    
  return false;

}

bool DMAnalysisTreeMaker::isInVector(std::vector<std::string> v, std::string s){
  for(size_t i = 0;i<v.size();++i){
    if(v.at(i)==s)return true;
  }
  return false;
}


//BTag weighter
bool DMAnalysisTreeMaker::BTagWeight::filter(int t)
{
    return (t >= minTags && t <= maxTags);
}

float DMAnalysisTreeMaker::BTagWeight::weight(vector<JetInfo> jetTags, int tags)
{
    if (!filter(tags))
    {
      //std::cout << "nThis event should not pass the selection, what is it doing here?" << std::endl;
        return 0;
    }
    //cout << "ciao" << endl;
    int njetTags = jetTags.size();
    //    cout<< " njettags "<< njetTags<<endl;
    int comb = 1 << njetTags;
    float pMC = 0;
    float pData = 0;
    for (int i = 0; i < comb; i++)
    {
        float mc = 1.;
        float data = 1.;
        int ntagged = 0;
        for (int j = 0; j < njetTags; j++)
        {
            bool tagged = ((i >> j) & 0x1) == 1;
            if (tagged)
            {
                ntagged++;
                mc *= jetTags[j].eff;
                data *= jetTags[j].eff * jetTags[j].sf;
            }
            else
            {
                mc *= (1. - jetTags[j].eff);
                data *= (1. - jetTags[j].eff * jetTags[j].sf);
            }
        }

        if (filter(ntagged))
        {
	  //	  std::cout << mc << " " << data << endl;
            pMC += mc;
            pData += data;
        }
    }

    if (pMC == 0) return 0;
    return pData / pMC;
}


double DMAnalysisTreeMaker::MCTagEfficiencySubjet(string algo, int flavor, double pt, double eta){
if ((abs(flavor) ==5)){ 
	 if (algo=="csvt") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.296489; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.283186; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.22293; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.345938; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.345226; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.334581; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.390734; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.392072; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.364133; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.416481; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.400909; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.380812; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.407575; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.409937; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.350389; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.389089; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.380354; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.312596; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.358209; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.346685; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.264417; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.270419; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.268743; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.206701; 
	 }
	 if (algo=="csvm") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.467274; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.478342; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.382431; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.509244; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.525654; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.514811; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.57471; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.585932; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.542646; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.611787; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.603202; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.554784; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.602366; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.608457; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.545904; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.592466; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.587298; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.513208; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.566559; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.558705; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.472549; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.491582; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.482574; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.415809; 
	 }
	 if (algo=="csvl") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.626355; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.6381; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.609076; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.635854; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.648998; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.684752; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.7; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.707276; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.708857; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.734967; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.731502; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.720738; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.734769; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.742253; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.745629; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.752044; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.752473; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.755872; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.744627; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.7484; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.76149; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.740626; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.736961; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.757653; 
	 }
}
 if ((abs(flavor) ==4)){ 
   if (algo=="csvt") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.0305344; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.0189983; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0175439; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.0308756; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.0314319; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.0321429; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.0402351; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.032538; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.0250804; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.0392806; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.0303384; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.032345; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.0380644; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.0400807; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.0308721; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.0408119; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.0315339; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.0374813; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.0329393; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.0369307; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.0279534; 
   }
	 if (algo=="csvm") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.115522; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.102476; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0845704; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.147465; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.142026; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.143452; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.159584; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.154013; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.13119; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.169901; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.163944; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.147574; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.164608; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.174439; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.138765; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.160192; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.142437; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.141829; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.142984; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.148338; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.130893; 
	}
	 if (algo=="csvl") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.295165; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.287277; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.278453; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.337327; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.335274; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.33631; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.350362; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.339479; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.358842; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.368197; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.373979; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.366577; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.358575; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.368036; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.372056; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.364688; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.363976; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.396102; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.380335; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.393311; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.42629; 
	}
}
 if ((abs(flavor)!= 5 && (abs(flavor)!= 4 ))){ 
   if (algo=="csvt") {
     if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.0101643; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.00741036; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.00463209; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.00674189; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.0056723; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.00512842; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.00625346; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.00420484; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.00454994; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.00560886; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.00504171; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.00524952; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.0061885; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.00530979; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.00545241; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.00512219; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.00430636; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.00391961; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.00515681; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.00461102; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.00430911; 
   }
   if (algo=="csvm") {
     if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.0352963; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.0317506; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0222538; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.0377191; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.0332173; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.0299433; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.0353044; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.0293838; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.0293827; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.0345339; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.0321109; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.0305427; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.0303497; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.0291583; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.0312996; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.0247239; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.0247349; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.0250143; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.0262878; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.0251433; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.0286476; 
   }
   if (algo=="csvl") {
     if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.199172; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.202331; 
     if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.222045; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.136186; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.135791; 
     if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.151454; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.128908; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.13025; 
     if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.156014; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.123812; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.131985; 
     if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.158645; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.116918; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.124735; 
     if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.162019; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.117095; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.129155; 
     if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.181763; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.144244; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.153178; 
     if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.225239; 
   }
 }
 return 1.0;
}

double DMAnalysisTreeMaker::MCTagEfficiency(string algo, int flavor, double pt, double eta){
if ((abs(flavor) ==5)){ 
	 if (algo=="csvt") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.352168; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.319601; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.23781; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.526553; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.497584; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.417645; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.541338; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.528332; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.440321; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.543049; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.521693; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.436618; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.528846; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.510575; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.413917; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.482124; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.462141; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.364861; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.434455; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.411495; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.301717; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.296532; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.305188; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.224686; 
	 }
	 if (algo=="csvm") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.529907; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.525015; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.407325; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.718092; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.70719; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.62003; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.742626; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.740707; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.642757; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.749783; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.740894; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.641026; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.743477; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.727679; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.631493; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.704655; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.691343; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.593679; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.672252; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.650178; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.537782; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.538615; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.538796; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.45657; 
	 }
	 if (algo=="csvl") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.722738; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.721697; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.642741; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.845194; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.838713; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.800013; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.861446; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.86854; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.820012; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.868555; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.869711; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.817321; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.870079; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.861107; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.831645; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.861654; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.854613; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.845222; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.854853; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.843449; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.830373; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  0) && (pt < 40))) return 0.807141; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  40) && (pt < 60))) return 0.802808; 
		if(((eta>=5.6) && (eta <6.4)) && ((pt >=  60) && (pt < 80))) return 0.801898; 
	 }
}
if ((abs(flavor) ==4)){ 
	 if (algo=="csvt") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.0279674; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.0217746; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0177389; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.0460526; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.0357602; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.0354151; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.0430039; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.0390738; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.0383872; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.0457265; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.0422265; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.037679; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.0455882; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.0398756; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.0437271; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.0445344; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.0435157; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.0494012; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.041203; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.0440075; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.0340389; 
	 }
	 if (algo=="csvm") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.123621; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.116604; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0851939; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.194646; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.173913; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.155941; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.185173; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.198263; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.167081; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.196581; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.208253; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.165411; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.198284; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.195136; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.165522; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.178506; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.192159; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.178518; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.161718; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.178047; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.152218; 
	}
	 if (algo=="csvl") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.338892; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.326619; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.317029; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.404719; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.414716; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.403656; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.41656; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.413169; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.41066; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.426068; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.425144; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.405049; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.414951; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.426753; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.424679; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.407067; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.420939; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.454341; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.432387; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.458657; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.489522; 
	}
}
 if ((abs(flavor)!= 5 && (abs(flavor)!= 4 ))){ 
	 if (algo=="csvt") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.00576088; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.00446663; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.00241878; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.00399483; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.00426005; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.00308828; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.00272815; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.00338576; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.00254755; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.00244485; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.00300245; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.00272979; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.00196052; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.00253911; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.00270194; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.00200782; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.00174762; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.00318708; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.00305354; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.00392678; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.00398442; 
	 }
	 if (algo=="csvm") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.0191537; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.018406; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.0141139; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.0169362; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.020136; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.0225764; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.0158694; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.0161145; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.0198304; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.0142888; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.0159926; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.0216754; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.0133914; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.0159441; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.0216155; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.0117827; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.0163915; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.02327; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.0212436; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.0252931; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.0322783; 
	 }
	 if (algo=="csvl") {
		if(((eta>=0) && (eta <0.6)) && ((pt >=  0) && (pt < 40))) return 0.240319; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  40) && (pt < 60))) return 0.255961; 
		if(((eta>=0) && (eta <0.6)) && ((pt >=  60) && (pt < 80))) return 0.279793; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  0) && (pt < 40))) return 0.0960913; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  40) && (pt < 60))) return 0.117138; 
		if(((eta>=0.6) && (eta <1.2)) && ((pt >=  60) && (pt < 80))) return 0.16714; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  0) && (pt < 40))) return 0.0837656; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  40) && (pt < 60))) return 0.10333; 
		if(((eta>=1.2) && (eta <2.4)) && ((pt >=  60) && (pt < 80))) return 0.153027; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  0) && (pt < 40))) return 0.0856786; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  40) && (pt < 60))) return 0.102328; 
		if(((eta>=2.4) && (eta <3.2)) && ((pt >=  60) && (pt < 80))) return 0.15678; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  0) && (pt < 40))) return 0.086396; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  40) && (pt < 60))) return 0.102759; 
		if(((eta>=3.2) && (eta <4)) && ((pt >=  60) && (pt < 80))) return 0.161679; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  0) && (pt < 40))) return 0.0964282; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  40) && (pt < 60))) return 0.114439; 
		if(((eta>=4) && (eta <4.8)) && ((pt >=  60) && (pt < 80))) return 0.194674; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  0) && (pt < 40))) return 0.139306; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  40) && (pt < 60))) return 0.160493; 
		if(((eta>=4.8) && (eta <5.6)) && ((pt >=  60) && (pt < 80))) return 0.245989; 
	 }
}
  return 1.0;
  }
    
double DMAnalysisTreeMaker::TagScaleFactorSubjet(string algo, int flavor, string syst, double pt){ 
  double x = pt;
  if(algo == "csvl"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 0.935202;
	if (pt >= 140  && pt < 180) return 0.975405;
	if (pt >= 180  && pt < 240) return 0.944102;
	if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 0.935202;
	if (pt >= 140  && pt < 180) return 0.975405;
	if (pt >= 180  && pt < 240) return 0.944102;
	if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if ( pt>=20 && pt < 1000) return ((0.917093+(0.00164417*x))+(-2.68598e-06*(x*x)))+(1.3727e-09*(x*(x*x)));
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 140) return 0.935202;
        if (pt >= 140  && pt < 180) return 0.975405;
        if (pt >= 180  && pt < 240) return 0.944102;
        if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 140) return 0.935202;
        if (pt >= 140  && pt < 180) return 0.975405;
        if (pt >= 180  && pt < 240) return 0.944102;
        if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt < 1000) return ((0.960807+(0.0018048*x))+(-3.09836e-06*(x*x)))+(1.64823e-09*(x*(x*x)));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 140) return 0.935202;
        if (pt >= 140  && pt < 180) return 0.975405;
        if (pt >= 180  && pt < 240) return 0.944102;
        if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 140) return 0.935202;
        if (pt >= 140  && pt < 180) return 0.975405;
        if (pt >= 180  && pt < 240) return 0.944102;
        if (pt >= 240  && pt < 450) return 0.965785;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt < 1000) return ((0.873402+(0.00148169*x))+(-2.27029e-06*(x*x)))+(1.09633e-09*(x*(x*x)));
      }
    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 1.07507337034;
	if (pt >= 140  && pt < 180) return 1.05887537845;
	if (pt >= 180  && pt < 240) return 1.02505525203;
	if (pt >= 240  && pt < 450) return 1.02947460786;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 1.21494474069;
	if (pt >= 140  && pt < 180) return 1.1423457569;
	if (pt >= 180  && pt < 240) return 1.10600850407;
	if (pt >= 240  && pt < 450) return 1.09316421573;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>= 20 && pt < 1000) return ((0.960807+(0.0018048*x))+(-3.09836e-06*(x*x)))+(1.64823e-09*(x*(x*x)));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 0.795325287574;
	if (pt >= 140  && pt < 180) return 0.89192934327;
	if (pt >= 180  && pt < 240) return 0.863145180357;
	if (pt >= 240  && pt < 450) return 0.902127534532;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 0.655448575149;
	if (pt >= 140  && pt < 180) return 0.808453686539;
	if (pt >= 180  && pt < 240) return 0.782188360714;
	if (pt >= 240  && pt < 450) return 0.838470069063;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>= 20 && pt < 1000) return ((0.873402+(0.00148169*x))+(-2.27029e-06*(x*x)))+(1.09633e-09*(x*(x*x)));
      }
    } 
  }
    if(algo == "csvm"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 0.889127;
	if (pt >= 140  && pt < 180) return 0.931134;
	if (pt >= 180  && pt < 240) return 0.909904;
	if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 0.889127;
	if (pt >= 140  && pt < 180) return 0.931134;
	if (pt >= 180  && pt < 240) return 0.909904;
	if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if ( pt>=20 && pt < 1000) return ((0.859405+(-0.0011568*x))+(4.62759e-06*(x*x)))+(-2.96878e-09*(x*(x*x)));
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 140) return 0.889127;
        if (pt >= 140  && pt < 180) return 0.931134;
        if (pt >= 180  && pt < 240) return 0.909904;
        if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 140) return 0.889127;
        if (pt >= 140  && pt < 180) return 0.931134;
        if (pt >= 180  && pt < 240) return 0.909904;
        if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt < 1000) return ((0.909612+(-0.000983575*x))+(4.37994e-06*(x*x)))+(-2.81911e-09*(x*(x*x)));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 140) return 0.889127;
        if (pt >= 140  && pt < 180) return 0.931134;
        if (pt >= 180  && pt < 240) return 0.909904;
        if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 140) return 0.889127;
        if (pt >= 140  && pt < 180) return 0.931134;
        if (pt >= 180  && pt < 240) return 0.909904;
        if (pt >= 240  && pt < 450) return 0.899026;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt < 1000) return ((0.809148+(-0.00132666*x))+(4.86772e-06*(x*x)))+(-3.11608e-09*(x*(x*x)));
      }
    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 1.08003212386;
	if (pt >= 140  && pt < 180) return 1.05432375412;
	if (pt >= 180  && pt < 240) return 1.01894082837;
	if (pt >= 240  && pt < 450) return 1.01737126216;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 1.27093724771;
	if (pt >= 140  && pt < 180) return 1.17751350825;
	if (pt >= 180  && pt < 240) return 1.12797765674;
	if (pt >= 240  && pt < 450) return 1.13571652433;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>= 20 && pt < 1000) return ((0.909612+(-0.000983575*x))+(4.37994e-06*(x*x)))+(-2.81911e-09*(x*(x*x)));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 140) return 0.698203917937;
	if (pt >= 140  && pt < 180) return 0.807996622827;
	if (pt >= 180  && pt < 240) return 0.80100703546;
	if (pt >= 240  && pt < 450) return 0.780740827981;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 140) return 0.507280835873;
	if (pt >= 140  && pt < 180) return 0.684859245653;
	if (pt >= 180  && pt < 240) return 0.69211007092;
	if (pt >= 240  && pt < 450) return 0.662455655963;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>= 20 && pt < 1000) return ((0.809148+(-0.00132666*x))+(4.86772e-06*(x*x)))+(-3.11608e-09*(x*(x*x)));
      }
    } 
  }
    return 1.0;
}


double DMAnalysisTreeMaker::TagScaleFactor(string algo, int flavor, string syst, double pt){ //DONE, to update with mistag
  double x = pt;
if(algo == "csvt"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 670) return (0.857294+(3.75846e-05*x));
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 670) return (0.857294+(3.75846e-05*x));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if( pt>=20 && pt < 1000) return 0.688619+260.84/(x*x);
      }
    }
     if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return (0.857294+(3.75846e-05*x));
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return (0.857294+(3.75846e-05*x));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.688619+260.84/(x*x))*(1+(0.144982+0.000116685*x+-1.0021e-07*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670)  return 0.857294+(3.75846e-05*x);

      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return (0.857294+(3.75846e-05*x));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.688619+260.84/(x*x))*(1-(0.144982+0.000116685*x+-1.0021e-07*x*x));
      }
    }
  if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 50 ) return ((0.857294+(3.75846e-05*x))+0.028772119432687759);
        if (pt >= 50  && pt < 70 ) return (0.857294+(3.75846e-05*x))+0.027094315737485886;
        if (pt >= 70  && pt < 100) return (0.857294+(3.75846e-05*x))+0.029823761433362961;
        if (pt >= 100 && pt < 140) return (0.857294+(3.75846e-05*x))+0.038315325975418091;
        if (pt >= 140 && pt < 200) return (0.857294+(3.75846e-05*x))+0.043433453887701035;
        if (pt >= 200 && pt < 300) return (0.857294+(3.75846e-05*x))+0.093691818416118622;
        if (pt >= 300 && pt < 670) return (0.857294+(3.75846e-05*x))+0.084598012268543243;
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 50 ) return (0.857294+(3.75846e-05*x))+0.01438605971634388;
        if (pt >= 50  && pt < 70 ) return (0.857294+(3.75846e-05*x))+0.013547157868742943;
        if (pt >= 70  && pt < 100) return (0.857294+(3.75846e-05*x))+0.01491188071668148;
        if (pt >= 100 && pt < 140) return (0.857294+(3.75846e-05*x))+0.019157662987709045;
        if (pt >= 140 && pt < 200) return (0.857294+(3.75846e-05*x))+0.021716726943850517;
        if (pt >= 200 && pt < 300) return (0.857294+(3.75846e-05*x))+0.046845909208059311;
        if (pt >= 300 && pt < 670) return (0.857294+(3.75846e-05*x))+0.042299006134271622;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if( pt>= 20 && pt < 1000) return (0.688619+260.84/(x*x))*(1+(0.144982+0.000116685*x+-1.0021e-07*x*x));
      }
    }
if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
        if (pt >= 30  && pt < 50 ) return 0.857294+((3.75846e-05*x)-0.028772119432687759);
        if (pt >= 50  && pt < 70 ) return 0.857294+((3.75846e-05*x)-0.027094315737485886);
        if (pt >= 70  && pt < 100) return 0.857294+((3.75846e-05*x)-0.029823761433362961);
        if (pt >= 100 && pt < 140) return 0.857294+((3.75846e-05*x)-0.038315325975418091);
        if (pt >= 140 && pt < 200) return 0.857294+((3.75846e-05*x)-0.043433453887701035);
        if (pt >= 200 && pt < 300) return 0.857294+((3.75846e-05*x)-0.093691818416118622);
        if (pt >= 300 && pt < 670) return 0.857294+((3.75846e-05*x)-0.084598012268543243);
      }
      if(abs(flavor)==4){
        if (pt >= 30  && pt < 50 ) return 0.857294+((3.75846e-05*x)-0.01438605971634388);
        if (pt >= 50  && pt < 70 ) return 0.857294+((3.75846e-05*x)-0.013547157868742943);
        if (pt >= 70  && pt < 100) return 0.857294+((3.75846e-05*x)-0.01491188071668148);
        if (pt >= 100 && pt < 140) return 0.857294+((3.75846e-05*x)-0.019157662987709045);
        if (pt >= 140 && pt < 200) return 0.857294+((3.75846e-05*x)-0.021716726943850517);
        if (pt >= 200 && pt < 300) return 0.857294+((3.75846e-05*x)-0.046845909208059311);
        if (pt >= 300 && pt < 670) return 0.857294+((3.75846e-05*x)-0.042299006134271622);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt<1000) return(0.688619+260.84/(x*x))*(1-(0.144982+0.000116685*x+-1.0021e-07*x*x));
      }
    }
 }                 
  //Medium WP
  if(algo == "csvm"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.901114+(1.32145e-05*x);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.901114+(1.32145e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return 0.980777+-0.00109334*x+4.2909e-06*x*x+-2.78512e-09*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.901114+(1.32145e-05*x);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.901114+(1.32145e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.980777+-0.00109334*x+4.2909e-06*x*x+-2.78512e-09*x*x*x)*(1+(0.0672836+0.000102309*x+-1.01558e-07*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670)  return 0.901114+(1.32145e-05*x);

      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.901114+(1.32145e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.980777+-0.00109334*x+4.2909e-06*x*x+-2.78512e-09*x*x*x)*(1-(0.0672836+0.000102309*x+-1.01558e-07*x*x));
      }
    }

    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 50 ) return (0.901114+(1.32145e-05*x))+0.023258373141288757;
	if (pt >= 50  && pt < 70 ) return (0.901114+(1.32145e-05*x))+0.023003481328487396;
	if (pt >= 70  && pt < 100) return (0.901114+(1.32145e-05*x))+0.022424774244427681;
	if (pt >= 100 && pt < 140) return (0.901114+(1.32145e-05*x))+0.032237973064184189;
	if (pt >= 140 && pt < 200) return (0.901114+(1.32145e-05*x))+0.033661492168903351;
	if (pt >= 200 && pt < 300) return (0.901114+(1.32145e-05*x))+0.064984261989593506;
	if (pt >= 300 && pt < 670) return (0.901114+(1.32145e-05*x))+0.072893746197223663;
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 50 ) return (0.901114+(1.32145e-05*x))+0.011629186570644379;
	if (pt >= 50  && pt < 70 ) return (0.901114+(1.32145e-05*x))+0.011501740664243698;
	if (pt >= 70  && pt < 100) return (0.901114+(1.32145e-05*x))+0.01121238712221384;
	if (pt >= 100 && pt < 140) return (0.901114+(1.32145e-05*x))+0.0161189865320920949;
	if (pt >= 140 && pt < 200) return (0.901114+(1.32145e-05*x))+0.016830746084451675;
	if (pt >= 200 && pt < 300) return (0.901114+(1.32145e-05*x))+0.032492130994796753;
	if (pt >= 300 && pt < 670) return (0.901114+(1.32145e-05*x))+0.036446873098611832;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>=20 && pt<1000) return (0.980777+-0.00109334*x+4.2909e-06*x*x+-2.78512e-09*x*x*x)*(1+(0.0672836+0.000102309*x+-1.01558e-07*x*x));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 50 ) return 0.901114+((1.32145e-05*x)-0.023258373141288757);
	if (pt >= 50  && pt < 70 ) return 0.901114+((1.32145e-05*x)-0.023003481328487396);
	if (pt >= 70  && pt < 100) return 0.901114+((1.32145e-05*x)-0.022424774244427681);
	if (pt >= 100 && pt < 140) return 0.901114+((1.32145e-05*x)-0.032237973064184189);
	if (pt >= 140 && pt < 200) return 0.901114+((1.32145e-05*x)-0.033661492168903351);
	if (pt >= 200 && pt < 300) return 0.901114+((1.32145e-05*x)-0.064984261989593506);
	if (pt >= 300 && pt < 670) return 0.901114+((1.32145e-05*x)-0.072893746197223663);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 50 ) return 0.901114+((1.32145e-05*x)-0.011629186570644379);
	if (pt >= 50  && pt < 70 ) return 0.901114+((1.32145e-05*x)-0.011501740664243698);
	if (pt >= 70  && pt < 100) return 0.901114+((1.32145e-05*x)-0.01121238712221384);
	if (pt >= 100 && pt < 140) return 0.901114+((1.32145e-05*x)-0.016118986532092094);
	if (pt >= 140 && pt < 200) return 0.901114+((1.32145e-05*x)-0.016830746084451675);
	if (pt >= 200 && pt < 300) return 0.901114+((1.32145e-05*x)-0.032492130994796753);
	if (pt >= 300 && pt < 670) return 0.901114+((1.32145e-05*x)-0.036446873098611832);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (0.980777+-0.00109334*x+4.2909e-06*x*x+-2.78512e-09*x*x*x)*(1-(0.0672836+0.000102309*x+-1.01558e-07*x*x));
      }
    }
  }

  
  
  //Loose WP
  if(algo == "csvl"){
    double x=pt;
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x)*(1+(0.0539991+-6.29073e-06*x+-3.39895e-09*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.931535+(1.40704e-05*x);
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x)*(1-(0.0539991+-6.29073e-06*x+-3.39895e-09*x*x));
      }
    }
    

    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 50 ) return (0.931535+(1.40704e-05*x))+0.01710829883813858;
	if (pt >= 50  && pt < 70 ) return (0.931535+(1.40704e-05*x))+0.020177565515041351;
	if (pt >= 70  && pt < 100) return (0.931535+(1.40704e-05*x))+0.019350524991750717;
	if (pt >= 100 && pt < 140) return (0.931535+(1.40704e-05*x))+0.027258865535259247;
	if (pt >= 140 && pt < 200) return (0.931535+(1.40704e-05*x))+0.027310512959957123;
	if (pt >= 200 && pt < 300) return (0.931535+(1.40704e-05*x))+0.039027795195579529;
	if (pt >= 300 && pt < 670) return (0.931535+(1.40704e-05*x))+0.089093379676342122;
      }
      if(abs(flavor)==4){
      	if (pt >= 30  && pt < 50 ) return (0.931535+(1.40704e-05*x))+0.0085541494190692902;
	if (pt >= 50  && pt < 70 ) return (0.931535+(1.40704e-05*x))+0.010088782757520676;
	if (pt >= 70  && pt < 100) return (0.931535+(1.40704e-05*x))+0.0096752624958753586;
	if (pt >= 100 && pt < 140) return (0.931535+(1.40704e-05*x))+0.013629432767629623;
	if (pt >= 140 && pt < 200) return (0.931535+(1.40704e-05*x))+0.013655256479978561;
	if (pt >= 200 && pt < 300) return (0.931535+(1.40704e-05*x))+0.019513897597789764;
	if (pt >= 300 && pt < 670) return (0.931535+(1.40704e-05*x))+0.044546689838171005;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x)*(1+(0.0539991+-6.29073e-06*x+-3.39895e-09*x*x));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 50 ) return 0.931535+((1.40704e-05*x)-0.01710829883813858);
	if (pt >= 50  && pt < 70 ) return 0.931535+((1.40704e-05*x)-0.020177565515041351);
	if (pt >= 70  && pt < 100) return 0.931535+((1.40704e-05*x)-0.019350524991750717);
	if (pt >= 100 && pt < 140) return 0.931535+((1.40704e-05*x)-0.027258865535259247);
	if (pt >= 140 && pt < 200) return 0.931535+((1.40704e-05*x)-0.027310512959957123);
	if (pt >= 200 && pt < 300) return 0.931535+((1.40704e-05*x)-0.039027795195579529);
	if (pt >= 300 && pt < 670) return 0.931535+((1.40704e-05*x)-0.089093379676342122);
      }
      if(abs(flavor)==4){
      	if (pt >= 30  && pt < 50 ) return 0.931535+((1.40704e-05*x)-0.0085541494190692902);
	if (pt >= 50  && pt < 70 ) return 0.931535+((1.40704e-05*x)-0.010088782757520676);
	if (pt >= 70  && pt < 100) return 0.931535+((1.40704e-05*x)-0.0096752624958753586);
	if (pt >= 100 && pt < 140) return 0.931535+((1.40704e-05*x)-0.013629432767629623);
	if (pt >= 140 && pt < 200) return 0.931535+((1.40704e-05*x)-0.013655256479978561);
	if (pt >= 200 && pt < 300) return 0.931535+((1.40704e-05*x)-0.019513897597789764);
	if (pt >= 300 && pt < 670) return 0.931535+((1.40704e-05*x)-0.044546689838171005);
      }      
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x)*(1-(0.0539991+-6.29073e-06*x+-3.39895e-09*x*x));
      }
    }
  }


  return 1.0;

}

float DMAnalysisTreeMaker::BTagWeight::weightWithVeto(vector<JetInfo> jetsTags, int tags, vector<JetInfo> jetsVetoes, int vetoes)
{//This function takes into account cases where you have n b-tags and m vetoes, but they have different thresholds. 
    if (!filter(tags))
    {
        //   std::cout << "nThis event should not pass the selection, what is it doing here?" << std::endl;
        return 0;
    }
    int njets = jetsTags.size();
    if(njets != (int)(jetsVetoes.size()))return 0;//jets tags and vetoes must have same size!
    int comb = 1 << njets;
    float pMC = 0;
    float pData = 0;
    for (int i = 0; i < comb; i++)
    {
        float mc = 1.;
        float data = 1.;
        int ntagged = 0;
        for (int j = 0; j < njets; j++)
        {
            bool tagged = ((i >> j) & 0x1) == 1;
            if (tagged)
            {
                ntagged++;
                mc *= jetsTags[j].eff;
                data *= jetsTags[j].eff * jetsTags[j].sf;
            }
            else
            {
                mc *= (1. - jetsVetoes[j].eff);
                data *= (1. - jetsVetoes[j].eff * jetsVetoes[j].sf);
            }
        }

        if (filter(ntagged))
        {
            //  std::cout << mc << " " << data << endl;
            pMC += mc;
            pData += data;
        }
    }

    if (pMC == 0) return 0;
    return pData / pMC;
}

bool DMAnalysisTreeMaker::isEWKID(int id){
  bool isewk=false;
  int aid = abs(id);
  if((aid>10 && aid <17) || (aid==23 || aid==24)){isewk=true;}

  return isewk;
}

double DMAnalysisTreeMaker::pileUpSF(string syst)
{
  // if (syst == "PUUp" )return LumiWeightsUp_.weight( n0);
  // if (syst == "PUDown" )return LumiWeightsDown_.weight( n0);
  // return LumiWeights_.weight( n0);

  if (syst == "PUUp" )return LumiWeightsUp_.weight(n0);
  if (syst == "PUDown" )return LumiWeightsDown_.weight(n0);
  return LumiWeights_.weight( n0);

}

void DMAnalysisTreeMaker::endJob(){
  //  for(size_t s=0;s< systematics.size();++s){
  //    std::string syst  = systematics.at(s);
  /*  cout <<" init events are "<< nInitEvents <<endl;
  trees["EventHistory"]->SetBranchAddress("initialEvents",&nInitEvents);
  for(size_t i = 0; i < (size_t)nInitEvents;++i){
      //      trees["EventHistory"]->GetBranch("initialEvents")->Fill();
      
      //      trees["EventHistory"]->GetBranch("initialEvents")->Fill();
      trees["EventHistory"]->Fill();
      //      cout <<" i is "<< i << " entry is now "<< trees["EventHistory"]->GetBranch("initialEvents")->GetEntry()<<endl;
      
      }*/
  ;
}



//DMAnalysisTreeMaker::~DMAnalysisTreeMaker(const edm::ParameterSet& iConfig)
// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(DMAnalysisTreeMaker);

//  LocalWords:  firstidx
