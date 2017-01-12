import commands, os
#1;95;0c## *****************************************************************************************
### Usage:
###
### cmsRun topplusdmanaEDMntuples_cfg.py maxEvts=N sample="mySample/sample.root" version="1"7 outputLabel="myoutput"
###
### Default values for the options are set:
### maxEvts     = -1
### sample      = 'file:/scratch/decosa/ttDM/testSample/tlbsm_53x_v3_mc_10_1_qPV.root'
### outputLabel = 'analysisTTDM.root'
### *****************************************************************************************
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('maxEvts',
                 1000,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

SingleElTriggers = []
SingleMuTriggers = []
hadronTriggers = []

chan = "MET_Prompt"

chan = "TTbarDMJets_scalar_Mchi-50_Mphi-50"

chan = "DY"
chan = "WJ"
filedir= "/tmp/oiorio/"
cmd = "ls "+filedir+"/"+chan+"/"

status,ls_la = commands.getstatusoutput( cmd )
listFiles = ls_la.split(os.linesep)
files = []
#files = ["file:re-MiniAOD_17Jul/"+l for l in listFiles]
#files = ["file:"+filedir+"MET_Prompt/"+l for l in listFiles]
files = ["file:"+filedir+"/"+chan+"/"+l for l in listFiles]
options.register('sample',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p3/JetHT/Run2016B/JetHT/Run2016B-23Sep2016-v3_B2GAnaFW_80X_V2p3/161216_214635/0000/B2GEDMNtuple_1.root',
                 #'/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170104_182059/0000/B2GEDMNtuple_1.root',
                 '/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016B-23Sep2016-v3_B2GAnaFW_80X_v2p4/161221_152050/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p1/SingleMuon/Run2016B/SingleMuon/Run2016B-PromptReco-v2_B2GAnaFW_80X_V2p1/161027_151657/0000/B2GEDMNtuple_1.root/',
                 #                 files,
                 #'/SingleMuon/vorobiev-Run2016H-PromptReco-v1_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER',
                 #"/store/user/grauco/B2GAnaFW_80X_V2p4/BprimeBToHB_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/BprimeBToHB_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161219_152207/0000/B2GEDMNtuple_1.root",
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161222_110143/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016H-PromptReco-v1_B2GAnaFW_80X_v2p4/161221_180445/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p1/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1_B2GAnaFW_80X_V2p1/161018_195320/0000/B2GEDMNtuple_1.root',
                 #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p1/161018_070104/0000/B2GEDMNtuple_1.root',
                 #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p1/161018_070023/0000/B2GEDMNtuple_4.root',
                 #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_111914/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p1/JetHT/Run2016B/JetHT/Run2016B-PromptReco-v2_B2GAnaFW_80X_V2p1/161013_132109/0000/B2GEDMNtuple_1.root',
                 #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_80X_V2p1/161021_085128/0000/B2GEDMNtuple_1.root',
                 #'/store/user/ggiannin/B2GAnaFW_80X_V2p1/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISpring16MiniAODv2/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1/161022_101446/0000/B2GEDMNtuple_meow_1.root',
                 #'file:B2GEDMNtuple.root',
                 #'root://xrootd.ba.infn.it//store/user/grauco/NtuplesJETMET_Fwk80v1/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/b2ganafw80xv1_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_27Jul/160727_091755/0000/B2GEDMNtuple_1.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')

options.register('version',
                 #'53',
                 '71',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'ntuple version (53 or 71)')

options.register('outputLabel',
                 'treesTTJets.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('isData',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Is data?')

options.register('applyRes',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'ApplyResiduals?')

options.register('addPartonInfo',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Add parton info??')

options.register('changeJECs',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Apply new JECs?')

options.register('LHE',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Keep LHEProducts')

options.register('lhesource',
                 'externalLHEProducer',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'LHEProducts source')

options.register('channel',
                 '',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'channel for weight evaluation'
                 )

options.parseArguments()

if(options.isData):options.LHE = False

#options.LHE = False 
if(not options.isData): options.applyRes = False

l = ["singleTrigger"+str(s) for s in xrange(15)]
l = l + ["trigger2"+str(s) for s in xrange(15)]

SingleElTriggers = ["HLT_Ele27_eta2p1_WP75_Gsf"]
SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WP75_Gsf_v"+str(s) for s in xrange(15)]

PhotonTriggers = ["HLT_Photon155"]
PhotonTriggers = PhotonTriggers + ["HLT_Photon155_v"+str(s) for s in xrange(15)]

SingleMuTriggers = ["HLT_IsoMu20","HLT_IsoTkMu20"]
SingleMuTriggers = SingleMuTriggers + ["HLT_IsoMu20_v"+str(s) for s in xrange(15)]
SingleMuTriggers = SingleMuTriggers + ["HLT_IsoTkMu20_v"+str(s) for s in xrange(15)]

hadronTriggers = ["HLT_PFHT200_v2", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v7","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7", "HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v3", "HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v2", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v2", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v2", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v2", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v2",  "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight", "HLT_PFHT300_v3", "HLT_PFHT350_v4", "HLT_PFHT400_v3", "HLT_PFHT475_v3", "HLT_PFHT600_v4", "HLT_PFHT650_v4", "HLT_PFHT800_v3", "HLT_PFHT900_v2", "HLT_AK8PFJet360_TrimMass30_v4", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v4", "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v3", "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v2", "HLT_AK8DiPFJet280_200_TrimMass30_v2", "HLT_AK8DiPFJet250_200_TrimMass30_v2", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v2", "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v2", "HLT_PFHT300_v2", "HLT_PFHT350_v3", "HLT_PFHT400_v2", "HLT_PFHT475_v2", "HLT_PFHT600_v3", "HLT_PFHT650_v3", "HLT_PFHT800_v2", "HLT_AK8PFJet360_TrimMass30_v3", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v3", "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v2", "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v1", "HLT_AK8DiPFJet280_200_TrimMass30_v1", "HLT_AK8DiPFJet250_200_TrimMass30_v1", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v1", "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v1", "HLT_AK8PFJet200_v4", "HLT_AK8PFJet260_v5", "HLT_AK8PFJet320_v5", "HLT_AK8PFJet400_v5", "HLT_AK8PFJet450_v5", "HLT_AK8PFJet500_v5", "HLT_PFJet260_v9", "HLT_PFJet320_v9",  "HLT_PFJet400_v9",  "HLT_PFJet450_v9",  "HLT_PFJet500_v9", "HLT_PFHT400_v7", "HLT_PFHT475_v7", "HLT_PFHT600_v8", "HLT_PFHT650_v8", "HLT_PFHT800_v7", "HLT_PFHT900_v6", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v7", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v7", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7",  "HLT_AK8PFJet200_v1", "HLT_AK8PFJet260_v1", "HLT_AK8PFJet320_v1", "HLT_AK8PFJet400_v1", "HLT_AK8PFJet450_v1", "HLT_AK8PFJet500_v1", "HLT_PFJet200_v5", "HLT_PFJet260_v5", "HLT_PFJet320_v5", "HLT_PFJet400_v5", "HLT_PFJet450_v5", "HLT_PFJet500_v5", "HLT_PFHT300_v3", "HLT_PFHT350_v4", "HLT_PFHT400_v3", "HLT_PFHT475_v3", "HLT_PFHT600_v4", "HLT_PFHT650_v4", "HLT_PFHT800_v3", "HLT_PFHT900_v2", "HLT_PFHT600_v3", "HLT_PFHT650_v3", "HLT_PFHT800_v2", "HLT_PFHT900_v1", "HLT_PFJet320_v4",  "HLT_PFJet400_v4",  "HLT_PFJet450_v4",  "HLT_PFJet500_v4" ]

hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]

if(options.isData):

    SingleElTriggers = ["HLT_Ele27_eta2p1_WPLoose_Gsf", "HLT_Ele23_WPLoose_Gsf"]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele23_WPLoose_Gsf_v"+str(s) for s in xrange(15)]

    SingleMuTriggers = ["HLT_IsoMu20","HLT_IsoTkMu20"]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoMu20_v"+str(s) for s in xrange(15)]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoTkMu20_v"+str(s) for s in xrange(15)]

    hadronTriggers = ["HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v1", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v1", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v1", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v1","HLT_PFHT300_v2", "HLT_PFHT350_v3", "HLT_PFHT400_v2", "HLT_PFHT475_v2", "HLT_PFHT600_v3", "HLT_PFHT650_v3", "HLT_PFHT800_v2", "HLT_PFHT900_v1", "HLT_AK8PFJet360_TrimMass30_v3", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v3", "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v2", "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v1", "HLT_AK8DiPFJet280_200_TrimMass30_v2", "HLT_AK8DiPFJet250_200_TrimMass30_v2", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v2", "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v2", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_NoID", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_Mu17_Mu8_v1", "HLT_ZeroBias_v2", "HLT_Photon22_v2 ","HLT_Photon30_v3 ","HLT_Photon36_v3 ","HLT_Photon50_v3 ","HLT_Photon75_v3 ","HLT_Photon90_v3 ","HLT_Photon120_v3", "HLT_Photon165_HE10_v3", "HLT_AK8PFJet200_v4", "HLT_AK8PFJet260_v5", "HLT_AK8PFJet320_v5", "HLT_AK8PFJet400_v5", "HLT_AK8PFJet450_v5", "HLT_AK8PFJet500_v5", "HLT_AK8PFJet260_v9", "HLT_AK8PFJet320_v9",  "HLT_AK8PFJet400_v9",  "HLT_AK8PFJet450_v9",  "HLT_AK8PFJet500_v9", "HLT_PFHT400_v7", "HLT_PFHT475_v7", "HLT_PFHT600_v8", "HLT_PFHT650_v8", "HLT_PFHT800_v7", "HLT_PFHT900_v6", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v7", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v7", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v2", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v2", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v2", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v2", "HLT_AK8PFJet200_v1", "HLT_AK8PFJet260_v1", "HLT_AK8PFJet320_v1", "HLT_AK8PFJet400_v1", "HLT_AK8PFJet450_v1", "HLT_AK8PFJet500_v1", "HLT_PFJet200_v5", "HLT_PFJet260_v5", "HLT_PFJet320_v5", "HLT_PFJet400_v5", "HLT_PFJet450_v5", "HLT_PFJet500_v5", "HLT_PFHT300_v3", "HLT_PFHT350_v4", "HLT_PFHT400_v3", "HLT_PFHT475_v3", "HLT_PFHT600_v4", "HLT_PFHT650_v4", "HLT_PFHT800_v3", "HLT_PFHT900_v2", "HLT_PFHT600_v3", "HLT_PFHT650_v3", "HLT_PFHT800_v2", "HLT_PFHT900_v1", "HLT_PFJet320_v4",  "HLT_PFJet400_v4",  "HLT_PFJet450_v4",  "HLT_PFJet500_v4"]

    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_NoID_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]

process = cms.Process("ttDManalysisTrees")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )
### Source file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        options.sample
        )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag

process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2')  
#process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

### Rootplizer

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputLabel))

process.load("ttDM.treeDumper.topplusdmedmRootTreeMaker_cff")
#process.load("ttDM.treeDumper.topplusdmedmRootTreeMaker_with_cat_cff")
#process.load("B2GAnaFW.B2GAnaFW.topplusdmedmRootTreeMaker_cff")

process.DMTreesDumper.channelInfo.SingleElTriggers=cms.vstring(SingleElTriggers)
process.DMTreesDumper.channelInfo.SingleMuTriggers=cms.vstring(SingleMuTriggers)
process.DMTreesDumper.channelInfo.hadronicTriggers=cms.vstring(hadronTriggers)

if options.addPartonInfo:
    #if options.isData: #G
#        process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
#        process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
#        process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
#        process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)

    if not options.isData: #G                                                                                                                            
        process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)
#process.DMTreesDumper.lhes =cms.InputTag("externalLHEProducer")
process.DMTreesDumper.lhes =cms.InputTag(options.lhesource)
#process.DMTreesDumper.lhes =cms.InputTag(options.lhesource)
process.DMTreesDumper.changeJECs = cms.untracked.bool(options.changeJECs)
process.DMTreesDumper.isData = cms.untracked.bool(options.isData)#This adds the L2L3Residuals
process.DMTreesDumper.applyRes = cms.untracked.bool(options.applyRes)#This adds the L2L3Residuals
process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(True)
process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(True)

if options.channel == "ttbar":
    process.DMTreesDumper.getPartonTop  = cms.untracked.bool(True)
    #process.DMTreesDumper.channelInfo.getParticleWZ  = cms.untracked.bool(True) 
if options.channel == "wzjets":
    print "channel is " + options.channel 
    process.DMTreesDumper.channelInfo.getPartonW  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getParticleWZ  = cms.untracked.bool(True)

if options.isData:
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(False)
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)

if not options.isData:
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(True) #G
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(True)

#17Jul if(options.isData): del process.DMTreesDumper.physicsObjects[10000]
process.analysisPath = cms.Path(
    process.DMTreesDumper
    )


#if(options.isData):
for p in process.DMTreesDumper.physicsObjects:
    if(p.prefix == cms.string("el")):
        p.variablesF.append(cms.InputTag("electrons","elvidVeto"))
        p.variablesF.append(cms.InputTag("electrons","elvidLoose"))
        p.variablesF.append(cms.InputTag("electrons","elvidMedium"))
        p.variablesF.append(cms.InputTag("electrons","elvidTight"))
            #print "yes"
        p.toSave.append("elvidVeto")
        p.toSave.append( "elvidLoose")
        p.toSave.append("elvidMedium")
        p.toSave.append("elvidTight")
            
