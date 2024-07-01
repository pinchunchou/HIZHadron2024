OutputBase = '/store/group/phys_heavyions/pchou/PbPb2023/Forest'
# DatasetName = '/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbSpring21MiniAOD-FixL1CaloGT_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM'
DatasetName = '/HIPhysicsRawPrime31/HIRun2023A-PromptReco-v2/MINIAOD'
Tag = '20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime31'

from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = Tag

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'forest_miniAOD_run3_DATA.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5000
config.JobType.maxJobRuntimeMin = 360
# config.JobType.inputFiles = ['']

config.section_("Data")
config.Data.inputDataset = DatasetName
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = Tag
config.Data.outLFNDirBase = OutputBase + DatasetName + '/'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.blacklist = ['T2_IT_Pisa']
config.Site.whitelist = ['T2_US_*', 'T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'

