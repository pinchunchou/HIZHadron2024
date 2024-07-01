#OutputBase = '/store/group/phys_heavyions/pchou/PbPb2023/Forest'
OutputBase = '/store/user/pchou/Zhadron2024'
# DatasetName = '/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbSpring21MiniAOD-FixL1CaloGT_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM'
Tag = 'HINPbPbSpring23wmLHEGSHIMix'

from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = Tag

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'HIN-HINPbPbSpring23wmLHEGSHIMix-00018_1_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5000
config.JobType.maxJobRuntimeMin = 2750
# config.JobType.inputFiles = ['']

config.section_("Data")
#config.Data.inputDataset = DatasetName
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.splitting = 'EventBased'
config.Data.totalUnits = 250000
config.Data.unitsPerJob = 50
config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = Tag
config.Data.outLFNDirBase = OutputBase + '/'
config.Data.outputPrimaryDataset = 'pchou_DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.blacklist = ['T2_IT_Pisa']
config.Site.whitelist = ['T2_US_*', 'T2_CH_CERN']
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_US_MIT'

