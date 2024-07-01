#OutputBase = '/store/group/phys_heavyions/pchou/PbPb2023/Forest'
OutputBase = '/store/user/pchou/Zhadron2024'
DatasetName = '/pchou_DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/pchou-HINPbPbSpring23Digi-cc118d95320c1c5d965f7d5e145c7fa5/USER'
Tag = 'HINPbPbSpring23Reco'

from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = Tag

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HIN-HINPbPbSpring23Reco-00036_1_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5000
# config.JobType.maxJobRuntimeMin = 2750
# config.JobType.inputFiles = ['']

config.section_("Data")
config.Data.inputDataset = DatasetName
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.splitting = 'Automatic'
#config.Data.totalUnits = 4890
config.Data.unitsPerJob = 180
config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = Tag
config.Data.outLFNDirBase = OutputBase + '/'
#config.Data.outputPrimaryDataset = 'pchou_DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.blacklist = ['T2_IT_Pisa']
config.Site.whitelist = ['T2_US_*', 'T2_CH_CERN']
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_US_MIT'

