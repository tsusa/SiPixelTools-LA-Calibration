import os
import glob

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = '<requestName>'

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '<name of the config file to run>'
#only need the following if output is non-EDM, like the LATrees
config.JobType.outputFiles = ['<file name>.root']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 6000
config.JobType.numCores = 8

config.section_('Data')
# for inputDBS use phys03 if running on privately produced GEN-SIM, otherwise comment out
config.Data.inputDBS = 'phys03' 
config.Data.inputDataset = '<input dataset name>'
config.Data.outLFNDirBase = '/store/user/<username>/<directory name>/'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3
config.Data.publication = True
config.Data.ignoreLocality = True

config.section_('Site')
#config.Site.storageSite = 'T2_HU_Budapest'
config.Site.storageSite = 'T3_US_FNALLPC'
#whitelist for the ignoreLocality option
config.Site.whitelist = ['T2_DE_DESY','T2_FR_IPHC','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*']