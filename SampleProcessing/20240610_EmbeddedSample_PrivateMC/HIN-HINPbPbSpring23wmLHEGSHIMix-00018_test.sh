#!/bin/bash

export SCRAM_ARCH=el8_amd64_gcc11

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_13_0_18_HeavyIon/src ] ; then
  echo release CMSSW_13_0_18_HeavyIon already exists
else
  scram p CMSSW CMSSW_13_0_18_HeavyIon
fi
cd CMSSW_13_0_18_HeavyIon/src
eval `scram runtime -sh`

mv ../../Configuration .
scram b
cd ../..

# Maximum validation duration: 28800s
# Margin for validation duration: 30%
# Validation duration with margin: 28800 * (1 - 0.30) = 20160s
# Time per event for each sequence: 4.1225s
# Threads for each sequence: 1
# Time per event for single thread for each sequence: 1 * 4.1225s = 4.1225s
# Which adds up to 4.1225s per event
# Single core events that fit in validation duration: 20160s / 4.1225s = 4890
# Produced events limit in McM is 10000
# According to 0.4000 efficiency, validation should run 10000 / 0.4000 = 25000 events to reach the limit of 10000
# Take the minimum of 4890 and 25000, but more than 0 -> 4890
# It is estimated that this validation will produce: 4890 * 0.4000 = 1956 events
EVENTS=4890

# Random seed between 1 and 100 for externalLHEProducer
SEED=$(($(date +%s) % 100 + 1))


# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/HIN-HINPbPbSpring23wmLHEGSHIMix-00018-fragment.py --python_filename HIN-HINPbPbSpring23wmLHEGSHIMix-00018_1_cfg.py --eventcontent RAWSIM,LHE --pileup HiMixGEN --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM,LHE --fileout file:HIN-HINPbPbSpring23wmLHEGSHIMix-00018.root --pileup_input "dbs:/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23GS-130X_mcRun3_2023_realistic_HI_v18-v2/GEN-SIM" --conditions 130X_mcRun3_2023_realistic_HI_v18 --beamspot MatchHI --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(${SEED})"\\nprocess.source.numberEventsInLuminosityBlock="cms.untracked.uint32(250)" --step LHE,GEN,SIM --scenario HeavyIons --geometry DB:Extended --era Run3_pp_on_PbPb --no_exec --mc -n $EVENTS || exit $? ;

# End of HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh file
