#!/bin/bash

# GEN Script begin
rm -f request_fragment_check.py
wget -q https://raw.githubusercontent.com/cms-sw/genproductions/master/bin/utils/request_fragment_check.py
chmod +x request_fragment_check.py
./request_fragment_check.py --bypass_status --prepid HIN-HINPbPbSpring23wmLHEGSHIMix-00018
GEN_ERR=$?
if [ $GEN_ERR -ne 0 ]; then
  echo "GEN Checking Script returned exit code $GEN_ERR which means there are $GEN_ERR errors"
  echo "Validation WILL NOT RUN"
  echo "Please correct errors in the request and run validation again"
  exit $GEN_ERR
fi
echo "Running VALIDATION. GEN Request Checking Script returned no errors"
# GEN Script end

# Download fragment from McM
curl -s -k https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/HIN-HINPbPbSpring23wmLHEGSHIMix-00018 --retry 3 --create-dirs -o Configuration/GenProduction/python/HIN-HINPbPbSpring23wmLHEGSHIMix-00018-fragment.py
[ -s Configuration/GenProduction/python/HIN-HINPbPbSpring23wmLHEGSHIMix-00018-fragment.py ] || exit $?;

# Check if fragment contais gridpack path ant that it is in cvmfs
if grep -q "gridpacks" Configuration/GenProduction/python/HIN-HINPbPbSpring23wmLHEGSHIMix-00018-fragment.py; then
  if ! grep -q "/cvmfs/cms.cern.ch/phys_generator/gridpacks" Configuration/GenProduction/python/HIN-HINPbPbSpring23wmLHEGSHIMix-00018-fragment.py; then
    echo "Gridpack inside fragment is not in cvmfs."
    exit -1
  fi
fi

# Dump actual test code to a HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh file that can be run in Singularity
cat <<'EndOfTestFile' > HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh
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

# Run generated config
REPORT_NAME=HIN-HINPbPbSpring23wmLHEGSHIMix-00018_report.xml
# Run the cmsRun
cmsRun -e -j $REPORT_NAME HIN-HINPbPbSpring23wmLHEGSHIMix-00018_1_cfg.py || exit $? ;

# Parse values from HIN-HINPbPbSpring23wmLHEGSHIMix-00018_report.xml report
processedEvents=$(grep -Po "(?<=<Metric Name=\"NumberEvents\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
producedEvents=$(grep -Po "(?<=<TotalEvents>)(\d*)(?=</TotalEvents>)" $REPORT_NAME | tail -n 1)
threads=$(grep -Po "(?<=<Metric Name=\"NumberOfThreads\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
peakValueRss=$(grep -Po "(?<=<Metric Name=\"PeakValueRss\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
peakValueVsize=$(grep -Po "(?<=<Metric Name=\"PeakValueVsize\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalSize=$(grep -Po "(?<=<Metric Name=\"Timing-tstoragefile-write-totalMegabytes\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalSizeAlt=$(grep -Po "(?<=<Metric Name=\"Timing-file-write-totalMegabytes\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalJobTime=$(grep -Po "(?<=<Metric Name=\"TotalJobTime\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalJobCPU=$(grep -Po "(?<=<Metric Name=\"TotalJobCPU\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
eventThroughput=$(grep -Po "(?<=<Metric Name=\"EventThroughput\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
avgEventTime=$(grep -Po "(?<=<Metric Name=\"AvgEventTime\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
if [ -z "$threads" ]; then
  echo "Could not find NumberOfThreads in report, defaulting to 1"
  threads=1
fi
if [ -z "$eventThroughput" ]; then
  eventThroughput=$(bc -l <<< "scale=4; 1 / ($avgEventTime / $threads)")
fi
if [ -z "$totalSize" ]; then
  totalSize=$totalSizeAlt
fi
if [ -z "$processedEvents" ]; then
  processedEvents=$EVENTS
fi
echo "Validation report of HIN-HINPbPbSpring23wmLHEGSHIMix-00018 sequence 1/1"
echo "Processed events: $processedEvents"
echo "Produced events: $producedEvents"
echo "Threads: $threads"
echo "Peak value RSS: $peakValueRss MB"
echo "Peak value Vsize: $peakValueVsize MB"
echo "Total size: $totalSize MB"
echo "Total job time: $totalJobTime s"
echo "Total CPU time: $totalJobCPU s"
echo "Event throughput: $eventThroughput"
echo "CPU efficiency: "$(bc -l <<< "scale=2; ($totalJobCPU * 100) / ($threads * $totalJobTime)")" %"
echo "Size per event: "$(bc -l <<< "scale=4; ($totalSize * 1024 / $producedEvents)")" kB"
echo "Time per event: "$(bc -l <<< "scale=4; (1 / $eventThroughput)")" s"
echo "Filter efficiency percent: "$(bc -l <<< "scale=8; ($producedEvents * 100) / $processedEvents")" %"
echo "Filter efficiency fraction: "$(bc -l <<< "scale=10; ($producedEvents) / $processedEvents")

# End of HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh file
EndOfTestFile

# Make file executable
chmod +x HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh

if [ -e "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:amd64" ]; then
  CONTAINER_NAME="el8:amd64"
elif [ -e "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64" ]; then
  CONTAINER_NAME="el8:x86_64"
else
  echo "Could not find amd64 or x86_64 for el8"
  exit 1
fi
# Run in singularity container
# Mount afs, eos, cvmfs
# Mount /etc/grid-security for xrootd
export SINGULARITY_CACHEDIR="/tmp/$(whoami)/singularity"
singularity run -B /afs -B /eos -B /cvmfs -B /etc/grid-security -B /etc/pki/ca-trust --home $PWD:$PWD /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/$CONTAINER_NAME $(echo $(pwd)/HIN-HINPbPbSpring23wmLHEGSHIMix-00018_test.sh)