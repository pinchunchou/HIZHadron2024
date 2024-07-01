
Common="--MinZPT 5 --MaxZPT 200 --MinTrackPT 0.5 --MaxTrackPT 1000 --DoSumET false"

MC="--DoGenLevel true --DoGenCorrelation false --IsData false"
MCGen="--DoGenLevel true --DoGenCorrelation true --GenCorrelationCharged true --IsData false"
Data="--DoGenLevel false --IsData true"

AODPF="--PFTree \"pfcandAnalyzer/pfTree\""
MiniAODPF="--PFTree \"particleFlowAnalyser/pftree\""

TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_0_20.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_20_60.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_60_100.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_100_200.root"
TrackResidualPathPP="${ProjectBase}/CommonCode/root/20230531_TrackResidualCorrection_V12_pp.root"

GenTrack="--DoTrackEfficiency false --TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ --DoTrackResidual false"
PPRecoTrack="--DoTrackEfficiency true --TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ --DoTrackResidual true --TrackResidualPath $TrackResidualPathPP"
PbPbRecoTrack="--DoTrackEfficiency true --TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ --DoTrackResidual true --TrackResidualPath $TrackResidualPathPbPb"

BackgroundMC="   --DoBackground true --HFShift 682  --Tolerance 90 --ToleranceFraction 0.005 --Oversample 10 --HFCeiling 134000"
BackgroundGenMC="--DoBackground true --HFShift 1083 --Tolerance 90 --ToleranceFraction 0.005 --Oversample 10 --HFCeiling 156000 --VZTolerance 10000"
BackgroundData=" --DoBackground true --HFShift 660  --Tolerance 75 --ToleranceFraction 0.005 --Oversample 25 --HFCeiling  70000"

BackgroundMCUEUp25="     --DoBackground true --HFShift 699   --Tolerance 90 --ToleranceFraction 0.05 --Oversample 10 --HFCeiling 134000"
BackgroundGenMCUEUp25="  --DoBackground true --HFShift 1110  --Tolerance 10 --ToleranceFraction 0.01 --Oversample 10 --HFCeiling 156000 --VZTolerance 10000"
BackgroundDataUEUp25="   --DoBackground true --HFShift 676.5 --Tolerance 0  --ToleranceFraction 0.01 --Oversample 10 --HFCeiling  70000"
BackgroundMCUEDown25="   --DoBackground true --HFShift 665   --Tolerance 20 --ToleranceFraction 0.05 --Oversample 10 --HFCeiling 134000"
BackgroundGenMCUEDown25="--DoBackground true --HFShift 1056  --Tolerance 10 --ToleranceFraction 0.01 --Oversample 10 --HFCeiling 156000 --VZTolerance 10000"
BackgroundDataUEDown25=" --DoBackground true --HFShift 643.5 --Tolerance 0  --ToleranceFraction 0.01 --Oversample 10 --HFCeiling  70000"

BackgroundMCUEUp30="     --DoBackground true --HFShift 702.5  --Tolerance 20 --ToleranceFraction 0.05 --Oversample 10 --HFCeiling 134000"
BackgroundGenMCUEUp30="  --DoBackground true --HFShift 1115.5 --Tolerance 10 --ToleranceFraction 0.01 --Oversample 10 --HFCeiling 156000 --VZTolerance 10000"
BackgroundDataUEUp30="   --DoBackground true --HFShift 680    --Tolerance 0  --ToleranceFraction 0.01 --Oversample 10 --HFCeiling  70000"
BackgroundMCUEDown30="   --DoBackground true --HFShift 661.5  --Tolerance 20 --ToleranceFraction 0.05 --Oversample 10 --HFCeiling 134000"
BackgroundGenMCUEDown30="--DoBackground true --HFShift 1050.5 --Tolerance 10 --ToleranceFraction 0.01 --Oversample 10 --HFCeiling 156000 --VZTolerance 10000"
BackgroundDataUEDown30=" --DoBackground true --HFShift 640.2  --Tolerance 0  --ToleranceFraction 0.01 --Oversample 10 --HFCeiling  70000"

DHSet Setting.dh PPSignalMC          Nominal  string "$Common $MC    --IsPP true  $AODPF     $PPRecoTrack"
DHSet Setting.dh PPSignalGenMC       Nominal  string "$Common $MCGen --IsPP true  $AODPF     $GenTrack"
DHSet Setting.dh PPSignalData        Nominal  string "$Common $Data  --IsPP true  $AODPF     $PPRecoTrack"
DHSet Setting.dh PbPbSignalMC        Nominal  string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack"
DHSet Setting.dh PbPbSignalGenMC     Nominal  string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack"
DHSet Setting.dh PbPbSignalData      Nominal  string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack"
DHSet Setting.dh PbPbBackgroundMC    Nominal  string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundMC"
DHSet Setting.dh PbPbBackgroundGenMC Nominal  string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack       $BackgroundGenMC"
DHSet Setting.dh PbPbBackgroundData  Nominal  string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundData"
DHSet Setting.dh PbPbBackgroundMC    UEUp25   string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundMCUEUp25"
DHSet Setting.dh PbPbBackgroundGenMC UEUp25   string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack       $BackgroundGenMCUEUp25"
DHSet Setting.dh PbPbBackgroundData  UEUp25   string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundDataUEUp25"
DHSet Setting.dh PbPbBackgroundMC    UEDown25 string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundMCUEDown25"
DHSet Setting.dh PbPbBackgroundGenMC UEDown25 string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack       $BackgroundGenMCUEDown25"
DHSet Setting.dh PbPbBackgroundData  UEDown25 string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundDataUEDown25"
DHSet Setting.dh PbPbBackgroundMC    UEUp30   string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundMCUEUp30"
DHSet Setting.dh PbPbBackgroundGenMC UEUp30   string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack       $BackgroundGenMCUEUp30"
DHSet Setting.dh PbPbBackgroundData  UEUp30   string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundDataUEUp30"
DHSet Setting.dh PbPbBackgroundMC    UEDown30 string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundMCUEDown30"
DHSet Setting.dh PbPbBackgroundGenMC UEDown30 string "$Common $MCGen --IsPP false $MiniAODPF $GenTrack       $BackgroundGenMCUEDown30"
DHSet Setting.dh PbPbBackgroundData  UEDown30 string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack  $BackgroundDataUEDown30"

LooseTrack="--DoAlternateTrackSelection true --AlternateTrackSelection 1"
DHSet Setting.dh PbPbSignalData      Loose string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack $LooseTrack"
DHSet Setting.dh PbPbSignalMC        Loose string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack $LooseTrack"

TightTrack="--DoAlternateTrackSelection true --AlternateTrackSelection 2"
DHSet Setting.dh PbPbSignalData      Tight string "$Common $Data  --IsPP false $MiniAODPF $PbPbRecoTrack $TightTrack"
DHSet Setting.dh PbPbSignalMC        Tight string "$Common $MC    --IsPP false $MiniAODPF $PbPbRecoTrack $TightTrack"

# Finally set sample locations

DHSet Setting.dh Sample PbPbSignalMCRun2 string /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20240605_ZHadronMCDY/240605_151412/0000/
DHSet Setting.dh Sample PbPbSignalDataRun2 string /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIHardProbes/HIRun2018A-PbPb18_MiniAODv1-v2/MINIAOD/HIHardProbes/20230531_ZHadronSingleElectronWithMuTree/240531_190157/
DHSet Setting.dh Sample PbPbBackgroundMCRun2 string /eos/cms/store/group/phys_heavyions/chenyi/Samples/store/user/mitaylor/PhotonJet/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20230127PbPbMCMB/230127_211118/
DHSet Setting.dh Sample PbPbBackgroundDataRun2 string /eos/cms/store/group/phys_heavyions/chenyi/Samples/store/user/mitaylor/PhotonJet/HIMinimumBias4/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIMinimumBias4/20230130PbPbMB/230130_223511/
DHSet Setting.dh Sample PPSignalMCRun2 string /eos/cms/store/group/phys_heavyions/chenyi/pp2017/Forest/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2/AODSIM/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/20231020_ZHadronMLLWithMuTreePP/231020_152224/0000/
DHSet Setting.dh Sample PPSignalDataRun2 string /eos/cms/store/group/phys_heavyions/chenyi/pp2017/Forest/HighEGJet/Run2017G-17Nov2017-v2/AOD/HighEGJet/20230403_ZHadronSingleElectronWithMuTree/230403_180009/
DHSet Setting.dh Sample PbPbSignalMCRun3 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/pchou_DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/pchou-HINPbPbSpring23MiniAOD-b75d3b5331198c44808663383a6aa555/USER/pchou_DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/20240613_ZHadronMCDY/240613_195634/0000/
DHSet Setting.dh Sample PbPbSignalDataRawPrime0  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime0/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime0/240607_222243/
DHSet Setting.dh Sample PbPbSignalDataRawPrime10 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime10/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime10/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime10/240608_031138/
DHSet Setting.dh Sample PbPbSignalDataRawPrime11 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime11/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime11/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime11/240608_031232/
DHSet Setting.dh Sample PbPbSignalDataRawPrime12 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime12/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime12/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime12/240608_031325/
DHSet Setting.dh Sample PbPbSignalDataRawPrime13 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime13/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime13/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime13/240608_031340/
DHSet Setting.dh Sample PbPbSignalDataRawPrime14 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime14/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime14/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime14/240608_031358/
DHSet Setting.dh Sample PbPbSignalDataRawPrime15 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime15/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime15/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime15/240608_031410/
DHSet Setting.dh Sample PbPbSignalDataRawPrime16 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime16/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime16/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime16/240608_031431/
DHSet Setting.dh Sample PbPbSignalDataRawPrime17 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime17/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime17/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime17/240608_031447/
DHSet Setting.dh Sample PbPbSignalDataRawPrime18 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime18/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime18/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime18/240608_031508/
DHSet Setting.dh Sample PbPbSignalDataRawPrime19 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime19/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime19/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime19/240608_031550/
DHSet Setting.dh Sample PbPbSignalDataRawPrime1  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime1/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime1/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime1/240608_024339/
DHSet Setting.dh Sample PbPbSignalDataRawPrime20 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime20/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime20/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime20/240608_031605/
DHSet Setting.dh Sample PbPbSignalDataRawPrime21 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime21/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime21/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime21/240608_031620/
DHSet Setting.dh Sample PbPbSignalDataRawPrime22 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime22/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime22/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime22/240608_031635/
DHSet Setting.dh Sample PbPbSignalDataRawPrime23 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime23/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime23/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime23/240608_031652/
DHSet Setting.dh Sample PbPbSignalDataRawPrime24 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime24/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime24/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime24/240608_031707/
DHSet Setting.dh Sample PbPbSignalDataRawPrime25 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime25/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime25/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime25/240608_031727/
DHSet Setting.dh Sample PbPbSignalDataRawPrime26 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime26/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime26/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime26/240608_031741/
DHSet Setting.dh Sample PbPbSignalDataRawPrime27 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime27/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime27/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime27/240608_031749/
DHSet Setting.dh Sample PbPbSignalDataRawPrime28 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime28/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime28/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime28/240608_031758/
DHSet Setting.dh Sample PbPbSignalDataRawPrime29 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime29/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime29/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime29/240608_031807/
DHSet Setting.dh Sample PbPbSignalDataRawPrime2  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime2/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime2/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime2/240608_030624/
DHSet Setting.dh Sample PbPbSignalDataRawPrime30 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime30/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime30/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime30/240608_031818/
DHSet Setting.dh Sample PbPbSignalDataRawPrime31 string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime31/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime31/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime31/240608_031827/
DHSet Setting.dh Sample PbPbSignalDataRawPrime3  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime3/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime3/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime3/240608_030655/
DHSet Setting.dh Sample PbPbSignalDataRawPrime4  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime4/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime4/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime4/240608_030717/
DHSet Setting.dh Sample PbPbSignalDataRawPrime5  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime5/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime5/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime5/240608_030746/
DHSet Setting.dh Sample PbPbSignalDataRawPrime6  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime6/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime6/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime6/240608_031036/
DHSet Setting.dh Sample PbPbSignalDataRawPrime7  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime7/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime7/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime7/240608_031050/
DHSet Setting.dh Sample PbPbSignalDataRawPrime8  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime8/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime8/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime8/240608_031103/
DHSet Setting.dh Sample PbPbSignalDataRawPrime9  string /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime9/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime9/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime9/240608_031117/

