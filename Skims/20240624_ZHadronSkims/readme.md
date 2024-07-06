
# Skimming code from forest to plain root file

## Introduction

## Update Log

* Version 1 - 2024 June 28
   * Simplified from https://github.com/pinchunchou/PhysicsHIZHadron2022 
   * Enable it to skim both dimuon and dielectron candidates as well as run 2 / run 3 samples.
   * Select the best Z mass candidate per Kaya, rather than the first Z candidate.
   * New weights not implemented yet. Will be in newer versions.
   * Background matching not implemented yet. Will be in newer versions.

* Version 1a - 2024 July 2
   * Implement basic background matching. Event selections and ForceGenMatch not yet ready for background.
   * New weights not implemented yet. Will be in newer versions.
   * Plan to include rapidity (Y) in newer version.

* Version 1b - 2024 July 5
   * Fix selections: remove Z rapidity selections, fix electron SC selections.

* Version 1c - 2024 July 5
   * Fix electron selections and save rapidity information (in which we assume track mass = pion mass).

## Arguments

| Argument | type | default | comments |
| -------- | ---- | ------- | -------- |
| Input | vector<string> | _required_ | List of input files separated by comma |
| Output | string | _required_ | Name of the output file |
| DoMuon | bool | true | If we skim dimuon candidates or not |
| DoElectron | bool | true | If we skim dielectron candidates or not |
| DoGenLevel | bool | true | If we store gen muon/Z information or not |
| Fraction | double | 1.00 | What fraction of root file to use |
| MinZPT | double | 20 | Minimum Z PT to store correlated tracks.  If Z PT is below this tracks are not stored |
| MaxZPT | double | 2000000 | Maximum Z PT to store correlated tracks.  If Z PT is above this tracks are not stored |
| MinTrackPT | double | 1.00 | Minimum track PT to store |
| MaxTrackPT | double | 10000.00 | Maximum track PT to store |
| MinGenTrackPT | double | 0.40 | Minimum Gen track PT to store and to be used in GetGenHFSum calculation |
| MinPFPT | double | 0 | Minimum PF PT to be used in GetHFSum calculation |
| IsData | bool | false | whether this is data or not |
| IsPP | bool | false | whether this is pp or not |
| DoGenCorrelation | bool | false | whether to use gen particles as the "tracks" or not |
| GenCorrelationCharged | bool | false | whether to use all gen particles or only charged |
| DoBackground | bool | false | whether to do background correlation or not |
| Background | vector<string> | _required if DoBackground = true_ | list of file names for background mixing |
| HFShift | double | _required if DoBackground is true_ | Amount of shift to subtract from the HF number |
| HFTolerance | double | _required if DoBackground is true_ | Tolerance to HF sum in GeV |
| HFToleranceFraction | double | _required if DoBackground is true_ | Tolerance to HF sum in fractions |
| HFCeiling | double | -1 | If lower bound of HF matching is larger than this, it is lowered to this value.  Prevents ultra-central matching failures.  Set to negative to disable this. |
| VZTolerance | double | 2 | Tolerance to VZ |
| Oversample | int | 1 | How many times we mix every signal event |
| DoSumET | bool | true | Whether we use SumET or SumE in HF for event matching |
| MuonVeto | double | 0.01 | window for track-muon rejection |
| DoAlternateTrackSelection | false | Whether to use alternate track selection or not |
| AlternateTrackSelection | _required if DoAlternateTrackSelection is true_ | 0 is default, 1 is loose, 2 is tight |
| DoTrackEfficiency | bool | true | If we want to store track efficiency correction factor |
| TrackEfficiencyPath | string | _required if DoTrackEfficiency is true_ | Base path for track correction files |
| DoTrackResidual | bool | true | If we want to store track residual efficiency correction factor |
| TrackResidualPath | string | _required if DoTrackResidual is true_ | Path for track residual correction files |
| PFTreeName | string | "pfcandAnalyzer/pfTree" (IsPP true) or "particleFlowAnalyser/pftree" (IsPP false) | Name of the particle flow tree.  Default value depends on whether it is pp mode or not |
| GGTreeName | string | "ggHiNtuplizerGED/EventTree" (IsPP true) or "ggHiNtuplizer/EventTree" (IsPP false) | Name of the ggHi tree.  Default value depends on whether it is pp mode or not |
| DoMCHiBinShift | bool | true | Whether to shift PbPb MC hiBin to match data better |
| MCHiBinShift | double | 3 | Amount of hiBin to shift PbPb MC.  Defaults to Kaya AN number (1.5%, or 3) |
| DoMCHiBinShift | bool | true | Whether to shift PbPb MC hiBin to match data better |
| WithProgressBar | bool | false | Whether to display progress bar or not |

