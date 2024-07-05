PbPbSignalMC         := $(shell DHQuery Setting.dh Sample PbPbSignalMCRun2)
PbPbSignalSingleMuon := $(shell DHQuery Setting.dh Sample PbPbSignalDataRun2)
PbPbBackgroundMC     := $(shell DHQuery Setting.dh Sample PbPbBackgroundMCRun2)
PbPbBackgroundData   := $(shell DHQuery Setting.dh Sample PbPbBackgroundDataRun2)
PPSignalMC           := $(shell DHQuery Setting.dh Sample PPSignalMCRun2)
PPSignalData         := $(shell DHQuery Setting.dh Sample PPSignalDataRun2)
PbPbSignalMCRun3     := $(shell DHQuery Setting.dh Sample PbPbSignalMCRun3)
PbPbSignalDataRun3   := $(shell DHQuery Setting.dh Sample PbPbSignalDataRawPrime0)
PbPbSignalDataRun3Multi := /eos/cms/store/group/phys_heavyions/pchou/PbPb2023/Forest/HIPhysicsRawPrime*/HIRun2023A-PromptReco-v2/MINIAOD/HIPhysicsRawPrime*/20240607_ZHadronPromptRecoV2_SingleMu12_RawPrime*/*/

default: Prepare TestRun

Prepare:
	bash DefineNominalSetting.sh
	chmod u+x Expand
	
Execute: DoZHadronSkim.cpp
	g++ DoZHadronSkim.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \
		# $(FJExtra) $(FJCExtra)

TestRun: Execute
	./Execute --Input  $(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_696.root \
		--Output /tmp/$(USER)/Test.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalData Nominal`

TestRun3: Execute
	./Execute --Input `ls $(PbPbSignalDataRun3)/0000/*root | head -n 100 | tr '\n' ','` \
		--Output /tmp/$(USER)/TestRun3.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalData Nominal`

TestRunMC: Execute
	time ./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_1.root \
		--Output /tmp/$(USER)/TestMC.root --Fraction 0.1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalMC Nominal`

TestMCRun3: Execute
	time ./Execute --Input $(PbPbSignalMCRun3)/HiForestMiniAOD_MC_1.root \
		--Output /tmp/$(USER)/TestMCRun3.root --Fraction 0.1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalMC Nominal`

TestRunGen: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_114.root \
		--Output /tmp/$(USER)/TestGenMC.root --Fraction 0.1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalGenMC Nominal`

TestRunGenPP: Execute
	./Execute --Input `ls $(PPSignalMC)/*root | head -n 2 | tr '\n' ','` \
		--Output /tmp/$(USER)/TestGenPP.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PPSignalGenMC Nominal`

TestRunPP: Execute
	./Execute --Input `ls $(PPSignalMC)/*root | head -n 2 | tr '\n' ','` \
		--Output /tmp/$(USER)/TestPP.root --Fraction 0.1 --WithProgressBar true \
		`./Expand Setting.dh PPSignalMC Nominal`

TestRunPPData: Execute
	./Execute --Input $(PPSignalData)/0000/HiForestAOD_12.root \
		--Output /tmp/$(USER)/TestPPData.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PPSignalData Nominal`

TestRunBackground: Execute
	./Execute --Input $(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_696.root \
		--Output /tmp/$(USER)/TestDataBkg.root --Fraction 1 --DoMuon false --WithProgressBar true \
		--Background `ls $(PbPbBackgroundData)/*root | head -n 5 | tr '\n' ','` \
		`./Expand Setting.dh PbPbBackgroundData Nominal`

TestRunMCBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/TestMCBkg.root --Fraction 0.1 --DoMuon false  --WithProgressBar true \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_10.root \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

TestRunMCGenBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/TestMCGenBkg.root --Fraction 0.1  --DoMuon false --WithProgressBar true \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_10.root \
		`./Expand Setting.dh PbPbBackgroundGenMC Nominal`

EOSBase = /eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024
AFSBase = /afs/cern.ch/user/p/pchou
CERNBOX = /eos/home-p/pchou

PrepareToRun:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		20 LogData_v1b_ee/ $(EOSBase)/OutputData_v1b_ee/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRun3:
	bash PrepareInput.sh \
		'$(PbPbSignalDataRun3Multi)/*/' \
		100 LogData3_v1b/ $(EOSBase)/OutputData3_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRunMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v1b_ee/ $(EOSBase)/OutputMC_v1b_ee/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunMC3:
	bash PrepareInput.sh \
		$(PbPbSignalMCRun3)/ \
		5 LogMC3_v1b/ $(EOSBase)/OutputMC3_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunGenMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogGenMC_v1b_ee/ $(EOSBase)/OutputMCGen_v1b_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalGenMC Nominal`

PrepareToRunPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogPP_v1b_ee/ $(EOSBase)/OutputPPMC_v1b_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalMC Nominal`

PrepareToRunGenPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogGenPP_v1b_ee/ $(EOSBase)/OutputPPMCGen_v1b_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalGenMC Nominal`

PrepareToRunPPData:
	bash PrepareInput.sh \
		'$(PPSignalData)/*/' \
		100 LogPPData_v1b_ee/ $(EOSBase)/OutputPPData_v1b_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalData Nominal`

PrepareToRunLoose:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataLoose_v1b/ $(EOSBase)/OutputDataLoose_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Loose`

PrepareToRunMCLoose:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCLoose_v1b/ $(EOSBase)/OutputMCLoose_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Loose`

PrepareToRunTight:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataTight_v1b/ $(EOSBase)/OutputDataTight_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Tight`

PrepareToRunMCTight:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCTight_v1b/ $(EOSBase)/OutputMCTight_v1b/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Tight`

PrepareToRunBackground:
	bash PrepareInputRandomBackground.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		10 LogDataBkg_v1b_ee/ $(EOSBase)/OutputDataBkg_v1b_ee/ \
		$(PbPbBackgroundData)/ 100  --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundData Nominal`

PrepareToRunMCBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		1 LogMCBkg_v1b_ee/ $(EOSBase)/OutputMCBkg_v1b_ee/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

PrepareToRunRandomMCGenBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		10 LogMCGenbkg_v1b_ee/ $(EOSBase)/OutputMCbkgGen_v1b_ee/ \
		$(PbPbBackgroundMC)/ 5 --DoMuon false \
 		`./Expand Setting.dh PbPbBackgroundGenMC Nominal`

