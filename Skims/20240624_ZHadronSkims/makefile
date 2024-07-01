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
	./Execute --Input $(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_14.root \
		--Output /tmp/$(USER)/Test.root --Fraction 0.1 \
		--Background $(PbPbBackgroundData)/HiForestMiniAOD_118.root,$(PbPbBackgroundData)/HiForestMiniAOD_123.root,$(PbPbBackgroundData)/HiForestMiniAOD_127.root,$(PbPbBackgroundData)/HiForestMiniAOD_138.root,$(PbPbBackgroundData)/HiForestMiniAOD_148.root,$(PbPbBackgroundData)/HiForestMiniAOD_16.root,$(PbPbBackgroundData)/HiForestMiniAOD_17.root,$(PbPbBackgroundData)/HiForestMiniAOD_184.root,$(PbPbBackgroundData)/HiForestMiniAOD_18.root,$(PbPbBackgroundData)/HiForestMiniAOD_194.root \
		`./Expand Setting.dh PbPbSignalData Nominal`

TestRunBackgroundSimple: Execute
	./Execute --Input $(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_14.root \
		--Output /tmp/$(USER)/Test.root --Fraction 1 \
		--Background $(PbPbBackgroundData)/HiForestMiniAOD_118.root \
		`./Expand Setting.dh PbPbSignalData Nominal`

TestRunMCBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/Test.root --Fraction 0.1 \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_10.root \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

TestRunMCGenBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/Test.root --Fraction 0.1 \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_10.root \
		`./Expand Setting.dh PbPbSignalMC Nominal`

TestRunRandomMCBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/Test.root --Fraction 0.1 \
		--Background `ls $(PbPbBackgroundMC)/*root | head -n 37 | tr '\n' ','` \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

TestRunRandomMCGenBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_10.root \
		--Output /tmp/$(USER)/Test.root --Fraction 0.1 \
		--Background `ls $(PbPbBackgroundMC)/*root | head -n 37 | tr '\n' ','` \
		`./Expand Setting.dh PbPbSignalMC Nominal`

EOSBase = /eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024
AFSBase = /afs/cern.ch/user/p/pchou
CERNBOX = /eos/home-p/pchou

PrepareToRun:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		20 LogData_v1_ee/ $(EOSBase)/OutputData_v1_ee/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRun3:
	bash PrepareInput.sh \
		'$(PbPbSignalDataRun3Multi)/*/' \
		100 LogData3_v1/ $(EOSBase)/OutputData3_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRunMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v1_ee/ $(EOSBase)/OutputMC_v1_ee/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunMC3:
	bash PrepareInput.sh \
		$(PbPbSignalMCRun3)/ \
		5 LogMC3_v1/ $(EOSBase)/OutputMC3_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunGenMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogGenMC_v1_ee/ $(EOSBase)/OutputMCGen_v1_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalGenMC Nominal`

PrepareToRunPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogPP_v1_ee/ $(EOSBase)/OutputPPMC_v1_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalMC Nominal`

PrepareToRunGenPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogGenPP_v1_ee/ $(EOSBase)/OutputPPMCGen_v1_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalGenMC Nominal`

PrepareToRunPPData:
	bash PrepareInput.sh \
		'$(PPSignalData)/*/' \
		100 LogPPData_v1_ee/ $(EOSBase)/OutputPPData_v1_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalData Nominal`

PrepareToRunLoose:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataLoose_v1/ $(EOSBase)/OutputDataLoose_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Loose`

PrepareToRunMCLoose:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCLoose_v1/ $(EOSBase)/OutputMCLoose_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Loose`

PrepareToRunTight:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataTight_v1/ $(EOSBase)/OutputDataTight_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Tight`

PrepareToRunMCTight:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCTight_v1/ $(EOSBase)/OutputMCTight_v1/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Tight`

PrepareToRunBackground:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		10 LogDatabkg_v1/ $(EOSBase)/OutputDatabkg_v1/ \
		--Background /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_10.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_13.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_17.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_20.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_22.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_34.root,/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_37.root \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

PrepareToRunRandomBackground:
	bash PrepareInputRandomBackground.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		10 LogDataBackground_v1/ $(EOSBase)/OutputDataBackground_v1/ \
		$(PbPbBackgroundData)/ 100 \
		`./Expand Setting.dh PbPbBackgroundData Nominal`

PrepareToRunMCBackground:
	bash PrepareInput.sh \
		$(PbPbSignalMC) \
		1 LogMCBackground_v1/ $(EOSBase)/OutputMCBackground_v1/ \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_10.root,$(PbPbBackgroundMC)/HiForestMiniAOD_11.root,$(PbPbBackgroundMC)/HiForestMiniAOD_12.root,$(PbPbBackgroundMC)/HiForestMiniAOD_13.root,$(PbPbBackgroundMC)/HiForestMiniAOD_14.root,$(PbPbBackgroundMC)/HiForestMiniAOD_15.root,$(PbPbBackgroundMC)/HiForestMiniAOD_16.root,$(PbPbBackgroundMC)/HiForestMiniAOD_17.root,$(PbPbBackgroundMC)/HiForestMiniAOD_18.root,$(PbPbBackgroundMC)/HiForestMiniAOD_19.root,$(PbPbBackgroundMC)/HiForestMiniAOD_1.root,$(PbPbBackgroundMC)/HiForestMiniAOD_20.root \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

PrepareToRunRandomMCBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		1 LogMCBackground_v1/ $(EOSBase)/OutputMCBackground_v1/ \
		$(PbPbBackgroundMC)/ 15 \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

PrepareToRunRandomMCGenBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		10 LogMCGenbkg_v1/ $(EOSBase)/OutputMCbkgGen_v1/ \
		$(PbPbBackgroundMC)/ 5 \
 		`./Expand Setting.dh PbPbBackgroundGenMC Nominal`

