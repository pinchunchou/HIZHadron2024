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
	./Execute --Input $(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_696.root,$(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_697.root,$(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_698.root,$(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_699.root,$(PbPbSignalSingleMuon)/0000/HiForestMiniAOD_700.root \
		--Output /tmp/$(USER)/Test.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalData Nominal` --DoElectron false

TestRun3: Execute
	./Execute --Input `ls $(PbPbSignalDataRun3)/0000/*root | head -n 100 | tr '\n' ','` \
		--Output /tmp/$(USER)/TestRun3.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalData Nominal`

TestRunMC: Execute
	time ./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_1.root \
		--Output /tmp/$(USER)/TestMC.root --Fraction 1 --WithProgressBar true \
		`./Expand Setting.dh PbPbSignalMC Nominal`  --DoElectron false

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
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_1.root \
		--Output TestMCBkg.root --Fraction 1 --DoMuon false  --WithProgressBar true \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_8.root \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

TestRunMCBackgroundResGenMatch: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_1.root \
		--Output TestMCBkgResGenMatch.root --Fraction 1 --DoMuon false  --WithProgressBar true \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_8.root \
		`./Expand Setting.dh PbPbBackgroundMCResGenMatch Nominal`

TestRunMCGenBackground: Execute
	./Execute --Input $(PbPbSignalMC)/HiForestMiniAOD_1.root \
		--Output TestMCGenBkg.root --Fraction 1  --DoMuon false --WithProgressBar true \
		--Background $(PbPbBackgroundMC)/HiForestMiniAOD_8.root \
		`./Expand Setting.dh PbPbBackgroundGenMC Nominal`

EOSBase = /eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024
AFSBase = /afs/cern.ch/user/p/pchou
CERNBOX = /eos/home-p/pchou/SkimZHadron2024

PrepareToRun:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		20 LogData_v4d_ee/ $(EOSBase)/OutputData_v4d_ee/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRun3:
	bash PrepareInput.sh \
		'$(PbPbSignalDataRun3Multi)/*/' \
		100 LogData3_v4d/ $(EOSBase)/OutputData3_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Nominal`

PrepareToRunMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v4d_ee_v2/ $(CERNBOX)/OutputMC_v4d_ee_v2/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunMCNoZ:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v4d_ee_NoZ/ $(CERNBOX)/OutputMC_v4d_ee_NoZ/ \
		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalMCNoZ Nominal`

PrepareToRunMCNoShift:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v4d_ee_NoShift/ $(CERNBOX)/OutputMC_v4d_ee_NoShift/ \
		--Fraction 1 --DoMuon false --DoMCHiBinShift false `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunMCPF04:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMC_v4d_ee_PF04/ $(EOSBase)/OutputMC_v4d_ee_PF04/ \
		--Fraction 1 --DoMuon false --MinPFPT 0.4 `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunMC3:
	bash PrepareInput.sh \
		$(PbPbSignalMCRun3)/ \
		5 LogMC3_v4d/ $(EOSBase)/OutputMC3_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Nominal`

PrepareToRunGenMC:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogGenMC_v4d_ee/ $(EOSBase)/OutputMCGen_v4d_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PbPbSignalGenMC Nominal`

PrepareToRunPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogPP_v4d_ee/ $(EOSBase)/OutputPPMC_v4d_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalMC Nominal`

PrepareToRunGenPP:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogGenPP_v4d_ee/ $(EOSBase)/OutputPPMCGen_v4d_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalGenMC Nominal`

PrepareToRunPPPF04:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogPP_v4d_ee_PF04/ $(EOSBase)/OutputPPMC_v4d_ee_PF04/ \
 		--Fraction 1 --DoMuon false  --MinPFPT 0.4  `./Expand Setting.dh PPSignalMC Nominal`

PrepareToRunGenPPPF04:
	bash PrepareInput.sh \
		$(PPSignalMC)/ \
		10 LogGenPP_v4d_ee_PF04/ $(EOSBase)/OutputPPMCGen_v4d_ee_PF04/ \
 		--Fraction 1 --DoMuon false  --MinPFPT 0.4  `./Expand Setting.dh PPSignalGenMC Nominal`

PrepareToRunPPData:
	bash PrepareInput.sh \
		'$(PPSignalData)/*/' \
		100 LogPPData_v4d_ee/ $(EOSBase)/OutputPPData_v4d_ee/ \
 		--Fraction 1 --DoMuon false `./Expand Setting.dh PPSignalData Nominal`

PrepareToRunLoose:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataLoose_v4d/ $(EOSBase)/OutputDataLoose_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Loose`

PrepareToRunMCLoose:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCLoose_v4d/ $(EOSBase)/OutputMCLoose_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Loose`

PrepareToRunTight:
	bash PrepareInput.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		5 LogDataTight_v4d/ $(EOSBase)/OutputDataTight_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalData Tight`

PrepareToRunMCTight:
	bash PrepareInput.sh \
		$(PbPbSignalMC)/ \
		5 LogMCTight_v4d/ $(EOSBase)/OutputMCTight_v4d/ \
		--Fraction 1 `./Expand Setting.dh PbPbSignalMC Tight`

PrepareToRunBackground:
	bash PrepareInputRandomBackground.sh \
		'$(PbPbSignalSingleMuon)/*/' \
		10 LogDataBkg_v4d_ee/ $(EOSBase)/OutputDataBkg_v4d_ee/ \
		$(PbPbBackgroundData)/ 100  --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundData Nominal`

PrepareToRunMCBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee/ $(CERNBOX)/OutputMCBkg_v4d_ee/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMC Nominal`

PrepareToRunMCBackgroundRes:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_Rres_v2/ $(CERNBOX)/OutputMCBkg_v4d_ee_Rres_v2/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMCRes Nominal`

PrepareToRunMCBackgroundResNoShift:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_RresNoShift/ $(CERNBOX)/OutputMCBkg_v4d_ee_RresNoShift/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false --DoMCHiBinShift false \
		`./Expand Setting.dh PbPbBackgroundMCRes Nominal`

PrepareToRunMCBackgroundNoZ:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_NoZ/ $(CERNBOX)/OutputMCBkg_v4d_ee_NoZ/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMCNoZ Nominal`

PrepareToRunMCBackgroundPF04:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_PF04/ $(CERNBOX)/OutputMCBkg_v4d_ee_PF04/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false --MinPFPT 0.4 \
		`./Expand Setting.dh PbPbBackgroundMCPF04 Nominal`

PrepareToRunMCBackgroundResPF04:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_Rres_PF04/ $(CERNBOX)/OutputMCBkg_v4d_ee_Rres_PF04/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false --MinPFPT 0.4 \
		`./Expand Setting.dh PbPbBackgroundMCResPF04 Nominal`

PrepareToRunMCBackgroundResGenMatch:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_Rres_GenMatch/ $(CERNBOX)/OutputMCBkg_v4d_ee_Rres_GenMatch/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMCResGenMatch Nominal`

PrepareToRunMCBackgroundNoZGenMatch:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_NoZ_GenMatch/ $(CERNBOX)/OutputMCBkg_v4d_ee_NoZ_GenMatch/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMCNoZGenMatch Nominal`

PrepareToRunMCBackgroundGenMatch:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCBkg_v4d_ee_GenMatch/ $(CERNBOX)/OutputMCBkg_v4d_ee_GenMatch/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
		`./Expand Setting.dh PbPbBackgroundMCGenMatch Nominal`

PrepareToRunRandomMCGenBackground:
	bash PrepareInputRandomBackground.sh \
		$(PbPbSignalMC) \
		5 LogMCGenbkg_v4d_ee/ $(CERNBOX)/OutputMCbkgGen_v4d_ee/ \
		$(PbPbBackgroundMC)/ 15 --DoMuon false \
 		`./Expand Setting.dh PbPbBackgroundGenMC Nominal`

