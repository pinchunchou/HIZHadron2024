default: Prepare TestRun

DYMC    = /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20240605_ZHadronMCDY/240605_151412/0000/
CERNBOX = /eos/home-p/pchou/figs/ZHadron2024/TrkEff/20240725/
EOSBase = /eos/cms/store/group/phys_heavyions/pchou/HIZHadron2024/Corrections/20240715_TrkEff/20240725/
AFSBase = /afs/cern.ch/user/p/pchou/HIZHadron2024/Corrections/20240715_TrkEff/
PPMC    = /eos/cms/store/group/phys_heavyions/chenyi/pp2017/Forest/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2/AODSIM/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/20231020_ZHadronMLLWithMuTreePP/231020_152224/0000/
MBMC    = /eos/cms/store/group/phys_heavyions/chenyi/Samples/store/user/mitaylor/PhotonJet/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20230127PbPbMCMB/230127_211118/

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_0_20.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_20_60.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_60_100.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_100_200.root"
#TrackResidualPathPbPbMB="${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_0_20.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_20_60.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_60_100.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_100_200.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20230531_TrackResidualCorrection_V12_pp.root"

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_0_20.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_20_60.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_60_100.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_100_200.root"
#TrackResidualPathPbPbMB="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_0_20_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_20_60_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_60_100_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_100_200_MB.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V1_pp.root"

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_0_20.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_20_60.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_60_100.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_100_200.root"
#TrackResidualPathPbPbMB="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_0_20_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_20_60_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_60_100_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_100_200_MB.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V2_pp.root"

TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_0_20.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_20_60.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_60_100.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_100_200.root"
TrackResidualPathPbPbMB="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_0_20_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_20_60_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_60_100_MB.root,${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_100_200_MB.root"
TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_V3_pp.root"

VERSION = V4

Prepare:
	ln -sf $(DYMC) DYMC
	ln -sf $(PPMC) PPMC
	ln -sf $(MBMC) MBMC

TestRun: Execute
	./Execute --Input `ls DYMC/* | head -n1 | tr '\n' ','` \
		--Output EffPlotsPbPb_0725_0_20_${VERSION}.pdf --RootOutput TrkEffPbPb_0725_0_20_${VERSION}.root --Fraction 0.2 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20 --DoZSelection true
	mkdir -p $(CERNBOX)
	mv EffPlotsPbPb_0725_0_20_${VERSION}.pdf $(CERNBOX)


TestRunPP: Execute
	./Execute --Input `ls PPMC/* | head -n1 | tr '\n' ','` \
		--Output EffPlotsPP_0725_${VERSION}.pdf --RootOutput TrkEffPP_0725_${VERSION}.root --Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true  --DoZSelection true
	mkdir -p $(CERNBOX)
	mv EffPlotsPP_0725_${VERSION}.pdf $(CERNBOX)

TestRunMB: Execute
	./Execute --Input `ls MBMC/* | head -n1 | tr '\n' ','` \
		--Output EffPlotsPbPbMB_0725_0_20_${VERSION}.pdf --RootOutput TrkEffPbPbMB_0725_0_20_${VERSION}.root --Fraction 0.2 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPbMB)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20 --DoZSelection false
	mkdir -p $(CERNBOX)
	mv EffPlotsPbPbMB_0725_0_20_${VERSION}.pdf $(CERNBOX)

TestRunAll: Execute
	./Execute --Input `ls DYMC/* | tr '\n' ','` \
		--Output AllPbPbEffPlots.pdf --RootOutput /tmp/pchou/TrkEffAll.root --Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb) --DoIteration true \
		--MinHiBin 0 --MaxHiBin 100  --DoZSelection true
	mv /tmp/pchou/TrkEffAll.root .
	mkdir -p $(CERNBOX)
	cp AllPbPbEffPlots.pdf $(CERNBOX)

TestRunAllPP: Execute
	./Execute --Input `ls PPMC/* | tr '\n' ','` \
		--Output AllPPEffPlots.pdf --RootOutput /tmp/pchou/TrkEffAllPP.root --Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true  --DoZSelection true
	mv /tmp/pchou/TrkEffAllPP.root .
	mkdir -p $(CERNBOX)
	cp AllPPEffPlots.pdf $(CERNBOX)

Execute: DeriveTrkEff.cpp
	g++ DeriveTrkEff.cpp  -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \
		# $(FJExtra) $(FJCExtra)

PrepareToRun1: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_0_20_${VERSION}/ $(EOSBase)/PbPb_0_20_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20  --DoZSelection true

PrepareToRun2: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_20_60_${VERSION}/ $(EOSBase)/PbPb_20_60_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 20 --MaxHiBin 60  --DoZSelection true

PrepareToRun3: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_60_100_${VERSION}/ $(EOSBase)/PbPb_60_100_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 60 --MaxHiBin 100  --DoZSelection true

PrepareToRun4: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_100_200_${VERSION}/ $(EOSBase)/PbPb_100_200_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 100 --MaxHiBin 200  --DoZSelection true

PrepareToRunPP: Execute
	bash PrepareInput.sh \
		$(PPMC) 2 LogPP_${VERSION}/ $(EOSBase)/pp_${VERSION}/ \
		--Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true  --DoZSelection true

PrepareToRunMB1: Execute
	bash PrepareInput.sh \
		$(MBMC) 2 LogPbPbMB_0_20_${VERSION}/ $(EOSBase)/PbPbMB_0_20_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPbMB)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20  --DoZSelection false

PrepareToRunMB2: Execute
	bash PrepareInput.sh \
		$(MBMC) 2 LogPbPbMB_20_60_${VERSION}/ $(EOSBase)/PbPbMB_20_60_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPbMB)  --DoIteration true \
		--MinHiBin 20 --MaxHiBin 60  --DoZSelection false

PrepareToRunMB3: Execute
	bash PrepareInput.sh \
		$(MBMC) 2 LogPbPbMB_60_100_${VERSION}/ $(EOSBase)/PbPbMB_60_100_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPbMB)  --DoIteration true \
		--MinHiBin 60 --MaxHiBin 100  --DoZSelection false

PrepareToRunMB4: Execute
	bash PrepareInput.sh \
		$(MBMC) 2 LogPbPbMB_100_200_${VERSION}/ $(EOSBase)/PbPbMB_100_200_${VERSION}/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPbMB)  --DoIteration true \
		--MinHiBin 100 --MaxHiBin 200  --DoZSelection false

Merge: MergePbPb MergePP MergeMB
MergePbPb: MergePbPb1 MergePbPb2 MergePbPb3 MergePbPb4
MergeMB: MergePbPbMB1 MergePbPbMB2 MergePbPbMB3 MergePbPbMB4

MergePbPb1: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_0_20_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_0_20_${VERSION}.pdf \
	--RootOutput  ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_0_20.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 0 \
	> log/20240725_TrackResidualCorrection_${VERSION}_0_20.out 2> log/20240725_TrackResidualCorrection_${VERSION}_0_20.err < /dev/null &

MergePbPb2: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_20_60_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_20_60_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_20_60.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 1 \
	> log/20240725_TrackResidualCorrection_${VERSION}_20_60.out 2> log/20240725_TrackResidualCorrection_${VERSION}_20_60.err < /dev/null &

MergePbPb3: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_60_100_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_60_100_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_60_100.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 2 \
	> log/20240725_TrackResidualCorrection_${VERSION}_60_100.out 2> log/20240725_TrackResidualCorrection_${VERSION}_60_100.err < /dev/null &

MergePbPb4: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_100_200_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_100_200_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_100_200.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 3 \
	> log/20240725_TrackResidualCorrection_${VERSION}_100_200.out 2> log/20240725_TrackResidualCorrection_${VERSION}_100_200.err < /dev/null &

MergePP: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/pp_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPPEffPlots_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_pp.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPP) --TrackResIdx 0 \
	> log/20240725_TrackResidualCorrection_${VERSION}_pp.out 2> log/20240725_TrackResidualCorrection_${VERSION}_pp.err < /dev/null &

MergePbPbMB1: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPbMB_0_20_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbMBEffPlots_0_20_${VERSION}.pdf \
	--RootOutput  ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_0_20_MB.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPbMB) --TrackResIdx 0 \
	> log/20240725_TrackResidualCorrection_${VERSION}_0_20_MB.out 2> log/20240725_TrackResidualCorrection_${VERSION}_0_20_MB.err < /dev/null &

MergePbPbMB2: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPbMB_20_60_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbMBEffPlots_20_60_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_20_60_MB.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPbMB) --TrackResIdx 1 \
	> log/20240725_TrackResidualCorrection_${VERSION}_20_60_MB.out 2> log/20240725_TrackResidualCorrection_${VERSION}_20_60_MB.err < /dev/null &

MergePbPbMB3: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPbMB_60_100_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbMBEffPlots_60_100_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_60_100_MB.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPbMB) --TrackResIdx 2 \
	> log/20240725_TrackResidualCorrection_${VERSION}_60_100_MB.out 2> log/20240725_TrackResidualCorrection_${VERSION}_60_100_MB.err < /dev/null &

MergePbPbMB4: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPbMB_100_200_${VERSION}/*.root" \
	--Output $(CERNBOX)/AllPbPbMBEffPlots_100_200_${VERSION}.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240725_TrackResidualCorrection_${VERSION}_100_200_MB.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPbMB) --TrackResIdx 3 \
	> log/20240725_TrackResidualCorrection_${VERSION}_100_200_MB.out 2> log/20240725_TrackResidualCorrection_${VERSION}_100_200_MB.err < /dev/null &

ExecuteMerge: mergeHist.cpp
	g++ mergeHist.cpp  -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteMerge `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \
