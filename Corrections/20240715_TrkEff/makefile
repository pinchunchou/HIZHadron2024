default: Prepare TestRun

DYMC    = /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20240605_ZHadronMCDY/240605_151412/0000/
CERNBOX = /eos/home-p/pchou/figs/ZHadron2024/TrkEff/20240721/
EOSBase = /eos/cms/store/group/phys_heavyions/pchou/HIZHadron2024/Corrections/20240715_TrkEff/20240721/
AFSBase = /afs/cern.ch/user/p/pchou/HIZHadron2024/Corrections/20240715_TrkEff/
PPMC    = /eos/cms/store/group/phys_heavyions/chenyi/pp2017/Forest/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2/AODSIM/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/20231020_ZHadronMLLWithMuTreePP/231020_152224/0000/

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_0_20.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_20_60.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_60_100.root,${ProjectBase}/CommonCode/root/20230522_TrackResidualCorrection_V9_100_200.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20230531_TrackResidualCorrection_V12_pp.root"

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V11_0_20.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V11_20_60.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V11_60_100.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V11_100_200.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V11_pp.root"

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V12_0_20.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V12_20_60.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V12_60_100.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V12_100_200.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V12_pp.root"

TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V13_0_20.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V13_20_60.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V13_60_100.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V13_100_200.root"
TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V13_pp.root"

#TrackResidualPathPbPb="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_0_20.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_20_60.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_60_100.root,${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_100_200.root"
#TrackResidualPathPP="${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_pp.root"

Prepare:
	ln -sf $(DYMC) DYMC
	ln -sf $(PPMC) PPMC

TestRun: Execute
	./Execute --Input `ls DYMC/* | head -n1 | tr '\n' ','` \
		--Output EffPlotsPbPb_0721_0_20_V14.pdf --RootOutput TrkEffPbPb_0721_0_20_V14.root --Fraction 0.2 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20
	mkdir -p $(CERNBOX)
	mv EffPlotsPbPb_0721_0_20_V14.pdf $(CERNBOX)

TestRunPP: Execute
	./Execute --Input `ls PPMC/* | head -n1 | tr '\n' ','` \
		--Output EffPlotsPP_0721_V14.pdf --RootOutput TrkEffPP_0721_V14.root --Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true 
	mkdir -p $(CERNBOX)
	mv EffPlotsPP_0721_V14.pdf $(CERNBOX)

TestRunAll: Execute
	./Execute --Input `ls DYMC/* | tr '\n' ','` \
		--Output AllPbPbEffPlots.pdf --RootOutput /tmp/pchou/TrkEffAll.root --Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb) --DoIteration true \
		--MinHiBin 0 --MaxHiBin 100
	mv /tmp/pchou/TrkEffAll.root .
	mkdir -p $(CERNBOX)
	cp AllPbPbEffPlots.pdf $(CERNBOX)

TestRunAllPP: Execute
	./Execute --Input `ls PPMC/* | tr '\n' ','` \
		--Output AllPPEffPlots.pdf --RootOutput /tmp/pchou/TrkEffAllPP.root --Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true 
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
		$(DYMC) 2 LogPbPb_0_20_V14/ $(EOSBase)/PbPb_0_20_V14/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 0 --MaxHiBin 20

PrepareToRun2: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_20_60_V14/ $(EOSBase)/PbPb_20_60_V14/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 20 --MaxHiBin 60

PrepareToRun3: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_60_100_V14/ $(EOSBase)/PbPb_60_100_V14/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 60 --MaxHiBin 100

PrepareToRun4: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_100_200_V14/ $(EOSBase)/PbPb_100_200_V14/ \
		--Fraction 1.00 \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPbPb)  --DoIteration true \
		--MinHiBin 100 --MaxHiBin 200

PrepareToRunPP: Execute
	bash PrepareInput.sh \
		$(PPMC) 2 LogPP_V14/ $(EOSBase)/pp_V14/ \
		--Fraction 1.00 --IsPP true \
		--TrackEfficiencyPath ${ProjectBase}/CommonCode/root/ \
		--TrackResidualPath $(TrackResidualPathPP) --DoIteration true 

Merge: MergePbPb1 MergePbPb2 MergePbPb3 MergePbPb4 MergePP

MergePbPb1: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_0_20_V14/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_0_20_V14.pdf \
	--RootOutput  ${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_0_20.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 0 \
	> log/20240721_TrackResidualCorrection_V14_0_20.out 2> log/20240721_TrackResidualCorrection_V14_0_20.err < /dev/null &

MergePbPb2: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_20_60_V14/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_20_60_V14.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_20_60.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 1 \
	> log/20240721_TrackResidualCorrection_V14_20_60.out 2> log/20240721_TrackResidualCorrection_V14_20_60.err < /dev/null &

MergePbPb3: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_60_100_V14/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_60_100_V14.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_60_100.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 2 \
	> log/20240721_TrackResidualCorrection_V14_60_100.out 2> log/20240721_TrackResidualCorrection_V14_60_100.err < /dev/null &

MergePbPb4: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/PbPb_100_200_V14/*.root" \
	--Output $(CERNBOX)/AllPbPbEffPlots_100_200_V14.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_100_200.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPbPb) --TrackResIdx 3 \
	> log/20240721_TrackResidualCorrection_V14_100_200.out 2> log/20240721_TrackResidualCorrection_V14_100_200.err < /dev/null &

MergePP: ExecuteMerge
	mkdir -p $(CERNBOX) log
	./ExecuteMerge --Input "$(EOSBase)/pp_V14/*.root" \
	--Output $(CERNBOX)/AllPPEffPlots_V14.pdf \
	--RootOutput ${ProjectBase}/CommonCode/root/20240721_TrackResidualCorrection_V14_pp.root \
	--DoIteration true --TrackResidualPath $(TrackResidualPathPP) --TrackResIdx 0 \
	> log/20240721_TrackResidualCorrection_V14_pp.out 2> log/20240721_TrackResidualCorrection_V14_pp.err < /dev/null &

ExecuteMerge: mergeHist.cpp
	g++ mergeHist.cpp  -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteMerge `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \
