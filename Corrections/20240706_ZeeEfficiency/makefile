default: Prepare TestRun

DYMC    = /eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20240605_ZHadronMCDY/240605_151412/0000/
CERNBOX = /eos/home-p/pchou/figs/ZHadron2024/ZeeEfficiency/20240911/
EOSBase = /eos/cms/store/group/phys_heavyions/pchou/HIZHadron2024/Corrections/20240706_ZeeEfficiency/
AFSBase = /afs/cern.ch/user/p/pchou/HIZHadron2024/Corrections/20240706_ZeeEfficiency/
PPMC    = /eos/cms/store/group/phys_heavyions/chenyi/pp2017/Forest/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2/AODSIM/DYJetsToLL_MLL-50_TuneCP5_5020GeV-amcatnloFXFX-pythia8/20231020_ZHadronMLLWithMuTreePP/231020_152224/0000/

Prepare:
	ln -sf $(DYMC) DYMC
	ln -sf $(PPMC) PPMC

TestRun: Execute
	./Execute --Input `ls DYMC/* | head -n1 | tr '\n' ','` \
		--Output EfficiencyPlotsPbPb.pdf --RootOutput /tmp/pchou/ZEfficiencyPbPb.root --Fraction 1.00
	mv /tmp/pchou/ZEfficiencyPbPb.root .
	mkdir -p $(CERNBOX)
	cp EfficiencyPlotsPbPb.pdf $(CERNBOX)

TestRunPP: Execute
	./Execute --Input `ls PPMC/* | head -n5 | tr '\n' ','` \
		--Output EfficiencyPlotsPP.pdf --RootOutput /tmp/pchou/ZEfficiencyPP.root --Fraction 1.00 --IsPP true
	mv /tmp/pchou/ZEfficiencyPP.root .
	mkdir -p $(CERNBOX)
	cp EfficiencyPlotsPP.pdf $(CERNBOX)

TestRunAll: Execute
	./Execute --Input `ls DYMC/* | tr '\n' ','` \
		--Output AllEfficiencyPlots.pdf --RootOutput /tmp/pchou/ZEfficiencyAll.root --Fraction 1.00
	mv /tmp/pchou/ZEfficiencyAll.root .
	mkdir -p $(CERNBOX)
	cp AllEfficiencyPlots.pdf $(CERNBOX)

TestRunAllPP: Execute
	./Execute --Input `ls PPMC/* | tr '\n' ','` \
		--Output AllPPEfficiencyPlots.pdf --RootOutput /tmp/pchou/ZEfficiencyAllPP.root --Fraction 1.00
	mv /tmp/pchou/ZEfficiencyAllPP.root .
	mkdir -p $(CERNBOX)
	cp AllPPEfficiencyPlots.pdf $(CERNBOX)

Execute: DeriveZEfficiency.cpp
	g++ DeriveZEfficiency.cpp  -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \
		# $(FJExtra) $(FJCExtra)

PrepareToRun: Execute
	bash PrepareInput.sh \
		$(DYMC) 2 LogPbPb_0911_vF/ $(AFSBase)/PbPb_0911_vF/ \
		--Fraction 1.00

PrepareToRunPP: Execute
	bash PrepareInput.sh \
		$(PPMC) 2 LogPP_0911_v0/ $(AFSBase)/pp_0911_v0/ \
		--Fraction 1.00 --IsPP true
