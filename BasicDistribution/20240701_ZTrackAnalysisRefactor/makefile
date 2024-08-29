CERNBOX = /eos/home-p/pchou
EOSBase = /eos/cms/store/group/phys_heavyions/pchou
AFSBase = /afs/cern.ch/user/p/pchou

default: TestRun

TestRun: Execute
	mkdir -p log
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData_v2f_ee/ \
		--Output /tmp/$(USER)/GraphDataSig_v2f_ee.root --Fraction 0.1 

Run2All: Run2Sig Run2PP Run2Bkg
Run2Sig: Run2DataSig Run2MCSig Run2MCSigGen Run2MCSigGen0Sub 
Run2PP: Run2PPMCSig Run2PPMCSig0NPU Run2PPMCSig1Vtx Run2PPData Run2PPData1Vtx Run2PPMCSigGen Run2PPMCSigGen0NPU 
Run2Bkg: Run2DataBkg Run2MCBkg Run2MCBkgGen Run2MCBkgGen0Sub
Run2MC: Run2MCSig Run2MCSigGen Run2MCSigGen0Sub Run2MCBkg Run2MCBkgGen Run2PPMCSig Run2PPMCSig0NPU Run2PPMCSig1Vtx Run2PPMCSigGen Run2PPMCSigGen0NPU 

Run3All: Run3Sig Run3PP
Run3Sig: Run3DataSig Run3MCSig
Run3PP: Run3PPMCSig Run3PPMCSig0NPU Run3PPMCSig1Vtx Run3PPData Run3PPData1Vtx

Run2DataSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphDataSig_v2f_ee.root --Fraction 1 \
		 > log/GraphDataSig_v2f_ee.out 2> log/GraphDataSig_v2f_ee.err < /dev/null &

Run2MCSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMC_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSig_v2f_ee.root --Fraction 1 \
		 > log/GraphMCSig_v2f_ee.out 2> log/GraphMCSig_v2f_ee.err < /dev/null &

Run2MCSig1Vtx: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMC_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSig1Vtx_v2f_ee.root --Fraction 1  --OnlyOneNVertex true \
		 > log/GraphMCSig1Vtx_v2f_ee.out 2> log/GraphMCSig1Vtx_v2f_ee.err < /dev/null &

Run2PPMCSig: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC_v2f_ee.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPMC_v2f_ee.out 2> log/GraphPPMC_v2f_ee.err < /dev/null &

Run2PPMCSig0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC0NPU_v2f_ee.root --IgnoreCentrality true --OnlyZeroNPU true --Fraction 1 \
		 > log/GraphPPMC0NPU_v2f_ee.out 2> log/GraphPPMC0NPU_v2f_ee.err < /dev/null &

Run2PPMCSig1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC1Vtx_v2f_ee.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPMC1Vtx_v2f_ee.out 2> log/GraphPPMC1Vtx_v2f_ee.err < /dev/null &

Run2PPData: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData_v2f_ee.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPData_v2f_ee.out 2> log/GraphPPData_v2f_ee.err < /dev/null &

Run2PPData1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData1Vtx_v2f_ee.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPData1Vtx_v2f_ee.out 2> log/GraphPPData1Vtx_v2f_ee.err < /dev/null &

Run2MCSigGen: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSigGen_v2f_ee.root --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCSigGen_v2f_ee.out 2> log/GraphMCSigGen_v2f_ee.err < /dev/null &

Run2MCSigGen0Sub: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSigGen0Sub_v2f_ee.root --DoGenCorrelation true  --OnlyZeroSub true --Fraction 1 \
		 > log/GraphMCSigGen0Sub_v2f_ee.out 2> log/GraphMCSigGen0Sub_v2f_ee.err < /dev/null &

Run2PPMCSigGen: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMCGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMCGen_v2f_ee.root --IgnoreCentrality true --DoGenCorrelation true  --Fraction 1 \
		 > log/GraphPPMCGen_v2f_ee.out 2> log/GraphPPMCGen_v2f_ee.err < /dev/null &

Run2PPMCSigGen0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMCGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMCGen0NPU_v2f_ee.root --IgnoreCentrality true --OnlyZeroNPU true --DoGenCorrelation true --Fraction 1 \
		 > log/GraphPPMCGen0NPU_v2f_ee.out 2> log/GraphPPMCGen0NPU_v2f_ee.err < /dev/null &

Run2DataBkg: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputDataBkg_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphDataBkg_v2f_ee.root --Fraction 1 \
		 > log/GraphDataBkg_v2f_ee.out 2> log/GraphDataBkg_v2f_ee.err < /dev/null &

Run2MCBkg: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCBkg_v2e_ee_newresMB/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkg_v2e_ee_newresMB.root --Fraction 1 \
		 > log/GraphMCBkg_v2e_ee_newresMB.out 2> log/GraphMCBkg_v2e_ee_newresMB.err < /dev/null &

Run2MCBkg1Vtx: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCBkg_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkg1Vtx_v2f_ee.root --Fraction 1 --OnlyOneNVertex true \
		 > log/GraphMCBkg1Vtx_v2f_ee.out 2> log/GraphMCBkg1Vtx_v2f_ee.err < /dev/null &

Run2MCBkgGen: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCbkgGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkgGen_v2f_ee.root --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCBkgGen_v2f_ee.out 2> log/GraphMCBkgGen_v2f_ee.err < /dev/null &

Run2MCBkgGen0Sub: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCbkgGen_v2f_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkgGen0Sub_v2f_ee.root --OnlyZeroSub true --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCBkgGen0Sub_v2f_ee.out 2> log/GraphMCBkgGen0Sub_v2f_ee.err < /dev/null &

Run3DataSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphData3Sig_v2f.root --Fraction 1 \
		 > log/GraphData3Sig_v2f.out 2> log/GraphData3Sig_v2f.err < /dev/null &

Run3MCSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMC3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMC3Sig_v2f.root --Fraction 1 \
		 > log/GraphMC3Sig_v2f.out 2> log/GraphMC3Sig_v2f.err < /dev/null &

Run3PPMCSig: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC3_v2f.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPMC3_v2f.out 2> log/GraphPPMC3_v2f.err < /dev/null &

Run3PPMCSig0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC30NPU_v2f.root --IgnoreCentrality true --OnlyZeroNPU true --Fraction 1 \
		 > log/GraphPPMC30NPU_v2f.out 2> log/GraphPPMC30NPU_v2f.err < /dev/null &

Run3PPMCSig1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC31Vtx_v2f.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPMC31Vtx_v2f.out 2> log/GraphPPMC31Vtx_v2f.err < /dev/null &

Run3PPData: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData3_v2f.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPData3_v2f.out 2> log/GraphPPData3_v2f.err < /dev/null &

Run3PPData1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData3_v2f/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData31Vtx_v2f.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPData31Vtx_v2f.out 2> log/GraphPPData31Vtx_v2f.err < /dev/null &

Execute: ZTrackAnalysis.cpp
	mkdir -p log
	g++ ZTrackAnalysis.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17

Draw: ExecuteDraw
	./ExecuteDrawRefactor

ExecuteDraw: ZtrackDrawRefactor.cpp
	g++ ZtrackDrawRefactor.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteDrawRefactor `root-config --cflags --libs` -std=c++17

DrawBasicBkgSub: ExecuteBasicBkgSub
	./ExecuteBasicBkgSub

ExecuteBasicBkgSub: ZBasicBkgSub.cpp
	g++ ZBasicBkgSub.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteBasicBkgSub `root-config --cflags --libs` -std=c++17

DrawZMass: ExecuteZMassDraw
	./ExecuteDrawZMass

ExecuteZMassDraw: ZMassDraw.cpp
	g++ ZMassDraw.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteDrawZMass `root-config --cflags --libs` -std=c++17

DrawSimple: ExecuteSimple
	mkdir -p log
	time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 10 -TptL 1 -TptH 2  -isRres true \
	 > log/DrawSimple_40_200_0_10_1_2_Rres.out 2> log/DrawSimple_40_200_0_10_1_2_Rres.err < /dev/null &
	time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 10 -TptL 1 -TptH 2  -isRres false \
	> log/DrawSimple_40_200_0_10_1_2.out 2> log/DrawSimple_40_200_0_10_1_2.err < /dev/null &
	time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 10 -TptL 1 -TptH 2  -noZ true \
	> log/DrawSimple_40_200_0_10_1_2_noZ.out 2> log/DrawSimple_40_200_0_10_1_2_noZ.err < /dev/null &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 10 -TptL 2 -TptH 4  &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 10 -TptL 4 -TptH 10  &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 0 -centH 30 -TptL 1 -TptH 2  &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 10 -centH 30 -TptL 1 -TptH 2 &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 30 -centH 50 -TptL 1 -TptH 2 &
	#time ./ExecuteSimple -ptL 40 -ptH 200 -centL 50 -centH 90 -TptL 1 -TptH 2  &

ExecuteSimple: DrawSimple.cpp
	g++ DrawSimple.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteSimple `root-config --cflags --libs` -std=c++17


DrawNTrk: ExecuteNTrk
	time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 0  -centH 10 -TptL 1 -TptH 2  -MaxNtrk 4000 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 0  -centH 10 -TptL 2 -TptH 4  -MaxNtrk 1000 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 0  -centH 10 -TptL 4 -TptH 10 -MaxNtrk 100 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 0  -centH 30 -TptL 1 -TptH 2  -MaxNtrk 4000 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 10 -centH 30 -TptL 1 -TptH 2  -MaxNtrk 4000 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 30 -centH 50 -TptL 1 -TptH 2  -MaxNtrk 1500 &
	#time ./ExecuteNTrk -ptL 40 -ptH 200 -centL 50 -centH 90 -TptL 1 -TptH 2  -MaxNtrk 600 &

ExecuteNTrk: DrawNTrk_simple.cpp
	g++ DrawNTrk_simple.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteNTrk `root-config --cflags --libs` -std=c++17