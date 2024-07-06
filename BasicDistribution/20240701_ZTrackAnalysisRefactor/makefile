CERNBOX = /eos/home-p/pchou
EOSBase = /eos/cms/store/group/phys_heavyions/pchou
AFSBase = /afs/cern.ch/user/p/pchou

default: TestRun

TestRun: Execute
	mkdir -p log
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData_v1c_ee/ \
		--Output /tmp/$(USER)/GraphDataSig_v1c_ee.root --Fraction 0.1 

Run2All: Run2Sig Run2PP Run2Bkg
Run2Sig: Run2DataSig Run2MCSig Run2MCSigGen Run2MCSigGen0Sub 
Run2PP: Run2PPMCSig Run2PPMCSig0NPU Run2PPMCSig1Vtx Run2PPData Run2PPData1Vtx Run2PPMCSigGen Run2PPMCSigGen0NPU 
Run2Bkg: Run2DataBkg Run2MCBkg Run2MCBkgGen Run2MCBkgGen0Sub

Run3All: Run3Sig Run3PP
Run3Sig: Run3DataSig Run3MCSig
Run3PP: Run3PPMCSig Run3PPMCSig0NPU Run3PPMCSig1Vtx Run3PPData Run3PPData1Vtx

Run2DataSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphDataSig_v1c_ee_noZw.root --Fraction 1 --DoZWeight false \
		 > log/GraphDataSig_v1c_ee_noZw.out 2> log/GraphDataSig_v1c_ee_noZw.err < /dev/null &

Run2MCSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMC_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSig_v1c_ee_noZw.root --Fraction 1 --DoZWeight false \
		 > log/GraphMCSig_v1c_ee_noZw.out 2> log/GraphMCSig_v1c_ee_noZw.err < /dev/null &

Run2PPMCSig: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC_v1c_ee_noZw.root --IgnoreCentrality true --Fraction 1 --DoZWeight false \
		 > log/GraphPPMC_v1c_ee._noZwout 2> log/GraphPPMC_v1c_ee_noZw.err < /dev/null &

Run2PPMCSig0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC0NPU_v1c_ee_noZw.root --IgnoreCentrality true --OnlyZeroNPU true --Fraction 1 --DoZWeight false \
		 > log/GraphPPMC0NPU_v1c_ee_noZw.out 2> log/GraphPPMC0NPU_v1c_ee_noZw.err < /dev/null &

Run2PPMCSig1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC1Vtx_v1c_ee_noZw.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 --DoZWeight false \
		 > log/GraphPPMC1Vtx_v1c_ee_noZw.out 2> log/GraphPPMC1Vtx_v1c_ee_noZw.err < /dev/null &

Run2PPData: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData_v1c_ee_noZw.root --IgnoreCentrality true --Fraction 1 --DoZWeight false \
		 > log/GraphPPData_v1c_ee_noZw.out 2> log/GraphPPData_v1c_ee_noZw.err < /dev/null &

Run2PPData1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData1Vtx_v1c_ee_noZw.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 --DoZWeight false \
		 > log/GraphPPData1Vtx_v1c_ee_noZw.out 2> log/GraphPPData1Vtx_v1c_ee_noZw.err < /dev/null &

Run2MCSigGen: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSigGen_v1c_ee.root --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCSigGen_v1c_ee.out 2> log/GraphMCSigGen_v1c_ee.err < /dev/null &

Run2MCSigGen0Sub: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCSigGen0Sub_v1c_ee.root --DoGenCorrelation true  --OnlyZeroSub true --Fraction 1 \
		 > log/GraphMCSigGen0Sub_v1c_ee.out 2> log/GraphMCSigGen0Sub_v1c_ee.err < /dev/null &

Run2PPMCSigGen: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMCGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMCGen_v1c_ee.root --IgnoreCentrality true --DoGenCorrelation true  --Fraction 1 \
		 > log/GraphPPMCGen_v1c_ee.out 2> log/GraphPPMCGen_v1c_ee.err < /dev/null &

Run2PPMCSigGen0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMCGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMCGen0NPU_v1c_ee.root --IgnoreCentrality true --OnlyZeroNPU true --DoGenCorrelation true --Fraction 1 \
		 > log/GraphPPMCGen0NPU_v1c_ee.out 2> log/GraphPPMCGen0NPU_v1c_ee.err < /dev/null &

Run2DataBkg: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputDataBkg_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphDataBkg_v1c_ee_noZw.root --Fraction 1 --DoZWeight false \
		 > log/GraphDataBkg_v1c_ee_noZw.out 2> log/GraphDataBkg_v1c_ee_noZw.err < /dev/null &

Run2MCBkg: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCBkg_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkg_v1c_ee_noZw.root --Fraction 1 --DoZWeight false \
		 > log/GraphMCBkg_v1c_ee_noZw.out 2> log/GraphMCBkg_v1c_ee_noZw.err < /dev/null &

Run2MCBkgGen: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCbkgGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkgGen_v1c_ee.root --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCBkgGen_v1c_ee.out 2> log/GraphMCBkgGen_v1c_ee.err < /dev/null &

Run2MCBkgGen0Sub: Execute		 
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMCbkgGen_v1c_ee/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMCBkgGen0Sub_v1c_ee.root --OnlyZeroSub true --DoGenCorrelation true --Fraction 1 \
		 > log/GraphMCBkgGen0Sub_v1c_ee.out 2> log/GraphMCBkgGen0Sub_v1c_ee.err < /dev/null &

Run3DataSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputData3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphData3Sig_v1c.root --Fraction 1 \
		 > log/GraphData3Sig_v1c.out 2> log/GraphData3Sig_v1c.err < /dev/null &

Run3MCSig: Execute
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputMC3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphMC3Sig_v1c.root --Fraction 1 \
		 > log/GraphMC3Sig_v1c.out 2> log/GraphMC3Sig_v1c.err < /dev/null &

Run3PPMCSig: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC3_v1c.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPMC3_v1c.out 2> log/GraphPPMC3_v1c.err < /dev/null &

Run3PPMCSig0NPU: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC30NPU_v1c.root --IgnoreCentrality true --OnlyZeroNPU true --Fraction 1 \
		 > log/GraphPPMC30NPU_v1c.out 2> log/GraphPPMC30NPU_v1c.err < /dev/null &

Run3PPMCSig1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPMC3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPMC31Vtx_v1c.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPMC31Vtx_v1c.out 2> log/GraphPPMC31Vtx_v1c.err < /dev/null &

Run3PPData: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData3_v1c.root --IgnoreCentrality true --Fraction 1 \
		 > log/GraphPPData3_v1c.out 2> log/GraphPPData3_v1c.err < /dev/null &

Run3PPData1Vtx: Execute	
	./Execute --InputBase $(EOSBase)/SkimZHadron2024/OutputPPData3_v1c/ \
		--Output $(CERNBOX)/BasicPlots/ZHadron2024/GraphPPData31Vtx_v1c.root --IgnoreCentrality true --OnlyOneNVertex true --Fraction 1 \
		 > log/GraphPPData31Vtx_v1c.out 2> log/GraphPPData31Vtx_v1c.err < /dev/null &

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