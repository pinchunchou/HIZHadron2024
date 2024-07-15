default: TestRun

TestRun: Execute
	./Execute

Execute: RecoLevelSumHFExploration.cpp
	g++ RecoLevelSumHFExploration.cpp -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17