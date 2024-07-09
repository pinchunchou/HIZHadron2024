

default: TestRun

TestRun: Execute
	./Execute

TestRunSystematics: ExecuteSystematics
	./ExecuteSystematics

Execute: MakePlot.cpp
	g++ MakePlot.cpp  -I. -I$(ProjectBase)/CommonCode/include \
		-o Execute `root-config --cflags --libs` -std=c++17 \
		$(ProjectBase)/CommonCode/library/Messenger.o \

ExecuteSystematics: DeriveSystematics.cpp
	CompileRootMacro DeriveSystematics.cpp ExecuteSystematics
