ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin


all: mkdirBin bin/Z_EE_Channel.exe bin/Z_mumu_Channel.exe bin/massPeakPlots.exe bin/v2Plots.exe bin/Z_mumu_EfficienciesMC.exe bin/Z_EE_EfficienciesMC.exe bin/MCReweighting.exe bin/efficiencyPlots.exe bin/massPeakPlots_BkgSub.exe bin/prettyPlots_pTy.exe bin/prettyPlots_yields.exe

mkdirBin:
	$(MKDIR_BIN)


bin/Z_EE_Channel.exe: src/Z_EE_Channel.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/Z_EE_Channel.exe src/Z_EE_Channel.C

bin/Z_mumu_Channel.exe: src/Z_mumu_Channel.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/Z_mumu_Channel.exe src/Z_mumu_Channel.C

bin/Z_mumu_EfficienciesMC.exe: src/Z_mumu_EfficienciesMC.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/Z_mumu_EfficienciesMC.exe src/Z_mumu_EfficienciesMC.C

bin/Z_EE_EfficienciesMC.exe: src/Z_EE_EfficienciesMC.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/Z_EE_EfficienciesMC.exe src/Z_EE_EfficienciesMC.C

bin/massPeakPlots.exe: src/massPeakPlots.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/massPeakPlots.exe src/massPeakPlots.C

bin/massPeakPlots_BkgSub.exe: src/massPeakPlots_BkgSub.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/massPeakPlots_BkgSub.exe src/massPeakPlots_BkgSub.C

bin/v2Plots.exe: src/v2Plots.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/v2Plots.exe src/v2Plots.C

bin/MCReweighting.exe: src/MCReweighting.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/MCReweighting.exe src/MCReweighting.C

bin/efficiencyPlots.exe: src/efficiencyPlots.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/efficiencyPlots.exe src/efficiencyPlots.C

bin/prettyPlots_pTy.exe: src/prettyPlots_pTy.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/prettyPlots_pTy.exe src/prettyPlots_pTy.C

bin/prettyPlots_yields.exe: src/prettyPlots_yields.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/prettyPlots_yields.exe src/prettyPlots_yields.C

clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
