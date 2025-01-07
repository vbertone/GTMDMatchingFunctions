CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++20

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# PARTONS
PARTONSINCS = -I/usr/local/include/
PARTONSLIBS = -L/usr/local/lib/ -lElementaryUtils -lNumA++ -lPARTONS

# NangaParbat
NANGAPARBATINCS = $(shell NangaParbat-config --cppflags)
NANGAPARBATLIBS = $(shell NangaParbat-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(APFELPPINCS) $(PARTONSINCS) $(NANGAPARBATINCS)
LDFLAGS  += $(APFELPPINCS) $(PARTONSINCS) $(NANGAPARBATINCS)

CLIBS += $(APFELPPLIBS) $(PARTONSLIBS) $(NANGAPARBATLIBS)

install : all
all : EvolutionCheck GTMDMatchingkTxi GTMDMatchingkTQ GTMDMatchingxxi GTMDMatchingxxi_analytic
full: GTMDMatchingkTxi_full

EvolutionCheck: EvolutionCheck.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GTMDMatchingkTxi: GTMDMatchingkTxi.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GTMDMatchingkTxi_full: GTMDMatchingkTxi_full.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GTMDMatchingkTQ: GTMDMatchingkTQ.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GTMDMatchingxxi: GTMDMatchingxxi.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GTMDMatchingxxi_analytic: GTMDMatchingxxi_analytic.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la EvolutionCheck GTMDMatchingkTxi GTMDMatchingkTQ GTMDMatchingxxi GTMDMatchingxxi_analytic GTMDMatchingkTxi_full *~
