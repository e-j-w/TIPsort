CXX=g++
CXXFLAGS:= -O -Wall -fPIC $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
INC = -I$(ROOTINC)
LIBS:= -lm -ldl
OBJECTS:=main.cxx MainWindow.cxx RootManager.cxx GMMsgBox.cxx GMTextEntry.cxx
INCLUDES:=MainWindow.h RootManager.h GMMsgBox.h GMTextEntry.h
BIN = $(TIPsort)/bin

GateMaker: $(OBJECTS) $(INCLUDES)
	$(CXX) -o $(BIN)/GateMaker $(OBJECTS) $(CXXFLAGS) $(INC) $(LIBS) $(ROOTLIBS) 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< $(CXXFLAGS) $(INC)

clean: 
	rm -f *.o *Dict* *~
