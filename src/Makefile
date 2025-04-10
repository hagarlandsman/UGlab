# Compiler and flags
CXX = g++
CXXFLAGS = -pthread -std=c++17 -I./vcpkg/installed/x64-linux/include `root-config --cflags`
LDFLAGS = `root-config --libs`

# Main executable
SOURCES = cosmic_muons.c LinkDef_rdict.cxx UGburden.cpp
TARGET = muon

# Secondary test executable
UGTEST_SRC = test_ug.cpp UGburden.cpp
UGTEST_TARGET = test_ug

# Default build rule
all: $(TARGET)

all: muon UGburden main


$(TARGET): $(SOURCES)
	$(CXX) -o $@ $(SOURCES) $(CXXFLAGS) $(LDFLAGS)

$(UGTEST_TARGET): $(UGTEST_SRC)
	$(CXX) -o $@ $(UGTEST_SRC) $(CXXFLAGS) $(LDFLAGS)


UGburden: main.cpp LinkDef_rdict.cxx UGburden.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)


main: main.cpp LinkDef_rdict.cxx UGburden.cpp LinkDef_rdict.cxx
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

knn_search: read_geometry.c
	$(CXX) -o $@ $^ -Ivcpkg/installed/x64-linux/include -pthread `root-config --cflags --libs`


get_latlon: get_latlon.cpp latLon_To_NorthingEasting.cpp
	$(CXX) -o $@ $^ -Ivcpkg/installed/x64-linux/include  -I. -pthread `root-config --cflags --libs`


LinkDef_rdict.cxx: UGburden.h LinkDef.h
	rootcling -f $@ -c UGburden.h LinkDef.h


clean:

	rm -f $(TARGET) $(UGTEST_TARGET)
