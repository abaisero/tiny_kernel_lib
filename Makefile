#####################################################################################################
#                               Makefile for the Tiny Kernel Library                                #
#####################################################################################################
#                                                                                                   #
# Targets you should use:                                                                           #
#   <default>   - Builds bin/usage                                                                  #
#   debug       - Builds bin/usage with debug flags                                                 #
#   optim       - Builds bin/usage with optimization flags                                          #
#   doc         - Builds documentation, creates soft link to doc/html/index.html in main directory. #
#   zip         - Compresses current state of directory in TKL.zip                                  #
#                                                                                                   #
# To view the documentation, build it and open doc.html in the main directory of the project.       #
#####################################################################################################

CXX=g++
CXXFLAGS=-Wall

SRC=src
MD=md
BIN=bin
DOC=doc

DOXYFILE=Doxyfile
DOCFILE=doc.html
INDEX=$(DOC)/html/index.html
SOURCE=$(SRC)/usage.cpp
BINARY=$(BIN)/usage
TKL=TKL

all: $(BINARY)

debug: CXXFLAGS+= -g -O0
debug: $(BINARY)

optim: CXXFLAGS+= -O2
optim: $(BINARY)

$(BINARY): $(SRC)/*
	echo Building the binary..
	mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(BINARY)

doc: $(SRC)/* $(MD)/*
	echo Building the documentation..
	mkdir -p $(DOC)
	doxygen $(DOXYFILE)
	ln -fs $(INDEX) $(DOCFILE)
	notify-send "Tiny Kernel Library Documentation" Done!


.PHONY: clean zip
zip:
	echo Zipping directory in TKL.zip..
	-rm -rf $(TKL).zip
	zip -r $(TKL) .

clean:
	echo Cleaning the directory..
	-rm -rf $(BIN) $(DOC) $(DOCFILE) $(TKL).zip

