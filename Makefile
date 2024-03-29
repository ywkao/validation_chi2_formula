CC       := g++ # This is the main compiler
BUILDDIR := build
INC 	 := -I include
CFLAGS   := $(shell root-config --cflags) -g -O3 #-Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL #-Wall -Wextra
ROOTLIBS := $(shell root-config --libs) -lMinuit -lMLP -lXMLIO -lTMVA -lGenVector

TARGET   := bin/main

all: ${TARGET}

build/chi2_helper.o: chi2_helper.cc chi2_helper.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/main.o: main.cc chi2_helper.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

bin/main: build/main.o build/chi2_helper.o
	@echo "Linking..." && mkdir -p bin
	$(CC) -o $@ $^ $(CFLAGS) $(ROOTLIBS)

.PHONY: clean

clean:
	@echo " Cleaning..."; 
	$(RM) -r $(BUILDDIR) $(TARGET)
