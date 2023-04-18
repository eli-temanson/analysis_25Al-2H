CC=g++
EXE=gwm_test
CXX_FLAGS= -std=c++17 -g -Wall

CATIMA_PATH=../build/
LIB_PATH=$(CATIMA_PATH)lib/
INCLUDE_PATH=$(CATIMA_PATH)include/
CPP=gwm_test.cpp

.PHONY: all clean

all: $(EXE)

$(EXE): $(CPP)
	$(CC) $(CXX_FLAGS) -I$(INCLUDE_PATH) $^ $(LIB_PATH)libcatima.a -o $@

clean:
	$(RM) $(EXE)
