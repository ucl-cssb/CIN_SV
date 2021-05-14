# CCC = icc		# Using intel complier
CCC = g++
# CCC = /usr/local/opt/llvm/bin/clang
BOOST = /usr/local
FLAG = -O3
# FLAG =

OUT_DIR = ../bin

MKDIR_P = mkdir -p

.PHONY: directories


all: directories simsv


directories: ${OUT_DIR}


${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}


simsv:
	$(CCC) $(FLAG) -std=gnu++11 sim_sv.cpp -o ${OUT_DIR}/simsv -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -lgsl -lgslcblas -lm -lz -I$(BOOST)/include


clean:
	rm simsv
