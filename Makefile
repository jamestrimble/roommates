MINISAT_LIB=../minisat/minisat/core/lib.a

.PHONY: all
all: sr generate

sr: sr.cc generator/Generator.cc generator/Generator.h ranklookup.h
	g++ --std=c++11 -o sr -O3 -I ../minisat/minisat \
               sr.cc \
	           generator/Generator.cc \
	           -lboost_program_options \
               $(MINISAT_LIB) \
               -Wall

generate: generator/generator_main.cc generator/Generator.cc generator/Generator.h
	g++ --std=c++11 -o generate -O3 -I ../minisat/minisat \
	           -lboost_program_options \
	           generator/generator_main.cc \
	           generator/Generator.cc \
               -Wall

clean:
	rm sr generate
