MINISAT_LIB=../minisat/minisat/core/lib.a

sr: sr.cc generator/generator.cpp ranklookup.h
	g++ --std=c++11 -o sr -O3 -I ../minisat/minisat sr.cc \
	           generator/generator.cpp \
               $(MINISAT_LIB) \
	           -lboost_program_options \
               -Wall

sr_debug: sr.cc generator/generator.cpp ranklookup.h
	g++ --std=c++11 -o sr_debug -O0 -I ../minisat/minisat sr.cc \
	           generator/generator.cpp \
               $(MINISAT_LIB) \
	           -lboost_program_options \
			   -g -Wall

sr_gprof: sr.cc generator/generator.cpp ranklookup.h
	g++ --std=c++11 -o sr_gprof -O0 -I ../minisat/minisat sr.cc \
	           generator/generator.cpp \
               $(MINISAT_LIB) \
	           -lboost_program_options \
			   -pg

clean:
	rm sr sr_debug sr_gprof
