MINISAT_LIB=../minisat/minisat/core/lib.a

.PHONY: all
all: sr generate

sr: sr.cc generator/Generator.cc generator/Generator.h ranklookup.h
	clang++ --stdlib=libc++ --std=c++11 -o sr -O3 -I ../minisat/minisat \
		        -I /usr/local/Cellar/boost/1.60.0_1/include \
		        -I /usr/local/Cellar/boost/1.60.0_1/lib \
				sr.cc \
	           generator/Generator.cc \
               $(MINISAT_LIB) \
			   /usr/local/Cellar/boost/1.60.0_1/lib/libboost_program_options.a \
               -Wall

generate: generator/generator_main.cc generator/Generator.cc generator/Generator.h
	clang++ --stdlib=libc++ --std=c++11 -o generate -O3 -I ../minisat/minisat \
		        -I /usr/local/Cellar/boost/1.60.0_1/include \
		        -I /usr/local/Cellar/boost/1.60.0_1/lib \
	           generator/generator_main.cc \
	           generator/Generator.cc \
               -Wall

clean:
	rm sr generate
