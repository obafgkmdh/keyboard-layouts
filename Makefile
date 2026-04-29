CXX=g++

CXXFLAGS=-O3

all: generate

generate: generate.cpp
	${CXX} ${CXXFLAGS} -o $@ $<

clean:
	rm -f generate
