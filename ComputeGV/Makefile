#CXX = clang++
#CXXFLAGS = -std=c++17 -march=native -O3 -fPIC
CPPFLAGS += -I`brew --prefix`/include
CXXFLAGS = -std=c++17 -O3 -fPIC #-mcpu=apple-m1
LDFLAGS  = -L`brew --prefix`/lib -lmpfr -lgmp -pthread
LOADLIBS = -L`brew --prefix`/lib -lmpfr -lgmp -pthread
#CPPFLAGS = -L/usr/local/include

computeGV: computeGV.o
	c++ -o computeGV ${LOADLIBS} computeGV.o

clean:
	rm *.o computeGV
