#CXX = clang++
#CXXFLAGS = -std=c++17 -march=native -O3 -fPIC
CPPFLAGS += -I`brew --prefix`/include # -I/opt/homebrew/include
CXXFLAGS = -std=c++17 -O3 -fPIC #-mcpu=apple-m1
LDFLAGS = -L/opt/homebrew/lib -lmpfr -lgmp -pthread 
LOADLIBS = -L/opt/homebrew/lib -lmpfr -lgmp -pthread
#CPPFLAGS = -L/usr/local/include

computeGV: computeGV.o
	c++ -o computeGV ${LOADLIBS} computeGV.o

clean:
	rm *.o computeGV
