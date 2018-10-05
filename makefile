# Custom makefile for Interface
all: TestInterface

TestInterface: TestInterface.cpp interface.cpp interface.h
	g++ -std=c++0x -o TestInterface TestInterface.cpp interface.cpp
	# 	g++ -std=c++0x -lboost_iostreams -lboost_system -lboost_filesystem -o TestInterface TestInterface.cpp interface.cpp

clean:
	rm -f TestInterface core
