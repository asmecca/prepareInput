prepare: prepare.cpp Makefile prepare.hpp
	g++ -o prepare prepare.cpp -Wall -std=c++20 -ggdb3
