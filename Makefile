prepare: prepare.cpp Makefile prepare.hpp
	g++-15 -o prepare prepare.cpp -Wall -std=c++20 -ggdb3
