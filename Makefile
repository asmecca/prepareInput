# prepare: prepare.cpp Makefile prepare.hpp
# 	g++ -o prepare prepare.cpp -Wall -std=c++20 -ggdb3

pipi: pipi.cpp Makefile
	g++ -o pipi pipi.cpp -Wall -std=c++20 -O3
