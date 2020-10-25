all: hwbs

hwbs: hwbs.c
	gcc -Wall -o hwbs hwbs.c -lm -O3 -mcmodel=medium

clean: 
	rm -fr *~ hwbs