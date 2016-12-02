CC = gcc
CFLAGS = -O2 -Wall -ansi -pedantic


all: cnvscan

cnvscan: CNVScan.c
	$(CC) $(CFLAGS) -o CNVScan CNVScan.c ParseFile.c -lgbfp


