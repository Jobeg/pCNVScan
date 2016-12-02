#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "ParseFile.h"

/* Pattern size when not specified by user in options */
#define DEFAULT_PATTERN_SIZE 6

/*=========================================================================================================================
 CNVScan: Algorithm for searching Copy Number Variations in a sample file (*.pileup) among a list of genes in a ref file (*.gbk).
 REQUIRED PROGRAM(S): gbfp (http://code.google.com/p/gbfp/) a GenBank flatfile parser library with high speed
 REQUIRED OPTIONS:
 Pattern size. Size of the patterns for algorithms. Recommended number: 6. Too High number will decrease program's speed. Too low number reduce the sensitivity.
 Pileup sample file. Path to the input Pileup file.
 GeneBank reference file. Path to the input reference file.
 Version: 1.2.1
 Author: Johann Beghain
 Contact: johann.beghain@pasteur.fr
=========================================================================================================================*/

/*============================= Functions ================================*/
void printHelp(char*prog){
	fprintf(stderr, "CNVScan: Algorithm for searching Copy Number Variations in a sample file (*.pileup) among a list of genes in a ref file (A gff+ fasta file OR a gbk file)\n");
	fprintf(stderr, "OPTIONS:\n"
	"-p String [REQUIRED]. Path to a valid .pileup file. The alignment must have been done on the same reference as the reference file.\n"
	"-o String [REQUIRED]. Path to the output file. Output a list of copy numbers by genes.\n"
	"-i String [REQUIRED1]. Path to a valid input .gbk reference file. No fasta file required.\n"
	"-g String [REQUIRED2]. Path to a valid .gff reference file. If set you have to give also a fasta file in -f option.\n");
	fprintf(stderr, "-f String [REQUIRED2]. Path to a valid .fasta file. Should be the same reference as the .gff file.\n"
	"-s Integer [OPTIONAL]. Pattern size for iterative search. Default: 6.\n"
	"-h Print this help and exit.\n");
	fprintf(stderr, "Usage: %s -s [int Pattern size] -i [.gbk reference file] -p [Pileup sample file] -o [Output file]\n", prog);
}

/*============================= Main program =============================*/

int main(int argc, char *argv[]){
	int option = 0;
	char* gbkFile = NULL;
	char* pileupFile = NULL;
	char* gffFile = NULL;
	char* fastaFile = NULL;
	char* outFile = NULL;
	int patternSize = -1;
	C*reference=NULL;
	char*patternArg;
	while ((option = getopt(argc, argv,"hi:o:p:s:g:f:")) != -1) {
		switch (option) {
			case 'h' :
				printHelp(argv[0]);
				exit(0);
				break;
			case 'i' : gbkFile = optarg;/* Input .GBK file */
				break;
			case 'o' : outFile = optarg;/* Output file */
				break;
			case 'p' : pileupFile = optarg;/* Input .pileup file */
				break;
			case 's' : patternArg = optarg;/* Pattern size */
				patternSize = atoi(patternArg);
				break;
			case 'g' : gffFile = optarg;/* Input .gff file */
				break;
			case 'f' : fastaFile = optarg;/* Input .fasta file */
				break;
			default: fprintf(stderr, "Usage: %s -s [int Pattern size] -g [.gbk reference file] -p [Pileup sample file] -o [Output file]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	if(pileupFile==NULL || outFile==NULL){
		fprintf(stderr, "\nError: missing one or more required arguments.\n");
		printHelp(argv[0]);
		exit(EXIT_FAILURE);
	}
	if(patternSize==-1){/* Set default options */
		printf("No pattern size specified. Default: 6.\n");
		patternSize=DEFAULT_PATTERN_SIZE;
	}
	if(patternSize==0){
		fprintf(stderr, "Error: too short Patern size (%d) or wrong argument (%s): must be an integer\n", patternSize, patternArg);
		exit(EXIT_FAILURE);
	}
	/* Input format check */
	if(gbkFile==NULL){
		if(gffFile==NULL && fastaFile==NULL){
			fprintf(stderr, "\nError: missing or bad input.\n2 Options:\n"
			"-i <GBK file>"
			"OR -g <GFF file> -f <FASTA file> (GFF and FASTA should be in the same reference assembly.");
			printHelp(argv[0]);
			exit(EXIT_FAILURE);
		}
		else if(gffFile==NULL || fastaFile==NULL){
			fprintf(stderr, "\nError: missing one option: -f or -g\n");
			printHelp(argv[0]);
			exit(EXIT_FAILURE);
		}
		reference=ExtractAllFromGff(gffFile, fastaFile);
		/*int i;
		for(i=0; i<16; i++){
			C Chromosome=reference[i];
			printf("Chromosome name:%s - (%d)\n", getChrName(Chromosome), getChrLenght(Chromosome));
			G genes = getChrGenes(Chromosome);
			while(genes!=NULL){
				char*name=getGenName(genes);
				printf("Gene name: %s\t(%d - %d)\n", name, getGenStart(genes), getGenEnd(genes));
				genes=getNextGene(genes);
			}
		}*/
	}
	else{
		/*Extract reference informations from gbk file
		Each chromosome is on a special structure*/
		reference=ExtractAllFromGbk(gbkFile);
	}

	/*Extract the coverage for every position in each chromosome*/
	reference=ExtractCoverageFromPileup(reference, pileupFile);
	
	/*Compute the ratio for each genes. Output the result file*/
	ComputeRatios(reference, patternSize, outFile);

	return 0;
}
