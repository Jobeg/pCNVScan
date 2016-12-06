#define __EXTENSIONS__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ParseFile.h"
#include "uthash.h"

#include <sys/types.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>

/*Max line size for pileup file*/
#define MAX_LINE_SIZE 100000
#define MAX_CHR 250

/*============================== HASH functions =================================================*/
/*Structure of a pattern.*/
struct Pattern{
	char*pattern;		/* key : sequence of the pattern */
	float coverage;		/* coverage of the pattern (mean coverage among all genome) */
	int n;			/* number of patterns in the genome */
	UT_hash_handle hh;	/* makes this structure hashable */
}*P;

/*Struct declaration : Patterns*/
struct Pattern *patterns = NULL;    /* important! initialize to NULL */


void addPattern(char*sequence, float coverage){
	struct Pattern *p;
	const char*key=sequence;

	HASH_FIND_STR(patterns, key, p);  /* id already in the hash? */
	if (p==NULL) {
		p = (struct Pattern*)malloc(sizeof(struct Pattern));
		p->pattern = sequence;
		HASH_ADD_KEYPTR(hh, patterns, p->pattern, strlen(p->pattern), p );
		p->n = 1;
		p->coverage = coverage;
	}
	else{
		float Xn_1 = findPatternCoverage(sequence);
		p->n += 1;
		/* Average formula : _Xn= _Xn-1+(Xn-_Xn-1)/n */
		p->coverage = Xn_1 + (coverage-Xn_1)/p->n;
	}
}

float findPatternCoverage(char*sequence) {
	float coverage=0;
	const char*name=sequence;
	struct Pattern *p;
	HASH_FIND_STR( patterns, name, p );
	coverage=p->coverage;
	return coverage;
}

int findPatternN(char*sequence) {
	int n=0;
	const char*name=sequence;
	struct Pattern *p;
	HASH_FIND_STR( patterns, name, p );
	n=p->n;
	return n;
}

/*============================================ Getters and setters for structures ==================================*/
/*### Chromosome ######*/
char* getChrName(C c){
/* Pre-condition : c is non-empty */
	return c->name;
}

char* getChrSeq(C c){
/* Pre-condition : c is non-empty */
	return c->seq;
}

int getChrLenght(C c){
/* Pre-condition : c is non-empty */
	return c->length;
}
int*getChrCoverage(C c){
/* Pre-condition : c is non-empty */
	return c->coverage;
}

G getChrGenes(C c){
/* Pre-condition : c is non-empty */
	return c->genes;
}

void setChrGenes(C c, G genes){
	(*c).genes=genes;
}

C newChromosome(char*name, char*seq, int length){
	/*New element*/
	C newChromosome=NULL;
	int*coverage=calloc(length, sizeof(int));
	newChromosome = (C)malloc(sizeof(struct Chromosome));
	(*newChromosome).name=name;
	(*newChromosome).seq=seq;
	(*newChromosome).length=length;
	(*newChromosome).coverage=coverage;
	(*newChromosome).genes=NULL;
	return newChromosome;
}

/*##### Gene ######*/
char* getGenName(G g){
/* Pre-condition : g is non-empty */
	return g->name;
}

int getGenStart(G g){
/* Pre-condition : g is non-empty */
	return g->start;
}

int getGenEnd(G g){
/* Pre-condition : g is non-empty */
	return g->end;
}

G getNextGene(G g){
/* Pre-condition : g is non-empty */
	return g->next;
}

void setNextGene(G g, G next){
	(*g).next=next;
}

G newGene(char*name, int start, int end, G next){
	/*New element*/
	G newGene=NULL;
	newGene = (G)malloc(sizeof(struct Gene));
	(*newGene).name=name;
	(*newGene).start=start;
	(*newGene).end=end;
	setNextGene(newGene, next);
	return newGene;
}

/*======================================== Functions ====================================================*/
void chomp(const char *s){
	char *p;
	while (NULL != s && NULL != (p = strrchr(s, '\n'))){
		*p = '\0';
	}
}

/* Extract a substring in regExp */
char* regex(const char *regexp, char* string){
	char*subString=NULL;
	regex_t preg;
	int err;
			
	err = regcomp(&preg, regexp, REG_EXTENDED);
	if(err == 0){
		int match;
		size_t nmatch = 0;
		regmatch_t *pmatch=NULL;
		nmatch = 2;
		pmatch = malloc (sizeof (*pmatch) * nmatch);
		if (pmatch){
			match = regexec (&preg, string, nmatch, pmatch, 0);
			regfree (&preg);
			if (match == 0){
				char *site = NULL;
				int start = pmatch[1].rm_so;
				int end = pmatch[1].rm_eo;
				size_t size = end - start;
				site = malloc (sizeof (*site) * (size + 1));
				if (site){
					strncpy (site, &string[start], size);
					site[size] = '\0';
					subString=site;
				}
			}
		}
	}
	return subString;
}

/*static char *getQualValue(char *sQualifier, gb_feature *ptFeature) {
    static char *null = "";
    gb_qualifier *i;

    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++) 
        if (strcmp(sQualifier, i->sQualifier) == 0)
            return i->sValue;
    return null;
}*/

/*Extract informations from a gbk file.
Return chromosome sequence and all genes features.*/
/*C* ExtractAllFromGbk(char* gbk){
	C*reference=NULL;
*/	
	/* USE GBFP PROGRAM ON A GPL LICENSEv2 (http://code.google.com/p/gbfp/) */
/*	gb_data **pptSeqData, *ptSeqData;
	gb_feature *ptFeature;

	char *sQualifier = "locus_tag";
	char *sQualifier2 = "product";
	char *sFeature = "CDS";
	G currentGene=NULL;

	unsigned int i, j;
*/
/*	pptSeqData = parseGBFF(gbk);*/ /* parse a GBF file which contains more than one GBF sequence data */
	/*for (i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {*/ /* ptSeqData points a parsed data of a GBF sequence data */
		/* sDef point to the chromosome definition */
/*		C currentChromosome=newChromosome(ptSeqData->sDef, ptSeqData->sSequence, ptSeqData->lLength);
		currentGene=NULL;*/
		/* start of user process */
/*		for (j = 0; j < ptSeqData->iFeatureNum; j++) {
			ptFeature = (ptSeqData->ptFeatures + j);
			if (strcmp(sFeature, ptFeature->sFeature) == 0) {
				char*geneName;
				geneName=malloc(400*sizeof(char));*/
				/*printf(">%li-%li %s(%s)\n", \
				ptFeature->lStart, \
				ptFeature->lEnd, \
				getQualValue(sQualifier2, ptFeature), \
				getQualValue(sQualifier, ptFeature));*/
/*				sprintf(geneName, "%s(%s)", getQualValue(sQualifier2, ptFeature), getQualValue(sQualifier, ptFeature));
				currentGene=newGene(geneName, ptFeature->lStart, ptFeature->lEnd, currentGene);
				setChrGenes(currentChromosome,currentGene);
			}
		}
		reference=(C*)realloc(reference,(i+1)*sizeof(C));
		reference[i]=currentChromosome;*/
		/* end of user process */
	/*}
	reference=(C*)realloc(reference,(i+2)*sizeof(C));
	reference[i+1]=NULL;*/
	/* freeGBData(pptSeqData);  release memory space */
/*	return reference;
}*/

/*Extract informations from a gff AND a fasta file.
Return chromosome sequence and all genes features.*/
C* ExtractAllFromGff(char* gff, char* fasta){
	C*reference=NULL;
	char ligne[MAX_LINE_SIZE];
	char*chrName=malloc(400*sizeof(char));
	char*chrSeq=malloc(1*sizeof(char));
	int chrSize=0;
	int c;
	C currentChromosome;
	G currentGene;/* Initialize the list of genes for each chromosome */
	int nbChr=0;/* Number of chromosomes */
	FILE*fastaFile;
	FILE*gffFile;
	/* Read fasta file first to extract chromosomes */
	if((fastaFile=fopen(fasta, "r"))==NULL){
		fprintf(stderr, "Can not open file fasta!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(ligne, MAX_LINE_SIZE, fastaFile) != NULL){
		char name[100];
		chomp(ligne);
		sscanf(ligne, ">%s", name);
		if(name != NULL){
			if(strcmp(name,chrName)!=0){
				/* Record Chromosome */
				if(chrSize>0){
					currentChromosome=newChromosome(chrName, chrSeq, chrSize);
					reference=(C*)realloc(reference,(nbChr+1)*sizeof(C));
					reference[nbChr]=currentChromosome;
					nbChr++;
				}				
				/* Start a new one */
				chrName=malloc(400*sizeof(char));
				strcpy(chrName,name);
				chrSeq=malloc(1*sizeof(char));
				chrSize=0;
				continue;
			}
			else{
				int len=strlen(ligne);
				chrSize+=len;
				chrSeq =(char*)realloc(chrSeq, chrSize*(sizeof(char)+1));
				strcat(chrSeq, ligne);
			}
		}
	}
	/* Last chromosome */
	currentChromosome=newChromosome(chrName, chrSeq, chrSize);
	reference=(C*)realloc(reference,(nbChr+1)*sizeof(C));
	reference[nbChr]=currentChromosome;
	nbChr++;
	
	/* Set the last element to NULL */
	reference[nbChr]=NULL;


	/* Read gff file to extract genes */
	if((gffFile=fopen(gff, "r"))==NULL){
		fprintf(stderr, "Can not open file gff!\n");
		exit(EXIT_FAILURE);
	}
	/* Reset chromosome name */
	chrName=malloc(400*sizeof(char));
	strcpy(chrName, "NULL");
	currentGene=NULL;
	while(fgets(ligne, MAX_LINE_SIZE, gffFile) != NULL){
		/* Fields declaration */
		char chromosome[100];
		char source[100];
		char feature[100];
		int start;
		int end;
		char score[10];
		char strand[2];
		char phase[2];
		char attributes[1000];
		chomp(ligne);
		/* Search genes lines */
		sscanf(ligne, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s", chromosome, source, feature, &start, &end, score, strand, phase, attributes);
		if(strcmp(feature,"gene")==0){
			char *name;
			char *desc;
			char*geneName;
			/* Parse the attribute field */
			/* Extract Name and description */
			char delims[] = ";";
			char *result = NULL;
			geneName=malloc(400*sizeof(char));

			result = strtok( ligne, delims );
			while( result != NULL ){
				const char *str_regex = "Name=(.+)";
				char*subName=regex(str_regex, result);
				char*subDesc;
				if(subName!=NULL){
					name=subName;
				}
				str_regex = "description=(.+)";
				subDesc=regex(str_regex, result);
				if(subDesc!=NULL){
					desc=subDesc;
				}
				result = strtok( NULL, delims );
			}
			free(result);
			sprintf(geneName, "%s(%s)", desc, name);

			/* If this gene is on a new chromosome */
			if(strcmp(chromosome,chrName)!=0){
				if(strcmp(chrName,"NULL")!=0){
					/* Find the good chromosome */
					for(c=0; c<nbChr; c++){
						if(strcmp(getChrName(reference[c]),chrName)==0){
							setChrGenes(reference[c],currentGene);
							break;
						}
					}
				}
				strcpy(chrName, chromosome);
				/* Reset list of genes */
				currentGene=NULL;
				currentGene=newGene(geneName, start, end, currentGene);
			}
			else{
				currentGene=newGene(geneName, start, end, currentGene);
			}
		}
	}
	for(c=0; c<nbChr; c++){
		if(strcmp(getChrName(reference[c]),chrName)==0){
			setChrGenes(reference[c],currentGene);
			break;
		}
	}
	return reference;
}

/*Extract coverage from a pileup file.
Compute the value of int*coverage in structure Chromosome.*/
C* ExtractCoverageFromPileup(C*reference, char* pileup){
	C Chromosome=reference[0];
	FILE*file;
	char ligne[MAX_LINE_SIZE];
	char chromosome[400];/*Chromosome name*/
	int*Coverage=getChrCoverage(Chromosome);
	printf("Chromosome name: %s\n", getChrName(Chromosome));
	if((file=fopen(pileup, "r"))==NULL){
		fprintf(stderr, "Can not open file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(ligne, MAX_LINE_SIZE, file) != NULL){
		char delims[] = "\t";
		char *result = NULL;
		int i=0;/*iterator for lines*/
		int pos;
		result = strtok( ligne, delims );
		while( result != NULL ) {
			if(i==0){/*Chromosome*/
				int resultat = strcmp(chromosome, result);
				/*Change Chromosome*/
				if(resultat!=0){
					int c=0;
					int isFound=0;
					Chromosome=reference[c];
					strcpy(chromosome, result);
					printf("Chromosome name: %s\n", chromosome);
					while(Chromosome!=NULL){
						/*printf("Essai chromosome name: %s (%p)\n", getChrName(Chromosome), getChrName(Chromosome));*/
						if(strcmp(getChrName(Chromosome),chromosome)==0){
							Chromosome=reference[c];
							Coverage=getChrCoverage(Chromosome);
							isFound++;
							/*printf("Trouvé le bon chromosome: %s (%p)\n", getChrName(Chromosome), getChrName(Chromosome));*/
							break;
						}
						c++;
						Chromosome=reference[c];
					}
					if(isFound==0){
						fprintf(stderr, "Error: chromosome name not found in reference: %s\n", chromosome);
					}
				}
			}
			else if(i==1){/*Position*/
				sscanf( result, "%d", &pos );
				if(pos>getChrLenght(Chromosome)){
					printf("Pos %d outRange for chromosome %s! Perhaps you take the wrong reference file?\n", pos, getChrName(Chromosome));
					exit(1);
				}
			}
			else if(i==3){/*Coverage*/
				int cov;
				sscanf( result, "%d", &cov );
				Coverage[pos]=cov;
			}
			i++;
			result = strtok( NULL, delims );
		}
		free(result);
	}
	return reference;
}

void ComputeRatios(C*reference, int patternSize, char*outFile){
	/*Compute the Hash of patterns (mean coverage of patterns) among ALL genome/genes*/
	FILE*file;
	int i=0;
	C Chromosome=reference[i];
	/* Output file */
	if((file=fopen(outFile, "w"))==NULL){
		fprintf(stderr, "Error in creating output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(file,"Chromosome\tGene\tStart\tEnd\tRatio\n");/* Header */
	/* For each chromosome */
	while(Chromosome!=NULL){
		G genes=getChrGenes(Chromosome);
		char*seq=getChrSeq(Chromosome);
		int*cov=getChrCoverage(Chromosome);
		/* For each gene */
		while(genes!=NULL){
			int start=getGenStart(genes);
			int end=getGenEnd(genes);
			int c;
			/*printf("%s\t%d_%d\n", getGenName(genes), getGenStart(genes), getGenEnd(genes));*/
			/* For each nt on Gene */
			for(c=start-1; c<end; c++){
				char*pattern=malloc(patternSize*sizeof(char));
				float coverage=0;/* Mean coverage for one pattern*/
				int y;
				for(y=0; y<patternSize; y++){
					pattern[y]=seq[c+y];
					coverage+=cov[c+y];
				}
				coverage=coverage/patternSize;
				addPattern(pattern, coverage);
			}
			genes=getNextGene(genes);
		}
		i++;
		Chromosome=reference[i];
	}
	/* Compute the ratio among all genes */
	i=0;
	Chromosome=reference[i];/* For each chromosome */
	while(Chromosome!=NULL){
		G genes=getChrGenes(Chromosome);
		char*seq=getChrSeq(Chromosome);
		int*cov=getChrCoverage(Chromosome);
		/* For each gene */
		while(genes!=NULL){
			int start=getGenStart(genes);
			int end=getGenEnd(genes);
			float meanRatio=0;/* Average Ratio on this gene */
			int c;
			
			/* For each nt on Gene */
			for(c=start-1; c<end; c++){
				char*pattern=malloc(patternSize*sizeof(char));
				float coverage=0;/* Mean coverage for one pattern*/
				float theoricCoverage=0;
				float ratio=0;
				int y;
				for(y=0; y<patternSize; y++){
					pattern[y]=seq[c+y];
					coverage+=cov[c+y];
				}
				coverage=coverage/patternSize;
				theoricCoverage=findPatternCoverage(pattern);
				ratio=coverage/theoricCoverage;
				meanRatio+=ratio;
			}
			/* Average Ratio */
			meanRatio=meanRatio/(end-start);
			fprintf(file,"%s\t%s\t%d\t%d\t%f\n", getChrName(Chromosome), getGenName(genes), getGenStart(genes), getGenEnd(genes), meanRatio);
			genes=getNextGene(genes);
		}
		i++;
		Chromosome=reference[i];
	}
	fclose(file);
}


