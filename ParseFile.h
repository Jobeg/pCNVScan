typedef
struct Gene{
	char*name;/*Chromosome Name*/
	int start;/*Gene start position on Chromosome*/
	int end;/*Gene end position*/
	struct Gene*next;
}*G;
typedef
struct Chromosome{
	char*name;/*Name of Chromosome*/
	char*seq;/*Sequence of Chromosome*/
	int*coverage;/*Coverage among this chromosome*/
	int length;/*Size of Chromosome*/
	struct Gene*genes;/*List of genes among this chromosome, points to the first gene*/
}*C;

void addPattern(char*sequence, float coverage);
float findPatternCoverage(char*sequence);
int findPatternN(char*sequence);
char* getChrName(C c);
char* getChrSeq(C c);
int getChrLenght(C c);
int*getChrCoverage(C c);
G getChrGenes(C c);
void setChrGenes(C c, G genes);
C newChromosome(char*name, char*seq, int length);
char* getGenName(G g);
int getGenStart(G g);
int getGenEnd(G g);
G getNextGene(G g);
G newGene(char*name, int start, int end, G next);
void chomp(const char *s);
char* regex(const char *regexp, char* string);
C* ExtractAllFromGbk(char* gbk);
C* ExtractAllFromGff(char* gff, char* fasta);
C* ExtractCoverageFromPileup(C*reference, char* pileup);
void ComputeRatios(C*reference, int patternSize, char*outFile);
