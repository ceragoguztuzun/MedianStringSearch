#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <limits.h>

#define MAX_K 16
#define MAX_N 500
#define MAX_T 100

// global variables
uint64_t dna_matrix[MAX_T*MAX_N]; //stores dna sequences
int positions; //no of positions
uint64_t s; //best patten to be found

// get (x,y)th nucleotide of dna matrix
uint64_t getNucleotide_seq(int x, int y)
{
	return (dna_matrix[x*positions + (y >> 4)] >> (2*(y & 0x0F))) & 0x03;
}

// returns hamming distance between two numbers with bitwise operations
int hammingDist(int a, int b, int k)
{
    int hammingDistance = 0;
    for( int i = 0; i < k; i++)
    {
        if( ((s >> (2*i)) & 0x03) ^ (getNucleotide_seq(a, (b+i))) ) hammingDistance++;
    }
    return hammingDistance;
}

// returns sum of hamming distances for all DNA sequences
int getTotalDistance(int t, int k, int n)
{
	int hamming_total = 0;
	int min_distance_yet = INT_MAX;
	int curr_dist;

    for(int i = 0; i < t; i++)
    {
        min_distance_yet = INT_MAX;
        for(int j = 0; j < n - k + 1; j++)
        {
            curr_dist = hammingDist(i, j, k);
            if(curr_dist < min_distance_yet) min_distance_yet = curr_dist;
        }
        hamming_total += min_distance_yet;
    }

    return hamming_total;
}

// converts encoded form to nucleotides char by char
char getChar(int x) 
{
    char nucleotide;
    if (x == 0) nucleotide = 'A';
    else if (x == 1) nucleotide = 'T';
    else if (x == 2) nucleotide = 'G';
    else if (x == 3) nucleotide = 'C';

    return nucleotide;
}

int main(int argc, char *argv[])
{
 	printf ("MEDIAN STRING SEARCH\n\n"); 

 	// variables
 	int t = 0; // number of DNA sequences
 	int n = 0; // length of each DNA sequence
 	int k; // length of motif sequence

 	FILE *fp_inputfile; 
 	FILE *fp_outputfile;
 	char nucleotide;
 	int i, j;

  	//check argument number
	if( argc != 2)
	{
		printf("User must enter correct number of arguments\n");
		return 0;
	}

	// read file
	fp_inputfile = fopen(argv[1],"r");

	//check if file ptrs are null
	if( fp_inputfile == NULL)
	{
		printf("Input file is NULL\n");
		return 0;
	}

	// get K value
	fscanf(fp_inputfile, "%d\n", &k);
	printf("k: %d\n", k);
	if (k > MAX_K)
	{
		printf("K length exceeds max\n");
		return 0;
	}
	
	// get N value
	while((nucleotide = getc(fp_inputfile)) != '\n')
    {
        n++;
    }
    printf("n: %d\n", n);

    positions = (n / 16)+1; 

    // get T value
    while(nucleotide != EOF)
    {   
        while((nucleotide = getc(fp_inputfile)) != '\n' && nucleotide != EOF);
        t++;
    }
    t++;
    printf("t: %d\n", t); 
	
	// fill in dna_matix
	fclose(fp_inputfile);
    fp_inputfile = fopen(argv[1], "r");
    fscanf(fp_inputfile, "%d\n", &k);

	for(i = 0; i < t; i++)
	{
		for(j = 0; j < n + 1; j++)
		{
			int bit_encoding = 0;
			char base = getc(fp_inputfile);

	        /*
	         * A:=00 
	         * T:=01 
	         * G:=10 
	         * C:=11
	        */
	        if( base == 'A')
	        {
	        	bit_encoding = 0;
	        }
	        else if( base == 'T')
	        {
	        	bit_encoding = 1;
	        }
	        else if( base == 'G')
	        {
	        	bit_encoding = 2;
	        }
	        else if( base == 'C')
	        {
	        	bit_encoding = 3;
	        }

	        if( base != '\n')
	        {
	        	dna_matrix[i*positions + (j / 16)] += bit_encoding << ((j % 16) * 2);
	        }
		}
    } 
    fclose(fp_inputfile);

////////////////////////////////////////////////////////////

    int hamming_total;
    int min_dist = INT_MAX;

    uint32_t best_word = 0;

    for(s = 0; s < (1 << (k*2)); s++)
    {
        hamming_total = getTotalDistance( t, k, n);
    
        if( hamming_total < min_dist)
        {
            min_dist = hamming_total;
            best_word = s;
        }
    }

    // create output file
    fp_outputfile = fopen("output.txt", "w+");
    if (fp_outputfile == NULL)
    {
    	printf("Output file is NULL\n");
		return 0;
    }

    printf("BEST PATTERN: ");

    // form the output string and result
    for( i = 0; i < k; i++)
    {
    	nucleotide = (best_word >> (2*i)) & 0x03;
        printf( "%c", getChar( nucleotide));
        putc( getChar( nucleotide), fp_outputfile);
    }
    
    printf("\nTOTAL HAMMING DISTANCE: %d \n", min_dist);


    fclose(fp_outputfile);
    return 0;
}