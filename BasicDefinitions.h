// BasicDefinitions.h: constants and enumerations for GERMLINE

#ifndef BASICDEFINITIONS_H
#define BASICDEFINITIONS_H
// #define VERBOSE
#include <vector>
#include <fstream>
#include <iostream>
#include <string.h>


using namespace std;

// where we are in the sequence (markerset & physical)
extern unsigned int position_ms;
extern unsigned int num_sets;

// parameters
extern double MIN_MATCH_LEN;
extern int MARKER_SET_SIZE;
extern int MAX_ERR_HOM;
extern int MAX_ERR_HET;
extern bool MAX_ERRORS_FIXED;
extern bool PRINT_HAPS;
extern bool PRINT_MATCH_HAPS;
extern bool ROI;
extern bool POI;
extern bool HG;
extern bool HAP_EXT;
extern bool WIN_EXT;
extern bool ALLOW_HOM;
extern bool CONF;
extern double CUR_CONF;

class SNPs;
extern SNPs snps;
extern ofstream MATCH_FILE;

enum ErrorType{RECOMB=0,MI=1};
const int HET=2;
const int MIS=9;
const short PAR_M = 0;
const short PAR_F = 1;

// type for file format
enum FileFormat{HAPS,PED,HM};
// type for families
enum FamilyType{NOPARENT,SINGLEPARENT,TWOPARENT};
// type for sex
enum Sex{MALE,FEMALE,UNKNOWN};
// type for nucleotides
enum Nucleotide{A=1,C=2,G=3,T=4};
// type for chromosomes types
enum ChromosomeType{TRANS=0,UNTRANS=1,MISSING=-1};
// type for chromosome ids
enum ChromosomeID{ONE, TWO, THREE, FOUR, FIVE, SIX, SEVEN, EIGHT, NINE, TEN, ELEVEN, TWELVE,
                THIRTEEN, FOURTEEN, FIFTEEN, SIXTEEN, SEVENTEEN, EIGHTEEN, NINETEEN, TWENTY,
				TWENTYONE, TWENTYTWO, X, Y};

#endif

// end of BasicDefinitions.h
