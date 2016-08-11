// SNP.cpp: Stores a SNP's information

#include "SNP.h"
#include <iostream>
using namespace std;


// SNP(): default constructor
SNP::SNP()
{
	chr = '0';
	varSet = 0;
	count[0] = count[1] = 0;
	centimorgan = -1;
}

void SNP::setCentimorgan( float cm )
{
	centimorgan = cm;
}

float SNP::getCentimorgan()
{
	return centimorgan;
}

float SNP::getFrequencySq( bool var )
{
	double total = count[0] + count[1];
	float p = float(double(count[(int)var])/total);
	return p * p;
}

float SNP::getConfidence()
{
	double total = count[0] + count[1];
	float p = float(double(count[0])/total);
	float q = float(double(count[1])/total);
	return 1 - 2*p*p*q*q;
}

int SNP::mapNucleotide(char nt)
{
	int ret = -1;
	if (varSet == 2){
		// BOTH variants are set
		if (nt == getVariant(1))
			ret = 1;
		else if (nt == getVariant(0))
			ret = 0;
		else
		{
			cerr << endl << "WARNING:SNPs::setmapNucleotideToBinary():"
				<< nt << " is not one of variant alleles for SNP "
				<< SNPID
				<< endl;
		}
	} else if(varSet == 1){
		// 1ST variant is set
		if (nt == getVariant(0)){
			ret = 0;
		} else {
			setVariant(1,nt);
			ret = 1;
		}
	} else {
		// NO variants are set
		setVariant(0,nt);
		ret = 0;
	}
	if(ret != -1) count[ret]++;
	return ret;
}

// getSNPID(): accessor for SNPID
string SNP::getSNPID() const
{
	return SNPID;
}


// getPhysPos(): accessor for physPos
long SNP::getPhysPos() const
{
	return physPos;
}

string SNP::getChr() const 
{
	return chr;
}

// getVariant(): accessor for variants
char SNP::getVariant(int i) const
{
	return variant[i];
}

void SNP::setChr(const string& c)
{
	chr = c;
}

// setSNPID(): mutator for SNPID
void SNP::setSNPID(const string& sid)
{
	SNPID = sid;
}


// setPhysPos(): mutator for physPos
void SNP::setPhysPos(long pp)
{
	physPos = pp;
}

void SNP::setMarkerNumber( int i )
{
	num = i;
}

int SNP::getMarkerNumber()
{
	return num;
}

// setVariant(): mutator for variants
void SNP::setVariant(int i, char nt)
{
	if (i == 1){
		variant[1] = nt;
		varSet = 2;
	}
	else if (i == 0){
		variant[0] = nt;
		varSet = 1;
	}
	else
		cerr << "WARNING:SNP::setVariant():index out of bounds" << endl;
}


// end SNP.cpp
