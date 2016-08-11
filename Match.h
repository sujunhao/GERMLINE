#ifndef MATCH_H
#define MATCH_H

#include "BasicDefinitions.h"
#include "SNPs.h"
#include "MarkerSet.h"
#include "Individual.h"
#include "math.h"
#include <iomanip>

class Individual;

class Match
{
public:

	Match() { saved = false; reviewed = false; conf = CUR_CONF; }

	Individual * node[2];

	unsigned int start_ms, end_ms;
	bool saved , reviewed;

	void addConfidence() { conf += CUR_CONF; }
	long physicalStart(){	return snps.getSNP(start_ms*MARKER_SET_SIZE).getPhysPos(); }
	long physicalEnd(){		return snps.getSNP((end_ms+1)*MARKER_SET_SIZE-1).getPhysPos(); }
	string chromosome(){	return snps.getSNP(start_ms*MARKER_SET_SIZE).getChr(); }
	double length()
	{ 
		return snps.getDistance(start_ms*MARKER_SET_SIZE,(end_ms+1)*MARKER_SET_SIZE-1);
	}
	bool assertLength()
	{
		unsigned int buff_start;
		if ( start_ms > 0 ) buff_start = start_ms - 1; else buff_start = start_ms;
		return snps.getDistance(buff_start*MARKER_SET_SIZE,(end_ms+2)*MARKER_SET_SIZE-1) >= MIN_MATCH_LEN;
	}
	Individual * getTarget(Individual * me){ if( node[0] == me ) return node[1]; else return node[0]; }
	void print( ostream& );
	bool approxEqual();
	void erase();

	void extendBack();

private:
	int scanLeft( unsigned int ms );
	int scanRight( unsigned int ms );
	int diff( unsigned int ms );
	bool isHom( int n , unsigned int ms );
	double conf;
};

#endif

