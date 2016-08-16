// GERMLINE.cpp: GEnetic Relationship Miner for LINear Extension

#include "GERMLINE.h"
#include "math.h"
#include <iostream>

using namespace std;

ofstream MATCH_FILE;

// GERMLINE(): default constructor
GERMLINE::GERMLINE()
{
}


// mine(): main function for GERMLINE
void GERMLINE::mine()
{
	PolymorphicIndividualsExtractor * pie = inputManager.getPie();
	inputManager.getIndividuals();
	if ( ! pie->valid() ) return;
	string out = inputManager.getOutput();
	cout << "asd";

	pie->loadInput(individuals);
	MatchesBuilder mb( individuals , pie );

	ofstream fout( ( out + ".log" ).c_str() );
	MATCH_FILE.open( ( out + ".match" ).c_str() );
	fout << "Minimum match length: " << MIN_MATCH_LEN << " cM" << endl;
	if( MAX_ERRORS_FIXED ) fout << "Allowed mismatching bits: " << MAX_ERR_HOM << " " << MAX_ERR_HET << endl;
	else fout << "Allowed mismatching bits: auto" << endl;
	fout << "Word size: " << MARKER_SET_SIZE << endl;
	if ( ROI )
		fout << "Target Region: " << snps.getROIStart().getSNPID() << " - " << snps.getROIEnd().getSNPID() << endl;
	else
		fout << "Target Region: all" << endl;
	
	time_t timer[2]; time( &timer[0] );
	
	if ( ROI )
	{
		snps.beginChromosome();
		num_sets = (long)ceil((double)snps.currentSize()/(double)MARKER_SET_SIZE);
		mb.buildMatches();
		individuals.freeMatches();
		if(PRINT_HAPS) mb.printHaplotypes( out + "." + snps.getChromosome() + ".haps" );
		individuals.freeMarkers();
	}
	else
	{
		for ( snps.beginChromosome() ; snps.moreChromosome() ; snps.nextChromosome() )
		{
			num_sets = (long)ceil((double)snps.currentSize()/(double)MARKER_SET_SIZE);
			mb.buildMatches();
			individuals.freeMatches();
			if(PRINT_HAPS) mb.printHaplotypes( out + "." + snps.getChromosome() + ".haps" );
			individuals.freeMarkers();
		}
	}

	time( &timer[1] );

	fout << "Total runtime (sec): " << difftime( timer[1] , timer[0] ) << endl;
	fout.close(); MATCH_FILE.close();
}


// end GERMLINE.cpp
