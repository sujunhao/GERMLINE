// MatchesBuilder.cpp: builds matches from individuals

#include "MatchesBuilder.h"

unsigned int position_ms;
unsigned int num_sets;
int ALLOWED_MASKED;

SNPs snps;

// MatchesBuilder(): default constructor
MatchesBuilder::MatchesBuilder(Individuals &inds , PolymorphicIndividualsExtractor* pie )
{
	individualsP = &inds;
	pieP = pie;
}

// buildMatches(): builds matches from individuals
void MatchesBuilder::buildMatches()
{
	cerr << "Read Markers" << endl;
	ms_start = 0; ms_end = num_sets;
	readAllMarkers();

	cerr << "Match Markers" << endl;
	matchAllMarkers();
}

void MatchesBuilder::printHaplotypes(string fout_name)
{
	ofstream fout( fout_name.c_str() );

	for(individualsP->begin();individualsP->more();)
	{
		individualsP->next()->print(fout,0,num_sets);
	}
	fout.close();
}

// matchAllMarkers(): builds matches for individuals considering all markers
void MatchesBuilder::matchAllMarkers()
{
	double avg_errors = 0;
	for (position_ms = ms_start; position_ms < ms_end; position_ms++)
	{
		cerr << "\r" << position_ms*100/ms_end << "%" << flush;
		matchMarkerSet();
	}
	cerr << '\r' << "100%" << endl;
	if ( !MAX_ERRORS_FIXED ) cerr << avg_errors << " average allowed error rate" << endl;
}

void MatchesBuilder::readAllMarkers()
{
	for (position_ms = ms_start; position_ms < ms_end; position_ms++)
	{
		cerr << "\r" << position_ms*100/ms_end << "%" << flush;
		readMarkerSet();
	}
	cerr << '\r' << "100%" << endl;
}

void MatchesBuilder::readMarkerSet()
{
	Individual * i;
	for(individualsP->begin();individualsP->more();)
	{
		i = individualsP->next();
		pieP->getCompleteMarkerSet(i);
	}
}

// matchMarkerSet(): builds matches for individuals considering markers in marker set.
void MatchesBuilder::matchMarkerSet()
{
	// Match:
	for(individualsP->begin();individualsP->more();)
		matchFactory.hash( individualsP->next() );

	// Verify:
	matchFactory.assertShares();

	// Extend:
	for(individualsP->begin();individualsP->more();)
		individualsP->next()->assertShares();
	
	matchFactory.initialize();
}

// end MatchesBuilder.cpp
