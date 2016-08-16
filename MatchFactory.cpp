// MatchFactory.cpp: Generates matches from individuals

#include "MatchFactory.h"


// MatchFactory(): default constructor
MatchFactory::MatchFactory()
{}

void MatchFactory::assertShares()
{
	for ( iter = segments.begin() ; iter != segments.end() ; iter++ )
	{
		iter->second.assertMatches();
	}
}

int MatchFactory::size()
{
	return (int)segments.size();
}

// initialize(): initializes object
void MatchFactory::initialize()
{
	segments.clear();
}	
void MatchFactory::hash( Individual * i )
{
	int haps , het = i->numHet();
	
	i->initialize();
	if ( het == 0 ) haps = 1;
	else haps = 2;

	if ( ALLOW_HOM && het <= MAX_ERR_HOM + MAX_ERR_HET ) i->assertHomozygous();
	
	for ( int c = 0 ; c < haps ; c++ )
	{
		if (HG)
		{
			//if xmark empty
			if ( (i->getChromosome(c)->getMarkerSet()->xgetMarkerBits()).count() == 0 )
			{
				boost::dynamic_bitset<>& ms = i->getChromosome(c)->getMarkerSet()->getMarkerBits();
				if ( (iter = segments.find( ms )) == segments.end() )
					segments.insert( make_pair ( ms , Share( i , ms ) ) );
				else
					iter->second.add( i );
			}
			continue;
		}
		boost::dynamic_bitset<>& ms = i->getChromosome(c)->getMarkerSet()->getMarkerBits();
		if ( (iter = segments.find( ms )) == segments.end() )
			segments.insert( make_pair ( ms , Share( i , ms ) ) );
		else
			iter->second.add( i );
	}
}

// end MatchFactory.cpp
