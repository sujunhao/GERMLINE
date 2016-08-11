
#include "Individual.h"

Share::Share( Individual * cip , boost::dynamic_bitset<>& hash ) : ms( hash )
{
	add( cip );
}

void Share::assertMatches()
{
	set<Individual*>::iterator i , ii;
	Match * m;

	for ( i = matches.begin() ; i != matches.end() ; i++ )
	{
		ii = i;
		for ( ++ii ; ii != matches.end() ; ii++ )
		{
			// Check if this pair matched in previous word
			m = (*i)->getPreviousMatch( *ii );
			if ( m != NULL )
			{
				// This match can be incremented
				m->end_ms = position_ms;
				m->addConfidence();
				// Increment in the first individual
				(*i)->incrementMatch( *ii );
				// Increment in the second individual
				(*ii)->incrementMatch( *i );
			}
			else
			{
				// This match must be created
				m = createMatch( *i , *ii );
				// Extend the match backwards
				m->extendBack();
				// Mark asserted
				(*i)->addAsserted( *ii , m );
				(*ii)->addAsserted( *i , m );
			}

			// Check if match asserts
			if( !m->saved && m->assertLength() )
			{
				m->saved = true;
				(*i)->insertCompletedMatch( *ii , m );
				(*ii)->insertCompletedMatch( *i , m );
			}
		}
	}
}

Match * Share::createMatch(Individual * c1 , Individual * c2)
{
	Match * new_match = new Match();
	new_match->end_ms = new_match->start_ms = position_ms;

	new_match->node[0] = c1;
	new_match->node[1] = c2;

	return new_match;
}

void Share::add(Individual * cip)
{
	matches.insert( cip );
}
