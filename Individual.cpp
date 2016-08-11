// Individual.cpp: An individual with genetic data

#include "Individual.h"
using namespace std;

// Individual(): default constructor
Individual::Individual() : sex(MALE),affected(false)
{}

void Individual::freeMatches()
{
	for(map<Individual*, Match*>::iterator iter = previous_matches.begin(); iter != previous_matches.end();) 
	{
		deletePreviousMatch(iter);
	}
	previous_matches.clear();

	for(map<Individual * , Match * >::iterator iter = complete_matches.begin(); iter != complete_matches.end();)
	{
		deleteCompletedMatch(iter);
	}
	complete_matches.clear();
}

Match * Individual::getPreviousMatch( Individual * i )
{
	map<Individual*, Match*>::iterator it = previous_matches.find( i );
	if ( it == previous_matches.end() ) return NULL; else return it->second;
}

void Individual::assertHomozygous()
{
	map<Individual*, Match*>::iterator iter;
	Match * m;
	if( (iter = previous_matches.find(this)) != previous_matches.end() )
	{
		// increment this match
		m = iter->second;
		m->end_ms = position_ms;
		incrementMatch( iter );

	} else
	{	
		// this is a new match
		m = new Match();
		m->end_ms = m->start_ms = position_ms;
		m->node[0] = m->node[1] = this;
		m->extendBack();
		asserted_matches.insert( make_pair( this , m ) );
	}

	// Check if match asserts
	if( !m->saved && m->assertLength() )
	{
		m->saved = true;
		insertCompletedMatch( this , m );
	}
}

void Individual::assertShares()
{
	map<Individual*, Match*>::iterator iter;
	set<Individual*>::iterator cip;

	// try to extend previous matches that did not match currently
	for( iter = previous_matches.begin(); iter != previous_matches.end(); )
	{
		// Can we increment?
		if ( iter->second->approxEqual() )
		{
			iter->second->end_ms = position_ms;

			// check if it asserts
			if( !iter->second->saved && iter->second->assertLength() )
			{
				iter->second->saved = true;
				if ( iter->first != this ) iter->first->insertCompletedMatch( this , iter->second );
				insertCompletedMatch( iter->first , iter->second );
			}
			
			if ( iter->first != this ) iter->first->incrementMatch( this );
			incrementMatch( iter );
		}
		else
		{
			deletePreviousMatch( iter );
		}
	}
	previous_matches.swap(asserted_matches);
}

void Individual::deleteCompletedMatch( Individual * cip )
{
	complete_matches.erase( cip );
}

void Individual::deletePreviousMatch( Individual * cip )
{
	previous_matches.erase( cip );
}

void Individual::deleteCompletedMatch( map<Individual*, Match*>::iterator& iter )
{
	iter->second->print( MATCH_FILE );
	iter->first->deleteCompletedMatch( this );
	delete iter->second;
	complete_matches.erase( iter++ );
}

void Individual::deletePreviousMatch( map<Individual*, Match*>::iterator& iter )
{
	if ( iter->second->saved )
	{
		// print it
		iter->second->print( MATCH_FILE );
		// Delete it from the completed list
		complete_matches.erase( iter->first );
		if ( iter->first != this ) iter->first->deleteCompletedMatch( this );
	}

	delete iter->second;
	// erase from the list
	if ( iter->first != this ) iter->first->deletePreviousMatch( this );
	previous_matches.erase(iter++);
}

void Individual::insertCompletedMatch( Individual * cip , Match * m )
{
	complete_matches.insert( make_pair( cip , m ) );
}

void Individual::incrementMatch( map<Individual*, Match*>::iterator& iter)
{
	asserted_matches.insert( *iter );
	previous_matches.erase( iter++ );
}

void Individual::incrementMatch( Individual * i )
{
	map<Individual*, Match*>::iterator iter = previous_matches.find( i );
	if ( iter != previous_matches.end() )
	{
		asserted_matches.insert( *iter );
		previous_matches.erase( iter++ );
	}
}

void Individual::addAsserted(Individual * target, Match * m)
{
	asserted_matches.insert( make_pair(target,m) );
}

void Individual::initialize()
{
	shares.clear();
}


void Individual::print(ostream& out,long start,long end)
{
	for(int i=0;i<2;i++)
	{
		out << getFamily() << ' ' << getID() << '\t';
		h[i].print(out,start,end);
		out << endl;
	}
}

int Individual::numHet()
{
	return int(( h[0].getMarkerSet()->getMarkerBits() ^ h[1].getMarkerSet()->getMarkerBits() ).count());
}

bool Individual::isHeterozygous()
{
	return !( h[0].getMarkerSet()->equal( h[1].getMarkerSet() ) );
}

bool Individual::isHeterozygous(int i)
{
	return h[0].getMarkerSet()->getMarker(i) != h[1].getMarkerSet()->getMarker(i);
}

string Individual::getParentID(int i) const
{
	return p_id[i];
}

void Individual::setParent(int i,string s){p_id[i] = s;}
void Individual::setParent(string s0,string s1){p_id[0] = s0; p_id[1] = s1;}

void Individual::setOffset(streamoff o)
{
	offset = o;
}

streamoff Individual::getOffset()
{
	return offset;
}

string Individual::getFamily() const
{
	return family;
}

// getPedigree(): accessor for pedigree
string Individual::getPedigree() const
{
	return pedigree;
}


// getID(): accessor for ID
string Individual::getID() const
{
	return ID;
}

// getSex(): accessor for sex
Sex Individual::getSex() const
{
	return sex;
}


// isAffected(): answers whether individual is affected
bool Individual::isAffected() const
{
	return affected;
}

// getChild(): accessor for children
string Individual::getChild(int index) const
{
	if (index >= 0 && index < (int)children.size())
		return children[index];
	else
	{
		cerr << "WARNING:Individual::getChild():index out of range" << endl;
		return "";
	}
}

/*
// getMarkerSet(): gets MarkerSet from a chromosome
MarkerSet Individual::getMarkerSet(int ct)
{
	return ms[ct];
}
*/

Chromosome * Individual::getAlternateChromosome( Chromosome * c)
{
	if( &(h[0]) == c ) return &(h[1]); else return &(h[0]);
}

Chromosome * Individual::getChromosome(int ct)
{
	return &(h[ct]);
}

void Individual::setFamily(string f){
	family = f;
}

// setPedigree(): mutator for pedigree.
void Individual::setPedigree(string ped)
{
	pedigree = ped;
}


// setID(): mutator for ID.
void Individual::setID(string id)
{
	ID = id;
}


// setSex(): mutator for sex.
void Individual::setSex(Sex s)
{
	sex = s;
}


// setAffected(): mutator for affected.
void Individual::setAffected(bool aff)
{
	affected = aff;
}


// addChild(): mutator for children
void Individual::addChild(string child)
{
	children.push_back(child);
}

void Individual::clearMarkers()
{
	h[0].clear();
	h[1].clear();
}

// addMarkerSet(): adds MarkerSet to a chromosome
void Individual::addMarkerSet(int ct, const MarkerSet& ms)
{
	h[ct].addMarkerSet(ms);
}

// operator<<(): overloaded stream insertion operator
ostream& operator<<(ostream &fout, Individual& ind)
{
	fout << endl << ind.getFamily() << ":" << ind.getID() << endl;
	fout << ind.getChromosome(0) << endl;
	fout << ind.getChromosome(1) << endl;
	return fout;
}


// end Individual.cpp
