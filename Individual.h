// Individual.h: An individual with genetic data

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#pragma warning(disable: 4786)
#include "BasicDefinitions.h"
#include "Chromosome.h"
#include "Match.h"
#include <string>
#include <list>
#include <map>
#include <set>
#include <iostream>

using namespace std;

class Individual;
class Match;
class Share
{
public:
	Share( Individual * , boost::dynamic_bitset<>& );
	
	void add(Individual * );
	void assertMatches();

private:
	Match * createMatch(Individual * c1 , Individual * c2);
	set<Individual*> matches;
	boost::dynamic_bitset<>& ms;
};


class Individual
{
public:

	/** Match Tracking **/
	void initialize();

	void addShare(Share*);
	void assertShares();
	void assertHomozygous();

	list<Share*>& getShareList();

	void freeMatches();
	void insertCompletedMatch( Individual * , Match * );

	Match * getPreviousMatch( Individual * i );

	void deleteCompletedMatch( Individual* );
	void deleteCompletedMatch( map<Individual*, Match*>::iterator& );
	void deletePreviousMatch( Individual* );
	void deletePreviousMatch( map<Individual*, Match*>::iterator& );

	void incrementMatch( Individual* );
	void incrementMatch( map<Individual*, Match*>::iterator& );

	void addAsserted(Individual*,Match*);

	/** Match Tracking **/
	int relatedness( Individual * );
	void print(ostream&,long,long);
	
	bool isHeterozygous();
	bool isHeterozygous(int i);
	int numHet();
	string getParentID(int) const;

	void setParent(int i,string s);
	void setParent(string s0,string s1);
	
	// Individual(): default constructor
	// Precondition: None.
	// Postcondition: All strings and vectors are initialized
	//  to empty. sex is initialized to MALE.
	//  affected is initialized to false.
	Individual();

	// getPedigree(): accessor for pedigree
	// Precondition: None.
	// Postcondition: Returns value of pedigree.
	string getPedigree() const;

	// getID(): accessor for ID
	// Precondition: None.
	// Postcondition: Returns value for ID
	string getID() const;

	string getFamily() const;

	// getSex(): accessor for sex
	// Precondition: None.
	// Postcondition: Returns value for sex.
	Sex getSex() const;

	// isAffected(): answers whether individual is affected
	// Precondition: None.
	// Postcondition: Returns true if individual is affected;
	//  otherwise returns false.
	bool isAffected() const;

	void setOffset(streamoff);
	streamoff getOffset();

	// getChild(): accessor for children
	// Precondition: None.
	// Postcondition: if 0=<index<numberOfChildren(),
	//  returns value of child index; otherwise prints
	//  a warning and returns the empty string.
	string getChild(int index) const;

	// getMarkerSet(): gets MarkerSet from a chromosome
    // Precondition: None.
	// Postcondition: If memory has been allocated for the MarkerSet at position index,
	//   then the MarkerSet at position index in chromosome determined by ct has been
	//   returned; otherwise, a warning is printed and the default MarkerSet is returned.
	// MarkerSet getMarkerSet(int) const;

	// getChromosome(): accessor for Chromosome references
	// Precondition: None.
	// Postcondition: returns reference to Chromosome chrom.
	Chromosome* getChromosome(int);
	Chromosome* getAlternateChromosome(Chromosome * );

	// setPedigree(): mutator for pedigree.
	// Precondition: None.
	// Postcondition: pedigree has been set to ped.
    void setPedigree(string ped);

	void setFamily(string f);

	// setID(): mutator for ID.
	// Precondition: None.
	// Postcondition: ID has been set to id.
    void setID(string id);

	// setSex(): mutator for sex.
	// Precondition: None.
	// Postcondition: sex has been set to s.
    void setSex(Sex s);

    // setAffected(): mutator for affected.
	// Precondition: None.
	// Postcondition: affected has been set to aff.
    void setAffected(bool aff);

	// addChild(): mutator for children
	// Precondition: None.
	// Postcondition: child has been added to children.
	void addChild(string child);

	// addMarkerSet(): adds MarkerSet to a chromosome
    // Precondition: None.
	// Postcondition: If chromosome identified by ct has room, then ms has been
	//  added to the end of the chromosome; otherwise a warning message is printed..
	void addMarkerSet(int, const MarkerSet& ms);

	// clearMarkers(): clears all MarkerSets from this individual
	void clearMarkers();

private:

	// family ID of the individual if known
	string family;
	// pedigree of the individual
	string pedigree;
	// ID of the individual
	string ID;
	// parents
	string p_id[2];
	// children of the individual if known
	vector<string> children;
	// sex of the individual if known
	Sex sex;
	// stores whether or not the individual is affected
	// by the condition associated with the input file
	bool affected;

	Chromosome h[2];
	// sequence start in file
	streamoff offset;
	
	list<Share*> shares;
 	map<Individual*, Match*> previous_matches;
	map<Individual*, Match*> asserted_matches;
	map<Individual*, Match*> complete_matches;
};


// operator<<(): overloaded stream insertion operator
// Precondition: fout is a value ostream.
// Postcondition: ind has been sent to fout.
ostream &operator<<(ostream &fout, Individual& ind);

#endif

// end Individual.h
