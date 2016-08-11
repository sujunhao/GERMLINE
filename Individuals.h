// Individuals.h: A collection of individuals

#ifndef INDIVIDUALS_H
#define INDIVIDUALS_H

#pragma warning(disable: 4786)
#include "BasicDefinitions.h"
#include "Chromosome.h"
#include "Individual.h"
#include <map>
#include <ostream>
using namespace std;

class Individuals
{
public:

	// Individuals(): default constructor
	// Precondition: None.
	// Postcondition: individuals is empty.
	Individuals();
	~Individuals();

	list<Individual*>& getPedigree();
	list<unsigned int>& getGenerations();

	// numIndividuals(): returns number of individuals
	// Precondition: None.
	// Postcondition: Returns number of individuals.
	int numIndividuals() const;

	// addIndividual(): adds an Individual object
	// Precondition: None.
    // Postcondition: ind has been added to individuals
	void addIndividual(Individual * ind);

	bool more();
	Individual* next();
	void begin();
	
	void freeMatches();
	void freeMarkers();

private:

	void permuteMarkerSet(Chromosome *, int, MarkerSet);
	// stores the individuals
	map<string,Individual*> pedigree;
	map<string,Individual*>::iterator iter;

	long sets;
};

#endif

// end Individuals.h
