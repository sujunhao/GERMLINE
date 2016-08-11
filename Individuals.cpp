// Individuals.cpp: A collection of individuals

#include "Individuals.h"
#include <iostream>
using namespace std;


// Individuals(): default constructor
Individuals::Individuals()
{
}

Individuals::~Individuals()
{
	for(begin();more();)
	{
		iter->second->freeMatches();
		delete iter->second;
		pedigree.erase(iter++);
	}
}

void Individuals::freeMatches()
{
	for(begin();more();)
		next()->freeMatches();
}

void Individuals::freeMarkers()
{
	for(begin();more();)
		next()->clearMarkers();
}

// numIndividuals(): returns number of individuals
int Individuals::numIndividuals() const
{
	return (int)pedigree.size();
}

void Individuals::begin()
{
	iter = pedigree.begin();
}

bool Individuals::more()
{
	return iter != pedigree.end();
}

Individual * Individuals::next()
{
	return (iter++)->second;
}

// addIndividual(): adds an Individual object
void Individuals::addIndividual(Individual * ind)
{
	pedigree.insert(make_pair(ind->getFamily() + " " + ind->getID(), ind));
}

// end Individuals.cpp
