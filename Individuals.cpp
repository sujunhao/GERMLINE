// Individuals.cpp: A collection of individuals
#include <iostream>
#include "Individuals.h"

using namespace std;


Individual * poi;
void definepoi()
{
	cout << "define" << endl;
}
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
	if (POI && !poi)
	{
		poi = ind;
		// cout << "the targer in " + ind->getID() << endl;
	}
}

// end Individuals.cpp
