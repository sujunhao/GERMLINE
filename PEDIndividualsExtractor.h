#ifndef PEDIndividualsExtractor_H
#define PEDIndividualsExtractor_H

#include "BasicDefinitions.h"
#include "NucleotideMap.h"
#include "PolymorphicIndividualsExtractor.h"
#include <string>
using namespace std;

class PEDIndividualsExtractor : public PolymorphicIndividualsExtractor
{
public:

	// PEDIndividualsExtractor(): default constructor
	// Precondition: None.
	// Postcondition: three input files have been obtained from user
	PEDIndividualsExtractor();

	// getInput(): gets input from .haps file
	// Precondition: Input is a .ped file in a valid format
	// Postcondition: inds points to individuals from .ped file
    void getInput();
	void loadInput(Individuals &inds);
	void getCompleteMarkerSet(Individual * p);
private:

	// getIndividuals(FamilyType ft): runs specific Individuals parser based on ft
	// Precondition: "ft" is set to a valid FamilyType (matches file input)
	// Postcondition: correct get method is run for family type
    bool getIndividuals();

	// getCompleteMarkerSet(int markerSetNumber): read marker set from stream, add to chromosomes
	// Precondition: stream is open and queued to beginning of sequence,
	//				 sequence contains alternating haplotype values.
	void getCompleteMarkerSet();

	void stripWhitespace();

	// file name for input file
	string ped_file;
	// file name for info file
	string map_file;

	ifstream stream;

	// holds map of nucleotides to binary/quad values
	NucleotideMap inp;
};

#endif

// end PEDIndividualsExtractor.h
