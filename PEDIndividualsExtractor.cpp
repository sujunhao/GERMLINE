// KOSIndividualsExtractor.cpp: Manages input from Plink files

#include "PEDIndividualsExtractor.h"
#include <cctype>
#include <iostream>
using namespace std;

// PEDIndividualsExtractor(): default constructor
PEDIndividualsExtractor::PEDIndividualsExtractor()
{}

void PEDIndividualsExtractor::stripWhitespace()
{
	if (stream.is_open())
	{
		char c;
		while((c=stream.peek())!=EOF && isspace(c))
		stream.get();
	}
	else
	{
		cerr << "WARNING:PolymorphicIndividualsExtractor::stripWhiteSpace():stream is not open" << endl;
		valid_flag = false;
	}
}

void PEDIndividualsExtractor::loadInput(Individuals& inds)
{
	individualsP = &inds;
	snps.processMAPFile();
	snps.beginChromosome();
	numberOfMarkers = snps.size();

	while (! stream.eof() && getIndividuals())
	{
		// getIndividuals();
		stream.seekg(numberOfMarkers*4 + 1,ios::cur);
	}
	stream.clear();
}

// getInput(): gets individuals from .ped file
void PEDIndividualsExtractor::getInput()
{
	
	
	cerr << "Please enter the MAP file name" << endl;
	cin >> map_file;
	cerr << "Please enter the PED file name" << endl;
	cin >> ped_file;
	
	if ( !snps.setFile( map_file ) )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open map file" << endl;
		valid_flag = false;
		return;
	}

	stream.open( ped_file.c_str() );
	if ( !stream )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open ped file" << endl;
		valid_flag = false;
		return;
	}
}

// getIndividuals(): gets the next nuclear family from stream
bool PEDIndividualsExtractor::getIndividuals()
{
	string discard, ID, famID, parentID1, parentID2, gender;
	Sex sex;
	stream >> famID >> ID >> parentID1 >> parentID2 >> gender >> discard;
	if (famID=="") return false;
	if(!stream.good()) return false;
	sex = gender == "0"?MALE:FEMALE;
	Individual * new_ind = new Individual;
	new_ind->setOffset(stream.tellg());
	new_ind->setID(ID);
	new_ind->setParent(parentID1,parentID2);
	new_ind->setSex(sex);
	new_ind->setFamily(famID);
	individualsP->addIndividual(new_ind);
	return true;
}

void PEDIndividualsExtractor::getCompleteMarkerSet(Individual * p)
{
	ind = p;
	stream.seekg(ind->getOffset() + 4*snps.getROIStart().getMarkerNumber() + 4*position_ms*MARKER_SET_SIZE + 1);
	getCompleteMarkerSet();
}

// getCompleteMarkerSet(): extracts the markers for a MarkerSet from stream and
//  adds this MarkerSet to ind after mapping markers to binary using snps
void PEDIndividualsExtractor::getCompleteMarkerSet()
{
	MarkerSet markerSet[2];
	unsigned int maxsize = snps.currentSize();
	
	for (int position = 0; position < MARKER_SET_SIZE; position++)
	{
		if(position_ms*MARKER_SET_SIZE+position >= maxsize) break;
		for(int al=0;al<2;al++){
			stripWhitespace();

			char marker = stream.peek();
			if ( snps.mapNucleotideToBinary(marker,position_ms*MARKER_SET_SIZE+position) == 1 )
				markerSet[al].set(position , true );
			
			stream.get();
		}
	}

	ind->addMarkerSet(UNTRANS,markerSet[0]);
	ind->addMarkerSet(TRANS,markerSet[1]);
}


// end PEDIndividualsExtractor.cpp
