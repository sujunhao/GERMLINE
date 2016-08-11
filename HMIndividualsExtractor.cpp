// HAPSIndividualsExtractor.cpp: Manages input from .haps files

#include "HMIndividualsExtractor.h"
#include <cctype>
#include <cmath>
#include <iostream>
using namespace std;


// HMIndividualsExtractor(): default constructor
HMIndividualsExtractor::HMIndividualsExtractor()
{}

void HMIndividualsExtractor::stripWhitespace()
{
	if (stream_phased.is_open())
	{
		char c;
		while((c=stream_phased.peek())!=EOF && isspace(c))
		stream_phased.get();
	}
	else
	{
		cerr << "WARNING:PolymorphicIndividualsExtractor::stripWhiteSpace():stream is not open" << endl;
		valid_flag = false;
	}
}

void HMIndividualsExtractor::loadInput(Individuals& inds)
{
	individualsP = &inds;

	snps.processLegendFile();
	snps.beginChromosome();
	numberOfMarkers = snps.size();

	while (!stream_sample.eof() && !stream_phased.eof()) getIndividuals();
	
	stream_sample.close();
	stream_phased.clear();
}


// getInput(): gets individuals from .haps file.
void HMIndividualsExtractor::getInput()
{

	cerr << "Please enter the _legend file name" << endl;
	cin >> fileLegend;
	cerr << "Please enter the _phased file name" << endl;
	cin >> filePhased;
	cerr << "Please enter the _sample file name" << endl;
	cin >> fileSample;
	
	if ( !snps.setFile( fileLegend ) )
	{
		cerr << "WARNING:HMIndividualsExtractor::getInput():cannot open legend file" << endl;
		valid_flag = false;
		return;
	}

	stream_phased.open(filePhased.c_str());
	if (!stream_phased) 
	{
		valid_flag = false;
		cerr << "WARNING:HMIndividualsExtractor::getInput():cannot open phased file" << endl;
		return;
	}

	stream_sample.open(fileSample.c_str());
	if (!stream_sample) 
	{
		valid_flag = false;
		cerr << "WARNING: HMIndividualsExtractor::openFileStream(): sample stream could not be opened" << endl;
		return;
	}
}

// getIndividuals(): gets the next nuclear family from stream
void HMIndividualsExtractor::getIndividuals()
{
	string sex_string, ID ,discard;
	Sex s;
	
	streamoff offset = stream_phased.tellg(); if ( offset > 0 ) offset--;

	getline ( stream_phased , discard );
	if ( stream_phased.eof() ) return;
	offset_buffer = ( stream_phased.tellg() - offset);
	getline ( stream_phased , discard );

	Individual * new_ind = new Individual;
	stream_sample >> ID >> sex_string;
	if(ID == "") return;

	if(sex_string == "1") s = MALE;
	else if(sex_string == "2") s = FEMALE;
	else s = UNKNOWN;
	new_ind->setID(ID);
	new_ind->setSex(s);
	new_ind->setOffset( offset );
	new_ind->setFamily("0");
	individualsP->addIndividual(new_ind);
}

void HMIndividualsExtractor::getCompleteMarkerSet(Individual * p)
{
	ind = p;
	getCompleteMarkerSet();
	stream_phased.clear();
}

// getCompleteMarkerSet(): extracts the markers for a MarkerSet from stream_phased and
//  adds this MarkerSet to ind after mapping markers to binary using snps
void HMIndividualsExtractor::getCompleteMarkerSet()
{
	MarkerSet markerSet[2];
	unsigned int maxsize = snps.currentSize();
	for(int al=0;al<2;al++)
	{
			stream_phased.seekg(
				ind->getOffset()
				+ 2 * snps.getROIStart().getMarkerNumber()
				+ 2 * position_ms * MARKER_SET_SIZE 
				+ al * ( offset_buffer - 1 )
				);

		for (int position = 0; position < MARKER_SET_SIZE; position++)
		{
			if(position_ms*MARKER_SET_SIZE+position >= maxsize) break;
			
			stripWhitespace();
			char marker = stream_phased.peek();

			// if ( snps.mapNucleotideToBinary(marker,position_ms*MARKER_SET_SIZE+position) == 1 )
			if ( marker == '1' ) markerSet[al].set(position , true );

			stream_phased.get();
		}
	}

	ind->addMarkerSet(UNTRANS,markerSet[0]);
	ind->addMarkerSet(TRANS,markerSet[1]);
}


// end HMIndividualsExtractor.cpp
