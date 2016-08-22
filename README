#GERMLINE


##About

This version base on Germline 1.4.1 [detail](http://www1.cs.columbia.edu/~gusev/germline/)


##New feature


Flag		|	Default		|		Description
-----		|	---------		|		---------------
-poi		|	    -			|		person of interest, test between first individual with other
-hg		|	    -			| 		have gap, allow test with indivudual having missing snp(mark by 0 or '-')





##fix bug

 >in file PEDIndividualsExtractor.cpp, an infinite loop due to while statement


PEDIndividualsExtractor.cpp

```
	void PEDIndividualsExtractor::loadInput(Individuals& inds)
			.....
-			while (! stream.eof())
-			{
-				getIndividuals();
-				stream.seekg(numberOfMarkers*4 + 1,ios::cur);
-			}
+			while (! stream.eof() && getIndividuals())
+			{
+				// getIndividuals();
+				stream.seekg(numberOfMarkers*4 + 1,ios::cur);
+			}

...

-	void PEDIndividualsExtractor::getIndividuals()
		...
-		if(!stream.good()) return;
		...
		individualsP->addIndividual(new_ind);

+  	bool PEDIndividualsExtractor::getIndividuals()
		...
+		if (famID=="") return false;
+		if(!stream.good()) return false;
		...
		individualsP->addIndividual(new_ind);
+		return true;


```


PEDIndividualsExtractor.h

```
-		    void getIndividuals();
+		    bool getIndividuals();
```



##feature POI


### define flag -poi

BasicDefinitions.h

```
+	extern bool POI;
```

>init and accept flag

GERMLINE_0001.cpp
```
+	bool POI = false;
	....
+	else if( strncmp(argv[i], "-poi", strlen("-poi")) == 0 )	POI = true;

```

### init person of interest

> define and init the interest  individual

Individuals.h

```
+	#ifndef HAVEPOI
+	#define HAVEPOI
+	extern Individual * poi;
+	void definepoi();
+	#endif /* HAVEPOI */

```

Individuals.cpp

```
+		Individual * poi;

```

> save the poi

Individuals.cpp

```
	void Individuals::addIndividual(Individual * ind)
		.....
+		if (POI && !poi)
+		{
+			poi = ind;
+			// cout << "the targer in " + ind->getID() << endl;
+		}
```

###change the add Match condition

share.cpp
```
	void Share::assertMatches()
	...
+		if (POI)
+		{
+			if ((*i) != poi && (*ii)  != poi) 
+				break;
+		}

```




##HG

in order to allow missing snp, each individual should save their missing infomastion, so beside the hash table in each match, we add a bit to save missing snp

### define flag -hg

BasicDefinitions.h

```
+	extern bool HG;
```

>init and accept flag

GERMLINE_0001.cpp
```
+	bool POI = false;
	....
+	else if( strncmp(argv[i], "-hg", strlen("-hg")) == 0 )	HG = true;

```

### read missing snp postion
first define list to save missing postion, it will be named xmark in markset.h
and the define some function to apply at xmark, after that, we can read missing snp in below

PEDIndividualsExtractor.cpp
```
	void PEDIndividualsExtractor::getCompleteMarkerSet()
+		if (HG)
+			if (marker=='0' || marker=='-')	
+				markerSet[al].xset(position , true );
+			else if ( snps.mapNucleotideToBinary(marker,   \
				position_ms*MARKER_SET_SIZE+position) == 1 )
+				markerSet[al].set(position , true );
```
### change the extension condition
after add xmark data struct, we can use it in dynamic programming step of match function. Because if use xmark in init hashing step will involve hash table multipmeaning problem, so we add xmarker function in extension step. Like to check hom mismatch and het mismatch, we check add xmarker mismatch

Match.cpp
bool Match::approxEqual()

```
+	boost::dynamic_bitset<> xtmp;
+	xtmp.resize(MARKER_SET_SIZE);
+	xtmp.flip();
....
+		if (HG) xtmp = (node[0]->getChromosome( a )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( b )->getMarkerSet()->xgetMarkerBits()).flip();
-		if ( (int)( (node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits())).count() <= MAX_ERR_HOM )
+		if ( (int)( (node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()) & xtmp ).count() <= MAX_ERR_HOM )
...
+		if (HG)
+		{
+			xtmp = node[0]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[0]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits()
+				 | node[1]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits();
+				xtmp = xtmp.flip();
+		}
+		mask = mask & xtmp;
-		if ( ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits())).count() <= MAX_ERR_HET )
+		if ( ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= MAX_ERR_HET )
...
```
in `Match::scanRight`, `Match::scanLeft` and `Match::diff` function
change the check HET process and will finish the HG function


