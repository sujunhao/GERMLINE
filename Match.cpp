#include "Match.h"

void Match::extendBack()
{
	// save old position values
	unsigned int SAVE_pms = position_ms;
	
	// iterate backwards through genome
	while(position_ms > 0)
	{
		position_ms--;
		if( !approxEqual() )
		{
			position_ms++;
			break;
		}
	}
	start_ms = position_ms;
	// restore saved values
	position_ms = SAVE_pms;
}

bool Match::approxEqual()
{
	// homozygosity check
	boost::dynamic_bitset<> xtmp();
	xtmp.resize(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits().size());
	xtmp.flip();
	if ( node[0] == node[1] )
	{
		if ( ALLOW_HOM )
		{
			if ( ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).count() 
				 <= ( MAX_ERR_HOM + MAX_ERR_HET ) ) return true; else return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		// 1. Haplotype extension
		for ( int a = 0 ; a < 2 ; a++ ) {
			for ( int b = 0 ; b < 2 ; b++ ) { 
				if (HG) xtmp = (node[0]->getChromosome( a )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( b )->getMarkerSet()->xgetMarkerBits()).flip();
				if ( (int)( (node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()) & xtmp ).count() <= MAX_ERR_HOM )
				{
					if ( CONF )
					{
						// Calculate the likelihood that this occurs by chance
						boost::dynamic_bitset<> tmp_xor = node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits();
						for ( int i = 0 ; i < MARKER_SET_SIZE; i++ )
						{
							// If these two positions matched, calculate the likelihood by chance
							if ( !tmp_xor[i] )
								conf += log10(snps.getSNP( position_ms * MARKER_SET_SIZE + i ).getConfidence());
						}
					}
					return true;
				}
			}
		}

		if ( HAP_EXT ) return false;

		// 2. Genotype extension
		// identify common homozygous SNPs
		boost::dynamic_bitset<> mask
			= ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip()
			& ( node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip();

		if (HG)
		{
			xtmp = node[0]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[0]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits()
				 | node[1]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits();
			xtmp = xtmp.flip();
		}
		mask = mask & xtmp;
		// assert that homozygous SNPs are identical
		if ( ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= MAX_ERR_HET )
		{
			if ( CONF )
			{
				// Calculate the likelihood that this occurs by chance
				boost::dynamic_bitset<> tmp_xor = (node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask;
				for ( int i = 0 ; i < MARKER_SET_SIZE; i++ )
				{
					// If these two positions matched, calculate the likelihood by chance
					if ( !tmp_xor[i] ) conf += log10(snps.getSNP( position_ms * MARKER_SET_SIZE + i ).getConfidence());
				}
			}
			return true;
		}
		else return false;
	}
}

int Match::scanLeft( unsigned int ms )
{
	bool err = false;
	int marker;

	boost::dynamic_bitset<> xtmp();
	xtmp.resize(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits().size());
	xtmp.flip();

	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	if (HG)
	{
		xtmp = node[0]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[0]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits()
			 | node[1]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits();
		xtmp = xtmp.flip();
	}
	mask = mask & xtmp;
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	for( marker = MARKER_SET_SIZE - 1 ; marker >= 0 && !err ; marker-- )
		if ( mask[marker] ) err = true;
	
	return marker;
}

int Match::scanRight( unsigned int ms )
{
	bool err = false;
	int marker;
	boboost::dynamic_bitset<> xtmp();
	xtmp.resize(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits().size());
	xtmp.flip();

	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	if (HG)
	{
		xtmp = node[0]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[0]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits()
			 | node[1]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits();
		xtmp = xtmp.flip();
	}
	mask = mask & xtmp;
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	for( marker = 0 ; marker < MARKER_SET_SIZE && !err ; marker++ )
		if ( mask[marker] ) err = true;
	
	return marker;
}

int Match::diff( unsigned int ms )
{
	boost::dynamic_bitset<> xtmp();
	xtmp.resize(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits().size());
	xtmp.flip();

	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	if (HG)
	{
		xtmp = node[0]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[0]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits()
			 | node[1]->getChromosome( 0 )->getMarkerSet()->xgetMarkerBits() | node[1]->getChromosome( 1 )->getMarkerSet()->xgetMarkerBits();
		xtmp = xtmp.flip();
	}
	mask = mask & xtmp;
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	return int(mask.count());
}

bool Match::isHom( int n , unsigned int ms )
{
	return ( node[n]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[n]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).count() <= ( MAX_ERR_HOM + MAX_ERR_HET );
}

void Match::print( ostream& fout )
{
	// extend this match from both ends
	unsigned int snp_start = start_ms * MARKER_SET_SIZE;
	unsigned int snp_end = ( end_ms + 1 ) * MARKER_SET_SIZE - 1;
	int marker;

	
	if ( WIN_EXT )
	{
	// backwards
	if( start_ms > 0 )
	{
		marker = scanLeft( start_ms - 1 );
		snp_start -= (MARKER_SET_SIZE - marker - 2);
	}
	// forwards
	if( end_ms < num_sets - 1 )
	{
		marker = scanRight( end_ms + 1 );
		snp_end += marker - 1;
	}
	}
	

	bool genetic;
	float distance;
	if ( ( distance = snps.getDistance(snp_start,snp_end,genetic)) < MIN_MATCH_LEN ) return;
	// print

	fout << node[0]->getFamily() << " ";
	fout << node[0]->getID() << '\t';
	fout << node[1]->getFamily() << " ";
	fout << node[1]->getID() << '\t';

	fout << snps.getSNP(snp_start).getChr() << '\t';
	fout << snps.getSNP(snp_start).getPhysPos() << ' ';
	fout << snps.getSNP(snp_end).getPhysPos() << '\t';

	fout << snps.getSNP(snp_start).getSNPID() << ' ';
	fout << snps.getSNP(snp_end).getSNPID() << '\t';

	fout << ( snp_end - snp_start + 1) << '\t';
	fout << setiosflags(ios::fixed) << setprecision(2) << distance << '\t';
	if ( genetic ) fout << "cM" << '\t'; else fout << "MB" << '\t';

	// get hamming distance & ignored bit count
	int dif = 0;
	for( unsigned int i = start_ms; i <= end_ms ; i++) { dif += diff( i ); }
	fout << dif;

	// calculate if homozygous
	if ( node[0] == node[1] ) fout << '\t' << 1 << '\t' << 1;
	else
	{
		for ( int n = 0 ; n < 2 ; n++ )
		{
			bool hom = true;
			for ( unsigned int i = start_ms ; i<= end_ms && hom ; i++ )
			{
				hom = isHom( n , i );
			}
			if ( hom ) fout << '\t' << 1; else fout << '\t' << 0;
		}
	}

	if ( CONF ) fout << '\t' << conf;

	fout << endl;

	// print the match sequences
	if ( PRINT_MATCH_HAPS )
	{
//		node[0]->getChromP()->print_snps(fout,snp_start,snp_end+1); fout << endl;
//		node[1]->getChromP()->print_snps(fout,snp_start,snp_end+1); fout << endl;
	}
}

