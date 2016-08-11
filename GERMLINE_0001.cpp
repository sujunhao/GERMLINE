#include "GERMLINE.h"
#include <string.h>

double MIN_MATCH_LEN = 5;
bool MAX_ERRORS_FIXED = true;
int MARKER_SET_SIZE = 128;
bool PRINT_HAPS = false;
bool PRINT_MATCH_HAPS = false;
bool ROI = false;
bool HAP_EXT = false;
bool WIN_EXT = false;
bool ALLOW_HOM = false;
bool CONF = false;
double CUR_CONF = 0;
int MAX_ERR_HOM = 2;
int MAX_ERR_HET = 0;

// main(): runs GERMLINE
int main(int argc, char* argv[])
{
	// parse arguments
	string rs_range[2] , map; map = rs_range[0] = rs_range[1] = "";
	bool bad_param = false;
	for(int i=1;i<argc;i++){
		if( strncmp(argv[i], "-min_m", strlen("-min_m")) == 0 && i < argc-1)				MIN_MATCH_LEN = atof(argv[++i]);
		else if( strncmp(argv[i], "-err_hom", strlen("-max_err")) == 0 && i < argc-1)		{ MAX_ERR_HOM = atoi(argv[++i]); }
		else if( strncmp(argv[i], "-err_het", strlen("-max_err")) == 0 && i < argc-1)		{ MAX_ERR_HET = atoi(argv[++i]); }
		else if( strncmp(argv[i], "-from_snp", strlen("-from_snp")) == 0 && i < argc-1 )	rs_range[0] = argv[++i];
		else if( strncmp(argv[i], "-to_snp", strlen("-to_snp")) == 0 && i < argc-1 )		rs_range[1] = argv[++i];
		else if( strncmp(argv[i], "-haps", strlen("-haps")) == 0 )							PRINT_HAPS = true;
		else if( strncmp(argv[i], "-print", strlen("-print")) == 0 )						PRINT_MATCH_HAPS = true;
		else if( strncmp(argv[i], "-map", strlen("-map")) == 0 && i < argc-1)				map = argv[++i];
		else if( strncmp(argv[i], "-bits", strlen("-bits")) == 0 && i < argc-1)				MARKER_SET_SIZE = atoi(argv[++i]);
		else if( strncmp(argv[i], "-homoz", strlen("-homoz")) == 0 )						ALLOW_HOM = true;
		else if( strncmp(argv[i], "-h_extend", strlen("-h_extend")) == 0 )					HAP_EXT = true;
		else if( strncmp(argv[i], "-w_extend", strlen("-w_extend")) == 0 )					WIN_EXT = true;
		else bad_param = true;
	}

	if(MIN_MATCH_LEN < 0)
	{
		cerr << "-min_m must be non-negative" << endl << endl;
		bad_param = true;
	} else if(MAX_ERR_HOM < 0 || MAX_ERR_HET < 0 )
	{
		cerr << "-err_hom,-err_het must be non-negative" << endl << endl;
		bad_param = true;
	}

	if(bad_param)
	{
		cerr << "usage: " << argv[0] << "<flags (optional)>" << endl
		<< "flags:" << endl
		<< '\t' << "-min_m" << '\t' << "Minimum length for match to be used for imputation (in cM or MB)." << endl
		<< '\t' << "-err_hom" << '\t' << "Maximum number of mismatching homozygous markers (per slice)." << endl
		<< '\t' << "-err_het" << '\t' << "Maximum number of mismatching heterozygous markers (per slice)." << endl
		<< '\t' << "-from_snp" << '\t' << "Start SNP (rsID)." << endl
		<< '\t' << "-to_snp" << '\t' << "End SNP (rsID)." << endl
		<< '\t' << "-haps" << '\t' << "Print the resolved haplotypes in a seperate HAPS file." << endl
		<< '\t' << "-print" << '\t' << "Print the sequence information for matches." << endl
		<< '\t' << "-map" << '\t' << "Genetic distance map." << endl
		<< '\t' << "-bits" << '\t' << "Slice size." << endl
		<< '\t' << "-homoz" << '\t' << "Allow self matches (homozygosity)" << endl
		<< '\t' << "-h_extend" << '\t' << "Extend from seeds if *haplotypes* match" << endl
		<< '\t' << "-w_extend" << '\t' << "Extend, one marker at a time, beyong the boundaries of a found match" << endl;
		return 0;
	}

	if( rs_range[0] != "" && rs_range[1] != "" )
	{
		ROI = true;
		snps.setROI(rs_range);
	}

	if(map != "")
	{
		snps.loadGeneticDistanceMap( map );
	}

	GERMLINE germline;
    germline.mine();
    return 1;
}


// end GERMLINE_0001.cpp
