GERMLINE
--------
http://www1.cs.columbia.edu/~gusev/germline/


> version1.4.1
> fix bug and add new feature
change detail to check [detail](detail.md)

# new feature
allow only test between first individual to other | use  `-poi` (person of interest)

allow missing SNP in ped file | use `-hg` (have gap)


# About
GERMLINE is a program for discovering long shared segments of Identity by Descent (IBD) between pairs of individuals in a large population. It takes as input genotype or haplotype marker data for individuals (as well as an optional known pedigree) and generates a list of all pairwise segmental sharing.

GERMLINE uses a novel hashing & extension algorithm which allows for segment identification in haplotype data in time proportional to the number of individuals. With genotype data, GERMLINE implements several pedigree-based phasing techniques to impute data for related individuals, and then iteratively uses identified IBD segments to infer additional missing information. GERMLINE can identify shared segments of any specified length, as well as allow for any number of mismatching markers.

The program has been developed in Itsik Pe'er's Lab of Computational Genetics at Columbia University. It is built in C++ and tested in the Red Hat Linux environment; the source is distributed here in a tar.gz package under the GPL license. 

# Usage
Compile germline with make.
The executable is run as `germline <options>` which prompts the user for input/output file information and runs the algorithm.

# Input
GERMLINE accepts as input the following formats:

    * Plink / ped+map
    * PHASE / HapMap

GERMLINE also accepts an optional genetic map, formatted according to the HapMap standard described here ( http://ftp.hapmap.org/recombination/ ).

# Output
Upon completion, GERMLINE generates a .match file in the specified location. The first five lines contain meta-data detailing the run settings and executions time. The following rows detail the identified pairwise shared segments, one per row, with each row containing the following fields:

    * Family ID 1
    * Individual ID 1
    * Family ID 2
    * Individual ID 2
    * Segment start (bp)
    * Segment end (bp)
    * Mismatching SNPs in segment
    * Total SNPs in segment
    * Genetic Length of segment(cM)

