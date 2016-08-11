CC=g++
OPT=-O3 -I include
SRCS=GERMLINE_0001.cpp GERMLINE.cpp Share.cpp Chromosome.cpp ChromosomePair.cpp HMIndividualsExtractor.cpp MarkerSet.cpp Individual.cpp Individuals.cpp InputManager.cpp MatchFactory.cpp MatchesBuilder.cpp NucleotideMap.cpp PEDIndividualsExtractor.cpp Match.cpp PolymorphicIndividualsExtractor.cpp SNP.cpp SNPPositionMap.cpp SNPs.cpp
OBJS=GERMLINE_0001.o GERMLINE.o Chromosome.o Share.o ChromosomePair.o HMIndividualsExtractor.o MarkerSet.o Individual.o Individuals.o InputManager.o MatchFactory.o MatchesBuilder.o NucleotideMap.o PEDIndividualsExtractor.o Match.o PolymorphicIndividualsExtractor.o SNP.o SNPPositionMap.o SNPs.o
MAIN=germline

all: clean germline test_case

germline: $(OBJS)
	$(CC) $(OPT) -o $(MAIN) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(OPT) -c $*.cpp
clean:
	-rm -f *.o $(MAIN) test/generated.match test/generated.log test/generated.err
test_case:
	-@./$(MAIN) -bits 50 -min_m 1 -err_hom 2 -err_het 0 < test/test.run 2> test/generated.err | echo "Running Test Case"
	diff -q -s test/expected.match test/generated.match
