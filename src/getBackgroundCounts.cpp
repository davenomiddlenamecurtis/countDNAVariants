#include "getBackgroundCounts.hpp"
#include "getSequenceFromReference.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#if 0

Code to count number of times each pentanucleotide background sequence occurs in the human reference.
Reference sequences is in files in the working directory named CHR1.FA, CHR2.FA ... CHRX.FA. 
These are files in which all characters have been converted to upper case.
Takes one single argument for the chromosome but if this is omitted will use chromosome 22.

#endif

baseHasher hasher;

const char* chr = "22";

int backgroundCountTable[1024];

#define MAXGAPS 200

const char* centromeresFn = "centromeres.hg38.txt";
const char* gapsFn = "gaps.hg38.txt";
int centromere[2], gaps[MAXGAPS][2],nGaps;

// background sequence consists of 5 bases

char outputSequence[6];

int outputBackgroundCounts(FILE* fo, int level,int *counts) {
	int i;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[4 - level] = baseNames[i];
			fprintf(fo, "%s\t%d\n", outputSequence, counts[hasher.hashBases(outputSequence, 5)]);
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[4 - level] = baseNames[i];
			outputBackgroundCounts(fo, level - 1,counts);
		}
		return 1;
	}
	return 0;
}

int getCentromeresAndGaps(const char* chr)
{
	char chrName[20],line[1000],chrom[100],type[100];
	int st = 300000000, en = 0,s,e,i;
	FILE* fi;
	sprintf(chrName, "chr%s", chr);
	fi = fopen(centromeresFn, "r");
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %d %d", chrom, &s, &e);
		if (!strcmp(chrom, chrName)) {
			if (s < st)
				st = s;
			if (e > en)
				en = e;
		}
	}
	fclose(fi);
	centromere[0] = st;
	centromere[1] = en;
	fi = fopen(gapsFn, "r");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%*s %s %d %d %*s %*s %*s %s", chrom, &s, &e,type);
		if (!strcmp(chrom, chrName)) {
			if (strcmp(type, "telomere")) // can add others later
				continue;
			gaps[i][0] = s;
			gaps[i][1] = e;
			++i;
		}
	}
	fclose(fi);
	nGaps = i;
	return nGaps;
}

int main(int argc, char* argv[])
{
	int i, c, tot, ePos, g, gs, cs,lostPos;
	unsigned int hash;
	char fn[100];
	const char* ptr;
	FILE* fi, * fo;
	faSequenceFile fa;
	chr = argv[1];
	getCentromeresAndGaps(chr);
	sprintf(fn, "CHR%s.FA", chr);
	fa.init(fn);
	fi = fa.getFp();
	while ((c = fgetc(fi)) != EOF && (c == 'N' || isspace(c) || (ptr = strchr(baseNames, c)) == 0))
		;
	ePos = fa.getPos();
	if (ePos <= gaps[0][1]) {
		fa.getSequence(outputSequence, gaps[0][1], 5);
	}
	else {
		outputSequence[0] = c;
		for (i = 1; i < 5; ++i)
			outputSequence[i] = fgetc(fi);
	}
	g = 1;
	gs = gaps[g][0];
	cs = centromere[0];
	// now have first lot of five letters
	hash = hasher.hashBases(outputSequence, 5);
	++backgroundCountTable[hash];
	tot = 1;
	ePos = fa.getPos();
	lostPos = 0;
	while ((c = fgetc(fi))!=EOF) {
		if (isspace(c))
			continue;
		if (c == 'N' || (ptr = strchr(baseNames, c)) == 0) {
			lostPos = 1;
			continue;
		}
		if (lostPos) {
			ePos = fa.getPos(); // get back on track
			fa.getSequence(outputSequence, ePos, 5); // have to hope this is a kosher sequence
			hash = hasher.hashBases(outputSequence, 5);
			++backgroundCountTable[hash];
			ePos = fa.getPos();
			lostPos = 0;
		}
		++ePos;
		if (ePos == gs || ePos == cs) {
			if (ePos == gs) {
				fa.getSequence(outputSequence, gaps[g][1] + 1, 5);
				++g;
				gs = gaps[g][0];
			}
			else {
				fa.getSequence(outputSequence, centromere[1] + 1, 5);
				for (; gaps[g][0] < centromere[1]; ++g)
					;
				gs = gaps[g][0];

			}
			hash = hasher.hashBases(outputSequence, 5);
			++backgroundCountTable[hash];
			ePos = fa.getPos();
			continue;
		}
		hash &= 255;
		hash <<= 2;
		hash |= (ptr - baseNames);
		++backgroundCountTable[hash];
		++tot;
	}
	sprintf(fn, "backgroundCounts.%s.txt", chr);
	fo = fopen(fn, "w");
	outputBackgroundCounts(fo, 4,backgroundCountTable);
	fclose(fo);
	printf("Dealt with %d bases\n", tot);
	return 0;
}
