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
Takes first argument for chromosome, second for bed file. 

Default length is 5 but a different value can be provided as third argument(may need to be an odd nnumber).

#endif

#define MAXLEN 7
#define FOURTOTHEMAXLEN 16384

baseHasher hasher;

const char* chr;

int backgroundCountTable[FOURTOTHEMAXLEN],len=5;

// background sequence consists of len bases

char outputSequence[MAXLEN+1];

int outputBackgroundCounts(FILE* fo, int level,int *counts) {
	int i;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[len-1 - level] = baseNames[i];
			fprintf(fo, "%s\t%d\n", outputSequence, counts[hasher.hashBases(outputSequence, len)]);
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[len - 1 - level] = baseNames[i];
			outputBackgroundCounts(fo, level - 1,counts);
		}
		return 1;
	}
	return 0;
}

int main(int argc, char* argv[])
{
	int i, c, tot, g, gs, cs,st,en,nSeq,s,lastEn,count;
	unsigned int hash,mask;
	char fn[100],line[200],bedChr[100],*outputFileRoot;
	const char* ptr;
	FILE* fi, * fo, *fb;
	faSequenceFile fa;
	chr = argv[1];
	fb = fopen(argv[2], "r");
	outputFileRoot = argv[3];
	if (argc > 4)
		len = atoi(argv[4]);
	if (!strcmp(chr, "all"))
	{
		sprintf(fn, "%s.all.txt", outputFileRoot);
		fo = fopen(fn, "w");
		for (c = 0; c < 23; ++c) {
			if (c==22)
				sprintf(fn, "%s.X.txt", outputFileRoot);
			else
				sprintf(fn, "%s.%d.txt", outputFileRoot,c+1);
			fi = fopen(fn, "r");
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %d", outputSequence, &count);
				hash = hasher.hashBases(outputSequence, len);
				backgroundCountTable[hash] += count;
				if (c == 22)
					fprintf(fo, "%s\t%d\n", outputSequence, backgroundCountTable[hash]);
			}
			fclose(fi);
		}
		fclose(fo);
		return 0;
	}
	mask = 0;
	for (i = 0; i < len-1; ++i) {
		mask <<= 2;
		mask |= 3;
	}
	sprintf(fn, "CHR%s.FA", chr);
	fa.init(fn);
	fi = fa.getFp();
	while (fgets(line, 199, fb) && sscanf(line + 3, "%s", bedChr) && strcmp(bedChr, chr))
		;
	tot = 0;
	lastEn = 0;
	do {
		sscanf(line, "%*s %d %d", &st, &en);
		if (st < lastEn)
			continue; // ignore repeated and overlapping exons
		lastEn = en;
		nSeq = en - st - len;
		fa.getSequence(outputSequence, st+1, len);
		hash = hasher.hashBases(outputSequence, len);
		for (s = 0; s < nSeq; ++s) {
			++backgroundCountTable[hash];
			++tot;
			do {
				c = fgetc(fi);
				for (i = 0; i < 4; ++i) {
					if (baseNames[i] == c)
						break;
				}
			} while (i == 4); // skip over white space between lines
			hash &= mask;
			hash <<= 2;
			hash |= i;
		}
	} while (fgets(line, 199, fb) && sscanf(line + 3, "%s", bedChr) && ! strcmp(bedChr, chr));
	fclose(fb);
	sprintf(fn, "%s.%s.txt",outputFileRoot, chr);
	fo = fopen(fn, "w");
	outputBackgroundCounts(fo, len-1,backgroundCountTable);
	fclose(fo);
	printf("Dealt with %d bases\n", tot);
	return 0;
}
