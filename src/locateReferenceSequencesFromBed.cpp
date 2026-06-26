#include "getBackgroundCounts.hpp"
#include "getSequenceFromReference.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


#define MAXLEN 7
#define FOURTOTHEMAXLEN 16384
#define MAXSEQ 100

int len;
baseHasher hasher;

const char* chr;
char outputSequence[MAXLEN+2];
int backgroundCountTable[FOURTOTHEMAXLEN];

// chr bedFile seqsFile (len)
int main(int argc, char* argv[])
{
	int i,start,end,pos,AC,b,c,f,count,halfLen,allDone,st,en,lastEn,nSeqToFind,nSeq,s;
	FILE* fi, * fo[MAXSEQ], * ft[MAXSEQ], * ff, * fb;
	unsigned hash, hashesToFind[MAXSEQ],mask,nHashes;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100],acBinFileRoot[100],compSequence[MAXLEN+2],transcript[100],strand[100],gene[100];
	char* posPtr, * ptr, line[10000], sequence[MAXLEN+2], seqsToFind[MAXSEQ][MAXLEN+2], fn[100],bedLine[200],bedChr[10];
	chr = argv[1];
	fb = fopen(argv[2], "r");
	fi = fopen(argv[3], "r");
	nSeqToFind = 0;
	while (fgets(line, 99, fi) && sscanf(line, "%s", seqsToFind[nSeqToFind]) == 1) {
		len = strlen(seqsToFind[nSeqToFind]); // must be same for all
		sprintf(fn, "transcripts.baackground.%s.%s.txt", seqsToFind[nSeqToFind], chr);
		ft[nSeqToFind] = fopen(fn, "w");
		hashesToFind[nSeqToFind]=hasher.hashBases(seqsToFind[nSeqToFind], len);
		++nSeqToFind;
	}
	fclose(fi);
	mask = 0;
	nHashes = pow(4, len);
	for (i = 0; i < len - 1; ++i) {
		mask <<= 2;
		mask |= 3;
	}
	sprintf(fn, "CHR%s.FA", chr);
	fa.init(fn);
	fi = fa.getFp();
	gene[0] = '\0';
	lastEn = 0;
	while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && strcmp(bedChr, chr))
		;
	do {
		sscanf(bedLine, "%*s %d %d %s %s %s", &st, &en,transcript, strand, gene);
		if (st < lastEn)
			continue; // ignore repeated and overlapping exons
		lastEn = en;
		nSeq = en - st - len;
		fa.getSequence(outputSequence, st + 1, len);
		hash = hasher.hashBases(outputSequence, len);
		for (i = 0; i < nHashes; ++i)
			backgroundCountTable[i] = 0;
		for (s = 0; s < nSeq; ++s) {
			++backgroundCountTable[hash];
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
		for (i = 0; i < nSeqToFind; ++i)
			fprintf(ft[i], "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%.8f\n",
				chr, st, en, transcript, strand, gene, backgroundCountTable[hashesToFind[i]], backgroundCountTable[hashesToFind[i]] / (float)(en - st));


	} while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && !strcmp(bedChr, chr));
	fclose(fb);
	for (i = 0; i < nSeqToFind; ++i)
			fclose(ft[i]);

return 0;
}
