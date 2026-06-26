#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "locateVariantsInContextFromBed.hpp"


#define MAXLEN 7
#define FOURTOTHEMAXLEN 16384
#define MAXSEQ 100

int len = 5;

const char* chr;
char outputSequence[MAXLEN+2];

// chr bedFile seqsFile (len)
int main(int argc, char* argv[])
{
	int i,start,end,pos,AC,b,c,f,count,halfLen,allDone,st,en,lastEn,nSeqToFind,nVar[MAXSEQ];
	FILE* fi, * fo[MAXSEQ], * ft[MAXSEQ], * ff, * fb;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100],acBinFileRoot[100],compSequence[MAXLEN+2],transcript[100],strand[100],gene[100];
	char* posPtr, * ptr, line[10000], sequence[MAXLEN+2], seqsToFind[MAXSEQ][MAXLEN+2], fn[100],bedLine[200],bedChr[10];
	chr = argv[1];
	fb = fopen(argv[2], "r");
	fi = fopen(argv[3], "r");
	nSeqToFind = 0;
	while (fgets(line, 99, fi) && sscanf(line, "%s", seqsToFind[nSeqToFind]) == 1) {
		sprintf(fn, "locs.%s.%s.txt", seqsToFind[nSeqToFind], chr);
		fo[nSeqToFind] = fopen(fn, "w");
		sprintf(fn, "transcripts.%s.%s.txt", seqsToFind[nSeqToFind], chr);
		ft[nSeqToFind] = fopen(fn, "w");
		++nSeqToFind;
	}
	fclose(fi);
	if (argc > 4)
		len = atoi(argv[3]);
	sprintf(inFile, "ukb24308_c%s_b0_v1.head.pvar", chr);
	sprintf(refFile, "CHR%s.FA", chr);
	fa.init(refFile);
	fi = fopen(inFile, "r");
	fgets(line, 999, fi);
	halfLen = (len - 1)/2;
	sequence[len+1] = '\0';
	ptr = 0; posPtr = 0;//just to avoid error code from over-anxious compiler

	while (fgets(line, 9999, fi)) {
		ptr = line;
		while (!isspace(*ptr++));
		while (isspace(*ptr++));
		posPtr = ptr - 1;
		for (i = 0; i < 3; ++i)
			while (*ptr++ != ':');
		if (ptr[1] == ':'&& isspace(ptr[3]))
			break; // only use SNVs
	}
	start = atoi(posPtr) - halfLen;
	allDone = 0;
	lastEn = 0;
	gene[0] = '\0';
	while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && strcmp(bedChr, chr))
		;
	do {
		for (i = 0; i < nSeqToFind; ++i)
			nVar[i]= 0; 
		sscanf(bedLine, "%*s %d %d %s %s %s", &st, &en,transcript, strand, gene);
		if (st < lastEn)
			continue; // ignore repeated and overlapping exons
		lastEn = en;
		while (start < st + 1) {
			if (!fgets(line, 999, fi)) {
				allDone = 1;
				break;
			}
			ptr = line;
			while (!isspace(*ptr++));
			while (isspace(*ptr++));
			posPtr = ptr - 1;
			for (i = 0; i < 3; ++i)
				while (*ptr++ != ':');
			if (ptr[1] != ':' || !isspace(ptr[3]))
				continue; // only use SNVs
			pos = atoi(posPtr);
			start = pos- halfLen;
		}
		if (allDone)
			break;
		do {
			fa.getSequence(sequence, start, len);
			sequence[len] = ptr[2];

			if (sequence[halfLen] == sequence[len])
				printf("Problem with sequence %s at %d\n", sequence, start);
			
			ptr += 7;
			// 22 10510001 DRAGEN:chr22:10510001:G:T G T 0 LowGTR AC=2;AN=17868;AF=0.000111932;NS=490541;NS_GT=8934;NS_NOGT=7484;NS_NODATA=474123
			while (*ptr++ != 'A' || *ptr++ != 'C') // look for AC= string in line
				;
			AC = atoi(ptr + 1);
			for (i = 0; i < nSeqToFind; ++i) {
				if (matches(sequence, seqsToFind[i], len + 1)) {
					fprintf(fo[i], "%s\t%d\t%c\t%c\t%s\t%d\t%d\t%d\n", 
						chr, pos, sequence[halfLen], sequence[len], sequence, AC,pos-st,en-pos);
					++nVar[i];
					break;
				}
			}
			do {
				if (!fgets(line, 999, fi)) {
					allDone = 1;
					break;
				}
				ptr = line;
				while (!isspace(*ptr++));
				while (isspace(*ptr++));
				posPtr = ptr - 1;
				for (i = 0; i < 3; ++i)
					while (*ptr++ != ':');
			}
			while (ptr[1] != ':' || !isspace(ptr[3])); // only use SNVs
			pos = atoi(posPtr);
			start = pos - halfLen;
			end = start+len-1;
		} while (end <= en && !allDone);
		for (i = 0; i < nSeqToFind; ++i)
			fprintf(ft[i], "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%.8f\n",
				chr,st,en,transcript,strand,gene,nVar[i],nVar[i]/(float)(en-st));
	} while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && !strcmp(bedChr, chr) && ! allDone);
	fclose(fb);
	for (i = 0; i < nSeqToFind; ++i) {
		fclose(fo[i]);
		fclose(ft[i]);
	}

return 0;
}
