#include "countVariantsInContext.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

baseHasher hasher;

#define MAXGAPS 200

const char* centromeresFn = "centromeres.hg38.txt";
const char* gapsFn = "gaps.hg38.txt";
int centromere[2], gaps[MAXGAPS][2], nGaps;
int getCentromeresAndGaps(const char* chr)
{
	char chrName[20], line[1000], chrom[100], type[100];
	int st = 300000000, en = 0, s, e, i;
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
		sscanf(line, "%*s %s %d %d %*s %*s %*s %s", chrom, &s, &e, type);
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

char chr[10];
int countTable[3][4096];

char outputSequence[7];

int outputCounts(FILE* fo, int level,int *counts) {
	int i;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[5 - level] = baseNames[i];
			fprintf(fo, "%s\t%d\n", outputSequence, counts[hasher.hashBases(outputSequence, 6)]);
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[5 - level] = baseNames[i];
			outputCounts(fo, level - 1,counts);
		}
		return 1;
	}
	return 0;
}


int main(int argc, char* argv[])
{
	int i,start,AC;
	unsigned int hash;
	FILE* fi,*fo,*fr;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100],chrName[100],pos[100],ref[100],alt[100];
	sprintf(inFile, "aau1043_datas5_revision1.tsv");
	sprintf(outFileRoot, "counts.aau1043_datas5_revision1");
	fi = fopen(inFile, "r");
	char* posPtr, * ptr, line[1000],sequence[6],fn[100];
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %s %s %s", chrName, pos, ref, alt);
		if (strncmp(chrName, "chr", 3))
			continue;
		if (strcmp(chrName + 3, chr)) {
			strcpy(chr, chrName + 3);
			getCentromeresAndGaps(chr); 
			fa.close();
			sprintf(refFile, "CHR%s.FA", chr);
			fa.init(refFile);
		}
		if (ref[1] || alt[1])
			continue; // only use SNVs
		start = atoi(pos) - 2;
		if (start <= gaps[0][1]-2 || (start >= centromere[0]-2 && start <= centromere[1]-2))
			continue;
		if (start >= gaps[1][0]-2)
			break;
		fa.getSequence(sequence, start, 5);
		if (sequence[0] == 'N')
			continue;
		sequence[5] = alt[0];
		hash=hasher.hashBases(sequence,6);
		++countTable[0][hash];
	}
	fclose(fi);
		sprintf(fn, "%s.txt", outFileRoot);
		fo = fopen(fn, "w");
		outputCounts(fo, 5,countTable[0]);
		fclose(fo);
	return 0;
}
