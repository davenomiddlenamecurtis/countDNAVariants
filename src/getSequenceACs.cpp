#include "countVariantsInContext.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

const char* chr = "";
#define MAXAC 200
int acBinTable[4096][MAXAC + 2], acTriosBinTable[256][MAXAC + 2];
char outputSequence[7],trioNames[256][5];
int backgroundCountTable[1024],backgroundTrioCountTable[256];


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

int outputBinCounts(FILE* fb) {
	int f,b[6];
	unsigned int s;
	for (b[0] = 0; b[0] < 4; ++b[0]) {
		outputSequence[0] = baseNames[b[0]];
		for (b[1] = 0; b[1] < 4; ++b[1]) {
			outputSequence[1] = baseNames[b[1]];
			for (b[2] = 0; b[2] < 4; ++b[2]) {
				outputSequence[2] = baseNames[b[2]];
				for (b[3] = 0; b[3] < 4; ++b[3]) {
					outputSequence[3] = baseNames[b[3]];
					for (b[4] = 0; b[4] < 4; ++b[4]) {
						outputSequence[4] = baseNames[b[4]];
						for (b[5] = 0; b[5] < 4; ++b[5]) {
							outputSequence[5] = baseNames[b[5]];
							s = hasher.hashBases(outputSequence, 6);
							fprintf(fb, "%s", outputSequence);
							for (f = 1; f <= MAXAC + 1; ++f)
								fprintf(fb, "\t%d", acBinTable[s][f]);
							fprintf(fb, "\n");
						}
					}
				}
			}
		}
	}

	return 1;
}

int main(int argc, char* argv[])
{
	int i,start,AC,b,c,f,count;
	unsigned int hash,bgHash,trioHash;
	FILE* fi,*fo,*ff;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100],acBinFileRoot[100],compSequence[7];
	char* posPtr, * ptr, line[10000], sequence[7], fn[100],trioSeq[5];
	chr = argv[1];
	if (!strcmp(chr, "all")) {
		for (c = 0; c < 23; ++c) {
			if (c == 22) {
				sprintf(fn, "acBinCounts.X.txt");
			}
			else {
				sprintf(fn, "acBinCounts.%d.txt",c+1);
			}
			fi = fopen(fn, "r");
			while (fgets(line, 9999, fi) && sscanf(line, "%s", sequence) == 1) {
				hash = hasher.hashBases(sequence,6);
				ptr = line;
				for (f = 1; f <= MAXAC + 1; ++f) {
					while (!isspace(*ptr))
						++ptr;
					while (isspace(*ptr))
						++ptr;
					AC = atoi(ptr);
					acBinTable[hash][f] += AC;
				}
			}
		}
		fo = fopen("acBinCounts.total.txt", "w");
		outputBinCounts(fo);
		fclose(fo);
		fi = fopen("backgroundCounts.total.txt", "r");
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %d", sequence, &count);
			hash = hasher.hashBases(sequence, 5);
			backgroundCountTable[hash] = count;
			if (sequence[2] == 'C' || sequence[2] == 'T') {
				trioHash= hasher.hashBases(sequence+1, 3);
			}
			else {
				hasher.getComplement(compSequence, sequence + 1, 3);
				trioHash = hasher.hashBases(compSequence , 3);
			}
			backgroundTrioCountTable[trioHash] += count;
		}
		fclose(fi);
		fi = fopen("acBinCounts.total.txt", "r");
		fo = fopen("acBinFrequencies.overall.txt", "w");
		while (fgets(line, 9999, fi) && sscanf(line, "%s", sequence) == 1) {
			if (sequence[2] == sequence[5])
				continue;
			hash = hasher.hashBases(sequence, 6);
			bgHash = hasher.hashBases(sequence, 5);
			fprintf(fo, "%s", sequence);
			for (f=1;f<=MAXAC+1;++f)
				fprintf(fo, "\t%.8f", acBinTable[hash][f] / (float)backgroundCountTable[bgHash]);
			fprintf(fo, "\n");
			sprintf(trioSeq, "%c%c%c%c", sequence[1], sequence[2], sequence[3], sequence[5]);
			if (trioSeq[1] == 'C' || trioSeq[1] == 'T') {
				trioHash = hasher.hashBases(trioSeq, 4);
				strcpy(trioNames[trioHash], trioSeq);
			}
			else {
				hasher.getComplement(compSequence, trioSeq, 3);
				hasher.getComplement(compSequence+3, trioSeq+3, 1);
				trioHash = hasher.hashBases(compSequence, 4);
			}
			for (f = 1; f <= MAXAC + 1; ++f)
				acTriosBinTable[trioHash][f] += acBinTable[hash][f];
		}
		fclose(fo);
		fo = fopen("acTrioBinCounts.total.txt", "w");
		ff = fopen("acTrioBinFrequencies.overall.txt", "w");
		for (i = 0; i < 256; ++i) {
			if (trioNames[i][0] == '\0')
				continue;
			fprintf(fo, "%s", trioNames[i]);
			fprintf(ff, "%s", trioNames[i]);
			trioHash = hasher.hashBases(trioNames[i], 4);
			bgHash= hasher.hashBases(trioNames[i], 3);
			for (f = 1; f <= MAXAC + 1; ++f) {
				fprintf(fo, "\t%d", acTriosBinTable[trioHash][f]);
				fprintf(ff, "\t%.8f", acTriosBinTable[trioHash][f]/(float)backgroundTrioCountTable[bgHash]);
			}
			fprintf(fo, "\n");
			fprintf(ff, "\n");
		}
		fclose(fo);
		fclose(ff);

		return 0;
	}
	getCentromeresAndGaps(chr);
	sprintf(inFile, "ukb24308_c%s_b0_v1.head.pvar", chr);
	sprintf(acBinFileRoot, "acBinCounts.%s", chr);
	sprintf(refFile, "CHR%s.FA", chr);
	fa.init(refFile);
	fi = fopen(inFile, "r");
	fgets(line, 999, fi);
	while (fgets(line, 999, fi)) {
		ptr = line;
		while (!isspace(*ptr++));
		while (isspace(*ptr++));
		posPtr = ptr-1;
		for (i = 0; i < 3; ++i)
			while (*ptr++ != ':');
		if (ptr[1]!=':' || !isspace(ptr[3]))
			continue; // only use SNVs
		start = atoi(posPtr) - 2;
		if (start <= gaps[0][1]-2 || (start >= centromere[0]-2 && start <= centromere[1]-2))
			continue;
		if (start >= gaps[1][0]-2)
			break;
		fa.getSequence(sequence, start, 5);
		if (sequence[0] == 'N')
			continue;
		sequence[5] = ptr[2];
		hash=hasher.hashBases(sequence,6);
		ptr += 7;
		while (*ptr++ != 'A' || *ptr++ != 'C') // look for AC= string in line
			;
		AC = atoi(ptr + 1);
		if (AC>MAXAC)
			++acBinTable[hash][MAXAC+1];
		else
			++acBinTable[hash][AC];
// 22 10510001 DRAGEN:chr22:10510001:G:T G T 0 LowGTR AC=2;AN=17868;AF=0.000111932;NS=490541;NS_GT=8934;NS_NOGT=7484;NS_NODATA=474123
	}
	fclose(fi);
	sprintf(fn, "%s.txt", acBinFileRoot);
	fo = fopen(fn, "w");
	outputBinCounts(fo);
	fclose(fo);
	
return 0;
}
