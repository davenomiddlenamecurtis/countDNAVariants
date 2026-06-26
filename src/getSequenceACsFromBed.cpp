#include "countVariantsInContext.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>


#define MAXLEN 7
#define FOURTOTHEMAXLEN 16384

int len = 5;

const char* chr;
#define MAXAC 200
int acBinTable[FOURTOTHEMAXLEN*4][MAXAC + 2], acTriosBinTable[256][MAXAC + 2];
char outputSequence[MAXLEN+2],trioNames[256][5];
int backgroundCountTable[FOURTOTHEMAXLEN],backgroundTrioCountTable[256];


baseHasher hasher;

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
							if (outputSequence[5] == outputSequence[2])
								continue;
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

int outputLenBinCounts(FILE* fb, int level) {
// recursive so any length but harder to control order
	int i,f,s;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[len - level] = baseNames[i];
			if (outputSequence[len] == outputSequence[(len-1)/2]) // len must be odd
				continue;
			s = hasher.hashBases(outputSequence, len+1);
			fprintf(fb, "%s", outputSequence);
			for (f = 1; f <= MAXAC + 1; ++f)
				fprintf(fb, "\t%d", acBinTable[s][f]);
			fprintf(fb, "\n");
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[len - level] = baseNames[i]; // was len-1 for background counts, but variant has length len+1
			outputLenBinCounts(fb, level - 1);
		}
		return 1;
	}
	return 0;
}

int main(int argc, char* argv[])
{
	int i,start,end,AC,b,c,f,count,halfLen,allDone,st,en,lastEn;
	unsigned int hash,bgHash,trioHash;
	FILE* fi,*fo,*ff,*fb;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100],acBinFileRoot[100],compSequence[MAXLEN+2], outputFileRoot[100], * filenameSpec;
	char* posPtr, * ptr, line[10000], sequence[MAXLEN+2], fn[100],trioSeq[5],bedLine[200],bedChr[10];
	chr = argv[1];
	fb = fopen(argv[2], "r");
	filenameSpec = argv[3];
	sprintf(outputFileRoot, "acBinCounts.%s", filenameSpec);
	if (argc > 4)
		len = atoi(argv[4]);
	if (!strcmp(chr, "all")) {
		for (c = 0; c < 23; ++c) {
			if (c == 22)
				sprintf(fn, "%s.X.txt", outputFileRoot);
			else
				sprintf(fn, "%s.%d.txt", outputFileRoot, c + 1);
			fi = fopen(fn, "r");
			while (fgets(line, 9999, fi) && sscanf(line, "%s", sequence) == 1) {
				hash = hasher.hashBases(sequence,len+1);
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
		sprintf(fn, "%s.total.txt", outputFileRoot);
		fo = fopen(fn, "w");
		fprintf(fo, "Variant");
		for (f = 1; f <= MAXAC; ++f)
			fprintf(fo, "\tAC%d", f);
		fprintf(fo, "\tACrest\n");
		if (len == 5)		
			outputBinCounts(fo);
		else
			outputLenBinCounts(fo,len);
		fclose(fo);
		sprintf(fn, "backgroundCounts.%s.all.txt", filenameSpec);
		fi = fopen(fn, "r");
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %d", sequence, &count);
			hash = hasher.hashBases(sequence, len);
			backgroundCountTable[hash] = count;
			if (len == 5) {
				if (sequence[2] == 'C' || sequence[2] == 'T') {
					trioHash = hasher.hashBases(sequence + 1, 3);
				}
				else {
					hasher.getComplement(compSequence, sequence + 1, 3);
					trioHash = hasher.hashBases(compSequence, 3);
				}
				backgroundTrioCountTable[trioHash] += count;
			}
		}
		fclose(fi);
		sprintf(fn, "%s.total.txt", outputFileRoot);
		fi = fopen(fn, "r");
		fgets(line, 9999, fi);
		sprintf(fn, "acBinFrequencies.%s.overall.txt", filenameSpec);
		fo = fopen(fn, "w");
		fprintf(fo, "Variant");
		for (f = 1; f <= MAXAC; ++f)
			fprintf(fo, "\tAC%d", f);
		fprintf(fo, "\tACrest\n");
		while (fgets(line, 9999, fi) && sscanf(line, "%s", sequence) == 1) {
			if (sequence[(len-1)/2] == sequence[len])
				continue;
			hash = hasher.hashBases(sequence, len+1);
			bgHash = hasher.hashBases(sequence, len);
			fprintf(fo, "%s", sequence);
			for (f=1;f<=MAXAC+1;++f)
				fprintf(fo, "\t%.8f", acBinTable[hash][f] / (float)backgroundCountTable[bgHash]);
			fprintf(fo, "\n");
			if (len == 5) 
				sprintf(trioSeq, "%c%c%c%c", sequence[1], sequence[2], sequence[3], sequence[5]);
			else
				sprintf(trioSeq, "%c%c%c%c", sequence[(len-1)/2-1], sequence[(len - 1) / 2], sequence[(len - 1) / 2+1], sequence[len]);
			if (trioSeq[1] == 'C' || trioSeq[1] == 'T') {
					trioHash = hasher.hashBases(trioSeq, 4);
					strcpy(trioNames[trioHash], trioSeq);
				}
				else {
					hasher.getComplement(compSequence, trioSeq, 3);
					hasher.getComplement(compSequence + 3, trioSeq + 3, 1);
					trioHash = hasher.hashBases(compSequence, 4);
				}
			for (f = 1; f <= MAXAC + 1; ++f)
				acTriosBinTable[trioHash][f] += acBinTable[hash][f];
		}
		fclose(fo);
		if (len == 5) {
			sprintf(fn, "acTrioBinCounts.%s.total.txt", filenameSpec);
			fo = fopen(fn, "w");
			fprintf(fo, "Variant");
			for (f = 1; f <= MAXAC; ++f)
				fprintf(fo, "\tAC%d", f);
			fprintf(fo, "\tACrest\n");
			sprintf(fn, "acTrioBinFrequencies.%s.overall.txt", filenameSpec);
			ff = fopen(fn, "w");
			fprintf(ff, "Variant");
			for (f = 1; f <= MAXAC; ++f)
				fprintf(ff, "\tAC%d", f);
			fprintf(ff, "\tACrest\n");
			for (i = 0; i < 256; ++i) {
				if (trioNames[i][0] == '\0')
					continue;
				fprintf(fo, "%s", trioNames[i]);
				fprintf(ff, "%s", trioNames[i]);
				trioHash = hasher.hashBases(trioNames[i], 4);
				bgHash = hasher.hashBases(trioNames[i], 3);
				for (f = 1; f <= MAXAC + 1; ++f) {
					fprintf(fo, "\t%d", acTriosBinTable[trioHash][f]);
					fprintf(ff, "\t%.8f", acTriosBinTable[trioHash][f] / (float)backgroundTrioCountTable[bgHash]);
				}
				fprintf(fo, "\n");
				fprintf(ff, "\n");
			}
			fclose(fo);
			fclose(ff);
		}
		return 0;
	}
	sprintf(inFile, "ukb24308_c%s_b0_v1.head.pvar", chr);
		sprintf(acBinFileRoot, "%s.%s",outputFileRoot, chr);
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
	while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && strcmp(bedChr, chr))
		;
	do {
		sscanf(bedLine, "%*s %d %d", &st, &en);
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
			start = atoi(posPtr) - halfLen;
		}
		if (allDone)
			break;
		do {
			fa.getSequence(sequence, start, len);
			sequence[len] = ptr[2];
			hash = hasher.hashBases(sequence, len + 1);

			if (sequence[2] == sequence[len] || hash == 0)
				printf("Problem with sequence %s at %d, hash = %d\n", sequence, start, hash);
			
			ptr += 7;
			// 22 10510001 DRAGEN:chr22:10510001:G:T G T 0 LowGTR AC=2;AN=17868;AF=0.000111932;NS=490541;NS_GT=8934;NS_NOGT=7484;NS_NODATA=474123
			while (*ptr++ != 'A' || *ptr++ != 'C') // look for AC= string in line
				;
			AC = atoi(ptr + 1);
			if (AC > MAXAC)
				++acBinTable[hash][MAXAC + 1];
			else
				++acBinTable[hash][AC];
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
			start = atoi(posPtr) - halfLen;
			end = start+len-1;
		} while (end <= en && !allDone);

	} while (fgets(bedLine, 199, fb) && sscanf(bedLine + 3, "%s", bedChr) && !strcmp(bedChr, chr) && ! allDone);
	fclose(fb);
	fclose(fi);
	sprintf(fn, "%s.txt", acBinFileRoot);
	fo = fopen(fn, "w");
	if (len == 5) {
		outputBinCounts(fo);
	} else {
		outputLenBinCounts(fo,len);
	}
	fclose(fo);
	
return 0;
}
