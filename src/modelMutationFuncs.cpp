#include "modelMutationFuncs.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hashBases.hpp"
#include "runModels.hpp"

#if 0
in general, each row will have a sequence and a count 

4 REFs

4 ALTs

4x3 MUTs

For each MUT:

4 U2s
4 D2s

16 U3s
16 F3s
16 D3s

64 U4s
64 D4s

256 F5s

maybe first use F3s, then U3s, then D3s, then all 3s

#endif
char chrNames[23][3], sequenceTable[NSEQUENCES][7], componentNames[NCOMPONENTS][9], componentSequences[NCOMPONENTS][7], backgroundSequenceTable[NBACKGROUND][6];
float sequenceCountTable[NSEQUENCES][2],backgroundCountTable[NBACKGROUND],xy[2][NCOMPONENTS],meanBackgroundCount;
baseHasher hasher;

// produce variants in same order as is output by SigProfilerAssignment (FFS)
int generateSBSVariant(char* s, int start)
{
	static int i, j, k, l;
	static const char* CT = "CT", * target[2] = { "AGT","ACG" }, * ACGT = "ACGT";
	if (start == 1) {
		i = j = k = 0;
		l = -1;
		return 1;
	}

	++l;
	if (l > 3) {
		l = 0;
		++k;
		if (k > 3) {
			k = 0;
			++j;
			if (j > 2) {
				j = 0;
				++i;
				if (i > 1)
					return 0;
			}
		}
	}
	sprintf(s, "_%c%c%c_%c", ACGT[k], CT[i], ACGT[l], target[i][j]);
	return 1;
}

int matches(char* seq, char* comp, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		if (comp[i] != '_' && comp[i] != seq[i])
			return 0;
	return 1;
}

double correl(float* x, float* y, int count)
{
	double xbar, ybar, bxy, byx, sigmaxy = 0, sigmax2 = 0, sigmay2 = 0, sigmax = 0, sigmay = 0, r;
	int i;
	for (i = 0; i < count; ++i)
	{
		sigmax += x[i];
		sigmay += y[i];
		sigmaxy += x[i] * y[i];
	}
	xbar = sigmax / count;
	ybar = sigmay / count;
	byx = bxy = 0;
	for (i = 0; i < count; ++i)
	{
		bxy += (x[i] - xbar) * (y[i] - ybar);
		sigmax2 += x[i] * x[i];
		sigmay2 += y[i] * y[i];
	}
	byx = bxy;
	bxy /= sigmax2 - sigmax * sigmax / count;
	byx /= sigmay2 - sigmay * sigmay / count;
	r = sqrt(bxy * byx) * (bxy < 0 ? -1 : 1);
	return r;
}

int getPredictedFreqs(FILE* flog, const char* compBetaFileName, const char* observedFreqsName, const char* outputFileName) {
	FILE* fi, * fo;
	char line[1000], compName[100];
	int i, c, nComps;
	unsigned int hash;
	float f1, * betas, beta, SE, z, predFreq, intBeta;
	betas = (float*)calloc(NCOMPONENTS, sizeof(float));
	fi = fopen(compBetaFileName, "r");
	i = 0;
	intBeta = 0;
	while (fgets(line, 999, fi)) {
		if (sscanf(line, "%s %f %f %f", compName, &beta, &SE, &z) == 4) {
			if (strchr("MUFD", compName[0])) {
				strcpy(componentSequences[i], compName + 2);
				betas[i] = beta;
				++i;
			}
			else if (!strcmp(compName, "Intercept"))
				intBeta = beta;
		}
	}
	fclose(fi);
	nComps = i;
	fi = fopen(observedFreqsName, "r");
	fgets(line, 999, fi);
	fo = fopen(outputFileName, "w");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %f", sequenceTable[i], &f1); // no real need to keep this in a table, but leave it for now
		if (sequenceTable[i][5] == sequenceTable[i][2])
			continue;
		for (beta = 0, c = 0; c < nComps; ++c) {
			if (matches(sequenceTable[i], componentSequences[c], 6))
				beta += betas[c];
		}
		beta += intBeta;
		predFreq = exp(beta) / (exp(beta) + 1);
		fprintf(fo, "%s\t%.8f\t%.8f\n", sequenceTable[i], f1, predFreq);
		xy[0][i] = predFreq;
		xy[1][i] = f1;
		++i;
	}
	fclose(fi);
	fclose(fo);
	printf("Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], i));
	fprintf(flog, "Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], i));
	free(betas);
	return 1;
}

int getPredictedCounts(FILE* flog, const char* compBetaFileName, const char* observedSeqsName, const char* outputFileName) {
	FILE* fi, * fo;
	char line[1000], compName[100];
	int i, c, nComps;
	unsigned int hash;
	float c1, * betas, beta, SE, z, predCount, intBeta;
	betas = (float*)calloc(NCOMPONENTS, sizeof(float));
	fi = fopen(compBetaFileName, "r");
	i = 0;
	intBeta = 0;
	while (fgets(line, 999, fi)) {
		if (sscanf(line, "%s %f %f %f", compName, &beta, &SE, &z) == 4) {
			if (strchr("MUFD", compName[0])) {
				strcpy(componentSequences[i], compName + 2);
				betas[i] = beta;
				++i;
			}
			else if (!strcmp(compName, "Intercept"))
				intBeta = beta;
		}
	}
	fclose(fi);
	nComps = i;
	fi = fopen(observedSeqsName, "r");
	fo = fopen(outputFileName, "w");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
		if (sequenceTable[i][5] == sequenceTable[i][2])
			continue;
		for (beta = 0, c = 0; c < nComps; ++c) {
			if (matches(sequenceTable[i], componentSequences[c], 6))
				beta += betas[c];
		}
		beta += intBeta;
		hash = hasher.hashBases(sequenceTable[i], 5); // just for background
		predCount = exp(beta) / (exp(beta) + 1) * meanBackgroundCount;
		//		predCount = exp(beta) / (exp(beta) + 1) * backgroundCountTable[hash];
		fprintf(fo, "%s\t%d\t%.2f\n", sequenceTable[i], (int)c1, predCount);
		xy[0][i] = predCount;
		xy[1][i] = c1;
		++i;
	}
	fclose(fi);
	fclose(fo);
	printf("Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], i));
	fprintf(flog, "Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], i));
	free(betas);
	return 1;
}
