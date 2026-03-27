#include "modelStrandBias.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "hashBases.hpp"
#include "runModels.hpp"

#if 0
in general, each row will have a sequence, a count and a probability
when I analyse corrected counts, the total count will be the number of times background occurred
(same for every sequence)
the names of component sequences will be generated, then betas will be 0 or 1 depending on whether they match full sequence
each sequence represents a complementary pair
#endif

#define MAXLEVELS 10 // for pseudo-recursion
// level will mean the opposite of what it did in outputCounts() - it means closer to the end
int generateSequence(char* s, int len,int start)
{
	static int level, i[MAXLEVELS],j;
	static char bases[MAXLEVELS + 1];
	if (start) {
		for (j = 0; j < len; ++j) {
			bases[j] = baseNames[0];
			i[j] = 0;
		}
		i[j-1] = -1;
		return 1;
	}
	else {
		if (i[len-1] < 3) {
			++i[len - 1];
			bases[len - 1] = baseNames[i[len - 1]];
			strcpy(s, bases);
			return 1;
		}
		else {
			for (j = len - 1; j > 0; --j)
				if (i[j] < 3)
					break;
			if (j == 0 && i[j] == 3)
				return 0;
			else {
				++i[j];
				bases[j] = baseNames[i[j]];
				for (++j; j < len; ++j) {
					i[j] = 0;
					bases[j] = baseNames[i[j]];
				}
				strcpy(s, bases);
				return 1;
			}
		}
	}
}

int writeRFiles(glfModel* m, char* root)
{
	char fn[100],line[10000];
	FILE* fo;
	int b,r,first;
	sprintf(fn, "%s.R.tab.txt", root);
	fo = fopen(fn, "w");
	for (b = 0; b < m->nCol; ++b)
		fprintf(fo, "%s\t",m->name[b]);
	fprintf(fo, "Y\tN\n");
	for (r = 0; r < m->nRow; ++r) {
		for (b = 0; b < m->nCol; ++b)
			fprintf(fo, "%f\t", m->X[r][b]);
		fprintf(fo, "0\t%d\n", (int)(m->F[r]- m->Y[r]* m->F[r]));
		for (b = 0; b < m->nCol; ++b)
			fprintf(fo, "%f\t", m->X[r][b]);
		fprintf(fo, "1\t%d\n", (int)(m->Y[r] * m->F[r]));
	}
	fclose(fo);
	sprintf(fn, "%s.R.formula.txt", root);
	fo = fopen(fn, "w");
	sprintf(line, "Y ~ ");
	first = 1;
	for (b = 0; b < m->nCol; ++b)
		if (m->toFit[b]) {
			sprintf(strchr(line, '\0'), "%s%s ", first ? "" : "+ ", m->name[b]);
				first = 0;
		}
	if (m->toFit[m->nCol] == 0)
		sprintf(strchr(line, '\0'), " - 1");
	fprintf(fo, "%s\n",line);
	fclose(fo);
	return 1;
}

// just add up counts to get averages for each beta
void setStartingBetasFromCounts(glfModel* m)
{
	float *N, *sigmaX,p,intOdds,SEp;
	int b, r;
	assert((sigmaX = (float*)calloc(m->nCol, sizeof(float))) != 0);
	assert((N = (float*)calloc(m->nCol, sizeof(float))) != 0);

	for (r = 0; r < m->nRow; ++r) {
		for (b = 0; b < m->nCol; ++b) 
			if (m->toFit[b]) {
				sigmaX[b] += m->X[r][b] * m->F[r]*m->Y[r]; // Y is a probability
				N[b] += m->X[r][b] * m->F[r];
			}
	}
	for (b = 0; b < m->nCol; ++b)
		if (m->toFit[b]) {
			p=sigmaX[b]/N[b];
			m->beta[b] = log(p / (1 - p));
			if (m->toUse[m->nCol]) {
				m->beta[b] -= m->beta[m->nCol]; // adjust for intercept odds
				// use formula from here: https://www.ncbi.nlm.nih.gov/books/NBK431098/
				// SE(beta)=sqrt(1/a+1/b!1/c*1/d)
				// as intercept has been provided, assume a and b are very large
				m->SE[b] = sqrt(1 / sigmaX[b] + 1 / (N[b] - sigmaX[b]));
			}
			else {
				SEp = sqrt(p * (1 - p) / N[b]);
				m->SE[b] = log((p + SEp) / (1 - (p + SEp))) - log(p / (1 - p));
				// get SE of p and use it to derive SE of odds
			}
		}
	free(sigmaX);
	free(N);
}

// produce variants in same order as is used by SigProfilerAssignment for input matrix (FFS)
int generateSBSVariantForInputMatrix(char* s, int start)
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
		if (k <= 3 && ACGT[k] == CT[j])
			++k;
		if (k > 3) {
			k = 0;
			++j;
			if (j > 1) { // C or T
				j = 0;
				++i;
				if (i > 3)
					return 0;
			}
		}
	}
	sprintf(s, "_%c%c%c_%c", ACGT[i], CT[j], ACGT[l], ACGT[k]);
	return 1;
}

int getBetas(double * beta, glfModel* m)
{
	double incBeta0=0;
	int b;
	if (m->isNormalised)
	{
		incBeta0 = 0;
		for (b = 0; b < m->nCol; ++b)
		{
			if (m->SD[b] && m->toUse[b])
			{
				incBeta0 -= m->beta[b] * m->mean[b] / m->SD[b];
				beta[b] = m->beta[b] / m->SD[b];
			}
			else
				beta[b] = 0;
		}
		beta[m->nCol] = m->beta[m->nCol] + incBeta0;
	}
	else
		for (b = 0; b < m->nCol + 1; ++b) {
			if (m->toUse[b])
				beta[b] = m->beta[b];
			else
				beta[b] = 0;
		}
	return 1;
}

void printModel(FILE* fo, const char* LLstr, double LL, glfModel* m,int clean)
{
	// change this to be row-wise and use names
	int b, bb, c;
	double* realBeta, * realSE;
	double incBeta0, incSEBeta0;
	realBeta = (double*)calloc(m->nCol + 1, sizeof(double));
	realSE = (double*)calloc(m->nCol + 1, sizeof(double));
	if (!clean) {
		fprintf(fo, "%s = %.2f\n", LLstr, LL);
		fprintf(fo, "var\tbeta\tvalue\tSE\tz\n");
	}
	if (m->isNormalised)
	{
		incSEBeta0 = incBeta0 = 0;
		for (c = 0; c < m->nCol; ++c)
		{
			if (m->SD[c] && m->toUse[c])
			{
				incBeta0 -= m->beta[c] * m->mean[c] / m->SD[c];
				incSEBeta0 -= m->SE[c] * m->mean[c] / m->SD[c];
				realBeta[c] = m->beta[c] / m->SD[c];
				realSE[c] = m->SE[c] / m->SD[c];
			}
		}
		realBeta[m->nCol] = m->beta[m->nCol] + incBeta0;
		realSE[m->nCol] = m->SE[m->nCol] + incSEBeta0;
	}
	else
		for (b = 0; b < m->nCol + 1; ++b)
		{
			realBeta[b] = m->beta[b];
			realSE[b] = m->SE[b];
		}
	for (b = 0; b < m->nCol + 1; ++b)
	{
		bb = (b + m->nCol) % (m->nCol + 1); // print last first
		if (m->toUse[bb])
			fprintf(fo, "%s\t%.5f\t%.5f\t%.5f\n", m->name[bb], realBeta[bb], realSE[bb], m->toFit[bb] ? realBeta[bb] / realSE[bb] : 0.0);
	}
	if (!clean) {
		fprintf(fo, "\n");
	}
	free(realBeta);
	free(realSE);
}

float evaluateModel(FILE* fo, glfModel* m, const char* name, int* useThese, int* fitThese, double *startingBetas,enum thingToMaximise useLinearRegression)
{
	float lnL;
	int b;
	m->setFunc(getMinusModelLnL);
	for (b = 0; b < m->nCol+1; ++b)
	{
		m->toUse[b] = useThese?useThese[b]:1;
		m->toFit[b] = fitThese?fitThese[b]:(useThese ? useThese[b] : 1);
	}
	if (startingBetas) {
		for (b = 0; b < m->nCol + 1; ++b)
			m->beta[b] = startingBetas[b];
	} // otherwise leave them as they were
	
	if (useLinearRegression)
		m->useLinearRegression(useLinearRegression);
	lnL = m->maximiseLnL();
	m->getSEs();
	if (fo)
		printModel(fo, name, lnL, m);
	return lnL;
}


