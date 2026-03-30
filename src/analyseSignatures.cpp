#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analyseSignatures.hpp"

#define nSBS 96
#define MAXSIGS 150

class SBSvariant {
public:
	char type[7], var[5];
	int count;
};

class SBSvariantSignature {
public:
	char type[7], var[5];
	float weight[MAXSIGS];
};

SBSvariant varCounts[nSBS];
SBSvariantSignature varSigs[4*4*4*4];  // because I do not know which strand will be used for hashing
char sigNames[MAXSIGS][100],inputNames[MAXSIGS][100]; // may be in different orders
const char* COSMICFILES[] = { "COSMIC_v3.4_SBS_GRCh38.txt", "human_sbs96_filtered_v1_0.txt" };

int main(int argc, char* argv[])
{
	FILE* fi, * fo;
	int i,nSig,onSpace,sig,useThese[MAXSIGS+1],ff,* forSBS,c1;
	unsigned int hash;
	float count,LL,zThreshold,scale;
	baseHasher hasher;
	glfModel model;
	char line[4001],buff[4000],source[100],type[8],var[5],*ptr, seq[7], component[6], fn[100],SBSsig[8];
	forSBS = (int*)calloc(4 * 4 * 4 * 4, sizeof(int)); // because we do not know numbering of strands

	fo = fopen("SBSVariantsForInputMatrix.txt", "w"); // may be useful e.g. for Excel
	generateSBSVariantForInputMatrix(seq, 1);
	while (generateSBSVariantForInputMatrix(seq)) {
		hasher.convertToSig(SBSsig, seq);
		fprintf(fo, "F3%s\t%s\n", seq, SBSsig);
	}
	fclose(fo);

	fo = fopen("SBSVariants.txt", "w"); // may be useful e.g. for Excel
	generateSBSVariant(seq, 1);
	while (generateSBSVariant(seq)) {
		hasher.convertToSig(SBSsig, seq);
		fprintf(fo, "F3%s\t%s\n", seq, SBSsig);
	}
	fclose(fo);

		for (i = 0; i < 4 * 4 * 4 * 4; ++i)
			forSBS[i] = 0;
		sprintf(fn, "counts.total.ac25.txt");
		fi = fopen(fn, "r");
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %d", seq, &c1);
			if (c1 == 0)
				continue;
			if (seq[2] != 'C' && seq[2] != 'T') {
				hasher.getVariantComplement(buff, seq);
				strcpy(seq, buff);
			}
			sprintf(component, "%c%c%c%c", seq[1], seq[2], seq[3], seq[5]);
			hash = hasher.hashBases(component, 4);
			forSBS[hash] += c1; // this means the counts for complementary pairs will be added together
		}
		fclose(fi);
		sprintf(fn, "SBSVariants.counts.txt");
		fo = fopen(fn, "w");
		fprintf(fo, "Variant\tcounts\n");
		generateSBSVariant(seq, 1);
		while (generateSBSVariant(seq)) {
			hasher.convertToSig(SBSsig, seq);
			sprintf(component, "%c%c%c%c", seq[1], seq[2], seq[3], seq[5]);
			hash = hasher.hashBases(component, 4);
			fprintf(fo, "%s\t%d\n", SBSsig, forSBS[hash]);
		}
		fclose(fo);
	
		for (ff = 0; ff < 2; ++ff) {
			fi = fopen(COSMICFILES[ff], "r");
			fgets(line, 4000, fi);
			ptr = line;
			nSig = 0;
			onSpace = isspace(*ptr);
			do {
				if (onSpace != isspace(*ptr))
				{
					if (!onSpace)
						++nSig;
					onSpace = isspace(*ptr);
				}
			} while (*++ptr);
			--nSig;
			model.init(nSBS, nSig);
			model.useLinearRegression(MLRNONNEGATIVELNLIKE); // fit square roots of coefficients, so they have to be positive


			sscanf(line, "%*s %[^\n]", buff);
			for (sig = 0; strcpy(line, buff), *buff = '\0', sscanf(line, "%s %[^\n]", sigNames[sig], buff) >= 1; ++sig)
				model.name[sig] = sigNames[sig];
			while (fgets(line, 4000, fi) && sscanf(line, "%s %[^\n]", type, buff) >= 1) {
				sprintf(var, "%c%c%c%c", type[0], type[2], type[6], type[4]);
				hash = hasher.hashBases(var, 4);
				strcpy(varSigs[hash].type, type);
				strcpy(varSigs[hash].var, var);
				for (sig = 0; strcpy(line, buff), *buff = '\0', sscanf(line, "%f %[^\n]", &varSigs[hash].weight[sig], buff) >= 1; ++sig)
					useThese[sig] = 1;
			}
			fclose(fi);
			useThese[sig] = 0; // no intercept
			sprintf(fn, "SBSVariants.counts.txt");
			fi = fopen(fn, "r");
			fgets(line, 1000, fi);
			sscanf(line, "%*s %s", source);
			i = 0;
			scale = -1;
			while (fgets(line, 1000, fi) && sscanf(line, "%s %f", type, &count) == 2)
			{
				strcpy(inputNames[i], type);
				sprintf(var, "%c%c%c%c", type[0], type[2], type[6], type[4]);
				hash = hasher.hashBases(var, 4);
				for (sig = 0; sig < nSig; ++sig)
					model.X[i][sig] = varSigs[hash].weight[sig];
				if (scale == -1) { // first line
					if (count > 10000) {
						scale = pow(10.0, (int)log10(count));
					}
					else {
						scale = 1.0;
					}
				}
				model.Y[i] = count / scale; // if counts are too large minimisation fails, because small starting step is squared
				model.F[i] = 1;
				++i;
			}
			fclose(fi);
			sprintf(line, "SBS.%s.model.%d.txt", source,ff);
			fo = fopen(line, "w");
			sprintf(line, "model %d from %s LL = ", ff, source);
			evaluateModel(fo, &model, line, useThese, 0, 0, MLRNONNEGATIVELNLIKE);
			fclose(fo);

			if (ff==1)
				zThreshold = 4;
			else
				zThreshold = 7;
			for (sig = 0; sig < nSig; ++sig)
				if (model.beta[sig] / model.SE[sig] > zThreshold)
					useThese[sig] = 1;
				else
					useThese[sig] = 0; // includes NaN 
			sprintf(line, "SBS.best.%s.model.%d.txt", source, ff);
			fo = fopen(line, "w");
			sprintf(line, "fitted model %d from %s LL = ", ff, source);
			evaluateModel(fo, &model, line, useThese, 0,0, MLRNONNEGATIVELNLIKE);
			fclose(fo);

			sprintf(line, "SBS.%s.%d.fitted.txt", source,ff);
			fo = fopen(line, "w");
			fprintf(fo, "Sig\t");
			for (sig = 0; sig < nSig; ++sig) {
				if (useThese[sig])
					fprintf(fo, "%s\t", model.name[sig]);
			}
			fprintf(fo, "Predicted\tObserved\n");

			for (i = 0; i < nSBS; ++i) {
				count = 0;
				fprintf(fo, "%s\t",inputNames[i]);
				for (sig = 0; sig < nSig; ++sig) {
					if (useThese[sig]) {
						fprintf(fo, "%d\t", (int)(model.X[i][sig] * model.beta[sig] * model.beta[sig] * scale));
						count += model.X[i][sig] * model.beta[sig] * model.beta[sig];
					}
				}
				fprintf(fo, "%d\t%d\n", (int)(count * scale), (int)(model.Y[i] * scale));
			}
			fclose(fo);
		}
	return 0;
}
