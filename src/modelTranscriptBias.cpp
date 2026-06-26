#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modelTranscriptBias.hpp"


#if 0
Do modelling similar to used for modelling strand bias but the odds to fit is for variant to be in a transcript or not
The intercept is the ratio of relevant background counts in one region or the other
The OR to model is whether a given variant is more likely to be in one region or the other

#endif

const char* regions[] = { "nonTranscribed","coding","nonCoding" };
enum region_t { NONTRANSCRIBED = 0,CODING,NONCODING };
// int totalCounts[3] = { 1457448, 848487,348195 };

float transcriptBackgroundCountTable[NBACKGROUND][2];
int nVariants;

int getPredictedTranscriptBiasBetas(FILE * flog, const char* compBetaFileName, const char* outputFileName) 
{
	FILE* fi, * fo;
	char line[1000], compName[100];
	int i, c, nComps,v,vv;
	unsigned int hash,bHash;
	float cFirst,cSecond, beta, SE, z, predCount, obsBeta;

	fo = fopen(outputFileName, "w");

// there will only be one component per model
	fi = fopen(compBetaFileName, "r");
	i = 0;
	vv = 0;
	while (fgets(line, 999, fi)) {
		if (sscanf(line, "%s %f %f %f", compName, &beta, &SE, &z) == 4) {
				if (strchr("MUFD", compName[0])) {
				strcpy(componentSequences[0], compName + 2);
				for (v = 0; v < nVariants; ++v) {
					if (matches(sequenceTable[v], componentSequences[0], 6)) {
						hash = hasher.hashBases(sequenceTable[v], 6);
						bHash= hasher.hashBases(sequenceTable[v], 5);
						if (sequenceCountTable[hash][0] == 0 || sequenceCountTable[hash][1] == 0)
							continue;
						obsBeta = log((sequenceCountTable[hash][1]/transcriptBackgroundCountTable[bHash][1]) / (sequenceCountTable[hash][0] / transcriptBackgroundCountTable[bHash][0]));
						// strictly speaking this is the ratio of frequencies, not the OR
						// but frequencies are so low that in practice there is no meaningful difference
						fprintf(fo, "%s\t%.5f\t%.5f\n", sequenceTable[v], obsBeta, beta);
						xy[0][vv] = obsBeta;
						xy[1][vv] = beta;
						++vv;
					}
				}
			}
		}
	}
	fclose(fi);
	fclose(fo);
	printf("Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], vv));
	fprintf(flog, "Correlation coefficient for %s: R = %f\n", outputFileName, correl(xy[0], xy[1], vv));
	if (vv != nVariants) {
		printf("(%d had zero count)\n", nVariants - vv);
		fprintf(flog,"(%d had zero count)\n", nVariants - vv);
	}

	return 1;
}

int main(int argc, char* argv[])
{
	int chr,i, c, f, l, b, bs, b1, b2, b3, b4, b5, comp, * useThis, * useThese, * fitThese, upToOne, upToTwo, upToThree, numTrios,v,r[2],rr,m;
	unsigned int hash,bHash;
	char target, source, first, second, third, fourth, fifth;
	float count, realBeta, countTotal, LL, * fittedLLs[4],totFirst,totSecond;
	double* fittedBetas[4], * fittedSEs[4], fitted5BetaThreshold = 0.4;
	FILE* fi, * fo, * flog, * fc, * ftab, * fplus;
	char chrName[3], line[1000], fn[100], fn2[100],fn3[100], seq[7], component[6], intercept[20], minusStrandSeq[7], sig[10], minusStrandSig[8];
	baseHasher hasher;
	glfModel model;
	useThis = (int*)calloc(NCOMPONENTS, sizeof(int)); // numbers of combinations for 1,2,3,4,5 bases
	useThese = (int*)calloc(NCOMPONENTS, sizeof(int));
	fitThese = (int*)calloc(NCOMPONENTS, sizeof(int));

	for (b = 0; b < 4; ++b) {
		fittedBetas[b] = (double*)calloc(NCOMPONENTS, sizeof(double));
		fittedSEs[b] = (double*)calloc(NCOMPONENTS, sizeof(double));
		fittedLLs[b] = (float*)calloc(NCOMPONENTS, sizeof(float));
	}


	flog = fopen("modelTranscriptBias.log.txt", "w");
	strcpy(intercept, "Intercept");
	for (m = 0; m < 3; ++m) {
		switch (m) {
		case 0:
			r[0] = NONTRANSCRIBED;
			r[1] = CODING;
			break;
		case 1:
			r[0] = NONTRANSCRIBED;
			r[1] = NONCODING;
			break;
		case 2:
			r[0] = NONCODING;
			r[1] = CODING;
			break;
		}
		
		for (rr = 0; rr < 2; ++rr) {
			sprintf(fn, "backgroundCounts.%s.txt", regions[r[rr]]);
			fi = fopen(fn, "r");
			fgets(line, 999, fi);
			while (fgets(line, 999, fi) && sscanf(line, "%s %f", seq, &count) == 2) {
				hash = hasher.hashBases(seq, 5);
				transcriptBackgroundCountTable[hash][rr] = count;
			}
			fclose(fi);
			sprintf(fn, "AC25.%s.txt",regions[r[rr]]);
			fi = fopen(fn, "r");
			i = 0;
			fgets(line, 999, fi);
			while (fgets(line, 999, fi) && sscanf(line, "%s %f", seq, &count) == 2) {
				strcpy(sequenceTable[i], seq);
				hash = hasher.hashBases(seq, 6);
				sequenceCountTable[hash][rr] = count;
				++i;
			}
			fclose(fi);
		}
		nVariants = i;

		sprintf(fn, "modelTranscriptBias.F1.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.F1.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "F1\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				sprintf(componentSequences[0], "__%c__%c", source, target);
				sprintf(componentNames[0], "F1__%c__%c", source, target);
				l = 1;
				model.init(4*4*4*4, l); // only those with correct target and source
				i = 0;
				totFirst = totSecond = 0;
				for (v = 0; v < nVariants; ++v) {
					if (!matches(sequenceTable[v], componentSequences[0], 6))
						continue;
					hash = hasher.hashBases(sequenceTable[v], 6);
					model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
					if (sequenceCountTable[hash][0] + sequenceCountTable[hash][1] == 0)
						model.Y[i]=0;
					else
						model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
					for (comp = 0; comp < l; ++comp)
						model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
					// this will be all of them but leave for now to show logic
					hash = hasher.hashBases(sequenceTable[v], 5);
					totFirst += transcriptBackgroundCountTable[hash][0];
					totSecond += transcriptBackgroundCountTable[hash][1];
					++i;
				}
				for (comp = 0; comp < l; ++comp)
					model.name[comp] = componentNames[comp];
				model.name[comp] = intercept;
				for (comp = 0; comp < l; ++comp)
					model.toUse[comp] = model.toFit[comp] = 1;
				model.toUse[comp] = 1; // intercept
				model.toFit[comp] = 0;
				model.beta[comp] = log(totSecond / totFirst);

				setStartingBetasFromCounts(&model);
				LL = model.getLnL();
				sprintf(line, "modelTranscriptBias.F1.%c.%c LL = ", source, target);
				printModel(fo, line, LL, &model);
				for (comp = 0; comp < l; ++comp) {
					hasher.convertToSig(sig, model.name[comp] + 2);
					fprintf(ftab, "%s\t%.3f\t(%.3f - %.3f)\n", sig, exp(model.beta[comp]),
						exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.F1.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.F1.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		// for each base outcome, model background which leads to it
		sprintf(fn, "modelTranscriptBias.U2.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.U2.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "U2\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			l = 0;
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					sprintf(componentSequences[0], "_%c%c__%c", first, source, target);
					hasher.getVariantComplement(minusStrandSeq, componentSequences[0]);
					sprintf(componentNames[0], "U2_%c%c__%c", first, source, target);
					l = 1;
					model.init(4*4*4, l); // only those with correct target, source and first
					i = 0;
					totFirst = totSecond = 0;
					for (v = 0; v < nVariants; ++v) {
						if (!matches(sequenceTable[v], componentSequences[0], 6))
							continue;
						hash = hasher.hashBases(sequenceTable[v], 6);
						model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
						model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
						for (comp = 0; comp < l; ++comp)
							model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
						// this will be all of them but leave for now to show logic
						hash = hasher.hashBases(sequenceTable[v], 5);
						totFirst += transcriptBackgroundCountTable[hash][0];
						totSecond += transcriptBackgroundCountTable[hash][1];
						++i;
					}
					for (comp = 0; comp < l; ++comp)
						model.name[comp] = componentNames[comp];
					model.name[comp] = intercept;
					for (comp = 0; comp < l; ++comp)
						model.toFit[comp] = model.toUse[comp] = 1;
					model.toFit[comp] = 0;
					model.toUse[comp] = 1; // intercept
					model.beta[comp] = log(totSecond / totFirst);

					setStartingBetasFromCounts(&model);
					LL = model.getLnL();
					sprintf(line, "modelTranscriptBias.U2.%c.%c.%c LL = ", first, source, target);
					printModel(fo, line, LL, &model);
					for (comp = 0; comp < l; ++comp) {
						hasher.convertToSig(sig, model.name[comp] + 2);
						fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
							exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.U2.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.U2.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.D2.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.D2.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "D2\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			l = 0;
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					sprintf(componentSequences[0], "__%c%c_%c", source, first, target);
					hasher.getVariantComplement(minusStrandSeq, componentSequences[0]);
					sprintf(componentNames[0], "D2__%c%c_%c", source, first, target);
					l = 1;
					model.init(NSEQUENCES / 48, l); // only those with correct target, source and first
					i = 0;
					totFirst = totSecond = 0;
					for (v = 0; v < nVariants; ++v) {
						if (!matches(sequenceTable[v], componentSequences[0], 6))
							continue;
						hash = hasher.hashBases(sequenceTable[v], 6);
						model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
						model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
						for (comp = 0; comp < l; ++comp)
							model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
						// this will be all of them but leave for now to show logic
						hash = hasher.hashBases(sequenceTable[v], 5);
						totFirst += transcriptBackgroundCountTable[hash][0];
						totSecond += transcriptBackgroundCountTable[hash][1];
						++i;
					}
					for (comp = 0; comp < l; ++comp)
						model.name[comp] = componentNames[comp];
					model.name[comp] = intercept;
					for (comp = 0; comp < l; ++comp)
						model.toFit[comp] = model.toUse[comp] = 1;
					model.toFit[comp] = 0;
					model.toUse[comp] = 1; // intercept
					model.beta[comp] = log(totSecond / totFirst);

					setStartingBetasFromCounts(&model);
					LL = model.getLnL();
					sprintf(line, "modelTranscriptBias.D2.%c.%c.%c LL = ", source, first, target);
					printModel(fo, line, LL, &model);
					for (comp = 0; comp < l; ++comp) {
						hasher.convertToSig(sig, model.name[comp] + 2);
						fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
							exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.D2.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.D2.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);


		// try adding F3s to F1s
		sprintf(fn, "modelTranscriptBias.eachF3WithF1s.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		fclose(fo);
		sprintf(fn, "modelTranscriptBias.fittedF3BetasWithF1s.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		fclose(fo);

		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			l = 0;
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				sprintf(componentSequences[0], "__%c__%c", source, target);
				sprintf(componentNames[0], "F1__%c__%c", source, target);
				l = 1;
				upToOne = 1;

				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						sprintf(componentSequences[l], "_%c%c%c_%c", first, source, second, target);
						sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
						++l;
					}
				}

				model.init(NSEQUENCES / 12, l); // only those with correct source and target
				i = 0;
				totFirst = totSecond = 0;
				for (v = 0; v < nVariants; ++v) {
					if (!matches(sequenceTable[v], componentSequences[0], 6))
						continue;
					hash = hasher.hashBases(sequenceTable[v], 6);
					model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
					model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
					for (comp = 0; comp < l; ++comp)
						model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
					hash = hasher.hashBases(sequenceTable[v], 5);
					totFirst += transcriptBackgroundCountTable[hash][0];
					totSecond += transcriptBackgroundCountTable[hash][1];
					++i;
				}
				fclose(fi);

				for (comp = 0; comp < l; ++comp)
					model.name[comp] = componentNames[comp];
				model.name[comp] = intercept;

				for (comp = 0; comp < upToOne; ++comp)
					fitThese[comp] = useThis[comp] = 1;
				for (comp = upToOne; comp < l; ++comp)
					fitThese[comp] = useThis[comp] = 0;
				for (c = 0; c < l; ++c)
					model.toFit[c] = model.toUse[c] = useThis[c];
				fitThese[comp] = model.toFit[comp] = 0;
				useThis[comp] = model.toUse[comp] = 1; // intercept
				model.beta[comp] = log(totSecond / totFirst);

				sprintf(fn, "modelTranscriptBias.eachF3WithF1s.%s.%s.txt",regions[r[0]],regions[r[1]]);
				fo = fopen(fn, "a");
				for (comp = upToOne; comp < l; ++comp) {
					setStartingBetasFromCounts(&model);
					useThis[comp] = fitThese[comp] = 1;
					evaluateModel(fo, &model, "modelTranscriptBias.eachF3withF1s LL = ", useThis, fitThese);
					realBeta = model.beta[comp];
					fittedBetas[b][comp] = model.beta[comp];
					fittedSEs[b][comp] = model.SE[comp];
					fittedLLs[b][comp] = model.getLnL();
					model.toFit[comp] = model.toUse[comp] = useThis[comp] = fitThese[comp] = 0;
					// need this for next setStartingBetasFromCounts()
				}
				fclose(fo);

				sprintf(fn, "modelTranscriptBias.fittedF3BetasWithF1s.%s.%s.txt",regions[r[0]],regions[r[1]]);
				fo = fopen(fn, "a");
				for (comp = upToOne; comp < l; ++comp) {
					fprintf(fo, "%s\t%f\t%f\t%f\t%f\n", componentNames[comp],
						fittedBetas[b][comp], fittedSEs[b][comp],
						fittedBetas[b][comp] / fittedSEs[b][comp], fittedLLs[b][comp]);
				}
				fclose(fo);
			}

		}

		sprintf(fn, "modelTranscriptBias.F3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.F3.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "F3\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						l = 0;
						sprintf(componentSequences[l], "_%c%c%c_%c", first, source, second, target);
						sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
						++l;
						model.init(16, l); // only those with correct trinucleotide variant
						i = 0;
						totFirst = totSecond = 0;
						for (v = 0; v < nVariants; ++v) {
							if (!matches(sequenceTable[v], componentSequences[0], 6))
								continue;
							hash = hasher.hashBases(sequenceTable[v], 6);
							model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
							model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
							for (comp = 0; comp < l; ++comp)
								model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
							// this will be all of them but leave for now to show logic
							hash = hasher.hashBases(sequenceTable[v], 5);
							totFirst += transcriptBackgroundCountTable[hash][0];
							totSecond += transcriptBackgroundCountTable[hash][1];
							++i;
						}

						for (comp = 0; comp < l; ++comp)
							model.name[comp] = componentNames[comp];
						model.name[comp] = intercept;
						for (comp = 0; comp < l; ++comp)
							model.toFit[comp] = model.toUse[comp] = 1;
						model.toFit[comp] = 0;
						model.toUse[comp] = 1; // intercept
						model.beta[comp] = log(totSecond / totFirst);

						setStartingBetasFromCounts(&model);
						LL = model.getLnL();
						sprintf(line, "modelTranscriptBias.F3.%c.%c.%c.%c LL = ", first, source, second, target);
						printModel(fo, line, LL, &model);
						for (comp = 0; comp < l; ++comp) {
							hasher.convertToSig(sig, model.name[comp] + 2);
							fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
								exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.F3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.F3.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.U3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.U3.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "U3\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						l = 0;
						sprintf(componentSequences[l], "%c%c%c__%c", second, first, source, target);
						sprintf(componentNames[l], "U3%c%c%c__%c", second, first, source, target);
						++l;
						model.init(16, l); // only those with correct trinucleotide variant
						i = 0;
						totFirst = totSecond = 0;
						for (v = 0; v < nVariants; ++v) {
							if (!matches(sequenceTable[v], componentSequences[0], 6))
								continue;
							hash = hasher.hashBases(sequenceTable[v], 6);
							model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
							model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
							for (comp = 0; comp < l; ++comp)
								model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
							// this will be all of them but leave for now to show logic
							hash = hasher.hashBases(sequenceTable[v], 5);
							totFirst += transcriptBackgroundCountTable[hash][0];
							totSecond += transcriptBackgroundCountTable[hash][1];
							++i;
						}

						for (comp = 0; comp < l; ++comp)
							model.name[comp] = componentNames[comp];
						model.name[comp] = intercept;
						for (comp = 0; comp < l; ++comp)
							model.toFit[comp] = model.toUse[comp] = 1;
						model.toFit[comp] = 0;
						model.toUse[comp] = 1; // intercept
						model.beta[comp] = log(totSecond / totFirst);

						setStartingBetasFromCounts(&model);
						LL = model.getLnL();
						sprintf(line, "modelTranscriptBias.U3.%c.%c.%c.%c LL = ", second, first, source, target);
						printModel(fo, line, LL, &model);
						for (comp = 0; comp < l; ++comp) {
							hasher.convertToSig(sig, model.name[comp] + 2);
							fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
								exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.U3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.U3.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.D3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.D3.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "D3\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						l = 0;
						sprintf(componentSequences[l], "__%c%c%c%c", source, first, second, target);
						sprintf(componentNames[l], "D3__%c%c%c%c", source, first, second, target);
						++l;
						model.init(16, l); // only those with correct trinucleotide variant
						i = 0;
						totFirst = totSecond = 0;
						for (v = 0; v < nVariants; ++v) {
							if (!matches(sequenceTable[v], componentSequences[0], 6))
								continue;
							hash = hasher.hashBases(sequenceTable[v], 6);
							model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
							model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
							for (comp = 0; comp < l; ++comp)
								model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
							// this will be all of them but leave for now to show logic
							hash = hasher.hashBases(sequenceTable[v], 5);
							totFirst += transcriptBackgroundCountTable[hash][0];
							totSecond += transcriptBackgroundCountTable[hash][1];
							++i;
						}

						for (comp = 0; comp < l; ++comp)
							model.name[comp] = componentNames[comp];
						model.name[comp] = intercept;
						for (comp = 0; comp < l; ++comp)
							model.toFit[comp] = model.toUse[comp] = 1;
						model.toFit[comp] = 0;
						model.toUse[comp] = 1; // intercept
						model.beta[comp] = log(totSecond / totFirst);

						setStartingBetasFromCounts(&model);
						LL = model.getLnL();
						sprintf(line, "modelTranscriptBias.D3.%c.%c.%c.%c LL = ", source, first, second, target);
						printModel(fo, line, LL, &model);
						for (comp = 0; comp < l; ++comp) {
							hasher.convertToSig(sig, model.name[comp] + 2);
							fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
								exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.D3.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.D3.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.U4.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.U4.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "U4\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						for (b3 = 0; b3 < 4; ++b3) {
							third = baseNames[b3];
							l = 0;
							sprintf(componentSequences[l], "%c%c%c%c_%c", third, second, source, first, target);
							sprintf(componentNames[l], "U4%s", componentSequences[l]);
							++l;
							model.init(4, l); // only those with correct five bases
							i = 0;
							totFirst = totSecond = 0;
							for (v = 0; v < nVariants; ++v) {
								if (!matches(sequenceTable[v], componentSequences[0], 6))
									continue;
								hash = hasher.hashBases(sequenceTable[v], 6);
								model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
								model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
								for (comp = 0; comp < l; ++comp)
									model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
								// this will be all of them but leave for now to show logic
								hash = hasher.hashBases(sequenceTable[v], 5);
								totFirst += transcriptBackgroundCountTable[hash][0];
								totSecond += transcriptBackgroundCountTable[hash][1];
								++i;
							}
							for (comp = 0; comp < l; ++comp)
								model.name[comp] = componentNames[comp];
							model.name[comp] = intercept;
							for (comp = 0; comp < l; ++comp)
								model.toFit[comp] = model.toUse[comp] = 1;
							model.toFit[comp] = 0;
							model.toUse[comp] = 1; // intercept
							model.beta[comp] = log(totSecond / totFirst);

							setStartingBetasFromCounts(&model);
							LL = model.getLnL();
							sprintf(line, "modelTranscriptBias.U4.%c.%c.%c.%c.%c LL = ", third, second, source, first, target);
							printModel(fo, line, LL, &model);
							for (comp = 0; comp < l; ++comp) {
								hasher.convertToSig(sig, model.name[comp] + 2);
								fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
									exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
							}
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.U4.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.U4.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.D4.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.D4.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "D4\tOR\t(CI)\n");
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						for (b3 = 0; b3 < 4; ++b3) {
							third = baseNames[b3];
							l = 0;
							sprintf(componentSequences[l], "_%c%c%c%c%c", first, source, second, third, target);
							sprintf(componentNames[l], "D4%s", componentSequences[l]);
							++l;
							model.init(4, l); // only those with correct five bases
							i = 0;
							totFirst = totSecond = 0;
							for (v = 0; v < nVariants; ++v) {
								if (!matches(sequenceTable[v], componentSequences[0], 6))
									continue;
								hash = hasher.hashBases(sequenceTable[v], 6);
								model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
								model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
								for (comp = 0; comp < l; ++comp)
									model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
								// this will be all of them but leave for now to show logic
								hash = hasher.hashBases(sequenceTable[v], 5);
								totFirst += transcriptBackgroundCountTable[hash][0];
								totSecond += transcriptBackgroundCountTable[hash][1];
								++i;
							}
							for (comp = 0; comp < l; ++comp)
								model.name[comp] = componentNames[comp];
							model.name[comp] = intercept;
							for (comp = 0; comp < l; ++comp)
								model.toFit[comp] = model.toUse[comp] = 1;
							model.toFit[comp] = 0;
							model.toUse[comp] = 1; // intercept
							model.beta[comp] = log(totSecond / totFirst);

							setStartingBetasFromCounts(&model);
							LL = model.getLnL();
							sprintf(line, "modelTranscriptBias.D4.%c.%c.%c.%c.%c LL = ", first, source, second, third, target);
							printModel(fo, line, LL, &model);
							for (comp = 0; comp < l; ++comp) {
								hasher.convertToSig(sig, model.name[comp] + 2);
								fprintf(ftab, "%s\t%.5f\t(%.5f - %.5f)\n", sig, exp(model.beta[comp]),
									exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]));
							}
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.D4.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.D4.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);

		sprintf(fn, "modelTranscriptBias.eachF5WithF3s.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		fclose(fo);
		sprintf(fn, "modelTranscriptBias.fittedF5BetasWithF3s.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		fclose(fo);
		sprintf(fn, "modelTranscriptBias.fittedF5BetasWithF3s.OR.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fplus = fopen(fn, "w");
		fprintf(fplus, "Variant\tOR\t(CI)\n");

		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						l = 0;
						sprintf(componentSequences[l], "_%c%c%c_%c", first, source, second, target);
						sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
						++l;
						upToThree = l;

						for (b3 = 0; b3 < 4; ++b3) {
							third = baseNames[b3];
							for (b4 = 0; b4 < 4; ++b4) {
								fourth = baseNames[b4];
								sprintf(componentSequences[l], "%c%c%c%c%c%c", third, first, source, second, fourth, target);
								sprintf(componentNames[l], "F5%s", componentSequences[l]);
								++l;
							}
						}

						model.init(16, l); // only those with correct trinucleotide variant
						i = 0;
						totFirst = totSecond = 0;
						for (v = 0; v < nVariants; ++v) {
							if (!matches(sequenceTable[v], componentSequences[0], 6))
								continue;
							hash = hasher.hashBases(sequenceTable[v], 6);
							model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
							model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
							for (comp = 0; comp < l; ++comp)
								model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
							hash = hasher.hashBases(sequenceTable[v], 5);
							totFirst += transcriptBackgroundCountTable[hash][0];
							totSecond += transcriptBackgroundCountTable[hash][1];
							++i;
						}

						for (comp = 0; comp < l; ++comp)
							model.name[comp] = componentNames[comp];
						model.name[comp] = intercept;

						for (comp = 0; comp < upToThree; ++comp)
							fitThese[comp] = useThis[comp] = 1;
						for (comp = upToThree; comp < l; ++comp)
							fitThese[comp] = useThis[comp] = 0;
						for (c = 0; c < l; ++c)
							model.toFit[c] = model.toUse[c] = useThis[c];
						fitThese[comp] = model.toFit[comp] = 0;
						useThis[comp] = model.toUse[comp] = 1; // intercept
						model.beta[comp] = log(totSecond / totFirst);

						sprintf(fn, "modelTranscriptBias.eachF5WithF3s.%s.%s.txt",regions[r[0]],regions[r[1]]);
						fo = fopen(fn, "a");
						for (comp = upToThree; comp < l; ++comp) {
							setStartingBetasFromCounts(&model);
							useThis[comp] = fitThese[comp] = 1;
							evaluateModel(fo, &model, "modelTranscriptBias.eachF5WithF3s LL = ", useThis, fitThese);
							fittedBetas[b][comp] = model.beta[comp];
							fittedSEs[b][comp] = model.SE[comp];
							fittedLLs[b][comp] = model.getLnL();
							model.toFit[comp] = model.toUse[comp] = useThis[comp] = fitThese[comp] = 0;
						}
						fclose(fo);
						sprintf(fn, "modelTranscriptBias.fittedF5BetasWithF3s.%s.%s.txt",regions[r[0]],regions[r[1]]);
						fo = fopen(fn, "a");
						for (comp = upToThree; comp < l; ++comp) {
							fprintf(fo, "%s\t%f\t%f\t%f\t%f\n", componentNames[comp],
								fittedBetas[b][comp], fittedSEs[b][comp],
								fittedBetas[b][comp] / fittedSEs[b][comp],
								fittedLLs[b][comp]);
							hasher.convertToSig(sig, componentNames[comp] + 2);
							fprintf(fplus, "%s\t%.5f\t(%.5f - %.5f)\n", sig,
								exp(fittedBetas[b][comp]),
								exp(fittedBetas[b][comp] - 2 * fittedSEs[b][comp]),
								exp(fittedBetas[b][comp] + 2 * fittedSEs[b][comp]));
						}
						fclose(fo);
					}
				}
			}

		}
		fclose(fplus);


		sprintf(fn, "modelTranscriptBias.F5.%s.%s.txt",regions[r[0]],regions[r[1]]);
		fo = fopen(fn, "w");
		sprintf(fn, "modelTranscriptBias.F5.ORs.%s.%s.txt",regions[r[0]],regions[r[1]]);
		ftab = fopen(fn, "w");
		fprintf(ftab, "F5\tOR\t(CI)\tz\t%sCount\t%sBackground\t%sCount\t%sBackground\n",
			regions[r[1]], regions[r[1]], regions[r[0]], regions[r[0]]);
		for (b = 0; b < 4; ++b) {
			target = baseNames[b];
			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						for (b3 = 0; b3 < 4; ++b3) {
							third = baseNames[b3];
							for (b4 = 0; b4 < 4; ++b4) {
								fourth = baseNames[b4];
								l = 0;
								sprintf(componentSequences[l], "%c%c%c%c%c%c", third, first, source, second, fourth, target);
								sprintf(componentNames[l], "F5%s", componentSequences[l]);
								++l;
								model.init(1, l); // only those with correct five bases
								i = 0;
								totFirst = totSecond = 0;
								for (v = 0; v < nVariants; ++v) {
									if (!matches(sequenceTable[v], componentSequences[0], 6))
										continue;
									hash = hasher.hashBases(sequenceTable[v], 6);
									model.F[i] = sequenceCountTable[hash][0] + sequenceCountTable[hash][1];
									model.Y[i] = sequenceCountTable[hash][1] / (sequenceCountTable[hash][0] + sequenceCountTable[hash][1]);
									for (comp = 0; comp < l; ++comp)
										model.X[i][comp] = matches(sequenceTable[v], componentSequences[comp], 6);
									// this will be all of them but leave for now to show logic
									hash = hasher.hashBases(sequenceTable[v], 5);
									totFirst += transcriptBackgroundCountTable[hash][0];
									totSecond += transcriptBackgroundCountTable[hash][1];
									++i;
								}

								for (comp = 0; comp < l; ++comp)
									model.name[comp] = componentNames[comp];
								model.name[comp] = intercept;
								for (comp = 0; comp < l; ++comp)
									model.toFit[comp] = model.toUse[comp] = 1;
								model.toFit[comp] = 0;
								model.toUse[comp] = 1; // intercept
								model.beta[comp] = log(totSecond / totFirst);

								setStartingBetasFromCounts(&model);
								LL = model.getLnL();
								sprintf(line, "modelTranscriptBias.F5.%c LL = ", target);
								printModel(fo, line, LL, &model);
								for (comp = 0; comp < l; ++comp) {
									hasher.convertToSig(sig, model.name[comp] + 2);
									hash = hasher.hashBases(model.name[comp] + 2, 6);
									bHash = hasher.hashBases(model.name[comp] + 2, 5);
									fprintf(ftab, "%s\t%.3f\t(%.3f - %.3f)\t%2f\t%d\t%d\t%d\t%d\n", sig, exp(model.beta[comp]),
										exp(model.beta[comp] - 2 * model.SE[comp]), exp(model.beta[comp] + 2 * model.SE[comp]),
										model.beta[comp] / model.SE[comp],
										(int)sequenceCountTable[hash][1], (int)transcriptBackgroundCountTable[bHash][1],
										(int)sequenceCountTable[hash][0], (int)transcriptBackgroundCountTable[bHash][0]);
								}
							}
						}
					}
				}
			}
		}
		fclose(fo);
		fclose(ftab);

		sprintf(fn, "modelTranscriptBias.F5.%s.%s.txt",regions[r[0]],regions[r[1]]);
		sprintf(fn2, "modelTranscriptBias.F5.predictions.%s.%s.txt",regions[r[0]],regions[r[1]]);
		getPredictedTranscriptBiasBetas(flog, fn, fn2);
	}



	free(useThis);
	free(useThese);
	free(fitThese);
	for (b = 0; b < 4; ++b) {
		free(fittedBetas[b]);
		free(fittedSEs[b]);
		free(fittedLLs[b]);
	}

	// fclose(fGlfLog);
	return 0;
}
