#include "modelMutation.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

get correlations with predictions for F1, U2, F3 and U4

then betas for each F1 then adding each F3 to relevant F1 and each F5 to relevant F3

#endif

const char* regions[] = {"nonTranscribed","coding","nonCoding"};

int main(int argc, char* argv[])
{
	int i, c, f, l, b, bs, b1, b2, b3, b4, b5, comp, * useThis, * useThese, upToOne, upToTwo, upToThree, numTrios, r;
	unsigned int hash;
	char target, source, first, second, third, fourth, fifth;
	float c1, realBeta, countTotal, LL, * fittedLLs[4];
	double* fittedBetas[4], * fittedSEs[4], fitted5BetaThreshold = 0.4;
	FILE* fi, * fo, * fplus, * fminus, * flog, * fpairs, * fc, * ftab;
	char* chr, line[1000], fn[100], fnn[3][100],seq[7], component[6], intercept[20], minusStrandSeq[7], sig[10], minusStrandSig[9];
	baseHasher hasher;
	glfModel model;
	useThis = (int*)calloc(NCOMPONENTS, sizeof(int)); // numbers of combinations for 1,2,3,4,5 bases
	useThese = (int*)calloc(NCOMPONENTS, sizeof(int));
	// fGlfLog = fopen("modelMutationGfl.log.txt", "w");
	flog = fopen("modelMutation.log.txt", "w");
	strcpy(intercept, "Intercept");
	for (b = 0; b < 4; ++b) {
		fittedBetas[b] = (double*)calloc(NCOMPONENTS, sizeof(double));
		fittedSEs[b] = (double*)calloc(NCOMPONENTS, sizeof(double));
		fittedLLs[b] = (float*)calloc(NCOMPONENTS, sizeof(float));
	}

	for (r = 0; r < 3; ++r) {
		sprintf(fn, "backgroundCounts.%s.txt", regions[r]);
		fi = fopen(fn, "r");
		fgets(line, 999, fi);
		i = 0;
		countTotal = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", backgroundSequenceTable[i], &c1);
			hash = hasher.hashBases(backgroundSequenceTable[i], 5);
			backgroundCountTable[hash] = c1;
			countTotal += c1;
			++i;
		}
	fclose(fi);
	meanBackgroundCount = countTotal / i;

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.F1.%s.txt",regions[r]);
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.F1.freq.%s.txt",regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "F1\tFrequency\tCI\n");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		for (bs = 0; bs < 4; ++bs) {
			if (bs == b)
				continue;
			source = baseNames[bs];
			sprintf(componentSequences[0], "__%c__%c", source, target);
			sprintf(componentNames[0], "F1__%c__%c", source, target);
			useThese[0] = 1;
			l = 1;
			useThese[l] = 0; // do not use intercept, just ends up close to 0 anyway
			model.init(NSEQUENCES / 12, l); // only those with correct target and source
			sprintf(fn, "AC25.%s.txt",regions[r]);
			fi = fopen(fn, "r");
			fgets(line, 999, fi);
			i = 0;
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
				if (!matches(sequenceTable[i], componentSequences[0], 6))
					continue;
				hash = hasher.hashBases(sequenceTable[i], 5); // just for background
				model.F[i] = backgroundCountTable[hash];
				model.Y[i] = c1 / backgroundCountTable[hash];
				for (comp = 0; comp < l; ++comp)
					model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
				// this will be all of them but leave for now to show logic
				++i;
			}
			fclose(fi);
			for (comp = 0; comp < l; ++comp)
				model.name[comp] = componentNames[comp];
			model.name[comp] = intercept;
			for (comp = 0; comp < l; ++comp)
				model.toFit[comp] = model.toUse[comp] = 1;
			setStartingBetasFromCounts(&model);
			LL = model.getLnL();
			sprintf(line, "modelMutation.F1.%c.%c LL = ", source, target);
			printModel(fo, line, LL, &model);
			for (comp = 0; comp < l; ++comp) {
				hasher.convertToSig(sig, model.name[comp] + 2);
				fprintf(ftab, "%s\t%.6f\t(%.6f - %.6f)\n", sig,
					exp(model.beta[comp]) / (exp(model.beta[comp]) + 1),
					exp(model.beta[comp] - 2 * model.SE[comp]) / (exp(model.beta[comp] - 2 * model.SE[comp]) + 1),
					exp(model.beta[comp] + 2 * model.SE[comp]) / (exp(model.beta[comp] + 2 * model.SE[comp]) + 1));
			}
		}
	}
	fclose(fo);
	fclose(ftab);
	sprintf(fnn[0], "modelMutation.F1.%s.txt", regions[r]);
	sprintf(fnn[1], "AC25.frequencies.%s.txt", regions[r]);
	sprintf(fnn[2], "modelMutation.F1.predictions.%s.txt", regions[r]);

	getPredictedFreqs(flog, fnn[0], fnn[1], fnn[2]);

	// must use frequencies for predictions because otherwise one is using the background counts to predict variant frequency

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.U2.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	sprintf(fn, "U2.complements.%s.txt", regions[r]);
	fpairs = fopen(fn, "w"); // helpful for Excel to find complementary sequence
	sprintf(fn, "modelMutation.U2.freq.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "U2\tFrequency\n");
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
				fprintf(fpairs, "U2%s\tD2%s\n", componentSequences[0], minusStrandSeq);
				sprintf(componentNames[0], "U2_%c%c__%c", first, source, target);
				useThese[0] = 1;
				l = 1;
				useThese[1] = 0;
				model.init(NSEQUENCES / 48, l); // only those with correct target, source and first
				sprintf(fn, "AC25.%s.txt", regions[r]);
				fi = fopen(fn, "r");
				fgets(line, 999, fi);
				i = 0;
				while (fgets(line, 999, fi)) {
					sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
					if (!matches(sequenceTable[i], componentSequences[0], 6))
						continue;
					hash = hasher.hashBases(sequenceTable[i], 5); // just for background
					model.F[i] = backgroundCountTable[hash];
					model.Y[i] = c1 / backgroundCountTable[hash];
					for (comp = 0; comp < l; ++comp)
						model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
					++i;
				}
				fclose(fi);

				for (comp = 0; comp < l; ++comp)
					model.name[comp] = componentNames[comp];
				model.name[comp] = intercept;
				for (comp = 0; comp < l; ++comp)
					model.toFit[comp] = model.toUse[comp] = 1;
				setStartingBetasFromCounts(&model);
				LL = model.getLnL();
				sprintf(line, "modelMutation.U2.%c.%c.%c LL = ", first, source, target);
				printModel(fo, line, LL, &model);
				for (comp = 0; comp < l; ++comp) {
					hasher.convertToSig(sig, model.name[comp] + 2);
					fprintf(ftab, "%s\t%.5f\n", sig, exp(model.beta[comp]) / (exp(model.beta[comp]) + 1));
				}
			}
		}
	}
	fclose(fo);
	fclose(ftab);
	fclose(fpairs);

	sprintf(fnn[0], "modelMutation.U2.%s.txt", regions[r]);
	sprintf(fnn[1], "AC25.frequencies.%s.txt", regions[r]);
	sprintf(fnn[2], "modelMutation.U2.predictions.%s.txt", regions[r]);

	getPredictedFreqs(flog, fnn[0], fnn[1], fnn[2]);
	
	// try adding F3s to F1s
	sprintf(fn, "modelMutation.eachF3WithF1s.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.fittedF3BetasWithF1s.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.fittedF3BetasWithF1s.plus.%s.txt", regions[r]);
	fplus = fopen(fn, "w");
	fprintf(fplus, "Variant\tOR\tCI\n");

	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			sprintf(componentSequences[0], "__%c__%c", source, target);
			sprintf(componentNames[0], "F1__%c__%c", source, target);
			useThese[l] = 1;
			l = 1;
			useThese[1] = 0;
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
			sprintf(fn, "AC25.%s.txt", regions[r]);
			fi = fopen(fn, "r");
			fgets(line, 999, fi);
			i = 0;
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
				if (!matches(sequenceTable[i], componentSequences[0], 6))
					continue;
				hash = hasher.hashBases(sequenceTable[i], 5); // just for background
				model.F[i] = backgroundCountTable[hash];
				model.Y[i] = c1 / backgroundCountTable[hash];
				for (comp = 0; comp < l; ++comp)
					model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
				++i;
			}
			fclose(fi);

			for (comp = 0; comp < l; ++comp)
				model.name[comp] = componentNames[comp];
			model.name[comp] = intercept;

			for (comp = 0; comp < upToOne; ++comp)
				useThis[comp] = 1;
			for (comp = upToOne; comp < l; ++comp)
				useThis[comp] = 0;
			for (c = 0; c < l; ++c)
				model.toFit[c] = model.toUse[c] = useThis[c];
			sprintf(fn, "modelMutation.eachF3WithF1s.%s.txt", regions[r]);
			fo = fopen(fn, "a");
			for (comp = upToOne; comp < l; ++comp) {
				setStartingBetasFromCounts(&model);
				model.toFit[comp] = model.toUse[comp] = useThis[comp] = 1;
				evaluateModel(fo, &model, "modelMutation.eachF3withF1s LL = ", useThis);
				realBeta = model.beta[comp];
				fittedBetas[b][comp] = model.beta[comp];
				fittedSEs[b][comp] = model.SE[comp];
				fittedLLs[b][comp] = model.getLnL();
				model.toFit[comp] = model.toUse[comp] = useThis[comp] = 0;
			}
			fclose(fo);

			sprintf(fn, "modelMutation.fittedF3BetasWithF1s.%s.txt", regions[r]);
			fo = fopen(fn, "a");

			for (comp = upToOne; comp < l; ++comp) {
				fprintf(fo, "%s\t%f\t%f\t%f\t%f\n", componentNames[comp],
					fittedBetas[b][comp], fittedSEs[b][comp],
					fittedBetas[b][comp] / fittedSEs[b][comp], fittedLLs[b][comp]);
				if (componentNames[comp][4] == 'C' || componentNames[comp][4] == 'T') {
					hasher.convertToSig(sig, componentNames[comp] + 2);
					fprintf(fplus, "%s\t%.3f\t(%.3f - %.3f)\n", sig,
						exp(fittedBetas[b][comp]),
						exp(fittedBetas[b][comp] - 2 * fittedSEs[b][comp]),
						exp(fittedBetas[b][comp] + 2 * fittedSEs[b][comp]));
				}
			}
			fclose(fo);
		}

	}
	fclose(fplus);


	sprintf(fn, "modelMutation.F3.plus.%s.txt", regions[r]);
	fplus = fopen(fn, "w");
	sprintf(fn, "modelMutation.F3.minus.%s.txt", regions[r]);
	fminus = fopen(fn, "w");
	fclose(fplus);
	fclose(fminus);
	sprintf(fn, "modelMutation.F3.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.F3.freq.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "F3\tFrequency\n");
	sprintf(fn, "modelMutation.F3.toPlot.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "F3.complements.%s.txt", regions[r]);
	fpairs = fopen(fn, "w");
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
					hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
					fprintf(fpairs, "F3%s\tF3%s\n", componentSequences[l], minusStrandSeq);
					sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
					useThese[l] = 1;
					++l;
					useThese[l] = 0;
					model.init(16, l); // only those with correct trinucleotide variant
					sprintf(fn, "AC25.%s.txt", regions[r]);
					fi = fopen(fn, "r");
					fgets(line, 999, fi);
					i = 0;
					while (fgets(line, 999, fi)) {
						sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
						if (!matches(sequenceTable[i], componentSequences[0], 6))
							continue;
						hash = hasher.hashBases(sequenceTable[i], 5); // just for background
						model.F[i] = backgroundCountTable[hash];
						model.Y[i] = c1 / backgroundCountTable[hash];
						for (comp = 0; comp < l; ++comp)
							model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
						++i;
					}
					for (comp = 0; comp < l; ++comp)
						model.name[comp] = componentNames[comp];
					model.name[comp] = intercept;
					fclose(fi);
					for (comp = 0; comp < l; ++comp)
						model.toFit[comp] = model.toUse[comp] = 1;
					setStartingBetasFromCounts(&model);
					LL = model.getLnL();
					sprintf(line, "modelMutation.F3.%c.%c.%c.%c LL = ", first, source, second, target);
					sprintf(fn, "modelMutation.F3.%s.txt", regions[r]);
					fo = fopen(fn, "a");
					printModel(fo, line, LL, &model);
					fclose(fo);
					for (comp = 0; comp < l; ++comp) {
						hasher.convertToSig(sig, model.name[comp] + 2);
						fprintf(ftab, "%s\t%.5f\n", sig, exp(model.beta[comp]) / (exp(model.beta[comp]) + 1));
					}
					sprintf(fn, "modelMutation.F3.toPlot.%s.txt", regions[r]);
					fo = fopen(fn, "a");
					printModel(fo, line, LL, &model, 1); // clean version with just coefficients
					fclose(fo);
					sprintf(fn, "modelMutation.F3.plus.%s.txt", regions[r]);
					fplus = fopen(fn, "a");
					sprintf(fn, "modelMutation.F3.minus.%s.txt", regions[r]);
					fminus = fopen(fn, "a");
					for (comp = 0; comp < l; ++comp) {
						strcpy(seq, model.name[comp] + 2);
						realBeta = model.beta[comp];
						if (seq[2] == 'C' || seq[2] == 'T') {
							hasher.convertToSig(sig, seq);
							fprintf(fplus, "%s\t%f\n", sig, realBeta);
							sprintf(component, "%c%c%c%c", seq[1], seq[2], seq[3], seq[5]);
						}
						else {
							hasher.getVariantComplement(minusStrandSeq, seq);
							hasher.convertToSig(sig, minusStrandSeq);
							fprintf(fminus, "%s\t%f\n", sig, realBeta);
						}
					}
					fclose(fplus);
					fclose(fminus);
				}
			}
		}
	}
	fclose(fpairs);
	fclose(ftab);

	sprintf(fnn[0], "modelMutation.F3.%s.txt", regions[r]);
	sprintf(fnn[1], "AC25.frequencies.%s.txt", regions[r]);
	sprintf(fnn[2], "modelMutation.F3.predictions.%s.txt", regions[r]);

	getPredictedFreqs(flog, fnn[0], fnn[1], fnn[2]);
	
	sprintf(fn, "modelMutation.U3.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.U3.freq.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "U3\tFrequency\n");
	sprintf(fn, "U3.complements.%s.txt", regions[r]);
	fpairs = fopen(fn, "w");
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
					hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
					fprintf(fpairs, "U3%s\tD3%s\n", componentSequences[l], minusStrandSeq);
					sprintf(componentNames[l], "U3%c%c%c__%c", second, first, source, target);
					useThese[l] = 1;
					++l;
					model.init(16, l); // only those with correct trinucleotide variant
					sprintf(fn, "AC25.%s.txt", regions[r]);
					fi = fopen(fn, "r");
					fgets(line, 999, fi);
					i = 0;
					while (fgets(line, 999, fi)) {
						sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
						if (!matches(sequenceTable[i], componentSequences[0], 6))
							continue;
						hash = hasher.hashBases(sequenceTable[i], 5); // just for background
						model.F[i] = backgroundCountTable[hash];
						model.Y[i] = c1 / backgroundCountTable[hash];
						for (comp = 0; comp < l; ++comp)
							model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
						++i;
					}
					for (comp = 0; comp < l; ++comp)
						model.name[comp] = componentNames[comp];
					model.name[comp] = intercept;
					fclose(fi);
					for (comp = 0; comp < l; ++comp)
						model.toFit[comp] = model.toUse[comp] = 1;
					setStartingBetasFromCounts(&model);
					LL = model.getLnL();
					sprintf(line, "modelMutation.U3.%c.%c.%c.%c LL = ", first, source, second, target);
					printModel(fo, line, LL, &model);
					for (comp = 0; comp < l; ++comp) {
						hasher.convertToSig(sig, model.name[comp] + 2);
						fprintf(ftab, "%s\t%.5f\n", sig, exp(model.beta[comp]) / (exp(model.beta[comp]) + 1));
					}
				}
			}
		}
	}
	fclose(fo);
	fclose(ftab);
	fclose(fpairs);

	sprintf(fnn[0], "modelMutation.U3.%s.txt", regions[r]);
	sprintf(fnn[1], "AC25.frequencies.%s.txt", regions[r]);
	sprintf(fnn[2], "modelMutation.U3.predictions.%s.txt", regions[r]);

	getPredictedFreqs(flog, fnn[0], fnn[1], fnn[2]);
	
	sprintf(fn, "modelMutation.U4.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.U4.freq.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "U4\tFrequency\n");
	sprintf(fn, "U4.complements.%s.txt", regions[r]);
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
						hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
						fprintf(fpairs, "U4%s\tD4%s\n", componentSequences[l], minusStrandSeq);
						sprintf(componentNames[l], "U4%s", componentSequences[l]);
						useThese[l] = 1;
						++l;
						model.init(4, l); // only those with correct five bases
			sprintf(fn, "AC25.%s.txt",regions[r]);
			fi = fopen(fn, "r");
			fgets(line, 999, fi);
						i = 0;
						while (fgets(line, 999, fi)) {
							sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
							if (!matches(sequenceTable[i], componentSequences[0], 6))
								continue;
							hash = hasher.hashBases(sequenceTable[i], 5); // just for background
							model.F[i] = backgroundCountTable[hash];
							model.Y[i] = c1 / backgroundCountTable[hash];
							for (comp = 0; comp < l; ++comp)
								model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
							++i;
						}
						for (comp = 0; comp < l; ++comp)
							model.name[comp] = componentNames[comp];
						model.name[comp] = intercept;
						fclose(fi);
						for (comp = 0; comp < l; ++comp)
							model.toFit[comp] = model.toUse[comp] = 1;
						setStartingBetasFromCounts(&model);
						LL = model.getLnL();
						sprintf(line, "modelMutation.U4.%c LL = ", target);
						printModel(fo, line, LL, &model);
						for (comp = 0; comp < l; ++comp) {
							hasher.convertToSig(sig, model.name[comp] + 2);
							fprintf(ftab, "%s\t%.5f\n", sig, exp(model.beta[comp]) / (exp(model.beta[comp]) + 1));
						}
					}
				}
			}
		}
	}
	fclose(fo);
	fclose(ftab);
	fclose(fpairs);

	sprintf(fnn[0], "modelMutation.U4.%s.txt", regions[r]);
	sprintf(fnn[1], "AC25.frequencies.%s.txt", regions[r]);
	sprintf(fnn[2], "modelMutation.U4.predictions.%s.txt", regions[r]);

	getPredictedFreqs(flog, fnn[0], fnn[1], fnn[2]);
	
	sprintf(fn, "modelMutation.eachF5WithF3s.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.fittedF5BetasWithF3s.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.fittedF5BetasWithF3s.plus.%s.txt", regions[r]);
	fplus = fopen(fn, "w");
	fprintf(fplus, "Variant\tOR\tCI\n");

	sprintf(fn, "F5.complements.%s.txt", regions[r]);
	fpairs = fopen(fn, "w");
	sprintf(fn, "F5.sigs.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
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
					useThese[l] = 1;
					++l;
					upToThree = l;

					for (b3 = 0; b3 < 4; ++b3) {
						third = baseNames[b3];
						for (b4 = 0; b4 < 4; ++b4) {
							fourth = baseNames[b4];
							sprintf(componentSequences[l], "%c%c%c%c%c%c", third, first, source, second, fourth, target);
							hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
							fprintf(fpairs, "F5%s\tF5%s\n", componentSequences[l], minusStrandSeq);
							hasher.convertToSig(sig, componentSequences[l]);
							fprintf(ftab, "%s\t%s\n", componentSequences[l], sig);
							sprintf(componentNames[l], "F5%s", componentSequences[l]);
							++l;
						}
					}

					model.init(16, l); // only those with correct trinucleotide variant
					sprintf(fn, "AC25.%s.txt", regions[r]);
					fi = fopen(fn, "r");
					fgets(line, 999, fi);
					i = 0;
					while (fgets(line, 999, fi)) {
						sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
						if (!matches(sequenceTable[i], componentSequences[0], 6))
							continue;
						hash = hasher.hashBases(sequenceTable[i], 5); // just for background
						model.F[i] = backgroundCountTable[hash];
						model.Y[i] = c1 / backgroundCountTable[hash];
						for (comp = 0; comp < l; ++comp)
							model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
						++i;
					}
					fclose(fi);

					for (comp = 0; comp < l; ++comp)
						model.name[comp] = componentNames[comp];
					model.name[comp] = intercept;

					for (comp = 0; comp < upToThree; ++comp)
						useThis[comp] = 1;
					for (comp = upToThree; comp < l; ++comp)
						useThis[comp] = 0;
					for (c = 0; c < l; ++c)
						model.toFit[c] = model.toUse[c] = useThis[c];
					sprintf(fn, "modelMutation.eachF5WithF3s.%s.txt", regions[r]);
					fo = fopen(fn, "a");
					for (comp = upToThree; comp < l; ++comp) {
						setStartingBetasFromCounts(&model);
						model.toFit[comp] = model.toUse[comp] = useThis[comp] = 1;
						sprintf(line, "modelMutation.eachF5WithF3s.%c LL = ", target);
						evaluateModel(fo, &model, line, useThis);
						fittedBetas[b][comp] = model.beta[comp];
						fittedSEs[b][comp] = model.SE[comp];
						fittedLLs[b][comp] = model.getLnL();
						model.toFit[comp] = model.toUse[comp] = useThis[comp] = 0;
					}
					fclose(fo);
					sprintf(fn, "modelMutation.fittedF5BetasWithF3s.%s.txt", regions[r]);
					fo = fopen(fn, "a");
					for (comp = upToThree; comp < l; ++comp) {
						fprintf(fo, "%s\t%f\t%f\t%f\t%f\n", componentNames[comp],
							fittedBetas[b][comp], fittedSEs[b][comp],
							fittedBetas[b][comp] / fittedSEs[b][comp],
							fittedLLs[b][comp]);
						if (componentNames[comp][4] == 'C' || componentNames[comp][4] == 'T') {
							hasher.convertToSig(sig, componentNames[comp] + 2);
							fprintf(fplus, "%s\t%.3f\t(%.3f - %.3f)\n", sig,
								exp(fittedBetas[b][comp]),
								exp(fittedBetas[b][comp] - 2 * fittedSEs[b][comp]),
								exp(fittedBetas[b][comp] + 2 * fittedSEs[b][comp]));
						}
					}
					fclose(fo);
				}
			}
		}

	}
	fclose(fpairs);
	fclose(fplus);
	fclose(ftab);

	sprintf(fn, "modelMutation.F5.%s.txt", regions[r]);
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.F5.freq.%s.txt", regions[r]);
	ftab = fopen(fn, "w");
	fprintf(ftab, "F5\tFrequency\n");
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
							hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
							sprintf(componentNames[l], "F5%s", componentSequences[l]);
							useThese[l] = 1;
							++l;
							model.init(1, l); // only those with correct five bases
							sprintf(fn, "AC25.%s.txt", regions[r]);
							fi = fopen(fn, "r");
							fgets(line, 999, fi);
							i = 0;
							while (fgets(line, 999, fi)) {
								sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
								if (!matches(sequenceTable[i], componentSequences[0], 6))
									continue;
								hash = hasher.hashBases(sequenceTable[i], 5); // just for background
								model.F[i] = backgroundCountTable[hash];
								model.Y[i] = c1 / backgroundCountTable[hash];
								for (comp = 0; comp < l; ++comp)
									model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
								++i;
							}
							for (comp = 0; comp < l; ++comp)
								model.name[comp] = componentNames[comp];
							model.name[comp] = intercept;
							fclose(fi);
							for (comp = 0; comp < l; ++comp)
								model.toFit[comp] = model.toUse[comp] = 1;
							setStartingBetasFromCounts(&model);
							LL = model.getLnL();
							sprintf(line, "modelMutation.F5.%c LL = ", target);
							printModel(fo, line, LL, &model);
							for (comp = 0; comp < l; ++comp) {
								hasher.convertToSig(sig, model.name[comp] + 2);
								fprintf(ftab, "%s\t%.5f\n", sig, exp(model.beta[comp]) / (exp(model.beta[comp]) + 1));
							}
						}
					}
				}
			}
		}
	}
	fclose(fo);
	fclose(ftab);


}

free(useThis);
free(useThese);
for (b = 0; b < 4; ++b) {
	free(fittedBetas[b]);
	free(fittedSEs[b]);
	free(fittedLLs[b]);
}

	// fclose(fGlfLog);
	return 0;
}
