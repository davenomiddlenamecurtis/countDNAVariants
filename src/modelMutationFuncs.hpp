#ifndef MODELMUTATIONFUNCSHPP
#define MODELMUTATIONFUNCSHPP
#include "dcerror.hpp"
#include "glfModel.hpp"
#include "hashBases.hpp"
#include "runModels.hpp"

#define NCOMPONENTS (4+4+4*3+4*3*(2*4+3*16+2*64+256))
#define NSEQUENCES (4*4*4*4*4*4) // not all real sequences but can get hashed to this number
#define NBACKGROUND (4*4*4*4*4)

extern char chrNames[23][3], sequenceTable[NSEQUENCES][7], componentNames[NCOMPONENTS][9], componentSequences[NCOMPONENTS][7], backgroundSequenceTable[NBACKGROUND][6];
extern float sequenceCountTable[NSEQUENCES][2], backgroundCountTable[NBACKGROUND], xy[2][NCOMPONENTS],meanBackgroundCount;
extern baseHasher hasher;

double correl(float* x, float* y, int count);
int getPredictedFreqs(FILE* flog, const char* compBetaFileName, const char* observedFreqsName, const char* outputFileName);
int getPredictedCounts(FILE* flog, const char* compBetaFileName, const char* observedSeqsName, const char* outputFileName);
int generateSBSVariant(char* s, int start = 0);
int matches(char* seq, char* comp, int len);

#endif