#ifndef LRFMODELHPP
#define LRFMODELHPP

// Like lrmodel.hpp but each observation also has a frequency entry

class lrfModel {
public:
	double **X, *Y,*F; // F is number of trials and Y is proportion of successes
	double *sigmaT,*t,*beta,*SE,*mean,*SD;
	int *toFit,*toUse;
	int nRow, nCol,gotMeans,isNormalised;
	char **name;
	void getMeans();
	void normalise();
	void deNormalise();
	double getLnL();
	void getSEs();
	lrfModel() { toFit = 0; X = 0; beta = 0; mean = 0; nRow = nCol = 0; name = 0; gotMeans = isNormalised=0;  }
	~lrfModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
	virtual double penaltyFunction() { return 0;  }
};

class lrfRidgePenaltyModel : public lrfModel {
public:
	lrfRidgePenaltyModel() { lamda = 0; }
	virtual double penaltyFunction();
	float lamda; // for ridge regression
};

#endif