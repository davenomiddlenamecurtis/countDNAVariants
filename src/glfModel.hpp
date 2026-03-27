#ifndef GLFMODELHPP
#define GLFMODELHPP 1

#if 0
Copyright 2018 David Curtis

This file is part of the scoreassoc package.

scoreassoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

scoreassoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with scoreassoc.If not, see <http://www.gnu.org/licenses/>.
#endif

#include <math.h>
#ifndef M_PI
// Source: http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406 
#endif

enum thingToMaximise {LRLNLIKE=0,MLRLNLIKE, MLRNONNEGATIVELNLIKE}; // MLMINUSSSQ

double getMinusModelLnL();
double getRegularisedMinusModelLnL();

class glfModel {
public:
	double **X, *Y,LN2PI,*F; // F is frequency of each observation
	double *sigmaT,*t,*beta,*SE,*mean,*SD;
	int *toFit,*toUse;
	int nRow, nCol,gotMeans,isNormalised;
	char **name;
	void getMeans();
	void normalise();
	void deNormalise();
	double getLnL();
	void getSEs();
	glfModel() { 
		toFit = 0; toUse = 0;  X = 0; Y = 0; t = 0; t = 0; sigmaT = 0; F = 0;
		beta = 0; mean = 0; SE = 0; SD = 0; nRow = nCol = 0; name = 0; gotMeans = isNormalised = 0;
		LN2PI = log(2 * M_PI); thing = LRLNLIKE; func = getRegularisedMinusModelLnL;
	}
	~glfModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
	virtual double penaltyFunction() { return 0;  }
	void useLinearRegression(enum thingToMaximise l) { thing = l ;  }
	thingToMaximise thing;
	double (*func)();
	void setFunc(double (*f)()) { func = f; }
	
};

class glfRidgePenaltyModel : public glfModel {
public:
	glfRidgePenaltyModel() { lamda = 0; }
	virtual double penaltyFunction();
	float lamda; // for ridge regression
};

#endif