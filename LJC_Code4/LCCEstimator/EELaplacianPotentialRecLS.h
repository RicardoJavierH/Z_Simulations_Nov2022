#ifndef ERRORESTIMATORLAPLACIANPOTENTIALRECONSTRUCTIONLSHPP
#define ERRORESTIMATORLAPLACIANPOTENTIALRECONSTRUCTIONLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEELaplacianPotentialRecLS : TErrorEstimator {

public:
	TEELaplacianPotentialRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEELaplacianPotentialRecLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
