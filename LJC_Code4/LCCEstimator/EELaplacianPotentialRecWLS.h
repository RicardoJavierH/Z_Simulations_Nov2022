#ifndef ERRORESTIMATORLAPLACIANPOTENTIALRECONSTRUCTIONWLSHPP
#define ERRORESTIMATORLAPLACIANPOTENTIALRECONSTRUCTIONWLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEELaplacianPotentialRecWLS : TErrorEstimator {

public:
	TEELaplacianPotentialRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEELaplacianPotentialRecWLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
