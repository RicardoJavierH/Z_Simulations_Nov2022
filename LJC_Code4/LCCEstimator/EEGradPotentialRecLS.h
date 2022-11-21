#ifndef ERRORESTIMATORGRADPOTENTIALRECLSHPP
#define ERRORESTIMATORGRADPOTENTIALRECLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradPotentialRecLS : TErrorEstimator {

public:
	TEEGradPotentialRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradPotentialRecLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
