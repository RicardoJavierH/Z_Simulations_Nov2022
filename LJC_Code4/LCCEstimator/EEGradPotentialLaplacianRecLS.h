#ifndef ERRORESTIMATORGRADPOTENTIALANDLAPLACIANRECONSTRUCTIONLSHPP
#define ERRORESTIMATORGRADPOTENTIALANDLAPLACIANRECONSTRUCTIONLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradPotentialLaplacianRecLS : TErrorEstimator {

public:
	TEEGradPotentialLaplacianRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradPotentialLaplacianRecLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
