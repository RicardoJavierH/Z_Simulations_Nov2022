#ifndef ERRORESTIMATORGRADIENTANDLAPLACIANRECONSTRUCTIONLSHPP
#define ERRORESTIMATORGRADIENTANDLAPLACIANRECONSTRUCTIONLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradLaplacianRecLS : TErrorEstimator {

public:
	TEEGradLaplacianRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradLaplacianRecLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
