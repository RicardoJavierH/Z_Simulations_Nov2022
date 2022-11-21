#ifndef ERRORESTIMATORGRADIENTANDLAPLACIANRECONSTRUCTIONWLSHPP
#define ERRORESTIMATORGRADIENTANDLAPLACIANRECONSTRUCTIONWLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradLaplacianRecWLS : TErrorEstimator {

public:
	TEEGradLaplacianRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradLaplacianRecWLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
