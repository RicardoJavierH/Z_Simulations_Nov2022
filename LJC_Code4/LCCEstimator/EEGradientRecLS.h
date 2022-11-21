#ifndef ERRORESTIMATORGRADIENTRECONSTRUCTIONLSHPP
#define ERRORESTIMATORGRADIENTRECONSTRUCTIONLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradientRecLS : TErrorEstimator {
	
public:
	TEEGradientRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradientRecLS();
	
	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
