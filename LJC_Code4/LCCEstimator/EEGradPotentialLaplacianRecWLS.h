#ifndef ERRORESTIMATORGRADPOTENTIALANDLAPLACIANRECONSTRUCTIONWLSHPP
#define ERRORESTIMATORGRADPOTENTIALANDLAPLACIANRECONSTRUCTIONWLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEEGradPotentialLaplacianRecWLS : TErrorEstimator {

public:
	TEEGradPotentialLaplacianRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradPotentialLaplacianRecWLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
