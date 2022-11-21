#ifndef ERRORESTIMATORGRADPOTENTIALRECWLSHPP
#define ERRORESTIMATORGRADPOTENTIALRECWLSHPP

#include <iostream>
#include "ErrorEstimator.h"


class TEEGradPotentialRecWLS : TErrorEstimator {

public:
	TEEGradPotentialRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradPotentialRecWLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);

};

#endif
