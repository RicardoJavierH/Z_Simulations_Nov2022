#ifndef ERRORESTIMATORGRADIENTRECONSTRUCTIONWLSHPP
#define ERRORESTIMATORGRADIENTRECONSTRUCTIONWLSHPP

#include <iostream>
#include "ErrorEstimator.h"


class TEEGradientRecWLS : TErrorEstimator {
	
public:
	TEEGradientRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEEGradientRecWLS();
	
	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);

};

#endif
