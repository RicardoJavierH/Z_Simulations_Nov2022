#include <stdlib.h>

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "TPZMaterial.h"

#include "EEGradLaplacianRecWLS.h"

TEEGradLaplacianRecWLS::TEEGradLaplacianRecWLS(TPZCompMesh *cmesh, SimulationData simul) : TErrorEstimator(cmesh, simul) {
}
TEEGradLaplacianRecWLS::~TEEGradLaplacianRecWLS() {
}

/* Function as error estimator for each element */
REAL TEEGradLaplacianRecWLS::RecoveredGradient(TPZCompEl *cel) {
	TPZGeoEl *father = cel->Reference()->Father();
	if (!father) return 0.;
	int i, nsons = father->NSubElements();
	if (!nsons) return 100.;
	TPZVec<TPZCompEl *> Sons(nsons, 0);
	TPZVec<REAL> MeansBySubElement(nsons, 0.);
	for (i = 0; i < nsons; i++) {
		Sons[i] = father->SubElement(i)->Reference();
		if (!Sons[i])
			DebugStop();
		MeansBySubElement[i] = ((TPZInterpolatedElement *)(Sons[i]))->MeanSolution(fVariable);
	}
	REAL mean = 0., diff = 0.;
	for (i = 0; i < nsons; i++) {
		mean += MeansBySubElement[i];
	}
	mean /= nsons;
	for (i = 0; i < nsons; i++) {
		REAL diffi = mean - MeansBySubElement[i];
		if (diff < diffi)
			diff = diffi;
	}
	return diff;
}

/* Compute gradient with gradiente reconstruction by least square using neighboards. Retorna el numero de estimativas preenchidas no vetor Estimatives */
int TEEGradLaplacianRecWLS::ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives) {
	// Nada sera realizado para elementos con dimension diferente de la dimension del problema
	Estimatives.Fill(0.);
	if (!cel || cel->IsInterface()) return 0;
	int dim = cel->Mesh()->Dimension();
	if (cel->Dimension() != dim) return 0;
	int ncols = 5;
	if (dim == 3) ncols = 9;
	int j, k;
	int nstates = cel->Material()->NSolutionVariables(fVariable);

	TPZFMatrix<REAL> vectorrec(ncols, 1, 0.0);
	TPZFMatrix<REAL> gradrec(dim, 1, 0.0);
	TPZManVector<REAL, 3> center(3, 0.0), centerbeta(3, 0.0);
	TPZManVector<STATE> solalfa(nstates, 0.0), solbeta(nstates, 0.0);

	TPZStack<int64_t> realneighs;
	int nneighs = GettingNeighboardsWithoutDuplicates(cel, realneighs);
	if (!nneighs) return 0;  // 

	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> H(nneighs, ncols, 1.0);
	TPZFMatrix<REAL> F(nneighs, 1, 0.0);

	// Encontramos el centro del elemento corriente cel
	TPZGeoEl* gelalfa = cel->Reference();
	TPZManVector<REAL> centerpsi(gelalfa->Dimension(), 0.0);
	gelalfa->CenterPoint(gelalfa->NSides() - 1, centerpsi);
	gelalfa->X(centerpsi, center);
	TPZFNMatrix<3, REAL> gradXY(dim, 1, 0.);
	ComputingUhAndGradUhOnCenterOfTheElement(cel, centerpsi, solalfa, gradXY);

	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
	for (int ineighs = 0; ineighs < nneighs; ineighs++) {
		TPZGeoEl * gelbeta = cel->Mesh()->ElementVec()[realneighs[ineighs]]->Reference();
		if (!gelbeta)
			DebugStop();
		centerpsi.Resize(gelbeta->Dimension(), 0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides() - 1, centerpsi);
		gelbeta->X(centerpsi, centerbeta);
		gelbeta->Reference()->Solution(centerpsi, fVariable, solbeta);

		// Calculo da distancia do centro centerbeta para o x0(center)
		REAL dist = 0.L, weight;
		for (k = 0; k < dim; k++) {
			dist += (centerbeta[k] - center[k])*(centerbeta[k] - center[k]);
		}
		weight = (1.L / sqrt(dist));
		if (IsZero(dist)) DebugStop();   // Isto nao deve acontecer, não serão elementos sobrepostos para o caso
		for (j = 0; j < dim; j++) {
			H(ineighs, j) = weight*(centerbeta[j] - center[j]);
			H(ineighs, j + dim) = 0.5*weight*H(ineighs, j)*H(ineighs, j);
		}
		H(ineighs, 2 * dim) = weight * H(ineighs, 0)*H(ineighs, 1);
		if (dim == 3) {
			H(ineighs, 8) = weight * H(ineighs, 0)*H(ineighs, 2);
			H(ineighs, 9) = weight * H(ineighs, 1)*H(ineighs, 2);
		}
		F(ineighs, 0) = weight * (solbeta[fVariable] - solalfa[fVariable]);
	}

	int iestim = 0;
	SolveUsingLeastSquareMethod(H, F, vectorrec);
	// Formato completo do vectorrec = gradx, grady, gradz, D2xx, D2yy, D2zz, D2xy, D2xz , D2yz, potential
	for (j = 0; j < dim; j++)
		gradrec(j, 0) = vectorrec(j, 0);

	// norma do gradiente recuperado - primeira estimativa
	Estimatives[iestim++] = VectorNorm(gradrec);
	Estimatives[iestim++] = VectorNorm(gradXY);
	// third estimator: Computing fabs ( 1 - cos angle Grad_u_h e Grad_recovery )
	Estimatives[iestim++] = ComputingOneMinusCosAngleBetweenNormalVectors(gradrec, gradXY);
	// Estimando o potencial 
	iestim++;
	// Estimando o laplaciano
	TPZVec<REAL> ff(1);
	REAL laplacianrec = 0.;
	for (k = dim; k < 2 * dim; k++)
		laplacianrec += vectorrec(k, 0);
	cel->Material()->ForcingFunction()->Execute(center, ff);
	Estimatives[iestim++] = fabs(laplacianrec - ff[0]); // ff[0] is forcingf->Execute(center, ff);

	return iestim;
}
