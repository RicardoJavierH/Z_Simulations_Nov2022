/// Poisson example, with exact solutions as SinXSin, ArcTg, High oscillatory.
/// On 1D, 2D, 3D
/// Using H1, GD and mixed formulations. To Mixed formulations it considering enriched internal version

#include <iostream>
#include <cmath>
#include <time.h>
#include <stdio.h>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZReadGIDGrid.h"
#include "tpzgeoelrefpattern.h"

#include "TPZCompMeshTools.h"

#include "pzanalysis.h"
#include "pzstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "PoissonEquations.h"
#include "IdMaterialsToGIDFiles.h"
#include "CreateGeoMesh.h"
#include "Refinements.h"
#include "commonJC.h"

#include "hadaptive.h"
#include "EEGradientRecLS.h"
#include "EEGradientRecWLS.h"
#include "EEGradPotentialRecLS.h"
#include "EEGradPotentialRecWLS.h"
#include "EELaplacianRecLS.h"
#include "EELaplacianRecWLS.h"
#include "EEGradLaplacianRecLS.h"
#include "EEGradLaplacianRecWLS.h"
#include "EEGradPotentialLaplacianRecLS.h"
#include "EEGradPotentialLaplacianRecWLS.h"
#include "EELaplacianPotentialRecLS.h"
#include "EELaplacianPotentialRecWLS.h"

#include "pzvisualmatrix.h"

// extern variables
std::string Partition;

int DimProblem;

// PROGRAM
int main(int argc, char* argv[]) {
#if COMPUTETIME
	clock_t texec1 = clock();

	clock_t t1 = clock();
#endif
	std::string sout;
	sout = argv[0];
	PutNameDirectoryForResults(sout);

	// To load refinement pattern
	gRefDBase.InitializeAllUniformRefPatterns();

	// All information to make simulation
	SimulationData simulationdata;
	// Filename with geometric information from GID file
	std::string fromgid;

	// Importing user data from filedata.
	Partition = argv[0];
	std::ifstream *fimport = UserDataImport(Partition,fromgid, simulationdata);
	if (!fimport) {
		cout << "\nUser data is not found.\nThe executable need of the \"LJC_InputOutput\\userdata.txt\" file.\n\n";
		return 100;
	}
	// To print information of processing times, system dimension and calculated errors (H1, GD and HdivH1)
	std::ofstream saida(sout + "InfoRun_Errors.txt", ios::app);
	std::ofstream logError(sout+"LogErrors.xls");
	std::string namexls = sout + "LogErrors.xls";
	if (!saida.is_open() || !logError.is_open()) { std::cout << "\nLogErrors.xls or InfoRun is not open.\n"; return -3; }
	logError << "Modelo " << simulationdata.Equation << " \t Tolerance " << simulationdata.Tolerance << std::endl;

	if (simulationdata.applyadaptivity) {
		logError << "M_" << simulationdata.ErrorEstimator << "_T" << simulationdata.ErrorType;
		if (simulationdata.ErrorTypeEE2 > -1)
			logError << "_" << simulationdata.ErrorTypeEE2;
	}
	else
		logError << "NoAdaptive";
	logError << " \t Log(NDofs)" << " \t Log(L2_Error)" << " \t L2_Error \t NDofs \t " << "Log(NDofs)" << " \t " << "H1_Error" << " \t \t  " << "Log(NDofs)" << "\t " << "SemiNorm_Error" << std::endl;

	if (simulationdata.printrunninginfo)
		PrintRunningInformation(*fimport, saida);
	fimport->close(); if (fimport) delete fimport;

#if COMPUTETIME
	// end of time dedicated to initial processes
	clock_t t2 = clock();
	saida << "TIME PROCESS Initialization: \tTime elapsed = " << t2-t1 << std::endl;

	// start time dedicated to geometric mesh processed
	t1 = clock();
#endif
	/** Creating geometric mesh from GID file or not */
	TPZGeoMesh *gmesh = 0;
	if (fromgid.empty()) {
		switch(simulationdata.Equation) {
			case 3:
			case 4:
				gmesh = ConstructingFicheraCorner(1.,simulationdata,DimProblem);
				break;
			default:
				gmesh = CreatingGeometricMesh(simulationdata, DimProblem);
				break;
		}
	}
	else
		gmesh = CreatingGeometricMeshFromGIDFile(fromgid, simulationdata, DimProblem);
	if (!gmesh || !gmesh->NElements()) return -2;
	saida << "\nMODEL " << simulationdata.Equation << "\t Dimension " << DimProblem << "\t Element Type " << simulationdata.typeel << "\t Level of Refinement " << simulationdata.nref << "\n";
	std::cout << "\nMODEL " << simulationdata.Equation << std::endl;
	std::cout << " Dimension " << DimProblem << "\t Element Type " << simulationdata.typeel << "\t Level of Refinement " << simulationdata.nref << std::endl;

	// Refining mesh (uniform)
	UniformRefinement(simulationdata.nref, gmesh);
	// Printing geomesh only
	if (simulationdata.print) {
		std::stringstream gmeshfile;
		gmeshfile << sout << "GMesh_E" << simulationdata.typeel << "_Dim" << DimProblem << "D.vtk";
		std::ofstream fGMesh(gmeshfile.str());
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fGMesh);
	}

#if COMPUTETIME
	// end of time dedicated to geometric mesh processed
	t2 = clock();
	saida << "TIME PROCESS to geometric Mesh: \tTime elapsed = " << t2-t1 << std::endl;
#endif

	// creating computational mesh
	TPZCompMesh *cmesh = NULL;
	// Analysis
	TPZAnalysis* an = 0;
	TPZStepSolver<STATE> *step = 0;
	TPZManVector<REAL> erro(3,0.0L);

#if COMPUTETIME
	clock_t tmodel = clock();
#endif
	// Determining Equation problem and relationated exact solution (if exists)
	ChooseEquation(simulationdata, DimProblem);
	if (DimProblem == 1 && simulationdata.NumericMethod == 1) {
		std::cout << "GD is not implemented for 1D problems.\n";
		return 110;
	}
	saida << "  Solving with Method " << simulationdata.NumericMethod << "\n";
	std::cout << "  Solving with Method " << simulationdata.NumericMethod << "\t Linear Solver " << simulationdata.stepsolver << std::endl;

	//rodar formulacao mista
	REAL OldError = 100., NewError = 10.;
	int counter = 0, stepgrahp = 0;
	int64_t nequations = 0;

	if (simulationdata.NumericMethod < 2) {
		// To errors
		saida << "\nSaida do erro para formulacão ";
		if(simulationdata.NumericMethod==1) saida << "GD,";
		else saida << "H1,";
		saida << " com ordem p  = " << simulationdata.pOrder << "\n";
		if (!fromgid.empty())
			saida << "Malha irregular from gid file.\n";

#if COMPUTETIME
		t1 = clock();
#endif
		/// PRIORITARY -> RELATIVE TO COMPUTATIONAL MESH
		///Indicacao de malha DGFem. Se false, vamos criar malha H1
		if(cmesh) delete cmesh;
		cmesh = new TPZCompMesh(gmesh);
		TPZMaterial* material = CreatingComputationalMesh(cmesh, simulationdata, DimProblem);
		do {
			int ndofs = cmesh->NEquations();
			saida << "\n   Computational mesh created:\n" << "   NDofs = " << ndofs << "\t NDofs before condensed = " << cmesh->NEquations() << std::endl;
			std::cout << "\n   Computational mesh created:\n" << "   Grau de Liberdade Total = " << cmesh->NEquations() << std::endl;
#if COMPUTETIME
			t2 = clock();
			saida << "\nTIME PROCESS Creating CMesh : \tTime elapsed = " << t2 - t1 << std::endl;
#endif
			// PRIORITARY -> RELATIVE TO SOLVING THE PROBLEM
			int numthreads = 0;
			if(an) delete an;
			if(simulationdata.NumericMethod) an = new TPZAnalysis(cmesh,0);
			else an = new TPZAnalysis(cmesh);

			///Criando structural matrix and linear equation solver in analysis
			step = SetStructuralMatrixAndLinearEquationSolver(an,cmesh,simulationdata, numthreads);

			///Assemble da matriz de rigidez e vetor de carga
#if COMPUTETIME
			t1 = clock();
#endif
			an->Assemble();
#if COMPUTETIME
			t2 = clock();
#endif
#if COMPUTETIME
			saida << "\nTIME PROCESS Assembling : \tTime elapsed = " << t2 - t1 << std::endl;
			t1 = clock();
#endif
			saida << "Dimension of the linear system equations\n" << "K: " << an->Solver().Matrix()->Rows() << " x " << an->Solver().Matrix()->Cols() << "\n";

			///Resolução do sistema
			an->Solve();
#if COMPUTETIME
			t2 = clock();
			saida << "\nTIME PROCESS Solving : \tTime elapsed = " << t2 - t1 << std::endl;
#endif
			std::cout << "   Approximated solution computed." << std::endl;
			///Calculando erro da aproximacao
			if(simulationdata.computeUFluxErrors) {
#if COMPUTETIME
				t1 = clock();
#endif
				if(counter) OldError = erro[1];
				SaveUAndFluxErrors(an,erro,simulationdata,saida);
				std::cout << "   Approximated error computed." << std::endl;
#if COMPUTETIME
				t2 = clock();
				saida << "\nTIME PROCESS Computing L2 error : \tTime elapsed = " << t2 - t1 << std::endl;
#endif
				logError << (counter+1) << " \t  " << log(ndofs) << " \t " << log(erro[1]) << " \t " << erro[1] << " \t " << ndofs << " \t " << log(ndofs) << " \t " << log(erro[0]) << " \t \t " << log(ndofs) << " \t " << log(erro[2]) << std::endl;
			}
			// PRIORITARY -> RELATIVE TO PRINT SOLUTION TO VISUALIZATION BY PARAVIEW
			///Exportando para Paraview
			if (an) {
				std::stringstream ssout;
				ssout << sout;
				ssout << "Model" << simulationdata.Equation << "_" << DimProblem << "D_MESH_E" << simulationdata.typeel << "H" << simulationdata.nref << "_p" << simulationdata.pOrder;
				if (simulationdata.NumericMethod == 2) ssout << "_Hdiv";
				else if(simulationdata.NumericMethod == 3) ssout << "_HdivPlus";
				else if (simulationdata.NumericMethod==1) ssout << "_GD";
				else ssout << "_H1";
				ssout << ".vtk";
				
#if COMPUTETIME
				t1 = clock();
#endif
				an->SetStep(stepgrahp);
				an->DefineGraphMesh(DimProblem, simulationdata.sV, simulationdata.vV, ssout.str());
				int resolution = 0;
				an->PostProcess(resolution);
#if COMPUTETIME
				t2 = clock();
				saida << "\nTIME PROCESS writing vtk files : \tTime elapsed = " << t2 - t1 << std::endl;
#endif
				std::cout << "   File to visualization of the approximated variables was saved." << std::endl << std::endl;
			}
			// Doing hp adaptivity
			if(simulationdata.applyadaptivity && counter+1 < simulationdata.MaxIterations) {
				saida << "\n\nADAPTING PROCESS - Step " << counter+1 << std::endl;
				//Carregando o índice dos arquivos de saida para próximas rodadas
				stepgrahp = an->GetStep();
				// creando el estimador de error
				TErrorEstimator *estimator = 0;
				TAdaptive *adaptive = new TAdaptive(cmesh,simulationdata);
				estimator = adaptive->AllocatingEstimator(cmesh, simulationdata);
				if (!estimator) return -100;
#if COMPUTETIME
				t1 = clock();
#endif
				adaptive->Adapting(estimator,saida,sout);
#if COMPUTETIME
				t2 = clock();
				saida << "\nTIME PROCESS Adapting : \tTime elapsed = " << t2 - t1 << "\n\n";
#endif
				NewError = estimator->MaximeEstimatedFounded();
				if(adaptive) delete adaptive;
				if(estimator) delete estimator;
			}
			counter++;
			nequations = cmesh->NEquations();
		}while(simulationdata.applyadaptivity && NewError > simulationdata.Tolerance && !IsZero(OldError) && counter < simulationdata.MaxIterations && nequations < simulationdata.MaxEquations);
	}
	
	if (step) { delete step; step = 0; }

#if COMPUTETIME
	clock_t tmodel2 = clock();
	saida << "\nTIME PROCESS Model " << simulationdata.Equation << " : \tTime elapsed = " << tmodel2 - tmodel << std::endl;
#endif
	if(an) { delete an; an = 0; }
	if(cmesh) { delete cmesh; cmesh = 0; }

	if (gmesh) delete gmesh; gmesh = 0;
#if COMPUTETIME
	clock_t texec2 = clock();
	saida << "\nTIME PROCESS Simulation : \tTime elapsed = " << texec2 - texec1 << std::endl;
#endif

	saida.close();
	logError.close();
	ChangeDecimalPointByDecimalComma(namexls);

	return 0;
}

void (*fExact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) = 0;
//void PutExact(int n);
void PutExact(int n) {
	if(n==1) fExact = &SolExactSeno;
	else if(n==2) fExact = &SolExactSenoSeno;
	else if(n==3) fExact = &SolExactSenoSenoSeno;
	else if(n==4) fExact = &SolExactArcTg1D;
	else if(n==5) fExact = &SolExactArcTg2D;
	else if(n==6) fExact = &SolExactArcTg3D;
	else if(n==7) fExact = &SolExactStrongOsc1D;
	else if(n==8) fExact = &SolExactStrongOsc2D;
	else if(n==9) fExact = &SolExactStrongOsc3D;
	else fExact = 0;
}

