/**
 * @file
 * @brief Contains implementations of the TPZMatPoisson3d methods.
 */

#include "pzpoisson3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include "LCC_poissontime.h"
#include "pzlog.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.poisson3d"));
#endif


using namespace std;

TLCCMatPoissonTime::TLCCMatPoissonTime(int nummat, int dim) : TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
TPZMatPoisson3d(nummat) {
	fXf = 0.;
	fK = 1.;
	fC = 0.;
	fConvDir[0] = fConvDir[1] = fConvDir[2] = 0.;
	fSymmetry = 1.;
	fShapeHdiv = false;
	fPenaltyConstant = 1000.;
	fDim = dim;
	fSD = 0.;
	fDeltaT = 1.;
	fNeumann = false;
}

TLCCMatPoissonTime::TLCCMatPoissonTime():TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
TPZMatPoisson3d() {
	fXf = 0.;
	fDim = 1;
	fSD = 0.;
	fDeltaT = 1.0;
	fNeumann = false;
}


TLCCMatPoissonTime::~TLCCMatPoissonTime() {
}

void TLCCMatPoissonTime::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	
//    if(data.numberdualfunctions)
//    {
//        ContributeHDiv(data , weight , ek, ef);
//
//        return;
//    }
    
    if(fNeumann){
        
        LocalNeumanContribute(data , weight , ek, ef);
        
        return;
    }

    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();

    STATE fXfLoc = fXf;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fXfLoc = res[0];
    }

    REAL delx = 0.;
    STATE ConvDirAx[3] = {0.};
    if(fC != 0.0) {
        int di,dj;
        delx = 0.;
        for(di=0; di<fDim; di++) {
            for(dj=0; dj<fDim; dj++) {
                delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
            }
        }
        delx = 2./delx;
        
        
        switch(fDim) {
            case 1:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                break;
            case 2:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
                break;
            case 3:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
                ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
                break;
            default:
                PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
        }
    }

//    std::cout << " fxfloc " << fXfLoc << " weight " << weight << " prod " << fXfLoc*weight << std::endl;
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        STATE dphiic = 0;
        for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*(STATE)dphi(kd,in);
        ef(in, 0) += (STATE)weight * (-(fDeltaT * fXfLoc) + data.sol[0][0] ) * ( (STATE)phi(in,0) + (STATE)(0.5*delx*fC)*fSD*dphiic );
        for( int jn = 0; jn < phr; jn++ ) {
			ek(in, jn) += (STATE)weight * (STATE)(phi(in, 0) * phi(jn, 0));
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight * (
                                       +fK * fDeltaT *(STATE)( dphi(kd,in) * dphi(kd,jn) ) 
                                       - (STATE)(fC* dphi(kd,in) * phi(jn)) * ConvDirAx[kd]
                                       + (STATE)(0.5 * delx * fC * dphi(kd,jn)) * fSD * dphiic * ConvDirAx[kd]
                                       );
            }
        }
    }

    if (this->IsSymetric()){    
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
}
/*
void TPZMatPoissonTime::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

	if (data.numberdualfunctions)
	{
		ContributeHDiv(data, weight, ek, ef);

		return;
	}

	if (fNeumann) {

		LocalNeumanContribute(data, weight, ek, ef);

		return;
	}

	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZVec<REAL>  &x = data.x;
	TPZFMatrix<REAL> &axes = data.axes;
	TPZFMatrix<REAL> &jacinv = data.jacinv;
	int phr = phi.Rows();

	STATE fXfLoc = fXf;

	if (fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<STATE, 1> res(1);
		TPZFMatrix<STATE> dres(Dimension(), 1);
		fForcingFunction->Execute(x, res, dres);       // dphi(i,j) = dphi_j/dxi
		fXfLoc = res[0];
	}

	REAL delx = 0.;
	STATE ConvDirAx[3] = { 0. };
	if (fC != 0.0) {
		int di, dj;
		delx = 0.;
		for (di = 0; di < fDim; di++) {
			for (dj = 0; dj < fDim; dj++) {
				delx = (delx < fabs(jacinv(di, dj))) ? fabs(jacinv(di, dj)) : delx;
			}
		}
		delx = 2. / delx;


		switch (fDim) {
		case 1:
			ConvDirAx[0] = axes(0, 0)*fConvDir[0] + axes(0, 1)*fConvDir[1] + axes(0, 2)*fConvDir[2];
			break;
		case 2:
			ConvDirAx[0] = axes(0, 0)*fConvDir[0] + axes(0, 1)*fConvDir[1] + axes(0, 2)*fConvDir[2];
			ConvDirAx[1] = axes(1, 0)*fConvDir[0] + axes(1, 1)*fConvDir[1] + axes(1, 2)*fConvDir[2];
			break;
		case 3:
			ConvDirAx[0] = axes(0, 0)*fConvDir[0] + axes(0, 1)*fConvDir[1] + axes(0, 2)*fConvDir[2];
			ConvDirAx[1] = axes(1, 0)*fConvDir[0] + axes(1, 1)*fConvDir[1] + axes(1, 2)*fConvDir[2];
			ConvDirAx[2] = axes(2, 0)*fConvDir[0] + axes(2, 1)*fConvDir[1] + axes(2, 2)*fConvDir[2];
			break;
		default:
			PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
		}
	}

	//    std::cout << " fxfloc " << fXfLoc << " weight " << weight << " prod " << fXfLoc*weight << std::endl;
		//Equacao de Poisson
	for (int in = 0; in < phr; in++) {
		int kd;
		STATE dphiic = 0;
		for (kd = 0; kd < fDim; kd++) dphiic += ConvDirAx[kd] * (STATE)dphi(kd, in);
		ef(in, 0) += (STATE)weight * (fXfLoc + data.sol[0][0]) * ((STATE)phi(in, 0) + (STATE)(0.5*delx*fC)*fSD*dphiic);
		for (kd = 0; kd < fDim; kd++) {
			ef(in, 0) += (STATE)weight * (0.5* fK * fDeltaT *(STATE)(data.dsol[0](kd, 0) * dphi(kd, 0)));
		}
		for (int jn = 0; jn < phr; jn++) {
			ek(in, jn) += (STATE)weight * (STATE)(phi(in, 0) * phi(jn, 0));
			for (kd = 0; kd < fDim; kd++) {
				ek(in, jn) += (STATE)weight * (
					+ 0.5* fK * fDeltaT *(STATE)(dphi(kd, in) * dphi(kd, jn))
					- (STATE)(fC* dphi(kd, in) * phi(jn)) * ConvDirAx[kd]
					+ (STATE)(0.5 * delx * fC * dphi(kd, jn)) * fSD * dphiic * ConvDirAx[kd]
					);
			}
		}
	}

	if (this->IsSymetric()) {
		if (!ek.VerifySymmetry()) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}*/

void TLCCMatPoissonTime::ContributeBC(TPZMaterialData &data,REAL weight,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
//	if(fShapeHdiv==true)//data.fShapeType == TPZMaterialData::EVecandShape || data.fShapeType == TPZMaterialData::EVecShape)
//	{
		
//		ContributeBCHDiv(data , weight , ek, ef, bc);
		
		
//		return;
//	}
	
	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &axes = data.axes;
	int phr = phi.Rows();
	short in,jn;
	STATE v2[1];
	v2[0] = bc.Val2()(0,0);
	
	if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
		TPZManVector<STATE,1> res(1);
		bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
//        if(fabs(res[0]) > 1.e-6)
//        {
//            bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
//            DebugStop();
//        }
		v2[0] = res[0];
	}

	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += (STATE)(gBigNumber* fDeltaT * phi(in,0) * weight) * v2[0];
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) += gBigNumber * fDeltaT * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		case 1 :			// Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in,0) += v2[0] * (STATE)(phi(in,0) * fDeltaT * weight);
			}
			break;
		case 2 :		// mixed condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in, 0) += v2[0] * (STATE)(phi(in, 0) * weight);
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi(in,0) * phi(jn,0) * weight);     // peso de contorno => integral de contorno
				}
			}
			break;
		case 3: // outflow condition
			int id, il, jl;
			REAL normal[3];
			if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
			if (fDim == 2){
				normal[0] = axes(0,1);
				//normal[1] = axes(1,1);
			}
			if (fDim == 3){
				normal[0] = axes(0,2);
				normal[1] = axes(1,2);
				normal[2] = axes(2,2);
			}
			REAL ConvNormal = 0.;    
			for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];  
			if(ConvNormal > 0.) {
				for(il=0; il<phr; il++) {
					for(jl=0; jl<phr; jl++) {
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
					}
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
			break;
	}
	
	if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
		if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}


void TLCCMatPoissonTime::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
	
	TPZVec<STATE> pressure(1);
	TPZVec<REAL> pto(3);
	TPZFMatrix<STATE> flux(3,1);
	
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
   // Solout.Resize(this->NSolutionVariables(var));
    
#ifndef STATE_COMPLEX
	
	switch (var) {
	/*	case 7:
		{
			//			{ //MinusKGradU
			int id;
			TPZManVector<STATE> dsolp(3,0);
			dsolp[0] = data.dsol[0](0,0)*data.axes(0,0)+data.dsol[0](1,0)*data.axes(1,0);
			dsolp[1] = data.dsol[0](0,0)*data.axes(0,1)+data.dsol[0](1,0)*data.axes(1,1);			
			dsolp[2] = data.dsol[0](0,0)*data.axes(0,2)+data.dsol[0](1,0)*data.axes(1,2);
			for(id=0 ; id<fDim; id++) 
			{
				Solout[id] = -1. * this->fK * dsolp[id];
			}
		}
			break;*/
		case 8:
			Solout[0] = data.p;
			break;
		case 10:
			if (data.numberdualfunctions) {
				
				Solout[0]=data.sol[0][0];
				Solout[1]=data.sol[0][1];
                Solout[2]=data.sol[0][2];
				
			}
			else {
				this->Solution(data.sol[0], data.dsol[0], data.axes, 2, Solout);
			}
			
			break;
            
        case 21:
            for(int k=0;k<fDim;k++){
                Solout[k]=data.sol[0][k];
            }
			break;
            
		case 11:
			if (data.numberdualfunctions) {
				Solout[0]=data.sol[0][2];
			}
			else{
				Solout[0]=data.sol[0][0];
			}
			break;
			
		case 12:
//				fExactSol->Execute(data.x,pressure,flux);
			fForcingFunctionExact->Execute(data.x,pressure,flux);

				Solout[0]=pressure[0];
			break;
		case 13:
//				fExactSol->Execute(data.x,pressure,flux);
			fForcingFunctionExact->Execute(data.x,pressure,flux);

				Solout[0]=flux(0,0);
				Solout[1]=flux(1,0);
            break;
            
        case 14:
        {
			if (data.numberdualfunctions){
				Solout[0]=data.sol[0][data.sol[0].NElements()-1];
			}else{
                //Solout[0]=data.dsol[0](0,0)+data.dsol[0](1,1)+data.dsol[0](2,2);
                STATE val = 0.;
                for(int i=0; i<fDim; i++){
                    val += data.dsol[0](i,i);
                }
                Solout[0] = val;
            }
        }
            break;
          
        case 15:
        {
//            fExactSol->Execute(data.x,pressure,flux);
			  fForcingFunctionExact->Execute(data.x,pressure, flux);
			  Solout[0]=flux(fDim,0);
        }
            break;

            
        case 16:
            if (data.numberdualfunctions) {
					Solout[0]=data.sol[0][2];
            }
            else {
                std::cout<<"Pressao somente em Omega1"<<std::endl;
                Solout[0]=0;//NULL;
            }
				
            break;
        
        case 17:
            if (!data.numberdualfunctions) {
                Solout[0]=data.sol[0][0];
            }
            else {
                std::cout<<"Pressao somente em omega2"<<std::endl;
                Solout[0]=0;//NULL;
            }
				
            break;
        case 18:
            if( data.numberdualfunctions){
                Solout[0]=data.sol[0][0];//fluxo de omega1
                Solout[1]=data.sol[0][1];
                //	Solout[2]=data.sol[2];
                return;
            }
        
        case 19:
            if(data.numberdualfunctions){
                Solout[0]=data.dsol[0](0,0);//fluxo de omega1
                Solout[1]=data.dsol[0](1,0);
                Solout[2]=data.dsol[0](2,0);
                return;
            }
        case 20:
            if( data.numberdualfunctions){
                Solout[0]=data.dsol[0](0,1);//fluxo de omega1
                Solout[1]=data.dsol[0](1,1);
                Solout[2]=data.dsol[0](2,1);
                return;
            }
            else {
                std::cout<<"Pressao somente em omega2"<<std::endl;
                Solout[0]=0;//NULL;
            }
            break;
        default:
           
            if (data.sol[0].size() == 4) {
                
                data.sol[0][0] = data.sol[0][2];
            }

            this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
            break;
    }
#endif
}

#include "pzaxestools.h"
void TLCCMatPoissonTime::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
#ifndef STATE_COMPLEX
	Solout.Resize( this->NSolutionVariables( var ) );
	
	if(var == 1){
		Solout[0] = Sol[0];//function
		return;
	}
	if(var == 2) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9,STATE> dsoldx;
			TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
			Solout[id] = dsoldx(id,0);//derivate
		}
		return;
	}//var == 2
	if (var == 3){ //KDuDx
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(0,0) * this->fK;
		return;
	}//var ==3
	if (var == 4){ //KDuDy
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(1,0) * this->fK;
		return;
	}//var == 4 
	if (var == 5){ //KDuDz
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(2,0) * this->fK;
		return;
	}//var == 5
	if (var == 6){ //NormKDu
		int id;
		REAL val = 0.;
		for(id=0 ; id<fDim; id++){
			val += (DSol(id,0) * this->fK) * (DSol(id,0) * this->fK);
		}
		Solout[0] = sqrt(val);
		return;
	}//var == 6
	if (var == 7){ //MinusKGradU
		int id;
		//REAL val = 0.;
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		for(id=0 ; id<fDim; id++) {
			Solout[id] = -1. * this->fK * dsoldx(id,0);
		}
		return;
	}//var == 7  
	if(var == 9){//Laplac
		Solout.Resize(1);
		Solout[0] = DSol(2,0);
		return;
	}//Laplac
	
#endif
	TPZMaterial::Solution(Sol, DSol, axes, var, Solout);
	
}//method


void TLCCMatPoissonTime::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                             REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	
	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	
	TPZFNMatrix<660> dphiL, dphiR;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	
	int &LeftPOrder=dataleft.p;
	int &RightPOrder=dataright.p;
	
	REAL &faceSize=data.HSize;
	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	
	//Convection term
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id] * normal[id];
	if(ConvNormal > 0.) {
		for(il=0; il<nrowl; il++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
			}
		}
		for(ir=0; ir<nrowr; ir++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl);
			}
		}
	} else {
		for(ir=0; ir<nrowr; ir++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr);
			}
		}
		for(il=0; il<nrowl; il++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr);
			}
		}
	}
	
	if(IsZero(fK)) return;
	//diffusion term
	STATE leftK, rightK;
	leftK  = this->fK;
	rightK = this->fK;
	
	// 1) phi_I_left, phi_J_left
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += (STATE)(weight * fDeltaT * ( this->fSymmetry * (0.5)*dphiLinormal*phiL(jl,0)-(0.5)*dphiLjnormal*phiL(il,0))) * leftK;
		}
	}
	
	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += (STATE)(weight * fDeltaT * (this->fSymmetry * ((-0.5) * dphiRinormal * phiR(jr) ) + (0.5) * dphiRjnormal * phiR(ir))) * rightK;
		}
	}
	
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += (STATE)weight * fDeltaT * ((STATE)fSymmetry * ((STATE)((-0.5) * dphiLinormal * phiR(jr)) * leftK ) + (STATE)((-0.5) * dphiRjnormal * phiL(il))* rightK );
		}
	}
	
	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += (STATE)weight * fDeltaT * (
										 (STATE)(fSymmetry * (0.5) * dphiRinormal * phiL(jl)) * rightK + (STATE)((0.5) * dphiLjnormal * phiR(ir)) * leftK
										 );
		}
	}
	
	if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
	
	if (this->fPenaltyConstant == 0.) return;
	
	leftK  = this->fK;
	rightK = this->fK;
	//penalty = <A p^2>/h 
	REAL penalty = fPenaltyConstant * (0.5 * (abs(leftK)*LeftPOrder*LeftPOrder + abs(rightK)*RightPOrder*RightPOrder)) / faceSize;
	
	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
		
		// 1) left i / left j
		for(il=0; il<nrowl; il++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr; ir++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl; il++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr; ir++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
			}
		}
		
	}
	
	if (this->fPenaltyType == EFluxPenalty || this->fPenaltyType == EBoth){
		
		REAL NormalFlux_i = 0.;
		REAL NormalFlux_j = 0.;
		
		// 1) left i / left j
		for(il=0; il<nrowl; il++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiL(id,il)*normal[id];
			}
			for(jl=0; jl<nrowl; jl++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiL(id,jl)*normal[id];
				}
				ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr; ir++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiR(id,ir)*normal[id];
			}
			for(jr=0; jr<nrowr; jr++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiR(id,jr)*normal[id];
				}      
				ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl; il++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiL(id,il)*normal[id];
			}
			for(jr=0; jr<nrowr; jr++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiR(id,jr)*normal[id];
				}      
				ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr; ir++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiR(id,ir)*normal[id];
			}
			for(jl=0; jl<nrowl; jl++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiL(id,jl)*normal[id];
				}
				ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
			}
		}
		
	}  
	
}

/** Termos de penalidade. */
void TLCCMatPoissonTime::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                            REAL weight,
                                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder= dataleft.p;
	REAL faceSize=data.HSize;
	
	//  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];
	switch(bc.Type()) {
		case 0: // DIRICHLET
			
			//Diffusion
			for(il=0; il<nrowl; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL(id,il)*normal[id];
				}
				ef(il,0) += (STATE)(weight*dphiLinormal*fSymmetry)*fK*bc.Val2()(0,0);
				for(jl=0; jl<nrowl; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL(id,jl)*normal[id];
					}
					ek(il,jl) += (STATE)(weight*(fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0)))*fK;
				}
			}
			
			//Convection
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
					}
				}
			} else {
				for(il=0; il<nrowl; il++) {
					ef(il,0) -= (STATE)(weight * ConvNormal * phiL(il)) * bc.Val2()(0,0);
				}
			}
			
			break;
			
		case 1: // Neumann
			for(il=0; il<nrowl; il++) {
				ef(il,0) += (STATE)(weight*phiL(il,0))*bc.Val2()(0,0);
			}
			break;
			
		case 3: // outflow condition
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
					}
				}
			}
			else {
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
			break;
			
		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}
    if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
	
	if (this->fPenaltyConstant == 0.) return;
	
	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){  
		nrowl = phiL.Rows(); 
		const REAL penalty = fPenaltyConstant * abs(fK) * POrder * POrder / faceSize; //Ap^2/h
		REAL outflow = 0.;
		for(il=0; il<fDim; il++) outflow += fC * fConvDir[il] * normal[il];
		
		
		switch(bc.Type()) {
			case 0: // DIRICHLET  
				for(il=0; il<nrowl; il++) {
					ef(il,0) += (STATE)(weight * penalty * phiL(il,0)) * bc.Val2()(0,0);
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
					}
				}
				
				break;
			case 1: // Neumann
				if(outflow > 0.)
				{
					for(il=0; il<nrowl; il++)
					{
						for(jl=0; jl<nrowl; jl++)
						{
							ek(il,jl) += weight * outflow * phiL(il,0) * phiL(jl,0);
						}
					}
				}
				//nothing to be done
				break;
			default:
				PZError << "TPZMatPoisson3d::Wrong boundary condition type\n";
				break;
		}
        
	}
	
}


void TLCCMatPoissonTime::Write(TPZStream &buf, int withclassid) const{
	TPZDiscontinuousGalerkin::Write(buf, withclassid);
	buf.Write(&fXf, 1);
	buf.Write(&fDim, 1);
	buf.Write(&fK, 1);
	buf.Write(&fC, 1);
	buf.Write(fConvDir, 3);
	buf.Write(&fSymmetry, 1);
	buf.Write(&fSD, 1);
	buf.Write(&fPenaltyConstant,1);
	buf.Write(&gAlfa, 1);
}

void TLCCMatPoissonTime::Read(TPZStream &buf, void *context){
	TPZDiscontinuousGalerkin::Read(buf, context);
	buf.Read(&fXf, 1);
	buf.Read(&fDim, 1);
	buf.Read(&fK, 1);
	buf.Read(&fC, 1);
	buf.Read(fConvDir, 3);
	buf.Read(&fSymmetry, 1);
	buf.Read(&fSD, 1);
	buf.Read(&fPenaltyConstant,1);
	buf.Read(&gAlfa, 1);
}

int TLCCMatPoissonTime::ClassId() const{
    return Hash("TPZMatPoissonTime") ^ TPZMatPoisson3d::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TLCCMatPoissonTime>;
#endif
