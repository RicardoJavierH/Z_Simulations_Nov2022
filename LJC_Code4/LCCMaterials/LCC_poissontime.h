/**
 * @file
 * @brief Contains the TPZMatPoisson3d class.
 */

#ifndef MATPOISSONTIMEH
#define MATPOISSONTIMEH

#include <iostream>
#include "pzpoisson3d.h"
#include "pzfmatrix.h"


#ifdef _AUTODIFF
#include "fadType.h"
#endif


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 * \f$ du/dt -fK Laplac(u) + fC * div(fConvDir*u) = - fXf  \f$
 */
class TLCCMatPoissonTime : public TPZMatPoisson3d {
	
	protected :
	
	/* Time interval - step to Euler method
	 */
	REAL fDeltaT;
	
	
public:
	
	/** @brief Sets the time step  */
	void SetStepTime(REAL deltat){ fDeltaT = deltat;}
	
	
	TLCCMatPoissonTime(int nummat, int dim);
    
	TLCCMatPoissonTime(int matid) : TPZRegisterClassId(&TLCCMatPoissonTime::ClassId),
    TPZMatPoisson3d(matid)
    {
    }
	
	TLCCMatPoissonTime();
	
	virtual ~TLCCMatPoissonTime();
	
    /** @brief Returns the number of state variables associated with the material */
	virtual int NStateVariables() const override
    {
        return 1;
    }
	
	virtual std::string Name() override { return "TLCCMatPoissonTime"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) override {
		TPZMatPoisson3d::Contribute(data,weight,ef);
	}
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
		TPZMatPoisson3d::Contribute(datavec,weight,ek,ef);
    }
	
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		TPZMatPoisson3d::ContributeBC(datavec,weight,ek,ef,bc);
    }

#ifdef _AUTODIFF
	/** @brief Computes contribution to the energy at an integration point */
	void ContributeEnergy(TPZVec<REAL> &x,
						  TPZVec<FADFADREAL> &sol,
						  TPZVec<FADFADREAL> &dsol,
						  FADFADREAL &U,
						  REAL weight);
#endif
	
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		TPZMatPoisson3d::ContributeBC(data,weight,ef,bc);
	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	

	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBCInterface(TPZMaterialData &data,TPZMaterialData &dataleft,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
	{
		TPZMatPoisson3d::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data,TPZMaterialData &dataleft,TPZMaterialData &dataright,REAL weight,
									 TPZFMatrix<STATE> &ef) override
	{
		TPZMatPoisson3d::ContributeInterface(data,dataleft,dataright,weight,ef);
	}

#ifdef _AUTODIFF
//	virtual void ContributeBCEnergy(TPZVec<REAL> &x,TPZVec<FADFADREAL> &sol, FADFADREAL &U,
//									REAL weight, TPZBndCond &bc);
#endif

	/** @} */
    
protected:
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override {
		TPZMatPoisson3d::Solution(datavec,var,Solout);
    }
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)  override {
		TPZMatPoisson3d::FillDataRequirements(datavec);
    }
	virtual void FillDataRequirements(TPZMaterialData &data) override
	{
		TPZMatPoisson3d::FillDataRequirements(data);
		data.fNeedsSol = true;
	}

	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
	
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override {
		TPZMatPoisson3d::Errors(data,u_exact,du_exact,errors);
    }
    
	virtual int NEvalErrors()  override {return 6;}
	
	virtual int IsInterfaceConservative() override { return 1;}
	
	virtual int ClassId() const override;

	
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	virtual void Read(TPZStream &buf, void *context) override;

};

#endif
