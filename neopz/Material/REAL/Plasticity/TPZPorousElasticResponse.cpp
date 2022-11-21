//
//  TPZPorousElasticResponse.cpp
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#include "TPZPorousElasticResponse.h"


int TPZPorousElasticResponse::ClassId() const {
    return Hash("TPZPorousElasticResponse");
}

TPZPorousElasticResponse::TPZPorousElasticResponse() {
    
    m_kappa = 0.0;
    m_pt_el = 0.0;
    m_e_0   = 0.0;
    m_p_0   = 0.0;
    m_nu    = 0.0;
    m_mu    = 0.0;
    m_is_G_constant_Q   = false;
    m_plane_stress_Q    = false;
}

TPZPorousElasticResponse::TPZPorousElasticResponse(const TPZPorousElasticResponse & other) {
    
    m_kappa = other.m_kappa;
    m_pt_el = other.m_pt_el;
    m_e_0   = other.m_e_0;
    m_p_0   = other.m_p_0;
    m_nu    = other.m_nu;
    m_mu    = other.m_mu;
    m_is_G_constant_Q   = other.m_is_G_constant_Q;
    m_plane_stress_Q    = other.m_plane_stress_Q;
}

TPZPorousElasticResponse & TPZPorousElasticResponse::operator=(const TPZPorousElasticResponse & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_kappa = other.m_kappa;
    m_pt_el = other.m_pt_el;
    m_e_0   = other.m_e_0;
    m_p_0   = other.m_p_0;
    m_nu    = other.m_nu;
    m_mu    = other.m_mu;
    m_is_G_constant_Q   = other.m_is_G_constant_Q;
    m_plane_stress_Q    = other.m_plane_stress_Q;
    return *this;
}

void TPZPorousElasticResponse::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&m_kappa);
    buf.Write(&m_pt_el);
    buf.Write(&m_e_0);
    buf.Write(&m_p_0);
    buf.Write(&m_nu);
    buf.Write(&m_mu);
    buf.Write(m_is_G_constant_Q);
    buf.Write(m_plane_stress_Q);
}

void TPZPorousElasticResponse::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&m_kappa);
    buf.Read(&m_pt_el);
    buf.Read(&m_e_0);
    buf.Read(&m_p_0);
    buf.Read(&m_nu);
    buf.Read(&m_mu);
    buf.Read(m_is_G_constant_Q);
    buf.Read(m_plane_stress_Q);
}

void TPZPorousElasticResponse::SetPorousElasticity(STATE kappa, STATE pt_el, STATE e_0, STATE p_0)
{
    m_kappa = kappa;
    m_pt_el = pt_el;
    m_e_0 = e_0;
    m_p_0 = p_0;
}

void TPZPorousElasticResponse::Setp_0(STATE p_0){
    m_p_0 = p_0;
}

void TPZPorousElasticResponse::Sete_0(STATE e_0){
    m_e_0 = e_0;
}

void TPZPorousElasticResponse::SetShearModulusConstant(STATE G){
    m_is_G_constant_Q = true;
    m_mu = G;
}

void TPZPorousElasticResponse::SetPoissonRatioConstant(STATE nu){
    m_is_G_constant_Q = false;
    m_nu = nu;
}

void TPZPorousElasticResponse::SetPlaneStrain(){
    m_plane_stress_Q = false;
}

void TPZPorousElasticResponse::SetPlaneStress(){
    m_plane_stress_Q = true;
}

const char * TPZPorousElasticResponse::Name() const {
    return "TPZPorousElasticResponse";
}

void TPZPorousElasticResponse::Print(std::ostream & out) const {
    out << this->Name();
    out << "\n Logarithmic bulk modulus = " << m_kappa;
    out << "\n Elastic tensile strength = " << m_pt_el;
    out << "\n Initial void ratio = " << m_e_0;
    out << "\n Initial hydrostatic stress = " << m_p_0;
    out << "\n Constant Shear modulus directive = " << m_is_G_constant_Q;
    if (m_is_G_constant_Q) {
        out << "\n Second Lamé parameter (Shear modulus) = " << m_mu;
    }else{
        out << "\n Poisson ratio = " << m_nu;
    }
    out << "\n Plane stress state directive = " << m_plane_stress_Q;
}

void TPZPorousElasticResponse::G(const TPZTensor<STATE> &epsilon, STATE & G, STATE & dG_desp_vol) const{

    STATE epsv = epsilon.I1();
    G = (3*(1 + m_e_0)*(1 + epsv)*(1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*m_kappa*(1 + m_nu));
    
    dG_desp_vol = (-3*pow(1 + m_e_0,2)*(1 + epsv)*
                  (1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*
     pow(m_kappa,2)*(1 + m_nu)) +
    (3*(1 + m_e_0)*(1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*
     m_kappa*(1 + m_nu));
}

void TPZPorousElasticResponse::Poisson(const TPZTensor<STATE> &epsilon, STATE & nu, STATE & dnu_desp_vol) const{
    
    if (m_is_G_constant_Q) {
        STATE lambda, K, dK_desp_vol, dnu_dK;
        this->K(epsilon, K, dK_desp_vol);
        nu = -1 + (9*K)/(2.*(3*K + m_mu));
        dnu_dK = (9*m_mu)/(2.*pow(3*K + m_mu,2));
        dnu_desp_vol = dnu_dK * dK_desp_vol;
    }else{
        nu = m_nu;
        dnu_desp_vol = 0.0;
    }
    
}

void TPZPorousElasticResponse::Poisson_linearized(const TPZTensor<STATE> &epsilon_ref,const TPZTensor<STATE> &epsilon, STATE & nu) const{
    
    STATE epsv = epsilon.I1();
    STATE epsv_ref = epsilon_ref.I1();
    nu = -1 + (9*(1 + m_e_0)*(m_p_0 + m_pt_el)*
          (m_mu*exp((epsv_ref*(1 + m_e_0))/m_kappa)*
           (-((epsv - epsv_ref)*(1 + epsv_ref)*
              (1 + m_e_0)) + (1 + epsv)*m_kappa) +
           3*pow(1 + epsv_ref,2)*(1 + m_e_0)*(m_p_0 + m_pt_el)))/
    (2.*pow(m_mu*exp((epsv_ref*(1 + m_e_0))/m_kappa)*m_kappa +
              3*(1 + epsv_ref)*(1 + m_e_0)*(m_p_0 + m_pt_el),2));
}

void TPZPorousElasticResponse::K(const TPZTensor<STATE> &epsilon, STATE & K, STATE & dK_desp_vol) const{
    
    STATE epsv = epsilon.I1();
    K = ((1 + epsv)*(1 + m_e_0)*(m_p_0 + m_pt_el))/(exp((epsv*(1 + m_e_0))/m_kappa)*m_kappa);
    
    dK_desp_vol = -(((1 + m_e_0)*(1 + epsv + m_e_0 + epsv*m_e_0 - m_kappa)*
                     (m_p_0 + m_pt_el))/(exp((epsv*(1 + m_e_0))/m_kappa)*pow(m_kappa,2)));
}


void TPZPorousElasticResponse::De_G_constant(const TPZTensor<STATE> & epsilon, TPZFMatrix<STATE> & De) const{
    
    STATE lambda, K, dK_desp_vol;
    this->K(epsilon, K, dK_desp_vol);
    lambda = K - (2.0/3.0)*m_mu;
    
    // Line 0
    De.PutVal(_XX_,_XX_, lambda + 2. * m_mu);
    De.PutVal(_XX_,_YY_, lambda);
    De.PutVal(_XX_,_ZZ_, lambda);
    
    // Line 1
    De.PutVal(_XY_,_XY_, 2. * m_mu);
    
    // Line 2
    De.PutVal(_XZ_,_XZ_, 2. * m_mu);
    
    // Line 3
    De.PutVal(_YY_,_XX_, lambda);
    De.PutVal(_YY_,_YY_, lambda + 2. * m_mu);
    De.PutVal(_YY_,_ZZ_, lambda);
    
    // Line 4
    De.PutVal(_YZ_,_YZ_, 2. * m_mu);
    
    // Line 5
    De.PutVal(_ZZ_,_XX_, lambda);
    De.PutVal(_ZZ_,_YY_, lambda);
    De.PutVal(_ZZ_,_ZZ_, lambda + 2. * m_mu);
    
    /// Nonlinear correction
    TPZFMatrix<STATE> De_nl(6,6,0.0);
    REAL constant = ( epsilon.XX() + epsilon.YY() + epsilon.ZZ() ) * dK_desp_vol;
    
    // Line 0
    De_nl.PutVal(_XX_, _XX_, 1.0);
    De_nl.PutVal(_XX_, _YY_, 1.0);
    De_nl.PutVal(_XX_, _ZZ_, 1.0);
    
    // Line 3
    De_nl.PutVal(_YY_, _XX_, 1.0);
    De_nl.PutVal(_YY_, _YY_, 1.0);
    De_nl.PutVal(_YY_, _ZZ_, 1.0);
    
    // Line 5
    De_nl.PutVal(_ZZ_, _XX_, 1.0);
    De_nl.PutVal(_ZZ_, _YY_, 1.0);
    De_nl.PutVal(_ZZ_, _ZZ_, 1.0);
    
    De+=constant*De_nl;
}

void TPZPorousElasticResponse::De_Poisson_constant(const TPZTensor<STATE> & epsilon, TPZFMatrix<STATE> & De) const{
    
    STATE lambda, G, dG_desp_vol;
    this->G(epsilon, G, dG_desp_vol);
    lambda = (2.0*G*m_nu)/(1.0-2.0*m_nu);
    
    // Line 0
    De.PutVal(_XX_,_XX_, lambda + 2. * G);
    De.PutVal(_XX_,_YY_, lambda);
    De.PutVal(_XX_,_ZZ_, lambda);
    
    // Line 1
    De.PutVal(_XY_,_XY_, 2. * G);
    
    // Line 2
    De.PutVal(_XZ_,_XZ_, 2. * G);
    
    // Line 3
    De.PutVal(_YY_,_XX_, lambda);
    De.PutVal(_YY_,_YY_, lambda + 2. * G);
    De.PutVal(_YY_,_ZZ_, lambda);
    
    // Line 4
    De.PutVal(_YZ_,_YZ_, 2. * G);
    
    // Line 5
    De.PutVal(_ZZ_,_XX_, lambda);
    De.PutVal(_ZZ_,_YY_, lambda);
    De.PutVal(_ZZ_,_ZZ_, lambda + 2. * G);
    
    /// Nonlinear correction
    TPZFMatrix<STATE> De_nl(6,6,0.0);
    REAL constant = (2.0/(-1.0+2.0*m_nu));
    
    // Line 0
    REAL l0_val = constant * ( epsilon.XX() * (m_nu - 1.0) - m_nu * (epsilon.YY() + epsilon.ZZ())) * dG_desp_vol;
    De_nl.PutVal(_XX_, _XX_, l0_val);
    De_nl.PutVal(_XX_, _YY_, l0_val);
    De_nl.PutVal(_XX_, _ZZ_, l0_val);
    
    // Line 1
    REAL l1_val = 2.0 * dG_desp_vol * epsilon.XY();
    De_nl.PutVal(_XY_, _XX_, l1_val);
    De_nl.PutVal(_XY_, _YY_, l1_val);
    De_nl.PutVal(_XY_, _ZZ_, l1_val);
    
    // Line 2
    REAL l2_val = 2.0 * dG_desp_vol * epsilon.XZ();
    De_nl.PutVal(_XZ_, _XX_, l2_val);
    De_nl.PutVal(_XZ_, _YY_, l2_val);
    De_nl.PutVal(_XZ_, _ZZ_, l2_val);
    
    // Line 3
    REAL l3_val = constant * ( epsilon.YY() * (m_nu - 1.0) - m_nu * (epsilon.XX() + epsilon.ZZ())) * dG_desp_vol;
    De_nl.PutVal(_YY_, _XX_, l3_val);
    De_nl.PutVal(_YY_, _YY_, l3_val);
    De_nl.PutVal(_YY_, _ZZ_, l3_val);
    
    // Line 4
    REAL l4_val = 2.0 * dG_desp_vol * epsilon.YZ();
    De_nl.PutVal(_YZ_, _XX_, l4_val);
    De_nl.PutVal(_YZ_, _YY_, l4_val);
    De_nl.PutVal(_YZ_, _ZZ_, l4_val);
    
    // Line 5
    REAL l5_val = constant * ( epsilon.ZZ() * (m_nu - 1.0) - m_nu * (epsilon.XX() + epsilon.YY())) * dG_desp_vol;
    De_nl.PutVal(_ZZ_, _XX_, l5_val);
    De_nl.PutVal(_ZZ_, _YY_, l5_val);
    De_nl.PutVal(_ZZ_, _ZZ_, l5_val);
    
    De+=De_nl;
    return;
}

void TPZPorousElasticResponse::De(const TPZTensor<STATE> & epsilon, TPZFMatrix<STATE> & De) const {
    if (m_is_G_constant_Q) {
        De_G_constant(epsilon, De);
    }else{
        De_Poisson_constant(epsilon, De);
    }
}

TPZElasticResponse TPZPorousElasticResponse::EvaluateElasticResponse(const TPZTensor<STATE> & epsilon) const{
    
    TPZElasticResponse LinearER;
    REAL Eyoung;
    /// The properties are computed as zero order approach, i.e. constant associated to epsilon
    STATE G;
    if (m_is_G_constant_Q) {
        STATE K,dK_desp_vol,nu;
        this->K(epsilon, K, dK_desp_vol);
        Eyoung = 9*K*m_mu/(3.0*K+m_mu);
        nu     = (3.0*K-2.0*m_mu)/(2.0*(3.0*K+m_mu));
        LinearER.SetEngineeringData(Eyoung, nu);
    }else{
        STATE dG_desp_vol;
        this->G(epsilon, G, dG_desp_vol);
        Eyoung = 2*G*(1.0+m_nu);
        LinearER.SetEngineeringData(Eyoung, m_nu);
    }
    
    /// Seeking for an equivalent residual strain
    {
        TPZTensor<REAL> linear_epsilon,sigma, eps_res;
        this->ComputeStress(epsilon, sigma);
        LinearER.ComputeStrain(sigma, linear_epsilon);
        eps_res = epsilon - linear_epsilon;
        LinearER.SetResidualStrainData(eps_res);
    }
    
    return LinearER;
}

TPZElasticResponse TPZPorousElasticResponse::LinearizedElasticResponse(const TPZTensor<STATE> & epsilon_ref, const TPZTensor<STATE> & epsilon) const{
    
    TPZElasticResponse LinearER;
    REAL Eyoung;
    /// The properties are computed as zero order approach, i.e. constant associated to epsilon
    STATE G;
    if (m_is_G_constant_Q) {
        STATE nu;
        this->Poisson_linearized(epsilon_ref,epsilon, nu);
        Eyoung = 2*m_mu*(1.0+nu);
        LinearER.SetEngineeringData(Eyoung, nu);
    }else{
        STATE dG_desp_vol;
        this->G(epsilon, G, dG_desp_vol);
        Eyoung = 2*G*(1.0+m_nu);
        LinearER.SetEngineeringData(Eyoung, m_nu);
        DebugStop(); // Implement linearization for G expression.
    }
    
    /// Seeking for an equivalent residual strain
    {
        TPZTensor<REAL> linear_epsilon,sigma, eps_res;
        this->ComputeStress(epsilon, sigma);
        LinearER.ComputeStrain(sigma, linear_epsilon);
        eps_res = epsilon - linear_epsilon;
        LinearER.SetResidualStrainData(eps_res);
    }
    
    return LinearER;
}
