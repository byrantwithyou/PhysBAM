/*
 *  DEFORMABLE_OBJECT_3D.h
 *
 *  Created by Joseph Teran on 12/27/10.
 *  Copyright 2010 UCLA. All rights reserved.
 *
 */
#ifndef _deformable_object_
#define _deformable_object_

#include "GEOMETRY.h"
#include "mesh_cutting_subd.h"

using namespace GEOMETRY;
using namespace ALGEBRA;

template <class T>
class HYPERELASTICITY_CONSTITUTIVE_MODEL_3D{
public:
    virtual T Psi(const MATRIX_3X3<T>& F,const int element){return (T)0;}
    virtual MATRIX_3X3<T> P(const MATRIX_3X3<T>& F,const int element){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> P32(const MATRIX_3X2<T>& F,const int element){return MATRIX_3X2<T>();}
    virtual MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F,const MATRIX_3X3<T>& dF,const int element){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> dP32(const MATRIX_3X2<T>& F,const MATRIX_3X2<T>& dF,const int element){return MATRIX_3X2<T>();}
    virtual T Psi(const int element){return (T)0;}
    //The following versions of the functions are used when we have stored F per element
    virtual MATRIX_3X3<T> P(const int element){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> P32(const int element){return MATRIX_3X2<T>();}
    virtual MATRIX_3X3<T> dP(const MATRIX_3X3<T>& dF,const int element){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> dP32(const MATRIX_3X2<T>& dF,const int element){return MATRIX_3X2<T>();}
    virtual void Element_Stifness_Matrix(const int element,const MATRIX_3X3<T>& Dm_Inverse,MATRIX_MXN<T>& element_stiffness){}
    virtual void Update_Position_Based_State(){}
    int Index(const int node, const int component){
        return 3*node+component;
    }
    VECTOR_3D<T> Natural_Interpolating_Function_Gradient(const int node_number){
        assert(node_number>=0);
        assert(node_number<4);
        
        if(node_number==0) return VECTOR_3D<T>(-1,-1,-1);
        else if(node_number==1) return VECTOR_3D<T>(1,0,0);
        else if(node_number==2) return VECTOR_3D<T>(0,1,0);
        else if(node_number==3) return VECTOR_3D<T>(0,0,1);
        else return VECTOR_3D<T>(0,0,0);        
    }
    virtual ~HYPERELASTICITY_CONSTITUTIVE_MODEL_3D(){};
};

template <class T>
class LINEAR_ELASTICITY_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
    T mu;
    VECTOR<MATRIX_3X3<T> >& F;
    
public:
    LINEAR_ELASTICITY_3D(const T youngs_modulus,const T poisson_ratio,VECTOR<MATRIX_3X3<T> >& F_input):F(F_input){
        lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
        mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
    }
    
    MATRIX_3X3<T> P(const MATRIX_3X3<T>& F,const int element){
        MATRIX_3X3<T> epsilon=(T).5*(F+F.Transposed())-MATRIX_3X3<T>::Identity();
        return (T)2*mu*epsilon+lambda*epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
    
    MATRIX_3X3<T> P(const int element){
        MATRIX_3X3<T> epsilon=(T).5*(F(element)+F(element).Transposed())-MATRIX_3X3<T>::Identity();
        return (T)2*mu*epsilon+lambda*epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
    
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F, const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T> d_epsilon=(T).5*(dF+dF.Transposed());
        return (T)2*mu*d_epsilon+lambda*d_epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
    
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T> d_epsilon=(T).5*(dF+dF.Transposed());
        return (T)2*mu*d_epsilon+lambda*d_epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
    
    void Element_Stifness_Matrix(const int element,const MATRIX_3X3<T>& Dm_Inverse,MATRIX_MXN<T>& element_stiffness){
        MATRIX_3X3<T> M=Dm_Inverse*Dm_Inverse.Transposed();
        
        T entry_ab;
        T jacobian_determinant=(T)1/Dm_Inverse.Determinant();
        
        for(int a=0;a<4;a++){
            VECTOR_3D<T> grad_Na=HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(a);
            for(int b=0;b<4;b++){
                VECTOR_3D<T> grad_Nb=HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(b);
                VECTOR_3D<T> image=M*grad_Nb;
                entry_ab=-(mu/(T)6)*jacobian_determinant*VECTOR_3D<T>::Dot_Product(grad_Na,image);
                MATRIX_3X3<T> diag=entry_ab*MATRIX_3X3<T>::Identity();
                for(int i=0;i<3;i++){
                    for(int j=0;j<3;j++){
                        element_stiffness(HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(a,i),HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(b,j))=diag(i,j);}}}}
    
        for(int a=0;a<4;a++){
            VECTOR_3D<T> grad_Na=Dm_Inverse.Transposed()*HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(a);
            for(int i=0;i<3;i++){
                for(int b=0;b<4;b++){
                    VECTOR_3D<T> grad_Nb=Dm_Inverse.Transposed()*HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(b);
                    for(int j=0;j<3;j++){
                        element_stiffness(HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(a,i),HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(b,j))-=((mu+lambda)*jacobian_determinant/(T)6)*grad_Na(i)*grad_Nb(j);}}}}
    }
};


template <class T>
class FIXED_COROTATED_ELASTICITY_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
    T mu;
    VECTOR<MATRIX_3X3<T> > U;
    VECTOR<MATRIX_3X3<T> > V;
    VECTOR<VECTOR_3D<T> > sigma;
    VECTOR<MATRIX_3X3<T> >& F;
    
public:
    FIXED_COROTATED_ELASTICITY_3D(const T youngs_modulus,const T poisson_ratio,VECTOR<MATRIX_3X3<T> >& F_input):
    U(F_input.Size()),V(F_input.Size()),sigma(F_input.Size()),F(F_input){
        lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
        mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
        //std::cout << "lambda: " << lambda << ", mu: " << mu << std::endl;
        //exit(1);
    }

    void Update_Position_Based_State(){//assumes F is up to date
        for(int t=0;t<F.Size();t++) {
            F(t).SVD(U(t),sigma(t),V(t));
            //std::cout << sigma(t)(0) << ", " << sigma(t)(1) << ", " << sigma(t)(2) << std::endl;
        }
    }
    
    T Psi(const MATRIX_3X3<T>& F,const int element){
        MATRIX_3X3<T> U_e,V_e;VECTOR_3D<T> sigma_e;
        F.SVD(U_e,sigma_e,V_e);
        T J=sigma_e.x()*sigma_e.y()*sigma_e.z();
        T I1=sigma_e.x()*sigma_e.x()+sigma_e.y()*sigma_e.y()+sigma_e.z()*sigma_e.z();
        T psi=mu*(I1-(T)2*(sigma_e.x()+sigma_e.y()+sigma_e.z())+(T)3)+(T).5*lambda*(J-(T)1)*(J-(T)1);
        return psi;
    }
    
    T Psi(const int element){
        const VECTOR_3D<T>& sigma_e = sigma(element);
        T J=sigma(element).x()*sigma(element).y()*sigma(element).z();
        T I1=sigma_e.x()*sigma_e.x()+sigma_e.y()*sigma_e.y()+sigma_e.z()*sigma_e.z();
        T psi=mu*(I1-(T)2*(sigma_e.x()+sigma_e.y()+sigma_e.z())+(T)3)+(T).5*lambda*(J-(T)1)*(J-(T)1);
        //T psi=mu*(sigma(element).x()*sigma(element).x()+sigma(element).y()*sigma(element).y()+sigma(element).z()*sigma(element).z())+(T).5*lambda*(J-(T)1)*(J-(T)1);
        return psi;
    }
    
    MATRIX_3X3<T> P(const MATRIX_3X3<T>& F,const int element){
        MATRIX_3X3<T> U_e,V_e;VECTOR_3D<T> sigma_e;
        F.SVD(U_e,sigma_e,V_e);
        MATRIX_3X3<T> JFinvT=U_e*MATRIX_3X3<T>(sigma_e.y()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.y())*V_e.Transposed();
        T J=sigma_e.x()*sigma_e.y()*sigma_e.z();
        MATRIX_3X3<T> P_e=(T)2*mu*(F-U_e*V_e.Transposed())+lambda*(J-1)*JFinvT;
        return P_e;
    }
    
    MATRIX_3X3<T> P(const int element){
        MATRIX_3X3<T> JFinvT=U(element)*MATRIX_3X3<T>(sigma(element).y()*sigma(element).z(),(T)0,(T)0,(T)0,sigma(element).x()*sigma(element).z(),(T)0,(T)0,(T)0,sigma(element).x()*sigma(element).y())*V(element).Transposed();
        T J=sigma(element).x()*sigma(element).y()*sigma(element).z();
        MATRIX_3X3<T> P_e=(T)2*mu*(F(element)-U(element)*V(element).Transposed())+lambda*(J-1)*JFinvT;        
        return P_e;
    }
    
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F,const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T> U_e,V_e;VECTOR_3D<T> sigma_e;
        F.SVD(U_e,sigma_e,V_e);
        MATRIX_3X3<T> R=U_e*V_e.Transposed();
        MATRIX_3X3<T> S=V_e*MATRIX_3X3<T>(sigma_e.x(),(T)0,(T)0,(T)0,sigma_e.y(),(T)0,(T)0,(T)0,sigma_e.z())*V_e.Transposed();
        MATRIX_3X3<T> dR=MATRIX_3X3<T>::Delta_R(R,S,dF);
        MATRIX_3X3<T> JFinvT=U_e*MATRIX_3X3<T>(sigma_e.y()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.y())*V_e.Transposed();
        T dJ=MATRIX_3X3<T>::Contract(JFinvT,dF);
        T J=sigma_e.x()*sigma_e.y()*sigma_e.z();
        MATRIX_3X3<T> dP_e=(T)2*mu*(dF-dR)+lambda*(dJ*JFinvT+(J-1)*MATRIX_3X3<T>::d_JF_Inverse(F,dF).Transposed());        
        return dP_e;
    }

    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T>& U_e=U(element);
        MATRIX_3X3<T>& V_e=V(element);
        VECTOR_3D<T>& sigma_e=sigma(element);
        MATRIX_3X3<T> R=U_e*V_e.Transposed();
        MATRIX_3X3<T> S=V_e*MATRIX_3X3<T>(sigma_e.x(),(T)0,(T)0,(T)0,sigma_e.y(),(T)0,(T)0,(T)0,sigma_e.z())*V_e.Transposed();
        MATRIX_3X3<T> dR=MATRIX_3X3<T>::Delta_R(R,S,dF);
        MATRIX_3X3<T> JFinvT=U_e*MATRIX_3X3<T>(sigma_e.y()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.z(),(T)0,(T)0,(T)0,sigma_e.x()*sigma_e.y())*V_e.Transposed();
        T dJ=MATRIX_3X3<T>::Contract(JFinvT,dF);
        T J=sigma_e.x()*sigma_e.y()*sigma_e.z();
        MATRIX_3X3<T> dP_e=(T)2*mu*(dF-dR)+lambda*(dJ*JFinvT+(J-1)*MATRIX_3X3<T>::d_JF_Inverse(F(element),dF).Transposed());
        return dP_e;
    }

};

template <class T>
class NEO_HOOKEAN_ELASTICITY_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
    T mu;
    VECTOR<MATRIX_3X3<T> >& F;
    
public:
    NEO_HOOKEAN_ELASTICITY_3D(const T youngs_modulus,const T poisson_ratio,VECTOR<MATRIX_3X3<T> >& F_input):F(F_input){
        lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
        mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
    }
    
    T r(const T x,const int kth_derivative){
        //this is a cubic approximation to log(x)
        
        T y=0;
        
        if(kth_derivative == 0){
            y=(x-1)-.5*(x-1)*(x-1)+((T)1/(T)3)*(x-1)*(x-1)*(x-1);
            //y=log(x);
        }
        else if(kth_derivative == 1){
            y=1-(x-1)+(x-1)*(x-1);
            //y=1/x;
        }
        else if(kth_derivative == 2){
            y=2*x-3;
            //y=-1/(x*x);
        }
        else if(kth_derivative == 3){
            y=2;
        }
        else{
            y=0;}
        
        return y;
    }
    
    T Psi(const int element){
        T J=F(element).Determinant();
        T tr = (F(element).Transposed() * F(element)).Trace();
        return mu*(tr-3)/2-mu*r(J,0)+lambda*r(J,0)*r(J,0)/2;
    }
    
    T Psi(const MATRIX_3X3<T>& F,const int element){
        return 0;
    }
    
    MATRIX_3X3<T> P(const int element){
        T J=F(element).Determinant();
        MATRIX_3X3<T> J_F_invT=F(element).Cofactor_Matrix_Transposed();
        return mu*F(element)+(lambda*r(J,0)-mu)*r(J,1)*J_F_invT;
    }
    
    MATRIX_3X3<T> P(const MATRIX_3X3<T>& F,const int element){
        return MATRIX_3X3<T>();
    }
    
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F,const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T> F_inv=F;
        F_inv.Invert();
        T J=F.Determinant();
        T log_J=r(J,0);
        T log_prime_J=r(J,1);
        T log_2_prime=r(J,2);
        T c=MATRIX_3X3<T>::Trace(F_inv*dF)*((lambda*log_J-mu)*J*log_prime_J+lambda*log_prime_J*log_prime_J*J*J+log_2_prime*(log_J*lambda-mu)*J);
        return mu*dF+c*F_inv.Transposed()-(log_J*lambda-mu)*J*log_prime_J*(F_inv*dF*F_inv).Transposed();
    }
    
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& dF,const int element){
        MATRIX_3X3<T> F_inv=F(element);
        F_inv.Invert();
        T J=F(element).Determinant();
        T log_J=r(J,0);
        T log_prime_J=r(J,1);
        T log_2_prime=r(J,2);
//        T c=MATRIX_3X3<T>::Trace(F_inv*dF)*((lambda*log_J-mu)*J*log_prime_J+lambda*log_prime_J*log_prime_J*J*J+log_2_prime*(log_J*lambda-mu)*J);
//        return mu*dF+c*F_inv.Transposed()-(log_J*lambda-mu)*J*log_prime_J*(F_inv*dF*F_inv).Transposed();
        
        T c=MATRIX_3X3<T>::Trace(F_inv*dF)*((lambda*log_J-mu)*J*log_prime_J+lambda*log_prime_J*log_prime_J*J*J+J*log_2_prime*(log_J*lambda-mu)*J);
        return mu*dF+c*F_inv.Transposed()-(log_J*lambda-mu)*J*log_prime_J*(F_inv*dF*F_inv).Transposed();
    }
//    MATRIX_3X2<T> P(const MATRIX_3X2<T>& F,const int element){
//        T J=F.Determinant_ABS();
//        VECTOR_3D<T> c0=F.Column(0);VECTOR_3D<T> c1=F.Column(1);
//        T a=c0.Magnitude_Squared();T b=c0.Dot(c1);T c=c1.Magnitude_Squared();
//        MATRIX_2X2<T> C(a,b,b,c);//F^TF
//        MATRIX_2X2<T> A(0,(T)1,-(T)1,0);
//        MATRIX_3X2<T> dJdF=((T)1/J)*F*A.Transposed()*C*A;        
//        return mu*F+(lambda*r(J,0)-mu)*r(J,1)*dJdF;
//    }
//    
//    MATRIX_3X2<T> dP(const MATRIX_3X2<T>& F, const MATRIX_3X2<T>& dF,const int element){
//        T J=F.Determinant_ABS();
//        VECTOR_3D<T> c0=F.Column(0);VECTOR_3D<T> c1=F.Column(1);
//        T a=c0.Magnitude_Squared();T b=c0.Dot(c1);T c=c1.Magnitude_Squared();
//        MATRIX_2X2<T> C(a,b,b,c);//F^TF
//        MATRIX_2X2<T> A(0,(T)1,-(T)1,0);
//        MATRIX_3X2<T> dJdF=((T)1/J)*F*A.Transposed()*C*A;
//        T d_one_over_J=((T)1/(J*J))*MATRIX_3X2<T>::Contract(dJdF,dF);
//        MATRIX_2X2<T> dF_transpose_times_F(dF.Column(0).Dot(F.Column(0)),dF.Column(1).Dot(F.Column(0)),dF.Column(0).Dot(F.Column(1)),dF.Column(1).Dot(F.Column(1)));
//        MATRIX_2X2<T> dC=dF_transpose_times_F.Symmetric_Part();
//        MATRIX_3X2<T> d_dJdF=d_one_over_J*dJdF+((T)1/J)*dF*A.Transposed()*C*A+((T)1/J)*F*A.Transposed()*dC*A;
//        return mu*dF+(lambda*r(J,0)-mu)*r(J,1)*d_dJdF+(lambda*r(J,1)*r(J,1)-mu*r(J,2))*MATRIX_3X2<T>::Contract(dJdF,dF)*dJdF;
//    }
    
};
template <class T>
class LAGRANGIAN_FORCES_3D{
    VECTOR<T> forces;
    SPARSE_MATRIX<T> stiffness_matrix;//this is -df/dx(x_k)
    bool stiffness_matrix_initialized;
    VECTOR<T> delta_x;//for linearizing the forces around the current configutaion: f_lin(x+dx)=df/dx(x)dx+f(x)
    int number_nodes;
public:
    LAGRANGIAN_FORCES_3D(const int number_nodes_input):forces(3*number_nodes_input),
    stiffness_matrix(3*number_nodes_input,3*number_nodes_input),
    stiffness_matrix_initialized(false),
    delta_x(3*number_nodes_input),number_nodes(number_nodes_input){}
    virtual void Update_Position_Based_State(){}
    virtual void Compute_Position_Based_Forces(){}
    
    virtual bool Stiffness_Matrix_Initialized() const {return stiffness_matrix_initialized;}
    virtual void Initialize_Stiffness_Matrix(){stiffness_matrix_initialized=true;}
    
    VECTOR<T>& Forces(){return forces;}
    SPARSE_MATRIX<T>& Stiffness_Matrix(){return stiffness_matrix;}
    VECTOR<T>& Delta_X(){return delta_x;}
    
    virtual int Index(const int node, const int component){
        return 3*node+component;
    }
    
    void Zero_Out_Forces(){
        forces.Set_To_Zero();
    }
    
    void Add_Force(const int particle_index,const VECTOR_3D<T>& force){
        forces(3*particle_index)+=force.x();
        forces(3*particle_index+1)+=force.y();
        forces(3*particle_index+2)+=force.z();
    }
    
    virtual ~LAGRANGIAN_FORCES_3D(){}
};

template <class T>
class FEM_HYPERELASTICITY_3D: public LAGRANGIAN_FORCES_3D<T>{
public:
    TETRAHEDRON_MESH* tet_mesh;//Currently you can only have either a tet mesh or a tri mesh. TODO: Add support for both.
    TRIANGLE_MESH* tri_mesh;
    VECTOR<T>& positions;
    MATRIX_MXN<T> element_stiffness_matrix;
    VECTOR<MATRIX_3X3<T> > Dm_inverse;
    VECTOR<MATRIX_2X2<T> > Dm_inverse_S3D;
    VECTOR<MATRIX_3X3<T> > F3X3;
    VECTOR<MATRIX_3X2<T> > F3X2;    
    HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>* constitutive_model;
    bool use_gravity;
public:
    using LAGRANGIAN_FORCES_3D<T>::Add_Force;
    FEM_HYPERELASTICITY_3D(TETRAHEDRON_MESH& mesh_input,VECTOR<T>& x_input):LAGRANGIAN_FORCES_3D<T>(x_input.Size()/3),
    tet_mesh(&mesh_input),tri_mesh(0),positions(x_input),
    element_stiffness_matrix(12,12),
    Dm_inverse(tet_mesh->Number_Of_Tetrahedra()),
    Dm_inverse_S3D(0),
    F3X3(tet_mesh->Number_Of_Tetrahedra()),
    F3X2(0),
    constitutive_model(0){}
    
    FEM_HYPERELASTICITY_3D(TRIANGLE_MESH& mesh_input,VECTOR<T>& x_input):LAGRANGIAN_FORCES_3D<T>(x_input.Size()/3),
    tet_mesh(0),tri_mesh(&mesh_input),positions(x_input),
    element_stiffness_matrix(9,9),
    Dm_inverse(0),
    Dm_inverse_S3D(tri_mesh->Number_Of_Triangles()),
    F3X3(0),
    F3X2(tri_mesh->Number_Of_Triangles()),
    constitutive_model(0){}
    
    VECTOR<MATRIX_3X3<T> >& F(){return F3X3;}
    
    void Nodal_Volume_Fractions(VECTOR<T>& fractions){
        fractions.Set_To_Zero();
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
                T jacobian_determinant=(T)1/Dm_inverse(t).Determinant();
                T volume=jacobian_determinant/(T)6;
                VECTOR_4D<int>& indices=tet_mesh->Nodes_Of_Element(t);
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),0))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),1))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),2))+=volume/(T)4;
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),0))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),1))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),2))+=volume/(T)4;
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),0))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),1))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),2))+=volume/(T)4;
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(3),0))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(3),1))+=volume/(T)4;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(3),2))+=volume/(T)4;}}
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                T jacobian_determinant=(T)1/Dm_inverse_S3D(t).Determinant();
                T area=jacobian_determinant/(T)2;
                VECTOR_3D<int>& indices=tri_mesh->Nodes_Of_Element(t);
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),0))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),1))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(0),2))+=area/(T)3;
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),0))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),1))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(1),2))+=area/(T)3;
                fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),0))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),1))+=area/(T)3;fractions(LAGRANGIAN_FORCES_3D<T>::Index(indices(2),2))+=area/(T)3;}}
    }
    
    void Initialize_Undeformed_Configuration(){
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
                VECTOR_4D<int>& indices=tet_mesh->Nodes_Of_Element(t);
                VECTOR_3D<T> x0=X(indices(0));VECTOR_3D<T> x1=X(indices(1));
                VECTOR_3D<T> x2=X(indices(2));VECTOR_3D<T> x3=X(indices(3));
                //for amadillo
                //Dm_inverse(t)=MATRIX_3X3<T>((x1-x0)*0.7,(x2-x0)*0.7,(x3-x0)*0.7);
                
                //for peeling
                Dm_inverse(t)=MATRIX_3X3<T>(x1-x0,x2-x0,x3-x0);
                
                //for letters
//                VECTOR_3D<T> x10=x1-x0;
//                VECTOR_3D<T> x20=x2-x0;
//                VECTOR_3D<T> x30=x3-x0;
//                x30(0)*=0.9; x30(1)*=0.9;
//                x20(0)*=0.9; x20(1)*=0.9;
//                x10(0)*=0.9; x10(1)*=0.9;
//                Dm_inverse(t)=MATRIX_3X3<T>(x10,x20,x30);
                
                
                Dm_inverse(t).Invert();}}
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                VECTOR_3D<int>& indices=tri_mesh->Nodes_Of_Element(t);
                VECTOR_3D<T> x0=X(indices(0));VECTOR_3D<T> x1=X(indices(1));
                VECTOR_3D<T> x2=X(indices(2));
                VECTOR_3D<T> v1=x1-x0;VECTOR_3D<T> v2=x2-x0;
                T component=v1.Magnitude();
                VECTOR_3D<T> n=VECTOR_3D<T>::Cross_Product(v1,v2);
                VECTOR_3D<T> e1=((T)1/component)*v1;
                VECTOR_3D<T> e2=VECTOR_3D<T>::Cross_Product(n,v1);
                e2.Normalize();
                Dm_inverse_S3D(t)=MATRIX_2X2<T>(component,(T)0,VECTOR_3D<T>::Dot_Product(e1,v2),VECTOR_3D<T>::Dot_Product(e2,v2));
                Dm_inverse_S3D(t).Invert();}}
    }
    
    void Initialize_Stiffness_Matrix(){
        assert(constitutive_model);
        LAGRANGIAN_FORCES_3D<T>::Initialize_Stiffness_Matrix();
        
        SPARSE_MATRIX<T>& stiffness_matrix=LAGRANGIAN_FORCES_3D<T>::Stiffness_Matrix();
        
        element_stiffness_matrix.Set_To_Zero();
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
            
                //Add entries from element stiffness matrix into global matrix
                VECTOR_4D<int>& nodes=tet_mesh->Nodes_Of_Element(t);
                for(int ie=0;ie<4;ie++){
                    int i=nodes(ie);
                    for(int ci=0;ci<3;ci++){
                        SPARSE_ROW<T>& row=stiffness_matrix.Row(LAGRANGIAN_FORCES_3D<T>::Index(i,ci));
                        for(int je=0;je<4;je++){
                            int j=nodes(je);
                            for(int cj=0;cj<3;cj++){
                                if(!row.Value_Exists_At_Entry(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))){
                                    row.Add_Entry(LAGRANGIAN_FORCES_3D<T>::Index(j,cj),element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj)));}
                                else{
                                    row(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))+=element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj));}}}}}}}
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                
                //Add entries from element stiffness matrix into global matrix
                VECTOR_3D<int>& nodes=tri_mesh->Nodes_Of_Element(t);
                for(int ie=0;ie<3;ie++){
                    int i=nodes(ie);
                    for(int ci=0;ci<3;ci++){
                        SPARSE_ROW<T>& row=stiffness_matrix.Row(LAGRANGIAN_FORCES_3D<T>::Index(i,ci));
                        for(int je=0;je<3;je++){
                            int j=nodes(je);
                            for(int cj=0;cj<3;cj++){
                                if(!row.Value_Exists_At_Entry(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))){
                                    row.Add_Entry(LAGRANGIAN_FORCES_3D<T>::Index(j,cj),element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj)));}
                                else{
                                    row(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))+=element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj));}}}}}}}
        
        stiffness_matrix.Zero_Out_Without_Changing_Sparsity();//We just set the sparsity *pattern*, will fill in actual values later.
    }
    
    void Set_Constitutive_Model(HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>& cons_input){
        constitutive_model=&cons_input;
    }
    
    MATRIX_3X3<T> Element_Nodewise_dF(const int tet,const int index,const VECTOR_3D<T>& dx){
        //This is for constructing the element stiffness matrix
        if(index==1){
            MATRIX_3X3<T> Ds(dx,VECTOR_3D<T>(0),VECTOR_3D<T>(0));
            return Ds*Dm_inverse(tet);}
        else if(index==2){
            MATRIX_3X3<T> Ds(VECTOR_3D<T>(0),dx,VECTOR_3D<T>(0));
            return Ds*Dm_inverse(tet);}
        else if(index==3){
            MATRIX_3X3<T> Ds(VECTOR_3D<T>(0),VECTOR_3D<T>(0),dx);
            return Ds*Dm_inverse(tet);}
        else if(index==0){
            MATRIX_3X3<T> Ds(-dx,-dx,-dx);
            return Ds*Dm_inverse(tet);}
        else{
            assert(false);
            return MATRIX_3X3<T>(0);}
    }
    
    MATRIX_3X2<T> Element_Nodewise_dF_S3D(const int tri,const int index,const VECTOR_3D<T>& dx){
        //This is for constructing the element stiffness matrix
        if(index==1){
            MATRIX_3X2<T> Ds(dx,VECTOR_3D<T>(0));
            return Ds*Dm_inverse_S3D(tri);}
        else if(index==2){
            MATRIX_3X2<T> Ds(VECTOR_3D<T>(0),dx);
            return Ds*Dm_inverse_S3D(tri);}
        else if(index==0){
            MATRIX_3X2<T> Ds(-dx,-dx);
            return Ds*Dm_inverse_S3D(tri);}
        else{
            assert(false);
            return MATRIX_3X2<T>(0);}
    }
    
    void Update_Stiffness_Matrix(){
        assert(constitutive_model);
        
        if(!LAGRANGIAN_FORCES_3D<T>::Stiffness_Matrix_Initialized()) Initialize_Stiffness_Matrix();
        
        SPARSE_MATRIX<T>& stiffness_matrix=LAGRANGIAN_FORCES_3D<T>::Stiffness_Matrix();
        
        stiffness_matrix.Zero_Out_Without_Changing_Sparsity();
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
                element_stiffness_matrix.Set_To_Zero();
                //MATRIX_3X3<T> F=F3X3(t);
                T jacobian_determinant=(T)1/Dm_inverse(t).Determinant();
                MATRIX_3X3<T> dP,G,dF;
                VECTOR_3D<T> df0,df1,df2,df3,dx;
                for(int node=0;node<4;node++){
                    for(int c=0;c<3;c++){
                        //cth column associated with node
                        dx=VECTOR_3D<T>::Standard_Basis_Vector(c);
                        dF=Element_Nodewise_dF(t,node,dx);
                        dP=constitutive_model->dP(dF,t);
                        G=-((T)1/(T)6)*jacobian_determinant*dP*(Dm_inverse(t).Transposed());
                        //G=-jacobian_determinant*dP*(Dm_inverse(t).Transposed());
                        df1=G.Column(0);df2=G.Column(1);df3=G.Column(2);
                        df0=-(df1+df2+df3);
                        Enter_Column_Of_Tet_Element_Stiffness_Matrix(3*node+c,df0,df1,df2,df3);}}
            
                //Add entries from element stiffness matrix into global matrix
                VECTOR_4D<int>& nodes=tet_mesh->Nodes_Of_Element(t);
                for(int ie=0;ie<4;ie++){
                    int i=nodes(ie);
                    for(int ci=0;ci<3;ci++){
                        SPARSE_ROW<T>& row=stiffness_matrix.Row(LAGRANGIAN_FORCES_3D<T>::Index(i,ci));
                        for(int je=0;je<4;je++){
                            int j=nodes(je);
                            for(int cj=0;cj<3;cj++){
                                row(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))+=element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj));}}}}}}
        
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                element_stiffness_matrix.Set_To_Zero();
                //MATRIX_3X2<T> F=F3X2(t);
                T jacobian_determinant=(T)1/Dm_inverse_S3D(t).Determinant();
                MATRIX_3X2<T> dP,G,dF;
                VECTOR_3D<T> df0,df1,df2,dx;
                for(int node=0;node<3;node++){
                    for(int c=0;c<3;c++){
                        //cth column associated with node
                        dx=VECTOR_3D<T>::Standard_Basis_Vector(c);
                        dF=Element_Nodewise_dF_S3D(t,node,dx);
                        dP=constitutive_model->dP32(dF,t);
                        G=-((T).5)*jacobian_determinant*dP*(Dm_inverse_S3D(t).Transposed());
                        df1=G.Column(0);df2=G.Column(1);
                        df0=-(df1+df2);
                        Enter_Column_Of_Tri_Element_Stiffness_Matrix(3*node+c,df0,df1,df2);}}
                
                //Add entries from element stiffness matrix into global matrix
                VECTOR_3D<int>& nodes=tri_mesh->Nodes_Of_Element(t);
                for(int ie=0;ie<3;ie++){
                    int i=nodes(ie);
                    for(int ci=0;ci<3;ci++){
                        SPARSE_ROW<T>& row=stiffness_matrix.Row(LAGRANGIAN_FORCES_3D<T>::Index(i,ci));
                        for(int je=0;je<3;je++){
                            int j=nodes(je);
                            for(int cj=0;cj<3;cj++){
                                row(LAGRANGIAN_FORCES_3D<T>::Index(j,cj))+=element_stiffness_matrix(LAGRANGIAN_FORCES_3D<T>::Index(ie,ci),LAGRANGIAN_FORCES_3D<T>::Index(je,cj));}}}}}}
    }
    
    void Enter_Column_Of_Tet_Element_Stiffness_Matrix(const int index,const VECTOR_3D<T>& v0,const VECTOR_3D<T>& v1,const VECTOR_3D<T>& v2,const VECTOR_3D<T>& v3){
        assert(tet_mesh);
        assert(element_stiffness_matrix.M()==element_stiffness_matrix.N());
        assert(element_stiffness_matrix.M()==12);
        element_stiffness_matrix(0,index)=v0.x();
        element_stiffness_matrix(1,index)=v0.y();
        element_stiffness_matrix(2,index)=v0.z();
        element_stiffness_matrix(3,index)=v1.x();
        element_stiffness_matrix(4,index)=v1.y();
        element_stiffness_matrix(5,index)=v1.z();
        element_stiffness_matrix(6,index)=v2.x();
        element_stiffness_matrix(7,index)=v2.y();
        element_stiffness_matrix(8,index)=v2.z();
        element_stiffness_matrix(9,index)=v3.x();
        element_stiffness_matrix(10,index)=v3.y();
        element_stiffness_matrix(11,index)=v3.z();
    }
    
    void Enter_Column_Of_Tri_Element_Stiffness_Matrix(const int index,const VECTOR_3D<T>& v0,const VECTOR_3D<T>& v1,const VECTOR_3D<T>& v2){
        assert(tri_mesh);
        assert(element_stiffness_matrix.M()==9 && element_stiffness_matrix.N()==9);
        element_stiffness_matrix(0,index)=v0.x();
        element_stiffness_matrix(1,index)=v0.y();
        element_stiffness_matrix(2,index)=v0.z();
        element_stiffness_matrix(3,index)=v1.x();
        element_stiffness_matrix(4,index)=v1.y();
        element_stiffness_matrix(5,index)=v1.z();
        element_stiffness_matrix(6,index)=v2.x();
        element_stiffness_matrix(7,index)=v2.y();
        element_stiffness_matrix(8,index)=v2.z();
    }
    
    void Update_Position_Based_State(){
        //save the deformation gradients
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
                element_stiffness_matrix.Set_To_Zero();
                F3X3(t)=Deformation_Gradient(t);}}
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                element_stiffness_matrix.Set_To_Zero();
                F3X2(t)=Deformation_Gradient_S3D(t);}}
        
        assert(constitutive_model);
        constitutive_model->Update_Position_Based_State();
        Update_Stiffness_Matrix();
        Compute_Position_Based_Forces();
    }
    
    void Compute_Position_Based_Forces(){
        assert(constitutive_model);
        LAGRANGIAN_FORCES_3D<T>::Zero_Out_Forces();
        if(tet_mesh){
            for(int t=0;t<tet_mesh->Number_Of_Tetrahedra();t++){
                MATRIX_3X3<T> P=constitutive_model->P(t);
                T jacobian_determinant=(T)1/Dm_inverse(t).Determinant();
                MATRIX_3X3<T> G=-((T)1/(T)6)*jacobian_determinant*P*(Dm_inverse(t).Transposed());
                //MATRIX_3X3<T> G=-jacobian_determinant*P*(Dm_inverse(t).Transposed());
                VECTOR_4D<int>& nodes=tet_mesh->Nodes_Of_Element(t);
                VECTOR_3D<T> force1=G.Column(0);VECTOR_3D<T> force2=G.Column(1);VECTOR_3D<T> force3=G.Column(2);
                Add_Force(nodes(1),force1);Add_Force(nodes(2),force2);Add_Force(nodes(3),force3);
                Add_Force(nodes(0),-(force1+force2+force3));}}
        
        if(tri_mesh){
            for(int t=0;t<tri_mesh->Number_Of_Triangles();t++){
                MATRIX_3X2<T> P=constitutive_model->P32(t);
                T jacobian_determinant=(T)1/Dm_inverse_S3D(t).Determinant();
                MATRIX_3X2<T> G=-((T).5)*jacobian_determinant*P*(Dm_inverse_S3D(t).Transposed());
                VECTOR_3D<int>& nodes=tri_mesh->Nodes_Of_Element(t);
                VECTOR_3D<T> force1=G.Column(0);VECTOR_3D<T> force2=G.Column(1);
                Add_Force(nodes(1),force1);Add_Force(nodes(2),force2);
                Add_Force(nodes(0),-(force1+force2));}}
    }
    
    MATRIX_3X3<T> Deformation_Gradient(const int t){
        VECTOR_4D<int>& indices=tet_mesh->Nodes_Of_Element(t);
        VECTOR_3D<T> x0=X(indices(0));VECTOR_3D<T> x1=X(indices(1));
        VECTOR_3D<T> x2=X(indices(2));VECTOR_3D<T> x3=X(indices(3));
        MATRIX_3X3<T> Ds(x1-x0,x2-x0,x3-x0);
        return Ds*Dm_inverse(t);
    }
    
    MATRIX_3X2<T> Deformation_Gradient_S3D(const int t){
        VECTOR_3D<int>& indices=tri_mesh->Nodes_Of_Element(t);
        VECTOR_3D<T> x0=X(indices(0));VECTOR_3D<T> x1=X(indices(1));
        VECTOR_3D<T> x2=X(indices(2));
        MATRIX_3X2<T> Ds(x1-x0,x2-x0);
        return Ds*Dm_inverse_S3D(t);
    }
    
    VECTOR_3D<T> X(const int particle_index){return VECTOR_3D<T>(positions(3*particle_index),positions(3*particle_index+1),positions(3*particle_index+2));}
    
};

template <class T>
class DEFORMABLE_OBJECT_3D{
    //State variables/geom
    //Currently can only have either a tet mesh or a tri mesh. TODO: add support for both.
    TETRAHEDRON_MESH tet_mesh;//mesh type = 0 (see constructor)
    TRIANGLE_MESH tri_mesh;//mesh type = 1 (see constructor)
    VECTOR<T> positions;
    VECTOR<T> velocities;
    int number_nodes;
    
public:
    DEFORMABLE_OBJECT_3D(const int num_elements_input,const int number_nodes_input,const int mesh_type_input):
    tet_mesh(mesh_type_input==0?num_elements_input:0,number_nodes_input),tri_mesh(mesh_type_input==1?num_elements_input:0,number_nodes_input),positions(3*number_nodes_input),velocities(3*number_nodes_input),number_nodes(number_nodes_input)
    {}
    
    DEFORMABLE_OBJECT_3D(GRID_3D<T>& grid):
    tet_mesh(grid.M(),grid.N(),grid.MN()),tri_mesh(0,grid.M()*grid.N()*grid.MN()),positions(3*grid.M()*grid.N()*grid.MN()),velocities(3*grid.M()*grid.N()*grid.MN()),number_nodes(grid.M()*grid.N()*grid.MN())
    {}
    
    VECTOR<T>& Positions(){return positions;}
    VECTOR<T>& Velocities(){return velocities;}
    
    TETRAHEDRON_MESH& Tetrahedron_Mesh(){return tet_mesh;}
    TRIANGLE_MESH& Triangle_Mesh(){return tri_mesh;}
    
    void Perturb_Positions(VECTOR<T>& x){
        for (int i = 0; i < 3*number_nodes; i++){
            positions(i)+=x(i);
        }
    }
    
    void Perturb_Positions_Negative(VECTOR<T>& x){
        for (int i = 0; i < 3*number_nodes; i++){
            positions(i)-=x(i);
        }
    }
    
    
    void Set_Position(const int particle_index,const VECTOR_3D<T>& x){
        positions(3*particle_index)=x.x();
        positions(3*particle_index+1)=x.y();
        positions(3*particle_index+2)=x.z();
    }
    
    void Set_Velocity(const int particle_index,const VECTOR_3D<T>& v){
        velocities(3*particle_index)=v.x();
        velocities(3*particle_index+1)=v.y();
        velocities(3*particle_index+2)=v.z();
    }
    
    VECTOR_3D<T> X(const int particle_index){return VECTOR_3D<T>(positions(3*particle_index),positions(3*particle_index+1),positions(3*particle_index+2));}
    VECTOR_3D<T> V(const int particle_index){return VECTOR_3D<T>(velocities(3*particle_index),velocities(3*particle_index+1),velocities(3*particle_index+2));}
    
};

template <class T>
class BACKWARD_EULER_TIME_STEPPING_3D{
public:
    T dt,final_time,start_time,time;
    DEFORMABLE_OBJECT_3D<T>& deformable_object;
    LAGRANGIAN_FORCES_3D<T>* elastic_forces;
    CONJUGATE_GRADIENT<T>* cg;
    MINRES<T>* minres;
    VECTOR<int>* bc_dofs;
    SPARSE_MATRIX<T> be_matrix;
    VECTOR<T> mass;
    VECTOR<T> be_rhs;
    VECTOR<T> old_positions;
    
    T rho;
    bool matrix_intialized;
    bool use_gravity;
    int linear_solver_type;//0=cg, 1=minres
    T alpha_damping;
    bool damping;

    BACKWARD_EULER_TIME_STEPPING_3D(const T dt_input,const T final_time_input, const T start_time_input,DEFORMABLE_OBJECT_3D<T>& object_input):
    dt(dt_input),final_time(final_time_input),start_time(start_time_input),deformable_object(object_input),elastic_forces(0),cg(0),minres(0),bc_dofs(0),
    be_matrix(deformable_object.Positions().Size(),deformable_object.Positions().Size()),
    mass(deformable_object.Positions().Size()),be_rhs(deformable_object.Positions().Size()),old_positions(deformable_object.Positions().Size()),
    rho((T)1000),matrix_intialized(false),use_gravity(true),linear_solver_type(-1), alpha_damping(0.1), damping(false){
        time=start_time;
    }
    
    virtual ~BACKWARD_EULER_TIME_STEPPING_3D(){if(cg) delete cg;if(minres) delete minres;}
    
    SPARSE_MATRIX<T>& BE_Matrix(){
        return be_matrix;
    }
    
    void Initialize_BE_Matrix(const VECTOR<T>& nodal_volume){
        //std::cout<< "Initialize_BE_Matrix" << std::endl;
        assert(elastic_forces);
        
        if(!elastic_forces->Stiffness_Matrix_Initialized()) elastic_forces->Initialize_Stiffness_Matrix();
        
        SPARSE_MATRIX<T>& elasticity_stiffness_matrix=elastic_forces->Stiffness_Matrix();
        //be_matrix.Write_DAT_File("be_matrix0.txt");
        
        for(int i=0;i<elasticity_stiffness_matrix.M();i++){
            SPARSE_ROW<T>& k_row=elasticity_stiffness_matrix.Row(i);
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            for(int j=0;j<k_row.Number_Nonzero();j++){
                T kij=k_row.Value_At_Sparse_Index(j);
                //std::cout << kij << " ";
                be_row.Add_Entry(k_row.Index(j),-dt*dt*kij);
                //be_row(k_row.Index(j)) -= dt*alpha_damping*kij;
            }
        }
        //be_matrix.Write_DAT_File("be_matrix1.txt");
        
        for(int i=0;i<nodal_volume.Size();i++) mass(i)=rho*nodal_volume(i);
        
        for(int i=0;i<be_matrix.M();i++){
            SPARSE_ROW<T>& row=be_matrix.Row(i);
            if(row.Number_Nonzero()) row(i)+=mass(i); }
//        be_matrix.Write_DAT_File("be_matrix2.txt");
//        exit(1);
        matrix_intialized=true;
    }
    
    void Initialize_CG(){
        cg=new CONJUGATE_GRADIENT<T>(be_matrix,elastic_forces->Delta_X(),be_rhs,200);
        linear_solver_type=0;
    }
    
    void Initialize_MINRES(){
        minres=new MINRES<T>(be_matrix,elastic_forces->Delta_X(),be_rhs,1000);
        linear_solver_type=1;
    }
    
    void Set_Elastic_Forces(LAGRANGIAN_FORCES_3D<T>& forces_input){
        elastic_forces=&forces_input;
    }
    
    void Update_BE_RHS_And_System_Matrix(const PhysBAM::MESH_CUTTING<double>& mcut, const T K){
        //std::cout << "Update_BE_RHS_And_System_Matrix" << std::endl;
        elastic_forces->Update_Position_Based_State();//this updates the elasticity stiffness matrix and computes elastic forces at the current configutation
//        for (int dragging_id = 0; dragging_id < mcut.dragging_elements.m; dragging_id++){
//            int node_id;
//            //std::cout << mcut.dragging_weights(dragging_id) << std::endl;
//            T w[4] = {mcut.dragging_weights(dragging_id)(0), mcut.dragging_weights(dragging_id)(1), mcut.dragging_weights(dragging_id)(2), 1-mcut.dragging_weights(dragging_id)(0)-mcut.dragging_weights(dragging_id)(1)-mcut.dragging_weights(dragging_id)(2)};
//            PhysBAM::VECTOR<T,3> p;
//            for (int j = 0; j < 4; j++){
//                node_id = mcut.sim_volume->mesh.elements(mcut.dragging_elements(dragging_id))(j);
//                for (int k = 0; k < 3; k++){
//                    p(k) += (deformable_object.Positions()(node_id*3+k)*w[j]);
//                }
//            }
//            PhysBAM::VECTOR<T, 3> sh = p - mcut.dragging_targets(dragging_id);
//            //std::cout << sh << std::endl;
//            for (int j = 0; j < 4; j++){
//                node_id = mcut.sim_volume->mesh.elements(mcut.dragging_elements(dragging_id))(j);
//                for(int k = 0; k < 3; k++){
////                        std::cout <<be_rhs(3*node_id+k)<<std::endl;
//                        //std::cout <<w[j]*K*sh(k)<<std::endl;
//                    elastic_forces->Forces()(3*node_id+k) -= (w[j]*K*sh(k));
////                        std::cout <<be_rhs(3*node_id+k)<<std::endl;
//                }  
//            }
//            
//            
//            for (int j = 0; j < 3; j++){
//                for (int k = 0; k < 4; k++){
//                    int jj = 3*mcut.sim_volume->mesh.elements(mcut.dragging_elements(dragging_id))(k)+j;
//                    for (int l = 0; l < 4; l++){
//                        int kk = mcut.sim_volume->mesh.elements(mcut.dragging_elements(dragging_id))(l)*3+j;
//                        //std::cout << (w[k]*w[l]*K) << std::endl;
//                        elastic_forces->Stiffness_Matrix()(jj, kk) -= (w[k]*w[l]*K);}}}
//        }
        
        
        VECTOR<T>& f_xn=elastic_forces->Forces();
        VECTOR<T>& vn=deformable_object.Velocities();
        
        SPARSE_MATRIX<T>& elasticity_stiffness_matrix=elastic_forces->Stiffness_Matrix();
        
        //first set the RHS
        for(int i=0;i<be_rhs.Size();i++) be_rhs(i)=dt*dt*f_xn(i)+dt*mass(i)*vn(i);
        
        //add damping terms
        if (damping) {
            for(int i=0;i<be_rhs.Size();i++) be_rhs(i)=dt*dt*f_xn(i)+dt*mass(i)*vn(i)-mass(i)*(deformable_object.Positions()(i)-old_positions(i));//with damping
            VECTOR<T> rrr(be_rhs.Size());
            for(int i = 0; i < be_rhs.Size(); i++){
                rrr(i) = deformable_object.Positions()(i)-old_positions(i);
            }
            VECTOR<T> rr(be_rhs.Size());
            elasticity_stiffness_matrix.Multiply(rrr, rr);
            rr.Scale(alpha_damping*dt);
            be_rhs += rr;
        }
        
        if(use_gravity) for(int i=0;i<be_rhs.Size()/3;i++) be_rhs(3*i+1)-=dt*dt*mass(3*i+1)*19.8;

        
        //now update the BE matrix
        assert(matrix_intialized);
        be_matrix.Zero_Out_Without_Changing_Sparsity();
        //be_matrix.Write_DAT_File("be_matrix00.txt");
        
        
        for(int i=0;i<elasticity_stiffness_matrix.M();i++){
            SPARSE_ROW<T>& k_row=elasticity_stiffness_matrix.Row(i);
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            for(int j=0;j<k_row.Number_Nonzero();j++){
                T kij=k_row.Value_At_Sparse_Index(j);
                be_row(k_row.Index(j))=-dt*dt*kij;
                if (damping) {
                    be_row(k_row.Index(j)) -= dt*alpha_damping*kij; //damping force
                }
            }
        }
        //be_matrix.Write_DAT_File("be_matrix111.txt", 1);
        
        for(int i=0;i<be_matrix.M();i++){
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            if(be_row.Number_Nonzero()) be_row(i)+=mass(i);}//std::cout << be_row(i) << std::endl;}
        //be_matrix.Write_DAT_File("be_matrix222.txt");
        //exit(1);
    }
    
    void Advance_One_Time_Step(const PhysBAM::MESH_CUTTING<double>& mcut, const T K, const int max_it=3, const bool check_convergence=false){
        //assuming the user set the boundary conditions already
        assert(elastic_forces);
        
        for(int p=0;p<old_positions.Size();p++) old_positions(p)=deformable_object.Positions()(p);
            
        for(int ns=0;ns<max_it;ns++){
            T step_size=(T)(ns+1)/(T)(max_it);
            //T step_size = 1;
            
            Update_BE_RHS_And_System_Matrix(mcut, K);
            //elastic_forces->Stiffness_Matrix().Write_DAT_File("matrix.txt");
            //exit(1);
            //the unknowns are the positions
            Solve_Linearized_System();
            for(int p=0;p<deformable_object.Positions().Size();p++)
                deformable_object.Positions()(p)+=step_size*elastic_forces->Delta_X()(p);
                
            if(check_convergence){
                //convergence check
                // norm of M*(xnp1-xn) - dt*M*vn -dt*dt*f(xnp1)
                VECTOR<T>& f_xn=elastic_forces->Forces();
                VECTOR<T>& vn=deformable_object.Velocities();
                VECTOR<T> residual(vn.Size());
                for(int i=0;i<be_rhs.Size();i++) residual(i)=mass(i)*(deformable_object.Positions()(i)-old_positions(i)-dt*vn(i))-dt*dt*f_xn(i);

                if(use_gravity) for(int i=0;i<be_rhs.Size()/3;i++) residual(3*i+1)+=dt*dt*mass(3*i+1)*(T)19.8;
            
                if(bc_dofs){for(int i=0;i<bc_dofs->Size();i++) residual((*bc_dofs)(i))=(T)0;}
            
                std::cout <<"Residual of Newton solve = " <<residual.Magnitude() << std::endl;
                
                if (ns == max_it-1) std::cout<<std::endl;    
            }
        }
        
        //update velocities
        for(int p=0;p<deformable_object.Positions().Size();p++)
            deformable_object.Velocities()(p)=(deformable_object.Positions()(p)-old_positions(p))/dt;
        
        time+=dt;
    }
    
//    void Update_BE_RHS_And_System_Matrix(const PhysBAM::MESH_CUTTING<T>& mcut, const T K){
//        elastic_forces->Update_Position_Based_State();//this updates the elasticity stiffness matrix and computes elastic forces at the current configutation
//        VECTOR<T>& f_xn=elastic_forces->Forces();
//        VECTOR<T>& vn=deformable_object.Velocities();
//        
//        //first set the RHS
//        for(int i=0;i<be_rhs.Size();i++) be_rhs(i)=dt*dt*f_xn(i)+dt*mass(i)*vn(i);
//        if(use_gravity) for(int i=0;i<be_rhs.Size()/3;i++) be_rhs(3*i+1)-=dt*dt*mass(3*i+1)*9.8;
//        
//        //now update the BE matrix
//        assert(matrix_intialized);
//        be_matrix.Zero_Out_Without_Changing_Sparsity();
//        
//        SPARSE_MATRIX<T>& elasticity_stiffness_matrix=elastic_forces->Stiffness_Matrix();
//        
//        for(int i=0;i<elasticity_stiffness_matrix.M();i++){
//            SPARSE_ROW<T>& k_row=elasticity_stiffness_matrix.Row(i);
//            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
//            for(int j=0;j<k_row.Number_Nonzero();j++){
//                T kij=k_row.Value_At_Sparse_Index(j);
//                be_row(k_row.Index(j))-=dt*dt*kij;}}
//        
//        for(int i=0;i<be_matrix.M();i++){
//            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
//            be_row(i)+=mass(i);}
//    }
//    
//    void Advance_One_Time_Step(const PhysBAM::MESH_CUTTING<T>& mcut, const T K, const int max_it=3, const bool check_convergence=false){
//        //assuming the user set the boundary conditions already
//        assert(elastic_forces);
//        
//        Update_BE_RHS_And_System_Matrix(mcut, K);
//        //the unknowns are the positions
//        Solve_Linearized_System();//this just does one step of Newton
//        for(int p=0;p<deformable_object.Positions().Size();p++){
//            deformable_object.Positions()(p)+=elastic_forces->Delta_X()(p);
//            deformable_object.Velocities()(p)=elastic_forces->Delta_X()(p)/dt;}
//        
//        time+=dt;
//    }
    
    T Time(){return time;}
    
    void Solve_Linearized_System(){
        
        if(linear_solver_type==0 && cg){
            //std::cout << "cg\n";
            cg->Solve();return;}
        else if(linear_solver_type==1 && minres){
            //std::cout << "minres\n";
            minres->Solve();return;}
        std::cout << "solver is not initialized!!!\n";
        exit(1);
    }
    
    void Set_Boundary_Conditions(VECTOR<int>& constrained_nodes,VECTOR<T>& constrained_node_locations){
        bc_dofs=&constrained_nodes;
        
        if(linear_solver_type==0 && cg){cg->Set_Dirichlet_Dofs(constrained_nodes);}
        else if(linear_solver_type==1 && minres){minres->Set_Dirichlet_Dofs(constrained_nodes);}
        
        for(int i=0;i<constrained_nodes.Size();i++){
            deformable_object.Positions()(constrained_nodes(i))=constrained_node_locations(i);
            elastic_forces->Delta_X()(constrained_nodes(i))=(T)0;}
    }
};

template <class T>
class QUASISTATIC_TIME_STEPPING_3D{
    T dt,final_time,start_time,time;
    DEFORMABLE_OBJECT_3D<T>& deformable_object;
    LAGRANGIAN_FORCES_3D<T>* elastic_forces;
    CONJUGATE_GRADIENT<T>* cg;
public:
    QUASISTATIC_TIME_STEPPING_3D(const T dt_input,const T final_time_input, const T start_time_input,DEFORMABLE_OBJECT_3D<T>& object_input):
    dt(dt_input),final_time(final_time_input),start_time(start_time_input),deformable_object(object_input),elastic_forces(0),cg(0){
        time=start_time;
    }
    
    virtual ~QUASISTATIC_TIME_STEPPING_3D(){if(cg) delete cg;}
    
    void Initialize_CG(SPARSE_MATRIX<T>& system_matrix,VECTOR<T>&x,VECTOR<T>&b){
        cg=new CONJUGATE_GRADIENT<T>(system_matrix,x,b,10*b.Size());
    }
    
    void Set_Elastic_Forces(LAGRANGIAN_FORCES_3D<T>& forces_input){
        elastic_forces=&forces_input;
    }
    
    void Advance_One_Time_Step(){
        //assuming the user set the boundary conditions already
        assert(elastic_forces);
        elastic_forces->Update_Position_Based_State();//this updates the matrix and computes forces at the current configutation
        VECTOR<T>& f_xn=elastic_forces->Forces();for(int i=0;i<f_xn.Size();i++) f_xn(i)=-f_xn(i);//need to negate these for the RHS in the Netwon solve: df/dx(x_old) dx = -f(x_old)
        Solve_For_Equilibrium_Configuration();//this just does one step of Newton
        for(int p=0;p<deformable_object.Positions().Size();p++) deformable_object.Positions()(p)+=elastic_forces->Delta_X()(p);
        time+=dt;
    }
    
    T Time(){return time;}
    
    void Solve_For_Equilibrium_Configuration(){if(cg) cg->Solve();}
    
    void Set_Boundary_Conditions(VECTOR<int>& constrained_nodes,VECTOR<T>& constrained_node_locations){
        if(cg) cg->Set_Dirichlet_Dofs(constrained_nodes);
        for(int i=0;i<constrained_nodes.Size();i++){
            deformable_object.Positions()(constrained_nodes(i))=constrained_node_locations(i);
            elastic_forces->Delta_X()(constrained_nodes(i))=(T)0;}
    }
};

#endif

