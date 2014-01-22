//#####################################################################
// Copyright 2011 Russell Howes
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################


/*#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
*/
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Deformables/Constitutive_Models/COROTATED.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <fstream>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
//    typedef float RW;
//    typedef VECTOR<T,2> TV;
    PARSE_ARGS parse_args(argc,argv);
    
    /*bool use_extended_neohookean;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_smooth;
    bool use_corotated;
    bool use_constant_ife;*/
    
    T stiffness_multiplier = 1.0;
    T damping_multiplier = 1.0;
    T stiffness= 1e6;
    T poissons_ratio = .45;
    T damping = .01;
    T energy,sv1,sv2;
    bool use_ext_neo=false,use_ext_neo2=false,use_corotated=false,use_constant_ife=false;
    VECTOR<T,2> singular_vals(1,1);

    parse_args.Add("-youngs_modulus",&stiffness,"value","parameter used by multiple tests to change the parameters of the test");
    parse_args.Add("-poissons_ratio",&poissons_ratio,"value","stiffness multiplier for various tests");
    parse_args.Add("-sv",&singular_vals,"sv sv","Singular Values");
    parse_args.Add("-use_ext_neo",&use_ext_neo,"use_ext_neo");
    parse_args.Add("-use_ext_neo2",&use_ext_neo2,"use_ext_neo2");
    parse_args.Add("-use_corotated",&use_corotated,"use_corotated");
    parse_args.Add("-use_constant_ife",&use_constant_ife,"use_constant_ife");   
    
    parse_args.Parse();
    sv1=singular_vals(0);
    sv2=singular_vals(1);
    
    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=0;
    if(use_ext_neo) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20);
    else if(use_ext_neo2) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20);
    else if(use_corotated) icm=new COROTATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else{
        NEO_HOOKEAN<T,2>* nh=new NEO_HOOKEAN<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
        icm=nh; std::cout << "Using regular Neo-Hookean" << std::endl;
        nh->use_constant_ife=use_constant_ife;}
    
    ///////////////////////////////////////////
    DIAGONAL_MATRIX<T,2> F(sv1,sv2);
    ///////////////////////////////////////////
    
    
    
    
    DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2> dPdF;
    energy = icm->Energy_Density(F,(int)-2);
    DIAGONAL_MATRIX<T,2> P((T)1,(T)1);
    P=icm->P_From_Strain(F,1.0,(int)1);
    icm->Isotropic_Stress_Derivative(F,dPdF,(int)-2);
    
    printf("Poisson's Ratio             : %10.5e \n",poissons_ratio);
    printf("Young's Modulus             : %10.5e \n\n",stiffness);
    printf("\n");
    printf("s1                          : %10.5e \n",sv1);
    printf("s2                          : %10.5e \n",sv2);
    printf("\n");    
    printf("Lambda                      : %10.5e \n",icm->constant_lambda);
    printf("Mu                          : %10.5e \n\n",icm->constant_mu);
    printf("Energy                      : %10.5e \n",energy);
    printf("\n");
    printf("P(1)                        : %10.5e \n",P(1));
    printf("P(2)                        : %10.5e \n",P(2));
    printf("\n");
    printf("dPdF/1111                   : %10.5e \n",dPdF.x0000);
    printf("dPdF/2222                   : %10.5e \n",dPdF.x1111);
    printf("dPdF/2211                   : %10.5e \n",dPdF.x1100);
    printf("\n");
    printf("dPdF/2121                   : %10.5e \n",dPdF.x1010);
    printf("dPdF/2112                   : %10.5e \n",dPdF.x1001);
    
    delete icm;
    return 0;
}
//#####################################################################
