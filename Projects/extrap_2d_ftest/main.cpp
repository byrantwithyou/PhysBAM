//#####################################################################
// Copyright 2011 Russell Howes
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################


/*#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
*/
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <fstream>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    PARSE_ARGS parse_args;
    
    /*bool use_extended_neohookean;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_extended_neohookean_smooth;
    bool use_corotated;
    bool use_corot_blend;
    bool use_constant_ife;*/
    
    T stiffness_multiplier = 1.0;
    T damping_multiplier = 1.0;
    T stiffness= 1e5;
    T poissons_ratio = .45;
    T damping = .01;
    T energy,sv1,sv2;
    
   // if(PARSE_ARGS::Find_And_Remove("-incomp",argc,argv)) example=new INCOMPRESSIBLE_TESTS<T>(stream_type);
   // else if(PARSE_ARGS::Find_And_Remove("-hair_sim_tests",argc,argv)) example=new HAIR_SIM_TESTS<T>(stream_type);
    
    parse_args.Add_Option_Argument("-dump_sv");
    parse_args.Add_Double_Argument("-youngs_modulus",1e6,"parameter used by multiple tests to change the parameters of the test");
    parse_args.Add_Double_Argument("-poissons_ratio",.45,"","stiffness multiplier for various tests");
    parse_args.Add_Vector_2D_Argument("-sv", VECTOR<double,2>((T)1.0,(T)1.0),"","Singular Values (default 1)");
    
    parse_args.Add_Option_Argument("-use_ext_neo");
    parse_args.Add_Option_Argument("-use_ext_neo2");
    parse_args.Add_Option_Argument("-use_ext_neo_ref");
    parse_args.Add_Option_Argument("-use_ext_neo_hyper");
    parse_args.Add_Option_Argument("-use_ext_neo_smooth");
    parse_args.Add_Option_Argument("-use_corotated");
    parse_args.Add_Option_Argument("-use_corot_blend");
    //parse_args.Add_Option_Argument("-use_ext_mooney");
    parse_args.Add_Option_Argument("-use_constant_ife");   
    
    parse_args.Parse(argc,argv);
    
    stiffness=(T)parse_args.Get_Double_Value("-youngs_modulus");
    poissons_ratio=(T)parse_args.Get_Double_Value("-poissons_ratio");
    VECTOR<double,2> singular_vals = parse_args.Get_Vector_2D_Value("-sv"); sv1=singular_vals(1); sv2=singular_vals(2);
    
    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=0;
    if(parse_args.Is_Value_Set("-use_ext_neo")) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20);
    else if(parse_args.Is_Value_Set("-use_ext_neo2")) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20);
    else if(parse_args.Is_Value_Set("-use_ext_neo_smooth")) icm=new NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.1);
    else if(parse_args.Is_Value_Set("-use_corotated")) icm=new COROTATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(parse_args.Is_Value_Set("-use_corot_blend")) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);

    else{
        NEO_HOOKEAN<T,2>* nh=new NEO_HOOKEAN<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
        icm=nh; std::cout << "Using regular Neo-Hookean" << std::endl;
        nh->use_constant_ife=parse_args.Is_Value_Set("-use_constant_ife");}
    
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
	printf("dPdF/1111                   : %10.5e \n",dPdF.x1111);
	printf("dPdF/2222                   : %10.5e \n",dPdF.x2222);
	printf("dPdF/2211                   : %10.5e \n",dPdF.x2211);
    printf("\n");
	printf("dPdF/2121                   : %10.5e \n",dPdF.x2121);
	printf("dPdF/2112                   : %10.5e \n",dPdF.x2112);
    
    delete icm;
    return 0;
}
//#####################################################################
