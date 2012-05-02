#include "PARAMETERS_COMMON.h"

template<class T> PARAMETERS_COMMON<T>::
PARAMETERS_COMMON()
    :time(0),dt(0),rho(0),mu(0),theta_threshold(0),cg_tolerance(0),print_matrix(false)
{
}

template<class T> void PARAMETERS_COMMON<T>::
Init_1(PARSE_ARGS& parse_args)
{
    parse_args.Add_Double_Argument("-time",(T)0,"starting time");
    parse_args.Add_Double_Argument("-dt",(T).01,"time step size");

    parse_args.Add_Double_Argument("-rho",rho,"Density");
    parse_args.Add_Double_Argument("-viscosity",0,"Viscosity");
    
    parse_args.Add_Option_Argument("-print_matrix","Print all matrices");

    parse_args.Add_Double_Argument("-theta_tol",1e-8,"zero tolerance for theta");
    parse_args.Add_Double_Argument("-cg_tol",1e-14,"zero tolerance for CG");
}
template<class TV> void PARAMETERS_COMMON<TV>::
Init_2(PARSE_ARGS& parse_args)
{
    time=parse_args.Get_Double_Value("-time");
    dt=parse_args.Get_Double_Value("-dt");

    rho=parse_args.Get_Double_Value("-rho");
    mu=parse_args.Get_Double_Value("-viscosity");

    print_matrix=parse_args.Is_Value_Set("-print_matrix");

    theta_threshold=parse_args.Get_Double_Value("-theta_tol");
    cg_tolerance=parse_args.Get_Double_Value("-cg_tol");
}

template struct PARAMETERS_COMMON<double>;
