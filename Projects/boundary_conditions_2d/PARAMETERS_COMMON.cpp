#include "PARAMETERS_COMMON.h"

template<class T> PARAMETERS_COMMON<T>::
PARAMETERS_COMMON()
    :time(0),dt((T)0.01),rho(0),mu(0),theta_threshold((T)1e-8),cg_tolerance((T)1e-14),print_matrix(false)
{
}

template<class T> void PARAMETERS_COMMON<T>::
Init_1(PARSE_ARGS& parse_args)
{
    parse_args.Add("-time",&time,"value","starting time");
    parse_args.Add("-dt",&dt,"value","time step size");

    parse_args.Add("-rho",&rho,"value","Density");
    parse_args.Add("-viscosity",&mu,"value","Viscosity");
    
    parse_args.Add("-print_matrix",&print_matrix,"Print all matrices");

    parse_args.Add("-theta_tol",&theta_threshold,"value","zero tolerance for theta");
    parse_args.Add("-cg_tol",&cg_tolerance,"value","zero tolerance for CG");
}
template<class TV> void PARAMETERS_COMMON<TV>::
Init_2(PARSE_ARGS& parse_args)
{
}

template struct PARAMETERS_COMMON<double>;
