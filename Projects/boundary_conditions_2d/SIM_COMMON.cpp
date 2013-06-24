#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <cmath>
#include <iomanip>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS.h"
#include "SIM_COMMON.h"

template<class TV> SIM_COMMON<TV>::
SIM_COMMON()
    :resolution(0),steps(1),base_resolution(4),proj_algo(proj_default),check_leaks(false),use_accuracy_samples(false),
    use_extrapolation(true),use_viscosity(true),use_projection(true),use_advection(true),use_proj_slip(false),no_gibou(false)
{
}
template<class TV> void SIM_COMMON<TV>::
Init_1(PARSE_ARGS& parse_args)
{
    param.Init_1(parse_args);

    parse_args.Add("-resolution",&resolution,"res","Resolution");
    parse_args.Add("-steps",&steps,"value","Perform multiple time steps");

    parse_args.Add("-leak",&check_leaks,"check for leaks");
    parse_args.Add("-no_gibou",&no_gibou,"Use Gibou for Neumann");
    parse_args.Add("-proj_slip",&use_proj_slip,"Use Slip for Neumann");

    parse_args.Add("-sample",&use_accuracy_samples,"Use sample points");
    parse_args.Add_Not("-no_extrap",&use_extrapolation,"Disable extrapolation");

    parse_args.Add_Not("-no_viscosity",&use_viscosity,"Do not apply viscosity");
    parse_args.Add_Not("-no_projection",&use_projection,"Do not perform pressure projection");
    parse_args.Add_Not("-no_advection",&use_advection,"Do not perform advection");
}
template<class TV> void SIM_COMMON<TV>::
Init_2(PARSE_ARGS& parse_args)
{
    param.Init_2(parse_args);
    param.rho=1;

    if(use_proj_slip) proj_algo=proj_slip;
    else if(!no_gibou) proj_algo=proj_gibou;
}
template<class TV> void SIM_COMMON<TV>::
Init_3()
{
    obj.bc->check_leaks=check_leaks;
    obj.ai.resolution=resolution;
    obj.ai.base_resolution=base_resolution;
    obj.ai.sample_domain.min_corner=VECTOR<int,TV::m>::All_Ones_Vector();
    obj.ai.sample_domain.max_corner=VECTOR<int,TV::m>::All_Ones_Vector()*resolution;
    Initialize_Grid_From_Domains(obj.grid,resolution,obj.bc->base_domain,obj.bc->bounding_box,obj.ai.sample_domain);

    obj.bc->Initialize_Phi(5);
    if(use_accuracy_samples) obj.ai.Compute();
    Prune_Outside_Sample_Points(obj.grid,*obj.bc,obj.ai);
    obj.ai.Print_Locations(obj.grid);
}

template struct SIM_COMMON<VECTOR<double,1> >;
template struct SIM_COMMON<VECTOR<double,2> >;
