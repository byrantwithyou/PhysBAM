#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <cmath>
#include <iomanip>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS.h"
#include "SIM_COMMON.h"

template<class TV> SIM_COMMON<TV>::
SIM_COMMON()
    :resolution(0),steps(0),base_resolution(4),proj_algo(proj_default),check_leaks(false),use_accuracy_samples(false),
    use_extrapolation(true),use_viscosity(false),use_projection(false),use_advection(false)
{
}

template<class TV> void SIM_COMMON<TV>::
Init_1(PARSE_ARGS& parse_args)
{
    param.Init_1(parse_args);

    parse_args.Add_Integer_Argument("-resolution",0,"Resolution");
    parse_args.Add_Integer_Argument("-steps",1,"Perform multiple time steps");

    parse_args.Add_Option_Argument("-leak","check for leaks");
    parse_args.Add_Option_Argument("-no_gibou","Use Gibou for Neumann");
    parse_args.Add_Option_Argument("-proj_slip","Use Slip for Neumann");

    parse_args.Add_Option_Argument("-sample","Use sample points");
    parse_args.Add_Option_Argument("-no_extrap","Disable extrapolation");

    parse_args.Add_Option_Argument("-no_viscosity","Do not apply viscosity");
    parse_args.Add_Option_Argument("-no_projection","Do not perform pressure projection");
    parse_args.Add_Option_Argument("-no_advection","Do not perform advection");
}
template<class TV> void SIM_COMMON<TV>::
Init_2(PARSE_ARGS& parse_args)
{
    param.Init_2(parse_args);
    param.rho=1;

    resolution=parse_args.Get_Integer_Value("-resolution");
    steps=parse_args.Get_Integer_Value("-steps");

    check_leaks=parse_args.Is_Value_Set("-leak");
    if(parse_args.Is_Value_Set("-proj_slip")) proj_algo=proj_slip;
    else if(!parse_args.Is_Value_Set("-no_gibou")) proj_algo=proj_gibou;

    use_accuracy_samples=parse_args.Is_Value_Set("-sample");
    use_extrapolation=!parse_args.Is_Value_Set("-no_extrap");

    use_viscosity=!parse_args.Is_Value_Set("-no_viscosity");
    use_advection=!parse_args.Is_Value_Set("-no_advection");
    use_projection=!parse_args.Is_Value_Set("-no_projection");
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
