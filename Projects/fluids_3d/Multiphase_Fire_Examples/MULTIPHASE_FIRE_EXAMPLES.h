//#####################################################################
// Copyright 2005, Frank Losasso
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIPHASE_FIRE_EXAMPLES 
//#####################################################################
#ifndef __MULTIPHASE_FIRE_EXAMPLES__
#define __MULTIPHASE_FIRE_EXAMPLES__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include "MULTIPHASE_FIRE_EXAMPLES_UNIFORM.h"
namespace PhysBAM{

template<class T_input>
class MULTIPHASE_FIRE_EXAMPLES:public MULTIPHASE_FIRE_EXAMPLES_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef MULTIPHASE_FIRE_EXAMPLES_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;using BASE::pseudo_dirichlet;using BASE::resolution;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::test_number;
    
    MULTIPHASE_FIRE_EXAMPLES(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :MULTIPHASE_FIRE_EXAMPLES_UNIFORM<TV>(stream_type_input,parse_args)
    {
        parse_args.Parse();

        LOG::cout<<"Running Multiphase Fire Example Number "<<test_number<<" at resolution "<<resolution<<std::endl;
        fluids_parameters.solve_neumann_regions=false;
        int cells=1*resolution;
        if(test_number==1) fluids_parameters.grid->Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(1,1,1)));
        if(test_number==2) fluids_parameters.grid->Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(1,1,1)));
        if(test_number==3) fluids_parameters.grid->Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(1,1,1)));
        if(!this->user_output_directory){
            if(!pseudo_dirichlet) output_directory=LOG::sprintf("Multiphase_Fire_Examples/Example_%d__Resolution_%d_%d_%d",test_number,
                (fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),(fluids_parameters.grid->counts.z-1));
            else output_directory=LOG::sprintf("Multiphase_Fire_Examples/Example_%d__Resolution_%d_%d_%d_pseudo_dirichlet",test_number,
                (fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),(fluids_parameters.grid->counts.z-1));}
    }

//#####################################################################
};
}
#endif
