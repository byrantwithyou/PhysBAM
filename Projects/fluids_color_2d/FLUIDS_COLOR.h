//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>

namespace PhysBAM{

template<class TV>
class FLUIDS_COLOR:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::particle_levelset_evolution;using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,const PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type)
    {
        int test_number=1;last_frame=200;
        int resolution=parse_args.Get_Integer_Value("-resolution");
        restart=parse_args.Get_Integer_Value("-restart");
        write_substeps_level=parse_args.Get_Integer_Value("-substep");
        output_directory=STRING_UTILITIES::string_sprintf("Water_Tests/Test_%d",test_number);
        grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Initialize()
    {
        dt=(T).1;
        mu.Append(1);
        face_velocities.Fill(0);
        levelset_color.phi.Fill(-1);
        levelset_color.color.Fill(0);
    }

//#####################################################################
};
}

#endif
