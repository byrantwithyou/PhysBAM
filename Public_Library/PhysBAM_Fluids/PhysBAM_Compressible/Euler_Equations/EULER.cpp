//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER
//##################################################################### 
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> EULER<T_GRID>::
EULER()
{
    cut_out_grid=false;
    use_solid_velocity_in_ghost_cells=false;
    use_force=false;
    boundary=&boundary_default;
    conservation=&conservation_default;
    if(use_solid_velocity_in_ghost_cells) conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_SOLID_VELOCITY<T_GRID>);
    else conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_EULER<T_GRID>);
    eos=&eos_default;
    Set_Max_Time_Step();
    Set_Gravity(0);
    Set_CFL_Number();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> EULER<T_GRID>::
~EULER()
{}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class T_GRID> void EULER<T_GRID>::
Log_Parameters() const
{
    LOG::SCOPE scope("EULER parameters");
    LOG::cout<<"cfl_number="<<cfl_number<<std::endl;
    LOG::cout<<"open_boundaries="<<open_boundaries<<std::endl;
    LOG::cout<<"cut_out_grid="<<cut_out_grid<<std::endl;
    LOG::cout<<"use_force="<<use_force<<std::endl;
    LOG::cout<<"max_time_step="<<max_time_step<<std::endl;
    conservation->Log_Parameters();
}
namespace PhysBAM{
template class EULER<GRID<VECTOR<float,1> > >;
template class EULER<GRID<VECTOR<float,2> > >;
template class EULER<GRID<VECTOR<float,3> > >;
template class EULER<GRID<VECTOR<double,1> > >;
template class EULER<GRID<VECTOR<double,2> > >;
template class EULER<GRID<VECTOR<double,3> > >;
}
