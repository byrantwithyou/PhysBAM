//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/NONCOPYABLE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Computations/GRADIENT_UNIFORM.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_AUXILIARY_DATA.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
#include <Compressible/Equations_Of_State/EOS.h>
#include <Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{
namespace COMPRESSIBLE_AUXILIARY_DATA{
template<class TV,class T_ARRAYS,class T_ARRAYS_BOOL_INPUT,class T_FACE_ARRAYS_DIMENSION_SCALAR>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const GRID<TV>& grid,
    const int number_of_ghost_cells,const T_ARRAYS& U,const T_ARRAYS_BOOL_INPUT& psi,const EOS<typename TV::SCALAR>& eos,const bool write_debug_data,const T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;

    STATIC_ASSERT((is_same<T_ARRAYS,T_ARRAYS_DIMENSION_SCALAR>::value));
    STATIC_ASSERT((is_same<T_ARRAYS_BOOL_INPUT,ARRAY<bool,TV_INT> >::value));

    RANGE<TV_INT> domain_indices=grid.Domain_Indices(number_of_ghost_cells);
    ARRAY<TV,TV_INT> velocity(domain_indices),velocity_plus_c(domain_indices),velocity_minus_c(domain_indices),momentum(domain_indices);
    ARRAY<T,TV_INT> density(domain_indices),energy(domain_indices),internal_energy(domain_indices),temperature(domain_indices),
        pressure(domain_indices),entropy(domain_indices),enthalpy(domain_indices),speedofsound(domain_indices),
        machnumber(domain_indices),density_gradient(domain_indices),pressure_gradient(domain_indices);
    for(CELL_ITERATOR<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(grid.Domain_Indices().Lazy_Inside_Half_Open(cell) && !psi(cell)) continue;
        density(cell)=EULER<TV>::Get_Density(U,cell);
        pressure(cell)=eos.p(density(cell),EULER<TV>::e(U(cell)));
        energy(cell)=EULER<TV>::Get_Total_Energy(U,cell);
        internal_energy(cell)=EULER<TV>::e(U(cell));
        temperature(cell)=eos.T(density(cell),internal_energy(cell));
        velocity(cell)=EULER<TV>::Get_Velocity(U,cell);
        if(write_debug_data){
            momentum(cell)=velocity(cell)*density(cell);
            entropy(cell)=eos.S(density(cell),eos.e_From_p_And_rho(pressure(cell),density(cell)));
            enthalpy(cell)=EULER<TV>::enthalpy(eos,U(cell));
            speedofsound(cell)=eos.c(density(cell),eos.e_From_p_And_rho(pressure(cell),density(cell)));
            if(speedofsound(cell)) machnumber(cell)=velocity(cell).Magnitude()/speedofsound(cell);
            else machnumber(cell)=(T)-1; // output non-physical negative value when machnumber is infinite
            velocity_plus_c(cell)=velocity(cell)+speedofsound(cell)*TV::All_Ones_Vector();
            velocity_minus_c(cell)=velocity(cell)-speedofsound(cell)*TV::All_Ones_Vector();}}

    ARRAY<T,FACE_INDEX<TV::m> > density_flux(grid),momentum_flux(grid),energy_flux(grid);
    if(fluxes){
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face_index=iterator.Face_Index();
            if(fluxes->Valid_Index(iterator.Full_Index())){int axis=iterator.Axis();
                density_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(0);
                momentum_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(1);
                energy_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(TV::m+1);}}}

    std::string f=LOG::sprintf("%d",frame);
    if(write_debug_data){
        GRADIENT::Compute_Magnitude(grid,number_of_ghost_cells,density,density_gradient);
        GRADIENT::Compute_Magnitude(grid,number_of_ghost_cells,pressure,pressure_gradient);
        
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density_flux",density_flux);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/momentum_flux",momentum_flux);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/energy_flux",energy_flux);

        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/momentum",momentum);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/energy",energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density_gradient",density_gradient);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure_gradient",pressure_gradient);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/entropy",entropy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/enthalpy",enthalpy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/speedofsound",speedofsound);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/machnumber",machnumber);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/internal_energy",internal_energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocity_plus_c",velocity_plus_c);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocity_minus_c",velocity_minus_c);}

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/temperature",temperature);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/centered_velocities",velocity);
}
template<class TV>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const COMPRESSIBLE_FLUID_COLLECTION<TV>& compressible_fluid_collection,const bool write_debug_data)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<TV_DIMENSION,FACE_INDEX<TV::m> > T_FACE_ARRAYS_DIMENSION_SCALAR;

    const GRID<TV>& grid=compressible_fluid_collection.grid;
    const T_ARRAYS_DIMENSION_SCALAR& U=compressible_fluid_collection.U;
    const EOS<T>* eos=compressible_fluid_collection.eos;
    const T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes = NULL;

    Write_Auxiliary_Files(stream_type,output_directory,frame,grid,0,U,compressible_fluid_collection.psi,*eos,write_debug_data,fluxes);
}
template void Write_Auxiliary_Files<VECTOR<float,1>,ARRAY<VECTOR<float,3>,VECTOR<int,1> >,ARRAY<bool,VECTOR<int,1> >,ARRAY<VECTOR<float,3>,FACE_INDEX<1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,1> > const&,int,ARRAY<VECTOR<float,3>,VECTOR<int,1> > const&,ARRAY<bool,VECTOR<int,1> > const&,EOS<VECTOR<float,1>::SCALAR> const&,const bool,const ARRAY<VECTOR<float,3>,FACE_INDEX<1> >*);
template void Write_Auxiliary_Files<VECTOR<float,2>,ARRAY<VECTOR<float,4>,VECTOR<int,2> >,ARRAY<bool,VECTOR<int,2> >,ARRAY<VECTOR<float,4>,FACE_INDEX<2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,2> > const&,int,ARRAY<VECTOR<float,4>,VECTOR<int,2> > const&,ARRAY<bool,VECTOR<int,2> > const&,EOS<VECTOR<float,2>::SCALAR> const&,const bool,const ARRAY<VECTOR<float,4>,FACE_INDEX<2> >*);
template void Write_Auxiliary_Files<VECTOR<float,3>,ARRAY<VECTOR<float,5>,VECTOR<int,3> >,ARRAY<bool,VECTOR<int,3> >,ARRAY<VECTOR<float,5>,FACE_INDEX<3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,3> > const&,int,ARRAY<VECTOR<float,5>,VECTOR<int,3> > const&,ARRAY<bool,VECTOR<int,3> > const&,EOS<VECTOR<float,3>::SCALAR> const&,const bool,const ARRAY<VECTOR<float,5>,FACE_INDEX<3> >*);
template void Write_Auxiliary_Files<VECTOR<float,1> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,1> > const&,const bool);
template void Write_Auxiliary_Files<VECTOR<float,2> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,2> > const&,const bool);
template void Write_Auxiliary_Files<VECTOR<float,3> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,3> > const&,const bool);
template void Write_Auxiliary_Files<VECTOR<double,1>,ARRAY<VECTOR<double,3>,VECTOR<int,1> >,ARRAY<bool,VECTOR<int,1> >,ARRAY<VECTOR<double,3>,FACE_INDEX<1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,1> > const&,int,ARRAY<VECTOR<double,3>,VECTOR<int,1> > const&,ARRAY<bool,VECTOR<int,1> > const&,EOS<VECTOR<double,1>::SCALAR> const&,const bool,const ARRAY<VECTOR<double,3>,FACE_INDEX<1> >*);
template void Write_Auxiliary_Files<VECTOR<double,2>,ARRAY<VECTOR<double,4>,VECTOR<int,2> >,ARRAY<bool,VECTOR<int,2> >,ARRAY<VECTOR<double,4>,FACE_INDEX<2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,2> > const&,int,ARRAY<VECTOR<double,4>,VECTOR<int,2> > const&,ARRAY<bool,VECTOR<int,2> > const&,EOS<VECTOR<double,2>::SCALAR> const&,const bool,const ARRAY<VECTOR<double,4>,FACE_INDEX<2> >*);
template void Write_Auxiliary_Files<VECTOR<double,3>,ARRAY<VECTOR<double,5>,VECTOR<int,3> >,ARRAY<bool,VECTOR<int,3> >,ARRAY<VECTOR<double,5>,FACE_INDEX<3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,3> > const&,int,ARRAY<VECTOR<double,5>,VECTOR<int,3> > const&,ARRAY<bool,VECTOR<int,3> > const&,EOS<VECTOR<double,3>::SCALAR> const&,const bool,const ARRAY<VECTOR<double,5>,FACE_INDEX<3> >*);
template void Write_Auxiliary_Files<VECTOR<double,1> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,1> > const&,const bool);
template void Write_Auxiliary_Files<VECTOR<double,2> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,2> > const&,const bool);
template void Write_Auxiliary_Files<VECTOR<double,3> >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,3> > const&,const bool);
}
}
