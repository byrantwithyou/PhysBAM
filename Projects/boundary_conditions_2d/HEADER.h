#ifndef __HEADER__
#define __HEADER__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Vectors/VECTOR.h>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS.h"
using namespace PhysBAM;

typedef float RW;

template<class T>
T Source_Velocity_Curve(bool d0,bool d1,T x);

template<class T,class TV,int d>
void Apply_Viscosity(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,T dt,T time,T viscosity,T density,int axis,T theta_threshold,
    T cg_tolerance,bool verbose);

enum proj_type {proj_default,proj_gibou,proj_slip};

template<class T,class TV,int d>
void Project_Incompressibility(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose);

template<class T,class TV,int d>
void Project_Incompressibility_Gibou(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose);

template<class T,class TV,int d>
void Project_Incompressibility_Slip(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& callback,const ACCURACY_INFO<d>& ai,T time,T density,
    T theta_threshold,T cg_tolerance,bool verbose);

template<class TV>
void Fill_Ghost_Cells(const GRID<TV>& grid,int ghost,int distance,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,const BOUNDARY_CONDITIONS<TV>& bc);

template<class TV>
void Initialize_Grid_From_Domains(GRID<TV>& grid,int resolution,const RANGE<TV>& sample_box,const RANGE<TV>& bounding_box,RANGE<VECTOR<int,TV::m> >& sample_domain);

extern std::string output_directory;

template<class TV>
void Prune_Outside_Sample_Points(const GRID<TV>& grid,const BOUNDARY_CONDITIONS<TV>& bc,ACCURACY_INFO<TV::m>& ai);

template<class TV> struct OBJECTS_COMMON;
template<class T> struct PARAMETERS_COMMON;
template<class T> struct SIM_COMMON;

template<class TV>
void Add_Viscosity(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation);

template<class TV>
void Add_Advection(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation);

template<class TV>
void Project_Pressure(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation,proj_type proj);

template<class TV>
void First_Order_Step(SIM_COMMON<TV>& sim,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u);

template<class TV>
void Second_Order_RE_Step(SIM_COMMON<TV>& sim,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u2);

template<class TV>
void Dump_Error(const SIM_COMMON<TV>& sim,const ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u2);

namespace PhysBAM{template<class T,class T2> class INTERPOLATION_CURVE;}

template<class T>
struct ERROR_COLOR_MAP
{
    T mn,mx;
    bool use_log,reverse_order,use_two_sets;
    INTERPOLATION_CURVE<T,VECTOR<T,3> >& colors;
    ERROR_COLOR_MAP(T min_value,T max_value,bool log_scale,bool reverse,bool two_sets);
    ~ERROR_COLOR_MAP();
    VECTOR<T,3> operator()(T x) const;
};

template<class TV>
void Dump_Error_Image(const SIM_COMMON<TV>& sim,const ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,VECTOR<int,TV::m> dim);

#endif

