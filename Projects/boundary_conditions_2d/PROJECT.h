#ifndef __PROJECT__
#define __PROJECT__
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <functional>
using namespace PhysBAM;

template<class T,class TV,class TV_INT>
void Project(const GRID<TV>& grid,int ghost,const ARRAY<T,TV_INT>& phi,std::function<TV(TV X)> u_star,
    std::function<TV(TV X)> u_projected,std::function<T(TV X)> p,T density,T theta_threshold,T cg_tolerance,
    bool use_p_null_mode,bool use_bc,bool print_matrix);
#endif
