#ifndef __BOUNDARY_CONDITIONS__
#define __BOUNDARY_CONDITIONS__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET.h>

using namespace PhysBAM;

template<class T> struct PARAMETERS_COMMON;

template<class TV>
struct BOUNDARY_CONDITIONS
{
    enum BOUNDARY_TYPE {none,free,noslip,slip};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    GRID<TV>& grid;
    RANGE<TV> base_domain,bounding_box;
    ARRAY<T,TV_INT> phi_array;
    bool check_leaks;
    bool use_analytic_solution;
    LEVELSET<TV>* phi;

    BOUNDARY_CONDITIONS(GRID<TV>& grid_input);
    virtual ~BOUNDARY_CONDITIONS();
    virtual bool Inside(const TV& X) const;
    virtual T Theta(const TV& X) const=0;
    virtual TV Gradient(const TV& X) const;
    virtual int Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const=0;
    virtual void Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const=0;
    virtual void Compute_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
    virtual TV Analytic_Velocity(const TV& X,T time) const=0;
    virtual void Update_Parameters(PARAMETERS_COMMON<T>& param);
    void Initialize_Phi(int ghost);
};

#endif
