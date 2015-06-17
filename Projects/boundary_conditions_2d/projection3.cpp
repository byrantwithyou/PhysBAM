#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <functional>
#include "PROJECT.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    PARSE_ARGS parse_args(argc,argv);
    int refine=1,resolution=32,interface=0,velocity_field=0;
    T rho=1,kg=1,m=1,s=1;
    bool test_analytic_diff=false,use_p_null_mode=false,dump_matrix=false,use_bc=true;
    parse_args.Extra(&interface,"number","interface to use");
    parse_args.Extra(&velocity_field,"number","velocity to use");
    parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
    parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
    parse_args.Add("-rho",&rho,"density","density for first fluid region");
    parse_args.Add("-m",&m,"scale","meter scale");
    parse_args.Add("-s",&s,"scale","second scale");
    parse_args.Add("-kg",&kg,"scale","kilogram scale");
    parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
    parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
    parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
    parse_args.Add_Not("-no_bc",&use_bc,"Disable u boundary condition term");
    parse_args.Parse();

    T unit_p=kg/(s*s)*pow(m,2-TV::m);

    std::function<T(TV X)> phi,p;
    std::function<TV(TV X)> u_star,u_projected;

    RANGE<TV> domain=RANGE<TV>::Unit_Box();

    switch(interface){
        case 0:phi=[](TV X){return (X-.5).Magnitude()-.3;};break;
        case 1:phi=[](TV X){return X.Magnitude()-.7;};domain=RANGE<TV>::Centered_Box();break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized interface");}

    GRID<TV> grid(TV_INT()+resolution,domain*m,true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    switch(velocity_field){
        case 0:
            u_star=[](TV X){return TV(0,0,0);};
            u_projected=[](TV X){return TV(0,0,0);};
            p=[](TV X){return 0;};
            break;
        case 1:
            u_star=[](TV X){return TV(2,1,3);};
            u_projected=[](TV X){return TV(2,1,3);};
            p=[](TV X){return 0;};
            break;
        case 2:
            u_star=[](TV X){return X*TV(2,-1,-1);};
            u_projected=[](TV X){return X*TV(2,-1,-1);};
            p=[](TV X){return 0;};
            break;
        case 3:
            u_star=[](TV X){return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y),0);};
            u_projected=[](TV X){return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y),0);};
            p=[](TV X){return 0;};
            break;
        case 4:
            u_star=[](TV X){return TV(2,1,3);};
            u_projected=[](TV X){return TV();};
            p=[=](TV X){return rho*X.Dot(TV(2,1,3));};
            break;
        case 5:
            u_star=[](TV X){return X*TV(2,-1,-1);};
            u_projected=[](TV X){return TV();};
            p=[=](TV X){return rho*X.Dot(X*TV(2,-1,-1))/2;};
            break;
        case 6:
            u_star=[](TV X){return X;};
            u_projected=[](TV X){return TV(1,0,0);};
            p=[=](TV X){return rho*(X.Magnitude_Squared()/2-X.x);};
            break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized velocity");}

    ARRAY<T,TV_INT> node_phi(grid.Node_Indices());
    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        node_phi(it.index)=phi(it.Location()/m)*m;

    Project<T,TV>(grid,3,node_phi,[=](TV X){return u_star(X/m)*m/s;},[=](TV X){return u_projected(X/m)*m/s;},[=](TV X){return p(X/m)*unit_p*s;},rho*kg*pow(m,-TV::m),
        (T)1e-14,(T)1e-14,use_p_null_mode,use_bc,dump_matrix);
    Flush_Frame<TV>("flush");
    return 0;
}


