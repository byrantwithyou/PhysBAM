#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Analytic_Tests/VORTEX_IMPLICIT_SURFACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "PROJECT.h"
#include <boost/function.hpp>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    PARSE_ARGS parse_args(argc,argv);
    int refine=1,resolution=32,interface=0,velocity_field=0;
    T rho=1,kg=1,m=1,s=1;
    bool test_analytic_diff=false,use_p_null_mode=false,dump_matrix=false;
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
    parse_args.Parse();

    T unit_p=kg/(s*s)*pow(m,2-TV::m);

    boost::function<T(TV X)> phi,p;
    boost::function<TV(TV X)> u_star,u_projected;

    RANGE<TV> domain=RANGE<TV>::Unit_Box();

    VORTEX_IMPLICIT_SURFACE<TV> vis;

    switch(interface){
        case 0:phi=[](TV X){return (X-.5).Magnitude()-.3;};break;
        case 1:phi=[](TV X){return X.Magnitude()-.7;};domain=RANGE<TV>::Centered_Box();break;
        case 2:vis.k=(T).2;phi=[=](TV X){return vis.Phi(X);};domain.max_corner*=(T)pi;break;
        case 3:phi=[](TV X){return X.Magnitude()-.8;};domain=RANGE<TV>::Centered_Box();break;
        case 4:vis.k=(T).2;phi=[=](TV X){return .2-sin(X.x)*sin(X.y);};domain.max_corner*=(T)pi;break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized interface");}

    GRID<TV> grid(TV_INT()+resolution,domain*m,true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    switch(velocity_field){
        case 0:
            u_star=[](TV X){return TV(0,0);};
            u_projected=[](TV X){return TV(0,0);};
            p=[](TV X){return 0;};
            break;
        case 1:
            u_star=[](TV X){return TV(2,1);};
            u_projected=[](TV X){return TV(2,1);};
            p=[](TV X){return 0;};
            break;
        case 2:
            u_star=[](TV X){return X*TV(1,-1);};
            u_projected=[](TV X){return X*TV(1,-1);};
            p=[](TV X){return 0;};
            break;
        case 3:
            u_star=[](TV X){return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y));};
            u_projected=[](TV X){return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y));};
            p=[](TV X){return 0;};
            break;
        case 4:
            u_star=[](TV X){return TV(2,1);};
            u_projected=[](TV X){return TV();};
            p=[=](TV X){return rho*X.Dot(TV(2,1));};
            break;
        case 5:
            u_star=[](TV X){return X*TV(1,-1);};
            u_projected=[](TV X){return TV();};
            p=[=](TV X){return rho*X.Dot(X*TV(1,-1))/2;};
            break;
        case 6:
            u_star=[](TV X){return X;};
            u_projected=[](TV X){return TV(1,0);};
            p=[=](TV X){return rho*(X.Magnitude_Squared()/2-X.x);};
            break;
        case 7:
            u_star=[](TV X){return TV(X.y,-X.x)*2;};
            u_projected=[](TV X){return TV(X.y,-X.x)*2;};
            p=[=](TV X){return 0;};
            break;
        case 8:
            u_star=[](TV X){return TV(sin(X.x)*cos(X.y)+(X.x*X.x-X.x)*(X.y*X.y*X.y/3-X.y*X.y/2),-cos(X.x)*sin(X.y)+(X.y*X.y-X.y)*(X.x*X.x*X.x/3-X.x*X.x/2));};
            u_projected=[](TV X){return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y));};
            p=[](TV X){return (X.x*X.x*X.x/3-X.x*X.x/2)*(X.y*X.y*X.y/3-X.y*X.y/2);};
            break;
        case 9:
            u_star=[](TV X){return TV((X.x*X.x-X.x)*(X.y*X.y*X.y/3-X.y*X.y/2),(X.y*X.y-X.y)*(X.x*X.x*X.x/3-X.x*X.x/2));};
            u_projected=[](TV X){return TV();};
            p=[](TV X){return (X.x*X.x*X.x/3-X.x*X.x/2)*(X.y*X.y*X.y/3-X.y*X.y/2);};
            break;
        case 10:
            u_star=[](TV X){return TV(X.y,X.x);};
            u_projected=[](TV X){return TV();};
            p=[](TV X){return X.x*X.y;};
            break;
        default: PHYSBAM_FATAL_ERROR("Unrecognized velocity");}

    ARRAY<T,TV_INT> node_phi(grid.Node_Indices());
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        node_phi(it.index)=phi(it.Location()/m)*m;

    Project<T,TV>(grid,3,node_phi,[=](TV X){return u_star(X/m)*m/s;},[=](TV X){return u_projected(X/m)*m/s;},[=](TV X){return p(X/m)*unit_p*s;},rho*kg*pow(m,-TV::m),(T)1e-14,(T)1e-14,use_p_null_mode);
    Flush_Frame<TV>("flush");
    return 0;
}
