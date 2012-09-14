#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include "LEVELSET_MESH_CUTTING_3D.h"

int main()
{
    typedef float T;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<T,3> TV;

    GRID<TV> grid(TV_INT()+2,RANGE<TV>::Unit_Box());

    TETRAHEDRALIZED_VOLUME<T>& tv=*TETRAHEDRALIZED_VOLUME<T>::Create();
    tv.Initialize_Cube_Mesh_And_Particles(grid);

    ARRAY<T> phi0(tv.particles.number),phi1(tv.particles.number);
    RANDOM_NUMBERS<T> random;
    random.Fill_Uniform(phi0,-1,1);
    random.Fill_Uniform(phi1,-1,1);


/*        struct TET
    {
        TV_INT4 parent,indices;
        VECTOR<TV4,4> weights;
    };
*/
    LOG::cout<<tv.mesh.elements<<std::endl;

    ARRAY<LEVELSET_MESH_CUTTING_3D::TET> cut_mesh;

    LEVELSET_MESH_CUTTING_3D::Subdivide(tv.mesh.elements,phi0,phi1,cut_mesh);

    for(int i=0;i<cut_mesh.m;i++){
        LOG::cout<<cut_mesh(i).parent<<"  "<<cut_mesh(i).indices<<"  "<<cut_mesh(i).weights<<std::endl;}

    delete &tv;
    return 0;
}
