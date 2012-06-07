#include <PhysBAM_Tools/Log/LOG.h>
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"
#include <PhysBAM_Dynamics/Meshing/COLLISION_AWARE_SUBDIVISION.h>

bool PHYSBAM_THREADED_RUN=false;

using namespace PhysBAM;

template<class T>
void Subdivide(int frame)
{
    LOG::Time("Reading files");
    DEFORMABLE_OBJECT_LIST_3D<T> deformable_list;
    deformable_list.template Read_Static_Variables<float>("./");
    deformable_list.template Read_Dynamic_Variables<float>("./",frame);
    RIGID_BODY_LIST_3D<T> rigid_list;
    rigid_list.template Read<T>(".",frame);

    LOG::Time("Writing original");
    FILE_UTILITIES::Write_To_File<T>("original.tri",*deformable_list(1).triangulated_surface);

    LOG::Time("Building subidvider");
    COLLISION_AWARE_SUBDIVISION<T> collision_aware_subdivision;collision_aware_subdivision.collision_tolerance=(T)1e-3;
    //collision_aware_subdivision.Add_Rigid_Bodies(rigid_list);
    collision_aware_subdivision.Add_Triangulated_Surface(*deformable_list(1).triangulated_surface);

    LOG::Time("Subidviding");
    collision_aware_subdivision.Subdivide(1);

    LOG::Time("Writing results");
    FILE_UTILITIES::Write_To_File<T>("subdivided.tri",*deformable_list(1).triangulated_surface);


}

int main(int argc,char *argv[])
{
    Subdivide<float>(atoi(argv[1]));
    return 0;
}
