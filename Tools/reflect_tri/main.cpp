//#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
//#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_IMPLICIT_SURFACE.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LIGHT.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
//#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
//#include <cstring>
#include <fstream>
//#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
//#include "../../Public_Library/Geometry/OCTREE_IMPLICIT_SURFACE.h"
//#include "../../Public_Library/Level_Sets/OCTREE_LEVELSET.h"
//#include "../../Public_Library/OpenGL/OPENGL_HEXAHEDRONS.h"
//#include "../../Public_Library/OpenGL/OPENGL_TETS.h"
#include "TRI_REFLECTOR.h"
using namespace PhysBAM;
using namespace std;

template<class T> void Add_Tri_File(const std::string& filename,OPENGL_WORLD& world,int number);
//#################################################################
// function main
//#################################################################
int main(int argc,char** argv)
{
     OPENGL_WORLD world;
    TRI_REFLECTOR<float> reflector(world);

    //setup lighting
    world.Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(1,1,1),float(.8)));
    world.Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(.3,-2,.1),float(.4)));
    world.Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(-2,.1,.5),float(.2)));
    world.Set_Ambient_Light(.2);
    world.Initialize("Reflect_Tri");

    reflector.Initialize_Tri_Reflector_Callbacks(argc,argv);

    world.Center_Camera_On_Scene();
    world.Start_Main_Loop();
    return 0;
}
//#################################################################
