/* Program to do point repulsion on a tri surface, by Jiayi Chong*/
/* To use the visualizer, you have to compile this program with fltk, opengl and opengl extensions.*/
/* That means you have to have preprocessor flags: WITH_FLTK, WITH_OPENGL_EXTENSIONS defined*/
/* Also, you have to link to the fltk, fltkgl and glux libraries*/

#define WITH_VISUALIZER 0
#include <PhysBAM_Geometry/Collisions/POINT_REPULSION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VBO_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <iostream>
#include "../../Public_Library/OpenGL/OPENGL_UNIVERSE.h"
#include <stdlib.h>
 
using namespace std;
using namespace PhysBAM;

template<class T> TRIANGULATED_SURFACE<T> *Add_Tri_File(const char* filename);

#if WITH_VISUALIZER
class myUniverse:public OPENGL_UNIVERSE {
public:
    int window1_id;
    TRIANGULATED_SURFACE<float> *triangulated_surface;
    POINT_REPULSION<float> *point_repulsion;
    OPENGL_VBO_TRIANGULATED_SURFACE<float> *ots;
    
    void Create() {
        window1_id = Add_Window("3-D Viewer", 0, 0, 640, 480);   
    }
    
    void Initialize_Window(int id) {
        Init();
    }   
    void Init() {
        ots=new OPENGL_VBO_TRIANGULATED_SURFACE<float>(*triangulated_surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.9),float(.9))));  
        ots->Create_VBO();
        Get_Window(window1_id)->Add_Object(ots,true,true);
        Get_Window(window1_id)->Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(1,1,1),float(.8)));
        Get_Window(window1_id)->Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(.3,-2,.1),float(.4)));
        Get_Window(window1_id)->Add_Light(new OPENGL_LIGHT(VECTOR_3D<double>(-2,.1,.5),float(.2)));
        Get_Window(window1_id)->Set_Ambient_Light(.2);
        Get_Window(window1_id)->Center_Camera_On_Scene();       
    }
    
    void Custom_Draw(int window_id) {
        glPointSize(3.0);
        glBegin(GL_POINTS);
        glColor3f(0.0,1.0,0.0);
        for(int i = 1; i <= point_repulsion->points.m; ++i) {
            glVertex3f(point_repulsion->points(i).position.x,point_repulsion->points(i).position.y,point_repulsion->points(i).position.z);
        }
        glEnd();
    }
};
#endif

int main(int argc,char *argv[])
{
    if(argc != 5 && argc != 6 && argc != 7) {
        std::cerr<<"Syntax: <input .tri file> <output point file> <iterations> <points> <scale factor> <Visualization 1/0>" << std::endl;
        exit(1);
    }
    
    cerr<<"Reading triangulated surface '"<<argv[1]<<"'..."<<endl;
    TRIANGULATED_SURFACE<float> *triangulated_surface=Add_Tri_File<float>(argv[1]);
    if (argc > 4)
        triangulated_surface->Rescale(atof(argv[5]));

    cerr<<endl<<"Starting point repulsion..."<<endl;
    POINT_REPULSION<float>* point_repulsion;
    
    if(atoi(argv[4])>0){
        point_repulsion=new POINT_REPULSION<float>(*triangulated_surface,abs(atoi(argv[4])));
        point_repulsion->Seed_Points();
    }
    else{
        point_repulsion=new POINT_REPULSION<float>(*triangulated_surface,100000);
        cerr<<"Doing recursive subdivide"<<endl;
        point_repulsion->Seed_Subdivided_Points(abs(atoi(argv[4])));
    }
    cerr<<endl<<"Seeding "<<point_repulsion->number_points<<" points..."<<endl;
    
    cerr<<endl<<"Moving points for "<<atoi(argv[3])<<" iterations..."<<endl;
    for(int i=0;i<atoi(argv[3]);++i){
        point_repulsion->Move_Points();
        cerr<<i+1<<" ";}
    cerr << endl;

    cerr<<endl<<"Writing data for "<<point_repulsion->points.m<<" points to '"<<argv[2]<<"'..."<<endl;
    std::fstream os;
    os.open(argv[2], std::ios::out|std::ios::binary);
//     assert(os.is_open());
    point_repulsion->template Write_Point_Data<float>(os);
    os.close();

    if((argc != 7) || !atoi(argv[6])) return 0;
#if WITH_VISUALIZER
    myUniverse displayUniverse;
    displayUniverse.triangulated_surface=triangulated_surface;
    displayUniverse.point_repulsion=point_repulsion;
    displayUniverse.Create();
    Fl::run();
#endif    
    return 0;
}

template<class T> TRIANGULATED_SURFACE<T> *Add_Tri_File(const char* filename)
{
    ifstream input(filename,ios::binary);
    if(!input) {cerr<<"Could not open "<<filename<<endl;exit(1);}
    TRIANGULATED_SURFACE<T> *surface=new TRIANGULATED_SURFACE<T>(*new TRIANGLE_MESH,*new SOLIDS_PARTICLES<T,VECTOR_3D<T> >);
    surface->template Read<T>(input);
    std::cerr << "Number of Vertices: " << surface->triangle_mesh.number_nodes << std::endl;
    return surface;
}
