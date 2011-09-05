#ifndef __TRI_REFLECTOR__
#define __TRI_REFLECTOR__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include "../../Public_Library/OpenGL/OPENGL_TETS.h"
#include "OPENGL_ORIENTED_BOX_3D.h"


#include <fstream>
#include <string>

using namespace PhysBAM;
using namespace std;
using namespace FILE_UTILITIES;

template<class T> class TRI_REFLECTOR;
//#################################################################
// Class OPENGL_CALLBACK_TOGGLE_ACTIVE_TRI
//#################################################################
template<class T>
class OPENGL_CALLBACK_TOGGLE_ACTIVE_TRI: public OPENGL_CALLBACK
{
public:
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& surfaces;
    int* active_surface;
    bool* highlight_active;
    MATRIX_4X4<T>*transform;
    OPENGL_CALLBACK_TOGGLE_ACTIVE_TRI(MATRIX_4X4<T>*transform_input,bool* highlight_active_input,ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& surfaces_input,int* active_surface_input):
    transform(transform_input),highlight_active(highlight_active_input),surfaces(surfaces_input),active_surface(active_surface_input)
    {
        if(surfaces.m>0) surfaces(*active_surface)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.5),float(.5))));
    }
    void operator()()
    {
        if(surfaces.m<1) return;
        if(*active_surface+1<=surfaces.m) *active_surface+=1;
        else *active_surface=1;
        for(int i=1;i<=surfaces.m;i++)
        {
            surfaces(i)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9), float(.9), float(.9))));
        }
        if (*highlight_active)
            surfaces(*active_surface)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.5),float(.5))));
        (*transform) = MATRIX_4X4<T>::Identity_Matrix();
    }
    void Print(std::ostream& out) { out << "Toggle active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_TOGGLE_HIGHLIGHT_ACTIVE
//#################################################################
template<class T>
class OPENGL_CALLBACK_TOGGLE_HIGHLIGHT_ACTIVE: public OPENGL_CALLBACK
{
public:
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& surfaces;
    int* active_surface;
    bool* highlight_active;
    OPENGL_CALLBACK_TOGGLE_HIGHLIGHT_ACTIVE(bool* highlight_active_input, ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& surfaces_input,int* active_surface_input):
    highlight_active(highlight_active_input),surfaces(surfaces_input),active_surface(active_surface_input)
    {
        if(surfaces.m>0&&*highlight_active) surfaces(*active_surface)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.5),float(.5))));
    }
    void operator()()
    {
        *highlight_active = !(*highlight_active);
        if(surfaces.m<1) return;
        if(*highlight_active)
            surfaces(*active_surface)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.5),float(.5))));
        else
            surfaces(*active_surface)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.9),float(.9))));
    }
    void Print(std::ostream& out) { out << "Toggle highlight active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_REFLECT_TRI
//#################################################################
template<class T>
class OPENGL_CALLBACK_REFLECT_TRI: public OPENGL_CALLBACK
{
public:
    ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces;
    ARRAY<VECTOR_3D<T>*>& centers;
    int* active_surface;
    MATRIX_4X4<T>*transform;
    OPENGL_CALLBACK_REFLECT_TRI(MATRIX_4X4<T>*transform_input,ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces_input,ARRAY<VECTOR_3D<T>*>& centers_input,int* active_surface_input):
    transform(transform_input),surfaces(surfaces_input),centers(centers_input),active_surface(active_surface_input)
    {}
    void operator()()
    {
        if(surfaces.m<1) return;
        TRIANGULATED_SURFACE<T>*surface = surfaces(*active_surface);
        MATRIX_4X4<T> cur_transform = MATRIX_4X4<T>::Translation_Matrix(-(*centers(*active_surface)));
        cur_transform = MATRIX_4X4<T>(-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1)*cur_transform;
        cur_transform = MATRIX_4X4<T>::Translation_Matrix(*centers(*active_surface))*cur_transform;
        for(int i=1;i<=surface->particles.number;i++)
        {
            surface->particles.X(i) = cur_transform * surface->particles.X(i);
        }
        for(int j=1;j<=surface->triangle_mesh.triangles.m;j++)
        {
            int temp = surface->triangle_mesh.triangles(1,j);
            surface->triangle_mesh.triangles(1,j) = surface->triangle_mesh.triangles(3,j);
            surface->triangle_mesh.triangles(3,j) = temp;
        }
        (*transform) = cur_transform * (*transform);
    }
    void Print(std::ostream& out) { out << "Reflect active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_TRANSLATE_TRI
//#################################################################
template<class T>
class OPENGL_CALLBACK_TRANSLATE_TRI: public OPENGL_CALLBACK
{
public:
    ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces;
    ARRAY<VECTOR_3D<T>*>& centers;
    int* active_surface;
    MATRIX_4X4<T>*transform;
    int orientation; //1=+x, 2=+y, 3=+z, 4=-x, 5=-y, 6=-z
    OPENGL_CALLBACK_TRANSLATE_TRI(MATRIX_4X4<T>*transform_input,ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces_input,ARRAY<VECTOR_3D<T>*>& centers_input,int* active_surface_input, int orientation_input):
    transform(transform_input),surfaces(surfaces_input),centers(centers_input),active_surface(active_surface_input),orientation(orientation_input)
    {}
    void operator()()
    {
        int i;
        if(surfaces.m<1) return;
        TRIANGULATED_SURFACE<T>*surface = surfaces(*active_surface);
        MATRIX_4X4<T> cur_transform;
        switch (orientation)
        {
        case 1:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(.0015, 0, 0));
            break;
        case 2:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, .0015, 0));
            break;
        case 3:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, 0, .0015));
            break;
        case 4:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(-.0015, 0, 0));
            break;
        case 5:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, -.0015, 0));
            break;
        case 6:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, 0, -.0015));
            break;
        default:
            cout << "Error; illegal translation orientation specified." << endl;
            return;
        }
        *(centers(*active_surface)) = VECTOR_3D<T>(0, 0, 0);
        for(i=1;i<=surface->particles.number;i++)
        {
            surface->particles.X(i) = cur_transform * surface->particles.X(i);
            *(centers(*active_surface)) += surface->particles.X(i);
        }
        *(centers(*active_surface)) /= surface->particles.number;
        (*transform) = cur_transform * (*transform);
    }
    void Print(std::ostream& out) { out << "Reflect active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_ROTATE_TRI
//#################################################################
template<class T>
class OPENGL_CALLBACK_ROTATE_TRI: public OPENGL_CALLBACK
{
public:
    ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces;
    ARRAY<VECTOR_3D<T>*>& centers;
    int* active_surface;
    int orientation; //1=+x, 2=+y, 3=+z, 4=-x, 5=-y, 6=-z
    MATRIX_4X4<T>*transform;
    OPENGL_CALLBACK_ROTATE_TRI(MATRIX_4X4<T>*transform_input,ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces_input,ARRAY<VECTOR_3D<T>*>& centers_input,int* active_surface_input, int orientation_input):
    transform(transform_input),surfaces(surfaces_input),centers(centers_input),active_surface(active_surface_input),orientation(orientation_input)
    {}
    void operator()()
    {
        int i;
        VECTOR_3D<T> temp, temp2;
        T theta;
        if(surfaces.m<1) return;
        TRIANGULATED_SURFACE<T>*surface = surfaces(*active_surface);
        VECTOR_3D<T>*center = centers(*active_surface);
        MATRIX_4X4<T> cur_transform;
        theta = 0.005;
        if (orientation >= 4)
            theta = -theta;
        cur_transform = MATRIX_4X4<T>::Translation_Matrix(-(*center));
        switch (orientation)
        {
        case 1:
        case 4:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_X_Axis(theta) * cur_transform;
            break;
        case 2:
        case 5:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_Y_Axis(theta) * cur_transform;
            break;
        case 3:
        case 6:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_Z_Axis(theta) * cur_transform;
            break;
        default:
            cout << "Error; illegal translation orientation specified." << endl;
            return;
        }
        cur_transform = MATRIX_4X4<T>::Translation_Matrix(*center) * cur_transform;
        for(i=1;i<=surface->particles.number;i++)
            surface->particles.X(i) = cur_transform * surface->particles.X(i);
        (*transform) = cur_transform * (*transform);
    }
    void Print(std::ostream& out) { out << "Reflect active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_COMBINE
//#################################################################
template<class T>
class OPENGL_CALLBACK_COMBINE: public OPENGL_CALLBACK
{
public:
    ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces;
    OPENGL_CALLBACK_COMBINE(ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces_input):
    surfaces(surfaces_input)
    {}
    void operator()()
    {
        int i, j, total_particles, total_triangles;
        total_particles = total_triangles = 0;
        for (i=1; i<=surfaces.m; i++) {
            total_particles += surfaces(i)->particles.number;
            total_triangles += surfaces(i)->triangle_mesh.triangles.m;
        }
        ARRAYS<int> *triangle_list = new ARRAYS<int>(3,total_triangles);
        SOLIDS_PARTICLES<T, VECTOR_3D<T> > *particles = new SOLIDS_PARTICLES<T, VECTOR_3D<T> >();
        particles->Add_Particles(total_particles);
        int cur_particles, cur_triangles;
        cur_particles = cur_triangles = 0;
        for (i=1; i<=surfaces.m; i++) {
            ARRAYS<int> *cur_list = &(surfaces(i)->triangle_mesh.triangles);
            for (j=1; j<=cur_list->m; j++)
                triangle_list->Set(cur_triangles + j,(*cur_list)(1,j)+cur_particles,(*cur_list)(2,j)+cur_particles,(*cur_list)(3,j)+cur_particles);
            for (j=1; j<=surfaces(i)->particles.number; j++)
                particles->Copy_Particle(surfaces(i)->particles, j, cur_particles + j);
            cur_particles += surfaces(i)->particles.number;
            cur_triangles += surfaces(i)->triangle_mesh.triangles.m;
        }
        TRIANGULATED_SURFACE<T> *surface = new TRIANGULATED_SURFACE<T>(*(new TRIANGLE_MESH(*triangle_list)),*particles);
        string filename = "combined.tri";
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename,true,false);
        if(!output) {cerr<<"Could not open "<<filename<<endl;return;}
        surface->template Write<T>(*output);delete output;
    }
    void Print(std::ostream& out) { out << "Reflect active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_WRITE_TRANSFORM
//#################################################################
template<class T>
class OPENGL_CALLBACK_WRITE_TRANSFORM: public OPENGL_CALLBACK
{
public:
    MATRIX_4X4<T> *transform;
    OPENGL_CALLBACK_WRITE_TRANSFORM(MATRIX_4X4<T>*transform_input):
    transform(transform_input)
    {}
    void operator()()
    {
        string filename = "VH_Bones/ribs_reflect_transform_new.txt";
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename,false,false);
        if(!output) {cerr<<"Could not open "<<filename<<endl;return;}
        (*output) << transform->x[0] << " " << transform->x[1] << " " << transform->x[2] << " " << transform->x[3] << " ";
        (*output) << transform->x[4] << " " << transform->x[5] << " " << transform->x[6] << " " << transform->x[7] << " ";
        (*output) << transform->x[8] << " " << transform->x[9] << " " << transform->x[10] << " " << transform->x[11] << " ";
        (*output) << transform->x[12] << " " << transform->x[13] << " " << transform->x[14] << " " << transform->x[15];
        delete output;        
    }
    void Print(std::ostream& out) { out << "Reflect active triangulated surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_READ_TRANSFORM
//#################################################################
template<class T>
class OPENGL_CALLBACK_READ_TRANSFORM: public OPENGL_CALLBACK
{
public:
    ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces, &original_surfaces;
    ARRAY<VECTOR_3D<T>*>& centers;
    int* active_surface;
    MATRIX_4X4<T>*transform;
    OPENGL_CALLBACK_READ_TRANSFORM(MATRIX_4X4<T>*transform_input,ARRAY<TRIANGULATED_SURFACE<T>*>& surfaces_input,ARRAY<VECTOR_3D<T>*>& centers_input,int* active_surface_input,
        ARRAY<TRIANGULATED_SURFACE<T>*>& original_surfaces_input):
    transform(transform_input),surfaces(surfaces_input),centers(centers_input),active_surface(active_surface_input),original_surfaces(original_surfaces_input)
    {}
    void operator()()
    {
        string filename = "VH_Bones/ribs_reflect_transform.txt";
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false,false);
        if(!input) {cerr<<"Could not open "<<filename<<endl;return;}
        (*input) >> transform->x[0] >> transform->x[1] >> transform->x[2] >> transform->x[3];
        (*input) >> transform->x[4] >> transform->x[5] >> transform->x[6] >> transform->x[7];
        (*input) >> transform->x[8] >> transform->x[9] >> transform->x[10] >> transform->x[11];
        (*input) >> transform->x[12] >> transform->x[13] >> transform->x[14] >> transform->x[15];
        delete input;        

        if(surfaces.m<1) return;
        TRIANGULATED_SURFACE<T>* surface, *original_surface;
        surface = surfaces(*active_surface);
        original_surface = original_surfaces(*active_surface);
        for(int i=1;i<=surface->particles.number;i++)
        {
            surface->particles.X(i) = (*transform) * original_surface->particles.X(i);
        }
        for(int j=1;j<=surface->triangle_mesh.triangles.m;j++)
        {
            surface->triangle_mesh.triangles(1,j) = original_surface->triangle_mesh.triangles(3,j);
            surface->triangle_mesh.triangles(3,j) = original_surface->triangle_mesh.triangles(1,j);
        }
    }
    void Print(std::ostream& out) { out << "Read and apply a transform from the file 'transform.txt'"; }
};
//#################################################################
// Class OPENGL_CALLBACK_SHOW_BONES
//#################################################################
template<class T>
class OPENGL_CALLBACK_SHOW_BONES: public OPENGL_CALLBACK
{
public:
    OPENGL_WORLD& world;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& gl_surfaces;
    OPENGL_CALLBACK_SHOW_BONES(OPENGL_WORLD& world_input,ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*>& gl_surfaces_input):
    world(world_input),gl_surfaces(gl_surfaces_input)
    {}
    void operator()()
    {
        for(int i=1;i<=gl_surfaces.m;i++)
            gl_surfaces(i)->visible = !(gl_surfaces(i)->visible);
    }
    void Print(std::ostream& out) { out << "Read and apply a transform from the file 'transform.txt'"; }
};
//#################################################################
// Class OPENGL_CALLBACK_SHOW_MUSCLE
//#################################################################
template<class T>
class OPENGL_CALLBACK_SHOW_MUSCLE: public OPENGL_CALLBACK
{
public:
    OPENGL_WORLD& world;
    ARRAY<OPENGL_TETS<T>*>& gl_vols;
    OPENGL_CALLBACK_SHOW_MUSCLE(OPENGL_WORLD& world_input,ARRAY<OPENGL_TETS<T>*>& gl_vols_input):
    world(world_input),gl_vols(gl_vols_input)
    {}
    void operator()()
    {
        for(int i=1;i<=gl_vols.m;i++)
            gl_vols(i)->visible = !(gl_vols(i)->visible);
    }
    void Print(std::ostream& out) { out << "Read and apply a transform from the file 'transform.txt'"; }
};
//#################################################################
// Class OPENGL_CALLBACK_TOGGLE_ACTIVE_ISURF
//#################################################################
template<class T>
class OPENGL_CALLBACK_TOGGLE_ACTIVE_ISURF: public OPENGL_CALLBACK
{
public:
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_isurfs;
    int* active_isurf;
    OPENGL_CALLBACK_TOGGLE_ACTIVE_ISURF(ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_isurfs_input,int* active_isurf_input):
    gl_isurfs(gl_isurfs_input),active_isurf(active_isurf_input)
    {        
        if(gl_isurfs.m>0) gl_isurfs(*active_isurf)->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.9),float(.7),float(.5))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
    }
    void operator()()
    {
        if(gl_isurfs.m<1) return;
        if(*active_isurf+1<=gl_isurfs.m) *active_isurf+=1;
        else *active_isurf=1;
        for(int i=1;i<=gl_isurfs.m;i++)
        {
            gl_isurfs(i)->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(.5),float(.3))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
        }
        gl_isurfs(*active_isurf)->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.9),float(.7),float(.5))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
    }
    void Print(std::ostream& out) { out << "Toggle active implicit surface"; }
};
//#################################################################
// Class OPENGL_CALLBACK_TRANSLATE_BOX
//#################################################################
template<class T>
class OPENGL_CALLBACK_TRANSLATE_BOX: public OPENGL_CALLBACK
{
public:
    ORIENTED_BOX_3D<T> *clip_box;
    int orientation; //1=+x, 2=+y, 3=+z, 4=-x, 5=-y, 6=-z
    OPENGL_CALLBACK_TRANSLATE_BOX(ORIENTED_BOX_3D<T> *clip_box_input, int orientation_input):
    clip_box(clip_box_input),orientation(orientation_input)
    {}
    void operator()()
    {
        MATRIX_4X4<T> cur_transform;
        switch (orientation)
        {
        case 1:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(.01, 0, 0));
            break;
        case 2:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, .01, 0));
            break;
        case 3:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, 0, .01));
            break;
        case 4:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(-.01, 0, 0));
            break;
        case 5:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, -.01, 0));
            break;
        case 6:
            cur_transform = MATRIX_4X4<T>::Translation_Matrix(VECTOR_3D<T>(0, 0, -.01));
            break;
        default:
            cout << "Error; illegal translation orientation specified." << endl;
            return;
        }
        clip_box->corner = cur_transform * clip_box->corner;
    }
    void Print(std::ostream& out) { out << "Translate clip box with orientation " << orientation; }
};
//#################################################################
// Class OPENGL_CALLBACK_ROTATE_BOX
//#################################################################
template<class T>
class OPENGL_CALLBACK_ROTATE_BOX: public OPENGL_CALLBACK
{
public:
    ORIENTED_BOX_3D<T> *clip_box;
    int orientation; //1=+x, 2=+y, 3=+z, 4=-x, 5=-y, 6=-z
    OPENGL_CALLBACK_ROTATE_BOX(ORIENTED_BOX_3D<T> *clip_box_input, int orientation_input):
    clip_box(clip_box_input),orientation(orientation_input)
    {}
    void operator()()
    {
        MATRIX_4X4<T> cur_transform;
        T theta = 0.01;
        if (orientation >= 4)
            theta = -theta;
        switch (orientation)
        {
        case 1:
        case 4:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_X_Axis(theta);
            break;
        case 2:
        case 5:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_Y_Axis(theta);
            break;
        case 3:
        case 6:
            cur_transform = MATRIX_4X4<T>::Rotation_Matrix_Z_Axis(theta);
            break;
        default:
            cout << "Error; illegal translation orientation specified." << endl;
            return;
        }
        clip_box->edge1 = cur_transform * clip_box->edge1;
        clip_box->edge2 = cur_transform * clip_box->edge2;
        clip_box->edge3 = cur_transform * clip_box->edge3;

    }
    void Print(std::ostream& out) { out << "Rotate clip box with orientation " << orientation; }
};
//#################################################################
// Class OPENGL_CALLBACK_RESIZE_BOX
//#################################################################
template<class T>
class OPENGL_CALLBACK_RESIZE_BOX: public OPENGL_CALLBACK
{
public:
    ORIENTED_BOX_3D<T> *clip_box;
    int orientation; //1=grow, 2=shrink
    OPENGL_CALLBACK_RESIZE_BOX(ORIENTED_BOX_3D<T> *clip_box_input, int orientation_input):
    clip_box(clip_box_input),orientation(orientation_input)
    {}
    void operator()()
    {
        T scale;
        switch (orientation)
        {
        case 1:
            scale = 2;
            break;
        case 2:
            scale = 0.5;
            break;
        default:
            cout << "Error; illegal translation orientation specified." << endl;
            return;
        }
        clip_box->edge1*=scale;
        clip_box->edge2*=scale;
        clip_box->edge3*=scale;
    }
    void Print(std::ostream& out) { out << "Scale clip box with orientation " << orientation; }
};
//#################################################################
// Class OPENGL_CALLBACK_CLIP_SKIN
//#################################################################
template<class T>
class OPENGL_CALLBACK_CLIP_SKIN: public OPENGL_CALLBACK
{
public:
    TRI_REFLECTOR<T> *tri_reflector;
    ORIENTED_BOX_3D<T> *clip_box;
    ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*>& isurfs;
    int *active_isurf;
    OPENGL_CALLBACK_CLIP_SKIN(TRI_REFLECTOR<T> *tri_reflector_input,ORIENTED_BOX_3D<T>* clip_box_input,ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*>& isurfs_input,int *active_isurf_input):
    tri_reflector(tri_reflector_input),clip_box(clip_box_input),isurfs(isurfs_input),active_isurf(active_isurf_input)
    {}
    void operator()()
    {
        if (isurfs.m < 1) return;
        LEVELSET_IMPLICIT_SURFACE<T> *isurf = isurfs(*active_isurf);
        if (!ORIENTED_BOX_3D<T>(isurf->box,QUATERNION<T>()).Intersection(*clip_box)) return;
        VECTOR_3D<T> location, local_location;
        T phi_value = isurf->levelset.grid.dx/2;
        for (int i = 1; i <= isurf->levelset.grid.m; i++)
            for (int j = 1; j <= isurf->levelset.grid.n; j++)
                for (int k = 1; k <= isurf->levelset.grid.mn; k++)
                {
                    location = isurf->levelset.grid.Node(i,j,k);
                    local_location = location - clip_box->corner;
                    if (VECTOR_3D<T>::Dot_Product(local_location,clip_box->edge1) < 0 ||
                        VECTOR_3D<T>::Dot_Product(local_location,clip_box->edge2) < 0 ||
                        VECTOR_3D<T>::Dot_Product(local_location,clip_box->edge3) < 0)
                    {
                        isurf->levelset.phi(i,j,k) = (isurf->levelset.phi(i,j,k) > 0) ? phi_value : -phi_value;
                        continue;
                    }
                    if (local_location.Projected(clip_box->edge1).Magnitude_Squared() > clip_box->edge1.Magnitude_Squared() ||
                        local_location.Projected(clip_box->edge2).Magnitude_Squared() > clip_box->edge2.Magnitude_Squared() ||
                        local_location.Projected(clip_box->edge3).Magnitude_Squared() > clip_box->edge3.Magnitude_Squared())
                    {
                        isurf->levelset.phi(i,j,k) = (isurf->levelset.phi(i,j,k) > 0) ? phi_value : -phi_value;
                        continue;
                    }
                    if (isurf->levelset.phi(i,j,k) < 0)
                        isurf->levelset.phi(i,j,k) = phi_value;
                }
        tri_reflector->gl_isurfs(*active_isurf)->Set_Levelset(isurf->levelset);
        tri_reflector->gl_isurfs(*active_isurf)->Generate_Triangulated_Surface(false,"");
        tri_reflector->gl_isurfs(*active_isurf)->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.9),float(.7),float(.5))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
        LEVELSET_IMPLICIT_SURFACE<T> *reflected_isurf = tri_reflector->Create_Reflected_Isurf(isurf);
        tri_reflector->reflected_isurfs.Remove_Index(*active_isurf);
        tri_reflector->reflected_isurfs.Insert(reflected_isurf,*active_isurf);
        tri_reflector->gl_reflected_isurfs(*active_isurf)->Set_Levelset(reflected_isurf->levelset);
        tri_reflector->gl_reflected_isurfs(*active_isurf)->Generate_Triangulated_Surface(false,"");
        tri_reflector->gl_reflected_isurfs(*active_isurf)->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(.5),float(.3))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
    }
    void Print(std::ostream& out) { out << "Write out the skin file 'new_skin_write.phi'"; }
};
//#################################################################
// Class OPENGL_CALLBACK_WRITE_SKIN
//#################################################################
template<class T>
class OPENGL_CALLBACK_WRITE_SKIN: public OPENGL_CALLBACK
{
public:
    ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*>& isurfs;
    int *active_isurf;
    OPENGL_CALLBACK_WRITE_SKIN(ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*>& isurfs_input,int *active_isurf_input):
    isurfs(isurfs_input),active_isurf(active_isurf_input)
    {}
    void operator()()
    {
        if (isurfs.m < 1) return;
        string filename = "new_skin_write.phi";
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename,true,false);
        if(!output) {cerr<<"Could not open "<<filename<<endl;return;}
        isurfs(*active_isurf)->template Write<T>(*output);delete output;
    }
    void Print(std::ostream& out) { out << "Write out the skin file 'new_skin_write.phi'"; }
};
//#################################################################
// Class OPENGL_CALLBACK_SHOW_SKIN
//#################################################################
template<class T>
class OPENGL_CALLBACK_SHOW_SKIN: public OPENGL_CALLBACK
{
public:
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_isurfs;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_reflected_isurfs;
    OPENGL_CALLBACK_SHOW_SKIN(ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_isurfs_input,ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*>& gl_reflected_isurfs_input):
    gl_isurfs(gl_isurfs_input),gl_reflected_isurfs(gl_reflected_isurfs_input)
    {}
    void operator()()
    {
        for(int i=1;i<=gl_isurfs.m;i++)
        {
            gl_isurfs(i)->visible = !(gl_isurfs(i)->visible);
            gl_reflected_isurfs(i)->visible = !(gl_reflected_isurfs(i)->visible);
        }
    }
    void Print(std::ostream& out) { out << "Read and apply a transform from the file 'transform.txt'"; }
};
//#################################################################
// Class TRI_REFLECTOR
//#################################################################
template<class T>
class TRI_REFLECTOR
{
public:
    OPENGL_WORLD& world;
    ARRAY<TRIANGULATED_SURFACE<T>*> surfaces;
    ARRAY<TRIANGULATED_SURFACE<T>*> original_surfaces;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> gl_surfaces;
    ARRAY<VECTOR_3D<T>*> surface_centers;
    ORIENTED_BOX_3D<T> clip_box;
    OPENGL_ORIENTED_BOX_3D<T> *gl_clip_box;
    int active_surface;
    bool highlight_active;
    MATRIX_4X4<T> transform, start_transform;
    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> vols;
    ARRAY<OPENGL_TETS<T>*> gl_vols;
    ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*> isurfs;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*> gl_isurfs;
    ARRAY<LEVELSET_IMPLICIT_SURFACE<T>*> reflected_isurfs;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*> gl_reflected_isurfs;
    int active_isurf;

    
    TRI_REFLECTOR(OPENGL_WORLD& world_input):world(world_input)
    {
        active_surface=active_isurf=0;
        highlight_active=true;
        transform=MATRIX_4X4<T>::Identity_Matrix();
        clip_box = ORIENTED_BOX_3D<T>();
        gl_clip_box = new OPENGL_ORIENTED_BOX_3D<T>(&clip_box);
        world.Add_Object(gl_clip_box);

        //read transform
        string transform_file = "VH_Bones/ribs_reflect_transform.txt";
        std::istream* input2=FILE_UTILITIES::Safe_Open_Input(transform_file,false,false);
        if(!input2) {cerr<<"Could not open "<<transform_file<<endl;return;}
        (*input2) >> start_transform.x[0] >> start_transform.x[1] >> start_transform.x[2] >> start_transform.x[3];
        (*input2) >> start_transform.x[4] >> start_transform.x[5] >> start_transform.x[6] >> start_transform.x[7];
        (*input2) >> start_transform.x[8] >> start_transform.x[9] >> start_transform.x[10] >> start_transform.x[11];
        (*input2) >> start_transform.x[12] >> start_transform.x[13] >> start_transform.x[14] >> start_transform.x[15];
        delete input2;
    }

    void Initialize_Tri_Reflector_Callbacks(int argc,char** argv){
        //parse args
        for(int i=1; i<argc; i++){
            int length=(int)strlen(argv[i]);
            if(length>=4 && strcmp(argv[i]+length-4,".tri")==0) Add_Tri_File(argv[i],world,i);
            else if(length>=4 && strcmp(argv[i]+length-4,".tet")==0) Add_Tet_File(argv[i],world,i);
            else if(length>=4 && strcmp(argv[i]+length-4,".phi")==0) Add_Phi_File(argv[i],world,i);
            else cerr<<"Unrecognized file "<<argv[i]<<endl;
        }
        //initilize callbacks
        active_surface=active_isurf=1;
        cout << "clip box corner initialized at " << clip_box.corner << endl;
        if (isurfs.m > 0)
        {
            cout << "setting clip box to " << isurfs(active_isurf)->box.Center() << endl;
            clip_box.corner = isurfs(active_isurf)->box.Center();
        }
        cout << "clip box corner now at " << clip_box.corner << endl;
        //triangulated surface callbacks
        OPENGL_CALLBACK_TOGGLE_ACTIVE_TRI<T>* surface_selector= new OPENGL_CALLBACK_TOGGLE_ACTIVE_TRI<T>(&transform,&highlight_active,gl_surfaces,&active_surface);
        OPENGL_CALLBACK_TOGGLE_HIGHLIGHT_ACTIVE<T>* highlight_active_toggler= new OPENGL_CALLBACK_TOGGLE_HIGHLIGHT_ACTIVE<T>(&highlight_active,gl_surfaces,&active_surface);
         OPENGL_CALLBACK_REFLECT_TRI<T>* surface_reflector= new OPENGL_CALLBACK_REFLECT_TRI<T>(&transform, surfaces, surface_centers, &active_surface);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_px= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 1);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_py= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 2);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_pz= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 3);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_nx= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 4);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_ny= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 5);
         OPENGL_CALLBACK_TRANSLATE_TRI<T>* surface_translator_nz= new OPENGL_CALLBACK_TRANSLATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 6);
        OPENGL_CALLBACK_COMBINE<T>* surface_combiner= new OPENGL_CALLBACK_COMBINE<T>(surfaces);
        OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_px= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 1);
         OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_py= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 2);
         OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_pz= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 3);
         OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_nx= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 4);
         OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_ny= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 5);
         OPENGL_CALLBACK_ROTATE_TRI<T>* surface_rotator_nz= new OPENGL_CALLBACK_ROTATE_TRI<T>(&transform, surfaces, surface_centers, &active_surface, 6);
        OPENGL_CALLBACK_WRITE_TRANSFORM<T>* transform_writer = new OPENGL_CALLBACK_WRITE_TRANSFORM<T>(&transform);
        OPENGL_CALLBACK_READ_TRANSFORM<T>* transform_reader = new OPENGL_CALLBACK_READ_TRANSFORM<T>(&transform, surfaces, surface_centers, &active_surface, original_surfaces);
        OPENGL_CALLBACK_SHOW_BONES<T>* bone_shower = new OPENGL_CALLBACK_SHOW_BONES<T>(world, gl_surfaces);
        //tetrahedralized volume callbacks
        OPENGL_CALLBACK_SHOW_MUSCLE<T>* muscle_shower = new OPENGL_CALLBACK_SHOW_MUSCLE<T>(world, gl_vols);
        //levelset implicit surface callbacks
        OPENGL_CALLBACK_TOGGLE_ACTIVE_ISURF<T>* isurf_selector= new OPENGL_CALLBACK_TOGGLE_ACTIVE_ISURF<T>(gl_isurfs,&active_isurf);
        OPENGL_CALLBACK_SHOW_SKIN<T>* skin_shower = new OPENGL_CALLBACK_SHOW_SKIN<T>(gl_isurfs, gl_reflected_isurfs);
        OPENGL_CALLBACK_WRITE_SKIN<T>* skin_writer= new OPENGL_CALLBACK_WRITE_SKIN<T>(isurfs,&active_isurf);
        OPENGL_CALLBACK_WRITE_SKIN<T>* reflected_skin_writer= new OPENGL_CALLBACK_WRITE_SKIN<T>(reflected_isurfs,&active_isurf);
         OPENGL_CALLBACK_CLIP_SKIN<T>* skin_clipper= new OPENGL_CALLBACK_CLIP_SKIN<T>(this, &clip_box, isurfs, &active_isurf);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_px= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 1);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_py= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 2);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_pz= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 3);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_nx= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 4);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_ny= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 5);
         OPENGL_CALLBACK_TRANSLATE_BOX<T>* box_translator_nz= new OPENGL_CALLBACK_TRANSLATE_BOX<T>(&clip_box, 6);
        OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_px= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 1);
         OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_py= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 2);
         OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_pz= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 3);
         OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_nx= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 4);
         OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_ny= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 5);
         OPENGL_CALLBACK_ROTATE_BOX<T>* box_rotator_nz= new OPENGL_CALLBACK_ROTATE_BOX<T>(&clip_box, 6);
         OPENGL_CALLBACK_RESIZE_BOX<T>* box_grower= new OPENGL_CALLBACK_RESIZE_BOX<T>(&clip_box, 1);
         OPENGL_CALLBACK_RESIZE_BOX<T>* box_shrinker= new OPENGL_CALLBACK_RESIZE_BOX<T>(&clip_box, 2);

        //ALL MODE CALLBACKS
        world.Bind_Key('t',surface_selector);  // Toggle Muscles
        world.Bind_Key('v',highlight_active_toggler);
        world.Bind_Key('R',surface_reflector);
        world.Bind_Key('1',surface_translator_px);
        world.Bind_Key('2',surface_translator_py);
        world.Bind_Key('3',surface_translator_pz);
        world.Bind_Key('!',surface_translator_nx);
        world.Bind_Key('@',surface_translator_ny);
        world.Bind_Key('#',surface_translator_nz);
        world.Bind_Key('c',surface_combiner);
        world.Bind_Key('4',surface_rotator_px);
        world.Bind_Key('5',surface_rotator_py);
        world.Bind_Key('6',surface_rotator_pz);
        world.Bind_Key('$',surface_rotator_nx);
        world.Bind_Key('%',surface_rotator_ny);
        world.Bind_Key('^',surface_rotator_nz);
        world.Bind_Key('w',transform_writer);
        world.Bind_Key('r',transform_reader);
        world.Bind_Key('b',bone_shower);
        
        world.Bind_Key('m',muscle_shower);

        world.Bind_Key('T',isurf_selector);
        world.Bind_Key('s',skin_shower);
        world.Bind_Key('W',skin_writer);
        world.Bind_Key('q',reflected_skin_writer);
        world.Bind_Key('C',skin_clipper);
        world.Bind_Key('[',box_translator_px);
        world.Bind_Key(']',box_translator_py);
        world.Bind_Key('\\',box_translator_pz);
        world.Bind_Key('{',box_translator_nx);
        world.Bind_Key('}',box_translator_ny);
        world.Bind_Key('|',box_translator_nz);
        world.Bind_Key(',',box_rotator_px);
        world.Bind_Key('.',box_rotator_py);
        world.Bind_Key('/',box_rotator_pz);
        world.Bind_Key('<',box_rotator_nx);
        world.Bind_Key('>',box_rotator_ny);
        world.Bind_Key('?',box_rotator_nz);
        world.Bind_Key('+',box_grower);
        world.Bind_Key('-',box_shrinker);
        world.Bind_Key('W',skin_writer);
    }

    void Add_Tri_File(const std::string& filename,OPENGL_WORLD& world,int number)
    {
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,true,false);
        if(!input) {cerr<<"Could not open "<<filename<<endl;return;}
        TRIANGULATED_SURFACE<T> *surface=new TRIANGULATED_SURFACE<T>(*new TRIANGLE_MESH,*new SOLIDS_PARTICLES<T,VECTOR_3D<T> >);
        surface->template Read<T>(*input);delete input;
        surfaces.Append(surface);
        TRIANGULATED_SURFACE<T> *surface2=new TRIANGULATED_SURFACE<T>(*new TRIANGLE_MESH(surface->triangle_mesh), *new SOLIDS_PARTICLES<T,VECTOR_3D<T> >(surface->particles));
        original_surfaces.Append(surface2);
        OPENGL_TRIANGULATED_SURFACE<T>*ots=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.9),float(.9))),OPENGL_MATERIAL::Plastic(OPENGL_COLOR::White()));
        gl_surfaces.Append(ots);
        world.Add_Object(ots,true,true);
        VECTOR_3D<T> *center = new VECTOR_3D<T>(0, 0, 0);
        for(int i=1;i<=surface->particles.number;i++)
        {
            (*center) += surface->particles.X(i);
        }
        (*center) /= surface->particles.number;
        cout << "Object " << filename.c_str() << " center is at " << center->x << "," << center->y << "," << center->z << endl;
        surface_centers.Append(center);
    }

    void Add_Tet_File(const std::string& filename,OPENGL_WORLD& world,int number)
    {
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,true,false);
        if(!input) {cerr<<"Could not open "<<filename<<endl;return;}
        TETRAHEDRALIZED_VOLUME<T> *vol=new TETRAHEDRALIZED_VOLUME<T>(*new TETRAHEDRON_MESH,*new SOLIDS_PARTICLES<T,VECTOR_3D<T> >);
        vol->template Read<T>(*input);delete input;vol->particles.Store_Velocity(false);vol->particles.Store_Mass(false);
        vols.Append(vol);
        cout << filename << " has " << vol->particles.array_size << " " << vol->particles.number << " tets." << endl;
        TETRAHEDRALIZED_VOLUME<T> *reflected_vol=new TETRAHEDRALIZED_VOLUME<T>(*new TETRAHEDRON_MESH(vol->tetrahedron_mesh), *new SOLIDS_PARTICLES<T,VECTOR_3D<T> >(vol->particles));
        cout << filename << " copy has " << reflected_vol->particles.X.array.m << " tets." << endl;
        vols.Append(reflected_vol);
        OPENGL_TETS<T>* tets=new OPENGL_TETS<T>(&(vol->tetrahedron_mesh),&(vol->particles),OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.1),float(.1))));
        world.Add_Object(tets,true,true);
        //transform reflected tets
        for(int i=1;i<=reflected_vol->particles.number;i++)reflected_vol->particles.X(i)=start_transform*vol->particles.X(i);
        for(int j=1;j<=reflected_vol->tetrahedron_mesh.tetrahedrons.m;j++){
            reflected_vol->tetrahedron_mesh.tetrahedrons(1,j)=vol->tetrahedron_mesh.tetrahedrons(3,j);
            reflected_vol->tetrahedron_mesh.tetrahedrons(3,j)=vol->tetrahedron_mesh.tetrahedrons(1,j);}
        OPENGL_TETS<T>* reflected_tets=new OPENGL_TETS<T>(&(reflected_vol->tetrahedron_mesh),&(reflected_vol->particles),OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.1),float(.1))));
        world.Add_Object(reflected_tets,true,true);
    }

    void Add_Phi_File(const std::string& filename,OPENGL_WORLD& world,int number)
    {
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,true,false);
        if(!input) {cerr<<"Could not open "<<filename<<endl;return;}
        LEVELSET_IMPLICIT_SURFACE<T> *surface=new LEVELSET_IMPLICIT_SURFACE<T>(*new GRID_3D<T>,*new ARRAYS<VECTOR<T,3> >);
        surface->template Read<T>(*input);delete input;
        cout<<"surface is defined on a "<<surface->levelset.phi.m<<" by "<<surface->levelset.phi.n<<" by "<<surface->levelset.phi.mn<<" grid"<<endl;
        int i,j,k;
        for(i=surface->levelset.phi.m_start; i<=surface->levelset.phi.m_end; ++i)
            for(j=surface->levelset.phi.n_start; j<=surface->levelset.phi.n_end; ++j)
                if(surface->levelset.phi(i,j,surface->levelset.phi.mn_start)<=0) {cout<<"phi<=0 on mn_start"<<endl;goto check2;}
        check2:
        for(i=surface->levelset.phi.m_start; i<=surface->levelset.phi.m_end; ++i)
            for(j=surface->levelset.phi.n_start; j<=surface->levelset.phi.n_end; ++j)
                if(surface->levelset.phi(i,j,surface->levelset.phi.mn_end)<=0) {cout<<"phi<=0 on mn_end"<<endl;goto check3;}
        check3:
        for(i=surface->levelset.phi.m_start; i<=surface->levelset.phi.m_end; ++i)
            for(k=surface->levelset.phi.mn_start; k<=surface->levelset.phi.mn_end; ++k)
                if(surface->levelset.phi(i,surface->levelset.phi.n_start,k)<=0) {cout<<"phi<=0 on n_start"<<endl;goto check4;}
        check4:
        for(i=surface->levelset.phi.m_start; i<=surface->levelset.phi.m_end; ++i)
            for(k=surface->levelset.phi.mn_start; k<=surface->levelset.phi.mn_end; ++k)
                if(surface->levelset.phi(i,surface->levelset.phi.n_end,k)<=0) {cout<<"phi<=0 on n_end"<<endl;goto check5;}
        check5:
        for(j=surface->levelset.phi.n_start; j<=surface->levelset.phi.n_end; ++j)
            for(k=surface->levelset.phi.mn_start; k<=surface->levelset.phi.mn_end; ++k)
                if(surface->levelset.phi(surface->levelset.phi.m_start,j,k)<=0) {cout<<"phi<=0 on m_start"<<endl;goto check6;}
        check6:
        for(j=surface->levelset.phi.n_start; j<=surface->levelset.phi.n_end; ++j)
            for(k=surface->levelset.phi.mn_start; k<=surface->levelset.phi.mn_end; ++k)
                if(surface->levelset.phi(surface->levelset.phi.m_end,j,k)<=0) {cout<<"phi<=0 on m_end"<<endl;goto end_of_checks;}
        end_of_checks:
        isurfs.Append(surface);
        OPENGL_LEVELSET_MULTIVIEW<T> *opengl_surface=new OPENGL_LEVELSET_MULTIVIEW<T>;
        opengl_surface->Set_Levelset(surface->levelset);
        opengl_surface->Generate_Triangulated_Surface(false,"");
        opengl_surface->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(.5),float(.3))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
        gl_isurfs.Append(opengl_surface);
        world.Add_Object(opengl_surface);
        LEVELSET_IMPLICIT_SURFACE<T> *reflected_isurf= Create_Reflected_Isurf(surface);
        reflected_isurfs.Append(reflected_isurf);
        OPENGL_LEVELSET_MULTIVIEW<T> *reflected_gl_isurf=new OPENGL_LEVELSET_MULTIVIEW<T>();
        reflected_gl_isurf->Set_Levelset(reflected_isurf->levelset);
        reflected_gl_isurf->Generate_Triangulated_Surface(false,"");
        reflected_gl_isurf->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(.5),float(.3))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
        gl_reflected_isurfs.Append(reflected_gl_isurf);
        world.Add_Object(reflected_gl_isurf);
    }

    LEVELSET_IMPLICIT_SURFACE<T> *Create_Reflected_Isurf(LEVELSET_IMPLICIT_SURFACE<T> *surface)
    {
        //transform grid and phi array
        bool flipx = false, flipy = false, flipz = false;
        VECTOR_3D<T> corners[8];
        corners[0] = VECTOR_3D<T>(surface->levelset.grid.xmin,surface->levelset.grid.ymin,surface->levelset.grid.zmin);
        corners[1] = VECTOR_3D<T>(surface->levelset.grid.xmin,surface->levelset.grid.ymin,surface->levelset.grid.zmax);
        corners[2] = VECTOR_3D<T>(surface->levelset.grid.xmin,surface->levelset.grid.ymax,surface->levelset.grid.zmin);
        corners[3] = VECTOR_3D<T>(surface->levelset.grid.xmin,surface->levelset.grid.ymax,surface->levelset.grid.zmax);
        corners[4] = VECTOR_3D<T>(surface->levelset.grid.xmax,surface->levelset.grid.ymin,surface->levelset.grid.zmin);
        corners[5] = VECTOR_3D<T>(surface->levelset.grid.xmax,surface->levelset.grid.ymin,surface->levelset.grid.zmax);
        corners[6] = VECTOR_3D<T>(surface->levelset.grid.xmax,surface->levelset.grid.ymax,surface->levelset.grid.zmin);
        corners[7] = VECTOR_3D<T>(surface->levelset.grid.xmax,surface->levelset.grid.ymax,surface->levelset.grid.zmax);
        VECTOR_3D<T> min_corner, max_corner, corner;
        corners[0] = min_corner = max_corner = start_transform*corners[0];
        for (int c = 1; c < 8; c++) {
            corner = start_transform*corners[c];
            for (int i = 1; i <= 3; i++) {
                if (corner[i] < min_corner[i])
                    min_corner[i] = corner[i];
                if (corner[i] > max_corner[i])
                    max_corner[i] = corner[i];
            }
        }
        BOX_3D<T>* bounding_box = new BOX_3D<T>(min_corner, max_corner);
        cout << "Reflected isurf bounding box " << *bounding_box << endl;
        GRID_3D<T> *grid=new GRID_3D<T>((bounding_box->xmax-bounding_box->xmin)/surface->levelset.grid.dx,
                                        (bounding_box->ymax-bounding_box->ymin)/surface->levelset.grid.dy,
                                        (bounding_box->zmax-bounding_box->zmin)/surface->levelset.grid.dz,*bounding_box);
        ARRAYS<VECTOR<T,3> > *phi=new ARRAYS<VECTOR<T,3> >(*grid);
        cout << "Calculating phi values for reflected implicit surface";
        MATRIX_4X4<T> inverse_transform = start_transform.Inverse();
        T phi_value = grid->dx/2;
        VECTOR_3D<T> reflected_node;
        for(int i=phi->m_start;i<=phi->m_end;i++)
            for(int j=phi->n_start;j<=phi->n_end;j++)
                for(int k=phi->mn_start;k<=phi->mn_end;k++)
                {
                    reflected_node = inverse_transform*grid->Node(i,j,k);
                    (*phi)(i,j,k) = (surface->box.Lazy_Inside(reflected_node)) ? (*surface)(reflected_node) : phi_value;
                }
        return new LEVELSET_IMPLICIT_SURFACE<T>(*grid,*phi);
    }
 };

#endif