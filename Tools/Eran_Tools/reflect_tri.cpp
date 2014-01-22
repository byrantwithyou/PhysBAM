#include <Tools/Matrices/MATRIX_4X4.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
typedef float T;

MATRIX<T,4> Reflect_About_Plane(const PLANE<T>& plane)
{
    MATRIX<T,3> reflection=MATRIX<T,3>::Identity_Matrix()-(T)2*MATRIX<T,3>::Outer_Product(plane.normal,plane.normal);
    VECTOR<T,3> translation=(T)2*plane.x0.Projected_On_Unit_Direction(plane.normal);
    return MATRIX<T,4>::Translation_Matrix(translation)*MATRIX<T,4>(reflection);
}

int main(int argc,char*argv[])
{
    if(argc<2) return 1;
    std::string filename=argv[1];
    TRIANGULATED_SURFACE<T>* tri_surf=TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File<T>(filename,*tri_surf);

    std::string rgd_filename=FILE_UTILITIES::Get_Basename(filename)+".rgd";
    if(FILE_UTILITIES::File_Exists(rgd_filename)){
        std::cout << "Putting in world frame using rgd file " << rgd_filename << std::endl;
        RIGID_BODY_3D<T> rigid_body;FILE_UTILITIES::Read_From_File<T>(rgd_filename,rigid_body);
        for(int i=0;i<tri_surf->particles.number;i++) tri_surf->particles.X(i)=rigid_body.frame*tri_surf->particles.X(i);    }

#if 0
    VECTOR<T,3> x0=(T).5*(VECTOR<T,3>(0.298154,0.147101,0.979182)+VECTOR<T,3>(0.298926,0.144395,0.979307))+VECTOR<T,3>(0.002,0,0);
    VECTOR<T,3> normal=VECTOR<T,3>(1,0.01,-0.03).Normalized();
    MATRIX<T,4> xform=Reflect_About_Plane(PLANE<T>(normal,x0));
#endif
#if 0
    VECTOR<T,3> x0=(T).5*(VECTOR<T,3>(0.298154,0.147101,0.979182)+VECTOR<T,3>(0.298926,0.144395,0.979307))+VECTOR<T,3>(0.002,0,0);
    VECTOR<T,3> normal=VECTOR<T,3>(1,-0.02,-0.01).Normalized();
    MATRIX<T,4> xform=Reflect_About_Plane(PLANE<T>(normal,x0));
#endif
    MATRIX<T,4> xform(-0.999,0.03998,0.01999,0,0.03998,0.9992,-0.0003998,0,0.01999,-0.0003998,0.9998,0,0.575378,-0.0115075,-0.00575377,1);

    for(int i=0;i<tri_surf->particles.number;i++) tri_surf->particles.X(i)=xform*tri_surf->particles.X(i);    
    // invert orientation
    for(int t=0;t<tri_surf->triangle_mesh.triangles.m;t++){
        int i,j,k;tri_surf->triangle_mesh.triangles.Get(t,i,j,k);tri_surf->triangle_mesh.triangles.Set(t,i,k,j);}
    FILE_UTILITIES::Write_To_File<T>(FILE_UTILITIES::Get_Basename(filename)+"_reflected.tri",*tri_surf);

    std::cout << "FINAL TRANSFORMATION:"<<std::endl;
    std::cout << MATRIX<T,4>::Transpose(xform)<<std::endl;
}
