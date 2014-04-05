//
//  mesh_cutting_test.cpp
//
//
//  Created by Yuting Wang on 5/24/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

//2. high res sim explosion: tune parameter or find bug in sim code...

#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Matrices/ROTATION.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>

#include <fstream>
#include <sstream>

#include "mesh_cutting_subd.h"
#include "DEFORMABLE_OBJECTS.h"

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <X11/Xlib.h>
#ifdef __APPLE__  // include Mac OS X verions of headers
//#  include <GL/glew.h>
#  include <OpenGL/OpenGL.h>
#  include <GLUT/glut.h>

#else // non-Mac OS X operating systems
#  include <GL/glew.h>
#  include <GL/freeglut.h>
#  include <GL/freeglut_ext.h>
#endif  // __APPLE__
namespace VS{
void read_tsc(){__asm__("rdtsc");}
struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}
}

using namespace PhysBAM;
using namespace std;

typedef double T;
typedef PhysBAM::VECTOR<T, 3> TV;
typedef PhysBAM::MATRIX<T, 3> TM;
typedef PhysBAM::VECTOR<int, 3> I3;

MESH_CUTTING<T> *mcut = NULL;
TETRAHEDRALIZED_VOLUME<T> *sim_volume = NULL;
TRIANGULATED_SURFACE<T> *cutting_tri_mesh = NULL;

T timestep = 40;
int ratio = 10;
T K = 10000;

const string dir = "/Users/yutingwang/PhysBAM";
const string volumeDir = dir + "/Tetrahedralized_Volumes";
void Initialize(const string& filename)
{
    TETRAHEDRALIZED_VOLUME<float> *sim_volume_float;
    sim_volume_float = TETRAHEDRALIZED_VOLUME<float>::Create();
    FILE_UTILITIES::Read_From_File<float>(filename, *sim_volume_float);
    
    sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
    sim_volume->particles.Add_Elements(sim_volume_float->particles.X.m);
    for (int i = 0; i < sim_volume_float->particles.X.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            sim_volume->particles.X(i)(j) = sim_volume_float->particles.X(i)(j);
        }
    }
    sim_volume->particles.Store_Velocity();
    sim_volume->mesh.elements = sim_volume_float->mesh.elements;
    
    Fit_In_Box<TV>(sim_volume->particles.X, RANGE<TV>(TV(-0.6,-0.6,-0.6),TV(0.6,0.6,0.6)));
    mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio, false);
    cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
}

template<class T>
void Fix_Orientation(TETRAHEDRALIZED_VOLUME<T>* volume)
{
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (TETRAHEDRON<T>(volume->particles.X(volume->mesh.elements(i)(0)),
                           volume->particles.X(volume->mesh.elements(i)(1)),
                           volume->particles.X(volume->mesh.elements(i)(2)),
                           volume->particles.X(volume->mesh.elements(i)(3))).Signed_Volume() < 0) {
            int temp = volume->mesh.elements(i)(0);
            volume->mesh.elements(i)(0) = volume->mesh.elements(i)(1);
            volume->mesh.elements(i)(1) = temp;
        }
    }
}

void initialize_cubes(int width, int height, int depth, T low, T high, T left, T right, T far, T near)
{
    sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
    sim_volume->particles.Store_Velocity();
    
    T lw = (right-left)/width;
    T lh = (high-low)/height;
    T ld = (far-near)/depth;
    int offset = (width+1)*(height+1);
    
    sim_volume->particles.Add_Elements(offset*(depth+1));
    sim_volume->Update_Number_Nodes();
    int node_index = 0;
    
    for (int i = 0; i < height+1; i++) {
        for (int j = 0; j < width+1; j++) {
            for (int k = 0; k < depth+1; k++) {
                sim_volume->particles.X(node_index+offset*k) = TV(left+j*lw, low+i*lh,near+k*ld);
            }
            node_index++;
        }
    }
    
    int tet_index = 0;
    for (int k = 0; k < depth; k++) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset, tet_index+1, tet_index, tet_index+width+1));
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+2, tet_index+width+1, tet_index+1, tet_index+offset));
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2, tet_index+offset+width+1, tet_index+offset));
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+1, tet_index+offset, tet_index+offset+1, tet_index+width+2));
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+width+1, tet_index+width+2, tet_index+offset+width+2, tet_index+offset+1));
                //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+1, tet_index+width+2, tet_index+offset, tet_index+offset+width+1));
                
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+1, tet_index+1, tet_index, tet_index+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index, tet_index+offset, tet_index+offset+1, tet_index+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+offset, tet_index+offset+1, tet_index+offset+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+1, tet_index+width+1, tet_index+offset+1, tet_index+width+2));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2, tet_index+offset+width+2, tet_index+offset+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2+offset, tet_index+offset+width+1, tet_index+offset+1));
                tet_index++;
            }
            tet_index++;
        }
        tet_index += (width+1);    
    }
    
    mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio, false);
    cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
}

template<class T>
void WriteToPovRaySmooth(TETRAHEDRALIZED_VOLUME<T>* volume, const string& outputDir, int frame)
{
    TRIANGULATED_SURFACE<T> *ts = TRIANGULATED_SURFACE<T>::Create();
    ts->Use_Vertex_Normals();
    int np = volume->particles.X.m;
    ts->particles.Add_Elements(np);
    for (int i = 0; i < np; ++i) {
        ts->particles.X(i) = volume->particles.X(i);
    }
    ts->Update_Number_Nodes();
    ts->mesh.elements = volume->mesh.boundary_mesh->elements;
    ts->avoid_normal_interpolation_across_sharp_edges = true;
    ts->Update_Vertex_Normals();
    
    stringstream ss;
    ss << frame;
    string filename = outputDir + "/povray" + ss.str() + ".pov";
    //cout << filename << endl;
    ofstream fs;
    fs.open(filename);
    fs << "#declare mesh = mesh2 {" << endl;
    fs << "vertex_vectors {" << endl;
    fs << ts->particles.X.m << ", " << endl;
    for (int i = 0; i < ts->particles.X.m; ++i) {
        fs << "<" << ts->particles.X(i)(0) << ", " << ts->particles.X(i)(1) << ", " << ts->particles.X(i)(2) << ">";
        if (i != ts->particles.X.m-1) {
            fs << ", ";
        }
    }
    fs << endl << "}" << endl;
    
//    ARRAY<TV> normals(volume->particles.X.m);
//    ARRAY<int> faceCount(volume->particles.X.m);
//    for (int i = 0; i < volume->mesh.boundary_mesh->elements.m; ++i) {
//        int n1 = volume->mesh.boundary_mesh->elements(i)(0);
//        int n2 = volume->mesh.boundary_mesh->elements(i)(1);
//        int n3 = volume->mesh.boundary_mesh->elements(i)(2);
//        TV e1 = volume->particles.X(n2) - volume->particles.X(n1);
//        TV e2 = volume->particles.X(n3) - volume->particles.X(n2);
//        TV normal = e1.Cross(e2);
//        normal.Normalize();
//        for (int j = 0; j < 3; ++j) {
//            normals(volume->mesh.boundary_mesh->elements(i)(j)) += normal;
//            ++faceCount(volume->mesh.boundary_mesh->elements(i)(j));
//        }
//    }
//    for (int i = 0; i < normals.m; ++i) {
//        if (faceCount(i) != 0) {
//            normals(i) /= faceCount(i);
//            normals(i).Normalize();
//        }
//    }
    
    fs << "normal_vectors {" << endl;
    fs << np << ", " << endl;
    for (int i = 0; i < np; ++i) {
        fs << "<" << (*(ts->vertex_normals))(i)(0) << ", " << (*(ts->vertex_normals))(i)(1) << ", " << (*(ts->vertex_normals))(i)(2) << ">";
        if (i != np-1) {
            fs << ", ";
        }
    }
    fs << endl << "}" << endl;
    
    fs << "face_indices {" << endl;
    fs << ts->mesh.elements.m << ", " << endl;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        fs << "<" << ts->mesh.elements(i)(0) << ", " << ts->mesh.elements(i)(1) << ", " << ts->mesh.elements(i)(2) << ">";
        if (i != ts->mesh.elements.m-1) {
            fs << ", ";
        }
    }
    fs << endl << "}" << endl;
    
    fs << "}" << endl;
    fs.close();
}

T laserZ1=-10;
T laserZ2=-9;
T laserX=0,laserY=0;

template<class T>
void WriteToPovRay(TETRAHEDRALIZED_VOLUME<T>* volume, const string& outputDir, int frame)
{
    TRIANGULATED_SURFACE<T> *ts = TRIANGULATED_SURFACE<T>::Create();
    ts->Use_Vertex_Normals();
    int np = volume->particles.X.m;
    ts->particles.Add_Elements(np);
    for (int i = 0; i < np; ++i) {
        ts->particles.X(i) = volume->particles.X(i);
    }
    ts->Update_Number_Nodes();
    ts->mesh.elements = volume->mesh.boundary_mesh->elements;
    //ts->Loop_Subdivide();//subd or go high res sim volume
    
    ts->avoid_normal_interpolation_across_sharp_edges = true;
    ts->normal_variance_threshold = 1;
    ts->Update_Vertex_Normals();
    
    stringstream ss;
    ss << frame;
    string filename = outputDir + "/povray" + ss.str() + ".pov";
    cout << filename << endl;
    ofstream fs;
    fs.open(filename);
    
    fs << "#declare LAZER_START=<" << laserX << ", " << laserY << ", " << laserZ1 << ">;" << endl;
    fs << "#declare LAZER_END=<" << laserX << ", " << laserY << ", " << laserZ2 << ">;" << endl;
    
    fs << "#declare mesh = mesh2 {" << endl;
    fs << "vertex_vectors {" << endl;
    fs << ts->mesh.elements.m*3 << ", " << endl;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            fs << "<" << ts->particles.X(ts->mesh.elements(i)(j))(0) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(1) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(2) << ">";
            if (i != ts->mesh.elements.m-1 && j != 2) {
                fs << ", ";
            }
        }
    }
    fs << endl << "}" << endl;
    
    //    ARRAY<TV> normals(volume->particles.X.m);
    //    ARRAY<int> faceCount(volume->particles.X.m);
    //    for (int i = 0; i < volume->mesh.boundary_mesh->elements.m; ++i) {
    //        int n1 = volume->mesh.boundary_mesh->elements(i)(0);
    //        int n2 = volume->mesh.boundary_mesh->elements(i)(1);
    //        int n3 = volume->mesh.boundary_mesh->elements(i)(2);
    //        TV e1 = volume->particles.X(n2) - volume->particles.X(n1);
    //        TV e2 = volume->particles.X(n3) - volume->particles.X(n2);
    //        TV normal = e1.Cross(e2);
    //        normal.Normalize();
    //        for (int j = 0; j < 3; ++j) {
    //            normals(volume->mesh.boundary_mesh->elements(i)(j)) += normal;
    //            ++faceCount(volume->mesh.boundary_mesh->elements(i)(j));
    //        }
    //    }
    //    for (int i = 0; i < normals.m; ++i) {
    //        if (faceCount(i) != 0) {
    //            normals(i) /= faceCount(i);
    //            normals(i).Normalize();
    //        }
    //    }
    
    fs << "normal_vectors {" << endl;
    fs << ts->mesh.elements.m*3 << ", " << endl;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            fs << "<" << (*(ts->face_vertex_normals))(i)(j)(0) << ", " << (*(ts->face_vertex_normals))(i)(j)(1) << ", " << (*(ts->face_vertex_normals))(i)(j)(2) << ">";
            if (i != ts->mesh.elements.m-1 && j != 2) {
                fs << ", ";
            }
        }
    }
    fs << endl << "}" << endl;
    
    fs << "face_indices {" << endl;
    fs << ts->mesh.elements.m << ", " << endl;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        fs << "<" << 3*i << ", " << 3*i+1 << ", " << 3*i+2 << ">";
        if (i != ts->mesh.elements.m-1) {
            fs << ", ";
        }
    }
    fs << endl << "}" << endl;
    
    fs << "}" << endl;
    fs.close();
}

template<class T>
void WriteToPovRay(TETRAHEDRALIZED_VOLUME<T>* volume, const string& outputDir, int frame, HASHTABLE<PhysBAM::VECTOR<int, 3> >& cuttingFaces)
{
    TRIANGULATED_SURFACE<T> *ts = TRIANGULATED_SURFACE<T>::Create();
    ts->Use_Vertex_Normals();
    int np = volume->particles.X.m;
    ts->particles.Add_Elements(np);
    for (int i = 0; i < np; ++i) {
        ts->particles.X(i) = volume->particles.X(i);
    }
    ts->Update_Number_Nodes();
    ts->mesh.elements = volume->mesh.boundary_mesh->elements;
    //ts->Loop_Subdivide();//subd or go high res sim volume
    
    ts->avoid_normal_interpolation_across_sharp_edges = true;
    ts->normal_variance_threshold = 1;
    ts->Update_Vertex_Normals();
    
    stringstream ss;
    ss << frame;
    string filename = outputDir + "/povray" + ss.str() + ".pov";
    cout << filename << endl;
    ofstream fs;
    fs.open(filename);
    
    const char* div = "";
    int total = ts->mesh.elements.m;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        if (cuttingFaces.Contains(ts->mesh.elements(i).Sorted())) {
            --total;
        }
    }
    
    //original faces
    fs << "#declare mesh = mesh2 {" << endl;
    fs << "vertex_vectors {" << endl;
    fs << ts->mesh.elements.m*3 << ", " << endl;
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            fs << div << "<" << ts->particles.X(ts->mesh.elements(i)(j))(0) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(1) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(2) << ">";
            div = ", ";
        }
    }
    fs << endl << "}" << endl;
    
    fs << "normal_vectors {" << endl;
    fs << ts->mesh.elements.m*3 << ", " << endl;
    div = "";
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            fs << div << "<" << (*(ts->face_vertex_normals))(i)(j)(0) << ", " << (*(ts->face_vertex_normals))(i)(j)(1) << ", " << (*(ts->face_vertex_normals))(i)(j)(2) << ">";
            div = ", ";
        }
    }
    fs << endl << "}" << endl;
    
    fs << "face_indices {" << endl;
    fs << total << ", " << endl;
    div = "";
    for (int i = 0; i < ts->mesh.elements.m; ++i) {
        if (cuttingFaces.Contains(ts->mesh.elements(i).Sorted())) {
            continue;
        }
        fs << div << "<" << 3*i << ", " << 3*i+1 << ", " << 3*i+2 << ">";
        div = ", ";
    }
    fs << endl << "}" << endl;
    fs << "}" << endl;
    
    if (total != ts->mesh.elements.m) {
        //cutting faces
        fs << "#declare cuttingFaces = mesh2 {" << endl;
        fs << "vertex_vectors {" << endl;
        fs << ts->mesh.elements.m*3 << ", " << endl;
        div = "";
        for (int i = 0; i < ts->mesh.elements.m; ++i) {
            for (int j = 0; j < 3; ++j) {
                fs << div << "<" << ts->particles.X(ts->mesh.elements(i)(j))(0) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(1) << ", " << ts->particles.X(ts->mesh.elements(i)(j))(2) << ">";
                div = ", ";
            }
        }
        fs << endl << "}" << endl;
        
        fs << "normal_vectors {" << endl;
        fs << ts->mesh.elements.m*3 << ", " << endl;
        div = "";
        for (int i = 0; i < ts->mesh.elements.m; ++i) {
            for (int j = 0; j < 3; ++j) {
                fs << div << "<" << (*(ts->face_vertex_normals))(i)(j)(0) << ", " << (*(ts->face_vertex_normals))(i)(j)(1) << ", " << (*(ts->face_vertex_normals))(i)(j)(2) << ">";
                div = ", ";
            }
        }
        fs << endl << "}" << endl;
        
        fs << "face_indices {" << endl;
        fs << ts->mesh.elements.m-total << ", " << endl;
        div = "";
        for (int i = 0; i < ts->mesh.elements.m; ++i) {
            if (cuttingFaces.Contains(ts->mesh.elements(i).Sorted())) {
                fs << div << "<" << 3*i << ", " << 3*i+1 << ", " << 3*i+2 << ">";
                div = ", ";
            }
        }
        fs << endl << "}" << endl;
        fs << "}" << endl;
    }
    
    
    
    fs.close();
}

template<class T>
void Write_Boundary_Mesh_To_File(const string& writing_directory, const string& filename, int frame, TETRAHEDRALIZED_VOLUME<T>* volume)
{
    stringstream ss;
    ss << frame;
    volume->mesh.Initialize_Boundary_Mesh();
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/"+filename+ss.str()+string(".tet.gz"), volume->particles, volume->mesh.boundary_mesh);
}


template<class T>
void Write_Volume_To_File(const string& writing_directory, const string& filename, int frame, TETRAHEDRALIZED_VOLUME<T>* volume)
{//assume volume boundary has been initialized
    stringstream ss;
    ss << frame;
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/"+filename+ss.str()+string(".tet.gz"), *volume);
}

//void energyTest() {
//    ofstream fs, fs1, fs2;
//    fs.open("results.txt");
//    fs1.open("results_hessian.txt");
//    fs2.open("scale.txt");
//    RANDOM_NUMBERS<T> rn;
//    ALGEBRA::VECTOR<T> perturbation1(mcut->sim_volume->particles.X.m*3);
//    for (int i = 0; i < perturbation1.Size(); i++){
//        perturbation1(i) = rn.Get_Number();
//        //cout << perturbation1(i) << endl;
//    }
//
//    for (int j = 1; j < 100; j++){
//        T eps = 1e-6*pow<T>(10, 0.03*j);//(T)j;
//        fs2 << -6 + j * 0.03 << endl;
//        
//        mcut->fem->Update_Position_Based_State();
//        ALGEBRA::VECTOR<T> f1(mcut->fem->Forces().Size());
//        for (int i = 0; i < f1.Size(); i++){f1(i) = mcut->fem->Forces()(i); }
//        
//        double e1 = 0;
//        for (int i = 0; i < mcut->deformable_object->Tetrahedron_Mesh().Number_Of_Tetrahedra(); i++){
//            e1 += mcut->le->Psi(i) / mcut->fem->Dm_inverse(i).Determinant() / 6;
//            //cout << mcut->le->Psi(i) << ", ";
//        }
//        //cout << "\n e1: " << e1 << endl;
//        
//        ALGEBRA::VECTOR<T> perturbation(mcut->sim_volume->particles.X.m*3);
//        for (int i = 0; i < perturbation.Size(); i++){
//            perturbation(i) = perturbation1(i) * eps;
//        }
//        ALGEBRA::VECTOR<T> df1(mcut->fem->Forces().Size());
//        mcut->fem->Stiffness_Matrix().Multiply(perturbation,df1);
//        
//        mcut->deformable_object->Perturb_Positions(perturbation);
//        mcut->fem->Update_Position_Based_State();
//        
//        ALGEBRA::VECTOR<T> df2(mcut->fem->Forces().Size());
//        mcut->fem->Stiffness_Matrix().Multiply(perturbation,df2);
//        
//        T e2 = 0;
//        for (int i = 0; i < mcut->deformable_object->Tetrahedron_Mesh().Number_Of_Tetrahedra(); i++){
//            e2 += mcut->le->Psi(i) / mcut->fem->Dm_inverse(i).Determinant() / 6.;
//        }
//        ALGEBRA::VECTOR<T> df(mcut->fem->Forces().Size());
//        for (int i = 0; i < df.Size(); i++){
//            //            cout << "df: " << mcut->fem->Forces()(i) - f1(i) << endl;
//            //            cout << "df12: " << (df1(i) + df2(i))/2 << endl;
//            df(i) = mcut->fem->Forces()(i) - f1(i) - (df1(i) + df2(i))/2;
//            
//        }
//        f1 += mcut->fem->Forces();
//        cout << (e2-e1) << ", " << f1.Dot(perturbation)/2. << endl;
//        //    cout << f1.Dot(perturbation)/2. << endl;
//        fs1 << log10(df.two_norm()/perturbation.two_norm()) << endl;
//        fs << log10(fabs((e2-e1+f1.Dot(perturbation)/2.)/perturbation.two_norm())) << endl;
//        
//        mcut->deformable_object->Perturb_Positions_Negative(perturbation);
//    }
//    fs.close();
//    fs1.close();
//    fs2.close();
//}

template<typename T>
void generateAndSaveRefinedVolume(TETRAHEDRALIZED_VOLUME<T>* refined_volume, int frame, const string& outputDir, const string& prefix) {
    mcut->Refine_And_Save_To(refined_volume);
    Fix_Orientation(refined_volume);
    
    TETRAHEDRALIZED_VOLUME<float> *f = TETRAHEDRALIZED_VOLUME<float>::Create();
    f->particles.Add_Elements(refined_volume->particles.X.m);
    f->Update_Number_Nodes();
    for (int i = 0; i < refined_volume->particles.X.m; ++i) {
        for (int j = 0; j < 3; ++j) {
            f->particles.X(i)(j) = refined_volume->particles.X(i)(j);
        }
    }
    f->mesh.elements = refined_volume->mesh.elements;
    Write_Boundary_Mesh_To_File(outputDir, prefix + "_boundary", frame, f);
    Write_Volume_To_File(outputDir, prefix, frame, f);
    
    ARRAY<int> l;
    refined_volume->mesh.Initialize_Boundary_Mesh();
    refined_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
    cout << "refined volume has " << l.Max() << " CCs\n";
}

int main(int argc, char** argv) {
    stringstream caseString(argv[1]);
    int caseN;
    caseString >> caseN;
    switch (caseN) {
        case 1:
        {
            string volumeFile(argv[2]);
            string outputDir(argv[3]);
            Initialize(volumeFile);
            cutting_tri_mesh->particles.Add_Elements(4);
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (mcut->sim_volume->particles.X(i).Magnitude_Squared() > 0.4 && !(fabs(mcut->sim_volume->particles.X(i)(0)) < 0.2 && mcut->sim_volume->particles.X(i)(1) > 0)) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
            //energyTest();
            //exit(1);
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int frame = 0;
            T xshift = 0.02;
            T yshift = 0;
            T x = 0.09;
            T y = x / 2;
//            for (int i = 0; i <= 30; ++i) {
//                ofstream fs;
//                stringstream ss;
//                ss << i;
//                fs.open("laser"+ss.str()+".inc");
//                if (i <= 15) {
//                    fs << "#declare LAZER_START=<" << cos(2*pi*i/20.)*x-xshift << ", " << y+sin(2*pi*i/20.)*y+yshift << ", " << -10 << ">;" << endl;
//                    fs << "#declare LAZER_END=<" << cos(2*pi*i/20.)*x-xshift << ", " << y+sin(2*pi*i/20.)*y+yshift << ", " << 0 << ">;" << endl;
//                }
//                else {
//                    T theta = pi/2-2*pi*(i-15)/20.;
//                    fs << "#declare LAZER_START=<" << cos(theta)*x-xshift << ", " << -y+sin(theta)*y+yshift << ", " << -10 << ">;" << endl;
//                    fs << "#declare LAZER_END=<" << cos(theta)*x-xshift << ", " << -y+sin(theta)*y+yshift << ", " << 0 << ">;" << endl;
//                }
//            }
//            exit(1);
            while (frame < 100) {
                ++frame;
//                if (frame == 20) {
//                    energyTest();
//                    exit(1);
//                }
                if (frame == 30) {
                    TRIANGULATED_SURFACE<T> *ts = new TRIANGULATED_SURFACE<T>();
//                    if (mcut->sim_volume->mesh.elements.m > 1e5) {
//                        mcut->ratio *= 40;
//                        ratio *= 40;
//                    }
                    for (int i = 0; i <= 30; ++i) {
                        if (i <= 15) {
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(2*pi*i/20.)*x-xshift, y+sin(2*pi*i/20.)*y+yshift, -1);
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(2*pi*i/20.)*x-xshift, y+sin(2*pi*i/20.)*y+yshift, 0);
                        }
                        else {
                            T theta = pi/2-2*pi*(i-15)/20.;
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, -y+sin(theta)*y+yshift, -1);
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, -y+sin(theta)*y+yshift, 0);
                        }
                        if (i > 0) {
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-2, 2*i-1, 2*i));
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-1, 2*i, 2*i+1));
                        }
                    }
                    mcut->Cut(*ts, 1);
                }
                
                T xx = 0.8;
                if (frame == 50) {
                    cutting_tri_mesh->particles.X(0) = TV(-xx, 0, -1);
                    cutting_tri_mesh->particles.X(1) = TV(-xx, 0, 1);
                    cutting_tri_mesh->particles.X(2) = TV(xx, 0, -1);
                    cutting_tri_mesh->particles.X(3) = TV(xx, 0, 1);
                    mcut->Cut(*cutting_tri_mesh, 1);
                }

//                if (frame == 50) {
//                    cutting_tri_mesh->particles.X(0) = TV(-0.5, -0.35, -1);
//                    cutting_tri_mesh->particles.X(1) = TV(-0.5, -0.35, 1);
//                    cutting_tri_mesh->particles.X(2) = TV(0.5, 0.65, -1);
//                    cutting_tri_mesh->particles.X(3) = TV(0.5, 0.65, 1);
//                    mcut->Cut(*cutting_tri_mesh);
//                }
//                
//                if (frame == 60) {
//                    cutting_tri_mesh->particles.X(0) = TV(-0.5, 0.65, -1);
//                    cutting_tri_mesh->particles.X(1) = TV(-0.5, 0.65, 1);
//                    cutting_tri_mesh->particles.X(2) = TV(0.5, -0.35, -1);
//                    cutting_tri_mesh->particles.X(3) = TV(0.5, -0.35, 1);
//                    mcut->Cut(*cutting_tri_mesh);
//                }
                
                if (frame == 80) {
                    cutting_tri_mesh->particles.X(0) = TV(0, -1, -1);
                    cutting_tri_mesh->particles.X(1) = TV(0, -1, 1);
                    cutting_tri_mesh->particles.X(2) = TV(0, 1, -1);
                    cutting_tri_mesh->particles.X(3) = TV(0, 1, 1);
                    mcut->Cut(*cutting_tri_mesh, 1);
                }

                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();

                generateAndSaveRefinedVolume(refined_volume, frame, outputDir, "refined_volume");
                
                WriteToPovRay(mcut->volume, outputDir, frame, mcut->cuttingFaces);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
                //Write_Volume_To_File(outputDir, "sim_volume", frame, mcut->sim_volume);
                Write_Volume_To_File(outputDir, "cutting_volume", frame, mcut->volume);
            }
            break;
        }

        //incrementally cut armadillo: change DEFORMABLE_OBJECT.h's Dm_inverse, use damping, no gravity
        case 2:
        {
            cout << "incremental\n";
            string volumeFile(argv[2]);
            string outputDir(argv[3]);
            Initialize(volumeFile);
            int n = 200;
            cutting_tri_mesh->particles.Add_Elements(2*(n+1));
            cutting_tri_mesh->Update_Number_Nodes();
            for (int i = 0; i < n; ++i) {
                cutting_tri_mesh->mesh.elements.Append(I3(2*i, 2*i+1, 2*i+2));
                cutting_tri_mesh->mesh.elements.Append(I3(2*i+1, 2*i+2, 2*i+3));
            }
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (mcut->sim_volume->particles.X(i).Magnitude_Squared() > 0.4 && !(fabs(mcut->sim_volume->particles.X(i)(0)) < 0.2 && mcut->sim_volume->particles.X(i)(1) > 0)) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int frame = 0;
            T xshift = 0.02;
            T yshift = 0.1;
            T x = 0.12;
            T y = x / 2;
            T xx = -0.3;
            T xxx = 0.01;
            T yy = 0.8;
            T yyy = 0.015;
            T x1 = 0, y1 = 0;
            int f1=30,f2=90,f3=105,f4=165,f5=200,f6=260;
            while (frame < 300) {
                ++frame;
//                if (frame == 2) {
//                    energyTest();
//                    exit(1);
//                }
                if (frame > f1 && frame <= f2) {
                    int i=frame-f1-1;
                    laserZ2=0.1;
                    T mid=(f2-f1)/2;
                    T m=mid*4./3;
                    T intv = 1.1 / n;
                    if (i < mid) {
                        x1=cos(2*pi*i/m)*x-xshift;
                        y1=y+sin(2*pi*i/m)*y+yshift;
                        for (int j = 0; j < n+1; ++j) {
                            cutting_tri_mesh->particles.X(2*j) = TV(x1, y1, -1+j*intv);
                        }
                        
                        ++i;
                        x1=cos(2*pi*i/m)*x-xshift;
                        y1=y+sin(2*pi*i/m)*y+yshift;
                        for (int j = 0; j < n+1; ++j) {
                            cutting_tri_mesh->particles.X(2*j+1) = TV(x1, y1, -1+j*intv);
                        }
                        
                    }
                    else {
                        T theta = pi/2-2*pi*(i-mid)/m;
                        x1=cos(theta)*x-xshift;
                        y1=-y+sin(theta)*y+yshift;
                        for (int j = 0; j < n+1; ++j) {
                            cutting_tri_mesh->particles.X(2*j) = TV(x1, y1, -1+j*intv);
                        }
                        
                        ++i;
                        theta -= (2*pi/m);
                        x1=cos(theta)*x-xshift;
                        y1=-y+sin(theta)*y+yshift;
                        for (int j = 0; j < n+1; ++j) {
                            cutting_tri_mesh->particles.X(2*j+1) = TV(x1, y1, -1+j*intv);
                        }
                    }
                    mcut->Cut(*cutting_tri_mesh, frame == f2, true);
                }
                else if (frame > f3 && frame <= f4) {
                    x1=xx+(frame-f3-1)*xxx;
                    y1=0;
                    laserZ2=10;
                    T intv = 2.0 / n;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(2*j) = TV(x1, y1, -1+j*intv);
                        cutting_tri_mesh->particles.X(2*j+1) = TV(x1+xxx, y1, -1+j*intv);
                    }
                    mcut->Cut(*cutting_tri_mesh, frame == f4, true);
                    if (frame == f4) {
                        Write_Volume_To_File(outputDir, "check", f4, mcut->volume);
                    }
                }
                else if (frame > f5 && frame <= f6) {
                    if (frame == f5 + 1) {
                        for (int i = 0; i < mcut->cutting_particle_material_space.m; ++i) {
                            mcut->cutting_particle_material_space(i) = mcut->volume->particles.X(i);
                        }
                    }
                    x1=0;
                    y1=yy-(frame-f5-1)*yyy;
                    laserZ2=10;
                    T intv = 2.0 / n;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(2*j) = TV(x1, y1, -1+j*intv);
                        cutting_tri_mesh->particles.X(2*j+1) = TV(x1, y1-yyy, -1+j*intv);
                    }
                    mcut->Cut(*cutting_tri_mesh, frame == f6, true);
                }
                else {
                    laserZ2=-9;
                }
                laserX=x1;
                laserY=y1;
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                generateAndSaveRefinedVolume(refined_volume, frame, outputDir, "cutting_volume");
                WriteToPovRay(refined_volume, outputDir, frame);
            }
            break;
        }
            
        case 3://cubes test
        {
            string outputDir(argv[2]);
            //initialize_cubes();

            cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
            cutting_tri_mesh->particles.Add_Elements(4);
            cutting_tri_mesh->particles.X(0) = TV(0.1, -0.5, -1.1);
            cutting_tri_mesh->particles.X(1) = TV(0.1, -0.5, 1);
            cutting_tri_mesh->particles.X(2) = TV(0.1, 0.6, -1.1);
            cutting_tri_mesh->particles.X(3) = TV(0.1, 0.6, 1);
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
            mcut->Cut(*cutting_tri_mesh, 1);

            Write_Volume_To_File(outputDir, "cutting_volume", 1, mcut->volume);
            break;
        }
            
        case 4://tets test
        {
            string outputDir(argv[2]);
            sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
            sim_volume->particles.Store_Velocity();
            
            sim_volume->particles.Add_Elements(6);
            sim_volume->particles.X(0) = TV(-0.5, 0.0, -0.5);
            sim_volume->particles.X(1) = TV(0.5, 0.0, -0.5);
            sim_volume->particles.X(2) = TV(0.5, 0.0, 0.5);
            sim_volume->particles.X(3) = TV(-0.5, 0.0, 0.5);
            sim_volume->particles.X(4) = TV(0.0, 0.0, 0.0);
            sim_volume->particles.X(5) = TV(0.0, 0.5, 0.0);
            for (int i = 0; i < 4; i++) {
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>((i+1)%4, i, 4, 5));
            }
            
            mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio, false);
            mcut->Initialize_Elasticity();
            
            cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
            cutting_tri_mesh->particles.Add_Elements(4);
            cutting_tri_mesh->particles.X(0) = TV(0.25, -1, -0.6);
            cutting_tri_mesh->particles.X(1) = TV(0.25, -1, 1);
            cutting_tri_mesh->particles.X(2) = TV(0.25, 1, -0.6);
            cutting_tri_mesh->particles.X(3) = TV(0.25, 1, 1);
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
            mcut->Cut(*cutting_tri_mesh, 1);
            
            Write_Volume_To_File(outputDir, "cutting_volume", 1, mcut->volume);
            break;
        }
            
        case 5://cubes, carve letters on it: change DEFORMABLE_OBJECT.h's Dm_inverse, change mesh_cutting_subd.cpp's dirichlet constraints
        {
            string outputDir(argv[2]);
            int width = 10;
            int height = 5;
            int depth = 2;
            T low = -0.5;
            T high = 0.5;
            T left = -1;
            T right = 1;
            T far = 0.05;
            T near = -0.05;
            initialize_cubes(width, height, depth, low, high, left, right, far, near);
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                TV p = mcut->sim_volume->particles.X(i);
                if (fabs(p(0)) > 0.99 || fabs(p(1)) > 0.49 || (PhysBAM::VECTOR<T,2>(p(0),p(1))-PhysBAM::VECTOR<T,2>(0.6534,0.15)).Magnitude_Squared() < 0.01 || (PhysBAM::VECTOR<T,2>(p(0),p(1))-PhysBAM::VECTOR<T,2>(-0.7,0.15)).Magnitude_Squared() < 0.0025 || (PhysBAM::VECTOR<T,2>(p(0),p(1))-PhysBAM::VECTOR<T,2>(-0.7,-0.15)).Magnitude_Squared() < 0.0025 || PhysBAM::VECTOR<T,2>(p(0),p(1)).Magnitude_Squared() < 0.01) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
            //energyTest();
            //exit(1);
            //TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int frame = 0;
            T xshift = 0;
            T yshift = 0.001;
            if (1) {
                ofstream fs;
                fs.open("sca_laser.txt");
                while (frame < 120) {
                    ++frame;
                    if (frame <= 15) {
                        T x = 0.2167;
                        T y = x / 1.4;
                        xshift = 0.7;
                        T theta = 2*pi*frame/20;
                        fs << cos(theta)*x-xshift << " " << y+sin(theta)*y+yshift << endl;
                    }
                    else if (frame <= 30) {
                        T x = 0.2167;
                        T y = x / 1.4;
                        xshift = 0.7;
                        T theta = pi / 2 - 2 * pi * (frame - 15) / 20;
                        fs << cos(theta)*x-xshift << " " << -y+sin(theta)*y+yshift << endl;
                    }
                    else if (frame <= 60) {
                        T x = 0.2567;
                        T y = 0.3;
                        xshift = 0;
                        T theta = 2 * pi * ((frame - 30) / 36. + 1. / 12.);
                        fs << cos(theta)*x-xshift << " " << sin(theta)*y+yshift << endl;
                    }
                    else if (frame <= 90) {
                        T x = 0.2567;
                        T y = 0.64;
                        yshift = -y / 2;
                        xshift = 0.6534;
                        T theta = pi * (frame - 60) / 30;
                        fs << cos(theta)*x+xshift << " " << sin(theta)*y+yshift << endl;;
                    }
                    else if (frame <= 105) {
                        T x = 0.249;
                        T dx = 2 * x / 15;
                        yshift = 0.003;
                        xshift = 0.6234;
                        fs << xshift - x + dx * (frame - 90) << " " << yshift << endl;;
                    }
                }
                fs.close();
                exit(1);
            }
            
            
            while (frame < 120) {
                ++frame;
                if (frame == 30) {
                    TRIANGULATED_SURFACE<T> *ts = new TRIANGULATED_SURFACE<T>();
                    T x = 0.2167;
                    T y = x / 1.4;
                    int f = 200;
                    xshift = 0.7;
                    for (int i = 0; i <= f; ++i) {
                        if (i <= f/2) {
                            ts->particles.Add_Elements(1);
                            T theta = 3*pi*i/f;
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, y+sin(theta)*y+yshift, -1);
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, y+sin(theta)*y+yshift, 1);
                        }
                        else {
                            T theta = 3*pi*(i-f/2)/f;
                            theta = pi/2 - theta;
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, -y+sin(theta)*y+yshift, -1);
                            ts->particles.Add_Elements(1);
                            ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, -y+sin(theta)*y+yshift, 1);
                        }
                        if (i > 0) {
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-2, 2*i-1, 2*i));
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-1, 2*i, 2*i+1));
                        }
                        
                    }
                    mcut->Cut(*ts, 1);
                }
                else if (frame == 60) {
                    TRIANGULATED_SURFACE<T> *ts = new TRIANGULATED_SURFACE<T>();
                    T x = 0.2567;
                    T y = 0.3;
                    xshift = 0;
                    int f = 240;
                    for (int i = 0; i <= f*10/12.; ++i) {
                        T theta = 2*pi*((T)i/f+1./12);
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, sin(theta)*y+yshift, -1);
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x-xshift, sin(theta)*y+yshift, 1);
                        if (i > 0) {
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-2, 2*i-1, 2*i));
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-1, 2*i, 2*i+1));
                        }
                    }
                    mcut->Cut(*ts, 1);
                }
                else if (frame == 90) {
                    TRIANGULATED_SURFACE<T> *ts = new TRIANGULATED_SURFACE<T>();
                    T x = 0.2567;
                    T y = 0.64;
                    yshift = -y / 2;
                    xshift = 0.6534;
                    int f = 200;
                    for (int i = 0; i <= f; ++i) {
                        T theta = pi * i / f;
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x+xshift, sin(theta)*y+yshift, -1);
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(cos(theta)*x+xshift, sin(theta)*y+yshift, 1);
                        if (i > 0) {
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-2, 2*i-1, 2*i));
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(2*i-1, 2*i, 2*i+1));
                        }
                    }
                    mcut->Cut(*ts, 1);
                }
                else if (frame == 105) {
                    int f = 30;
                    T x = 0.249;
                    T dx = 2 * x / f;
                    yshift = 0.003;
                    xshift -= 0.03;
                    TRIANGULATED_SURFACE<T> *ts = new TRIANGULATED_SURFACE<T>();
                    for (int i = 0; i <= f; ++i) {
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(xshift-x+dx*i, yshift, -1);
                        ts->particles.Add_Elements(1);
                        ts->particles.X(ts->particles.X.m-1) = TV(xshift-x+dx*(i+1), yshift, 1);
                        if (i > 0) {
                            int m = ts->particles.X.m;
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(m-4, m-3, m-2));
                            ts->mesh.elements.Append(PhysBAM::VECTOR<int,3>(m-2, m-3, m-1));
                        }
                    }
                    mcut->Cut(*ts, 1);
                }
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                WriteToPovRay(mcut->volume, outputDir, frame, mcut->cuttingFaces);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
                //Write_Volume_To_File(outputDir, "sim_volume", frame, mcut->sim_volume);
                Write_Volume_To_File(outputDir, "cutting_volume", frame, mcut->volume);
            }
            break;
        }

        case 6: //degeneracy cubes
        {
            string outputDir(argv[2]);
            int width = 10;
            int height = 10;
            int depth = 2;
            T low = -0.5;
            T high = 0.5;
            T left = -0.5;
            T right = 0.5;
            T far = 0.1;
            T near = -0.1;
            initialize_cubes(width, height, depth, low, high, left, right, far, near);
            sim_volume->particles.Store_Velocity();
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (fabs(mcut->sim_volume->particles.X(i)(0)) > 0.49 || fabs(mcut->sim_volume->particles.X(i)(1)) > 0.49) {
                    mcut->diri_nodes.Set(i);
                }
            }
            mcut->Initialize_Elasticity();
            
            int frame = 0;
            while (frame < 10) {
                ++frame;
                //                if (frame == 20) {
                //                    energyTest();
                //                    exit(1);
                //                }
                T x = 0.5;
                if (frame == 1) {
                    cutting_tri_mesh->particles.Add_Elements(4);
                    cutting_tri_mesh->Update_Number_Nodes();
                    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
                    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
                    cutting_tri_mesh->particles.X(0) = TV(x, 0.5, -1);
                    cutting_tri_mesh->particles.X(1) = TV(x, 0.5, 1);
                    cutting_tri_mesh->particles.X(2) = TV(-x, -0.5, -1);
                    cutting_tri_mesh->particles.X(3) = TV(-x, -0.5, 1);
                    mcut->Cut(*cutting_tri_mesh, 1);
                }
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();

                WriteToPovRay(mcut->volume, outputDir, frame, mcut->cuttingFaces);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
                //Write_Volume_To_File(outputDir, "sim_volume", frame, mcut->sim_volume);
                Write_Volume_To_File(outputDir, "cutting_volume", frame, mcut->volume);
            }
            break;
        }
        case 7://incremental cut cubes
        {
            string outputDir(argv[2]);
            int width = 10;
            int height = 10;
            int depth = 2;
            T low = -0.5;
            T high = 0.5;
            T left = -0.5;
            T right = 0.5;
            T far = 0.1;
            T near = -0.1;
            initialize_cubes(width, height, depth, low, high, left, right, far, near);
            sim_volume->particles.Store_Velocity();
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (fabs(mcut->sim_volume->particles.X(i)(0)) > 0.49 || fabs(mcut->sim_volume->particles.X(i)(1)) > 0.49) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
//            energyTest();
//            exit(1);
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            cutting_tri_mesh->particles.Add_Elements(4);
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
            int frame = 0;
            int f1 = 30;
            int f2 = 60;
            while (frame < 100) {
                ++frame;
                T y = 0.05;
                T x = -0.56;
                T l = -2 * x / (f2 - f1);
                if (frame >f1 && frame <= f2) {
                    cutting_tri_mesh->particles.X(0) = TV(x+l*(frame-f1-2), y, -1);
                    cutting_tri_mesh->particles.X(1) = TV(x+l*(frame-f1-2), y, 1);
                    cutting_tri_mesh->particles.X(2) = TV(x+l*(frame-f1), y, -1);
                    cutting_tri_mesh->particles.X(3) = TV(x+l*(frame-f1), y, 1);
                    mcut->Cut(*cutting_tri_mesh, frame == 60);
                }
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                mcut->Refine_And_Save_To(refined_volume);
                Fix_Orientation(refined_volume);
                refined_volume->Update_Number_Nodes();
                refined_volume->mesh.Initialize_Boundary_Mesh();
                Write_Boundary_Mesh_To_File(outputDir, "refined_volume_boundary", frame, refined_volume);
                Write_Volume_To_File(outputDir, "refined_volume", frame, refined_volume);
                
                WriteToPovRay(mcut->volume, outputDir, frame);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
            }
            break;
        }
        case 8: //cut on tet node or edge
        {
            string outputDir(argv[2]);
            sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
            sim_volume->particles.Store_Velocity();
            
            sim_volume->particles.Add_Elements(5);
            sim_volume->particles.X(0) = TV(0.5, 0.0, 0.0);
            sim_volume->particles.X(1) = TV(0.0, 0, -0.5);
            sim_volume->particles.X(2) = TV(-0.5, 0.0, 0.0);
            sim_volume->particles.X(3) = TV(0.0, 0.0, 0.5);
            sim_volume->particles.X(4) = TV(0.0, 0.5, 0);
            sim_volume->Update_Number_Nodes();
            
            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 3, 4));
            //sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 3, 4));

            mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio, false);
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (i == 0 || i == 2) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
            //energyTest();
            //exit(1);
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int frame = 0;
            //T xshift = 0;
            //T yshift = 0;
            while (frame < 10) {
                ++frame;
                //                if (frame == 20) {
                //                    energyTest();
                //                    exit(1);
                //                }
                //T xx = 0.8;
                if (frame == 1) {
                    cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
                    cutting_tri_mesh->particles.Add_Elements(4);
                    cutting_tri_mesh->Update_Number_Nodes();
                    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
                    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 2, 3));
                    cutting_tri_mesh->particles.X(0) = TV(0.1, 1, 1);
                    cutting_tri_mesh->particles.X(1) = TV(0.1, -1, 1);
                    cutting_tri_mesh->particles.X(2) = TV(0.1, -1, -1);
                    cutting_tri_mesh->particles.X(3) = TV(0.1, 1, -1);
                    mcut->Cut(*cutting_tri_mesh, 1);
                    
                    cutting_tri_mesh->particles.X(0) = TV(0.11, 1, 1);
                    cutting_tri_mesh->particles.X(1) = TV(0.11, -1, 1);
                    cutting_tri_mesh->particles.X(2) = TV(0.11, -1, -1);
                    cutting_tri_mesh->particles.X(3) = TV(0.11, 1, -1);
                    mcut->Cut(*cutting_tri_mesh, 1);
                }
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                   mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                generateAndSaveRefinedVolume(refined_volume, frame, outputDir, "refined_volume");
                
                WriteToPovRay(mcut->volume, outputDir, frame, mcut->cuttingFaces);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
                //Write_Volume_To_File(outputDir, "sim_volume", frame, mcut->sim_volume);
                Write_Volume_To_File(outputDir, "cutting_volume", frame, mcut->volume);
            }
            break;
        }
        case 9://material space procedural cut
        {
            string volumeFile(argv[2]);
            string outputDir(argv[3]);
            Initialize(volumeFile);
            cutting_tri_mesh->particles.Add_Elements(4);
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0, 1, 2));
            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1, 2, 3));
            
            //simulation setup: dirichlet, initial configuration
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if (mcut->sim_volume->particles.X(i).Magnitude_Squared() > 0.4 && !(fabs(mcut->sim_volume->particles.X(i)(0)) < 0.2 && mcut->sim_volume->particles.X(i)(1) > 0)) {
                    mcut->diri_nodes.Set(i);
                }
            }
            
            mcut->Initialize_Elasticity();
            
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            int frame = 0;
            while (frame < 100) {
                ++frame;
                T y = -0.0588;
                T l = 0.02;
                if (frame >30 && frame <= 60) {
                    cutting_tri_mesh->particles.X(0) = TV(-0.31+l*(frame-30), y, -1);
                    cutting_tri_mesh->particles.X(1) = TV(-0.31+l*(frame-30), y, 1);
                    cutting_tri_mesh->particles.X(2) = TV(-0.31+l*(frame-30+1), y, -1);
                    cutting_tri_mesh->particles.X(3) = TV(-0.31+l*(frame-30+1), y, 1);
                    mcut->Cut(*cutting_tri_mesh, frame == 60);
                }
                
                VS::start_timer();
                for (int i = 0; i < ratio; i++) {
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                mcut->Refine_And_Save_To(refined_volume);
                Fix_Orientation(refined_volume);
                refined_volume->Update_Number_Nodes();
                refined_volume->mesh.Initialize_Boundary_Mesh();
                Write_Boundary_Mesh_To_File(outputDir, "refined_volume_boundary", frame, refined_volume);
                Write_Volume_To_File(outputDir, "refined_volume", frame, refined_volume);
                
                WriteToPovRay(refined_volume, outputDir, frame);
                Write_Boundary_Mesh_To_File(outputDir, "cutting_volume_boundary", frame, mcut->volume);
            }
            break;
        }

        case 10://peel a ball: Dm_inverse = 1, use gravity, no damping
        {
            string volumeFile(argv[2]);
            string outputDir(argv[3]);
            Initialize(volumeFile);
            
            //dirichlet and peel nodes
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                // if in cylinder
                T myx=mcut->sim_volume->particles.X(i)(0);
                T myy=mcut->sim_volume->particles.X(i)(1);
                T myr=sqrt(sqr(myx)+sqr(myy));
                if(myr<0.5) mcut->diri_nodes.Set(i);
            }
            sim_volume->Update_Number_Nodes();
            sim_volume->mesh.Initialize_Boundary_Mesh();
            
            mcut->Initialize_Elasticity();
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int frame = 0;

            int f1=1,f2=61;
            int f = f2 - f1;
            T dtheta = 1.5 * pi / f;
            T r = 0.59;
            T w = sqrt(sqr(0.6)-sqr(r))*1.1;
            cout << "half width " << w / 1.1 << endl;
            T z1 = w;
            T z2 = -w;
            
            int n = 2;
            T dz = (z1 - z2) / n;
            cutting_tri_mesh->particles.Add_Elements(2*(n+1));
            cutting_tri_mesh->Update_Number_Nodes();
            for (int i = 0; i < n; ++i) {
                cutting_tri_mesh->mesh.elements.Append(I3(2*i, 2*i+1, 2*i+2));
                cutting_tri_mesh->mesh.elements.Append(I3(2*i+1, 2*i+2, 2*i+3));
            }
            
            while (frame < 70) {
                ++frame;
                if (1) {
                    if (frame == f1) {
                        for (int j = 0; j < n+1; ++j) {
                            T z = z1-dz*j;
                            cutting_tri_mesh->particles.X(2*j) = TV(-w, -r, z);
                            cutting_tri_mesh->particles.X(2*j+1) = TV(0, -r, z);
                        }
                        mcut->Cut(*cutting_tri_mesh, false, true);
                    }
                    else if (frame > f1 && frame <= f2) {
                        T theta = -pi / 2 + (frame - f1 - 1) * dtheta;
                        T x1 = r * cos(theta);
                        T y1 = r * sin(theta);
                        T x2 = r * cos(theta + dtheta);
                        T y2 = r * sin(theta + dtheta);
                        for (int j = 0; j < n+1; ++j) {
                            T z = z1-dz*j;
                            cutting_tri_mesh->particles.X(2*j) = TV(x1, y1, z);
                            cutting_tri_mesh->particles.X(2*j+1) = TV(x2, y2, z);
                        }
                        mcut->Cut(*cutting_tri_mesh, false, true);
                    }
                    else if (frame == f2 + 1) {
                        for (int j = 0; j < n+1; ++j) {
                            T z = z1-dz*j;
                            cutting_tri_mesh->particles.X(2*j) = TV(-r, 0, z);
                            cutting_tri_mesh->particles.X(2*j+1) = TV(-r, -w, z);
                        }
                        mcut->Cut(*cutting_tri_mesh, true, true);
                    }
                }
                
                if (1) {
                    VS::start_timer();
                    T theta = -pi/f/ratio;
                    PhysBAM::MATRIX<double,2> rotation(cos(theta),sin(theta),-sin(theta),cos(theta));
                    for (int r = 0; r < ratio; r++) {
                        int i = 0;
                        for (HASHTABLE_ITERATOR<int> it(mcut->diri_nodes); it.Valid(); it.Next()) {
                            int fixed_node = it.Key();
                            T oldx=mcut->deformable_object->X(fixed_node)(0);
                            T oldy=mcut->deformable_object->X(fixed_node)(1);
                            T oldz=mcut->deformable_object->X(fixed_node)(2);
                            PhysBAM::VECTOR<double,2> old(oldx,oldy);
                            PhysBAM::VECTOR<double,2> newp=rotation*old;
                            mcut->my_constrained->Set_Value(3*i,3*fixed_node);
                            mcut->my_constrained->Set_Value(3*i+1,3*fixed_node+1);
                            mcut->my_constrained->Set_Value(3*i+2,3*fixed_node+2);
                            mcut->my_constrained_locations->Set_Value(3*i,newp(0));
                            mcut->my_constrained_locations->Set_Value(3*i+1,newp(1));
                            mcut->my_constrained_locations->Set_Value(3*i+2,oldz);
                            i++;
                        }
                        mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
                        mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                    }
                    VS::stop_timer();
                    printf("frame %d time:    %f\n", frame, VS::get_time());
                    for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                        for (int k = 0; k<3; k++){
                            mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                            mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                        }
                    }
                    mcut->Update_Cutting_Particles();    
                }
                generateAndSaveRefinedVolume(refined_volume, frame, outputDir, "cutting_volume");
                WriteToPovRay(refined_volume, outputDir, frame);
            }
            break;
        }
            
        case 11://better peel a ball: Dm_inverse = 1, use gravity, no damping
        {
            string volumeFile(argv[2]);
            string outputDir(argv[3]);
            Initialize(volumeFile);
            
            //dirichlet and peel nodes
            for (int i = 0; i < mcut->sim_volume->particles.X.m; i++) {
                if(mcut->sim_volume->particles.X(i).Magnitude() < 0.55) mcut->diri_nodes.Set(i);
            }
            sim_volume->Update_Number_Nodes();
            sim_volume->mesh.Initialize_Boundary_Mesh();
            
            mcut->Initialize_Elasticity();
            TETRAHEDRALIZED_VOLUME<T> *refined_volume = new TETRAHEDRALIZED_VOLUME<T>();
            
            int f1 = 1, f2 = 181, f3 = 183, f4 = 273;
            int f = f2 - f1;
            int n = 5;//refine curve
            
            T r = 0.59;
            T w = sqrt(sqr(0.6)-sqr(r))*1.1;
            cout << "half width " << w / 1.1 << endl;
            
            T init_theta = -pi * 3 / 4;
            T dtheta = -4 * pi / f / ratio;
            T dtheta_cut = 4 * pi / f / n;
            T theta = init_theta;
            
            T init_phi = asin(-0.5 / 0.6);
            T final_phi = asin(0.5 / 0.6);
            T dphi_cut = (final_phi - init_phi) / f / n;
            T phi = init_phi;
            
            T pt = (final_phi - init_phi) / 4 / pi;
            
            cutting_tri_mesh->particles.Add_Elements((n+1)*3);
            cutting_tri_mesh->Update_Number_Nodes();
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < 2; ++k) {
                    int s = 3 * j + k;
                    cutting_tri_mesh->mesh.elements.Append(I3(s, s+3, s+1));
                    cutting_tri_mesh->mesh.elements.Append(I3(s+1, s+3, s+4));
                }
            }

            int frame = 0;
            while (frame < 280) {
                ++frame;
                if (frame == f1) {
                    TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                    TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                    TV dp = p.Cross(t).Normalized() * w;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(3*j) = p - dp - t * w * (n - j);
                        cutting_tri_mesh->particles.X(3*j+1) = p - t * w * (n - j);
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp - t * w * (n - j);
                    }
                    mcut->Cut(*cutting_tri_mesh, false, true);
                }
                else if (frame > f1 && frame <= f2) {
                    for (int j = 0; j < n+1; ++j) {
                        TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                        TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt);
                        TV dp = p.Cross(t).Normalized() * w;
                        cutting_tri_mesh->particles.X(3*j) = p - dp;
                        cutting_tri_mesh->particles.X(3*j+1) = p;
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp;
                        if (j != n) {
                            theta += dtheta_cut;
                            phi += dphi_cut;
                        }
                    }
                    mcut->Cut(*cutting_tri_mesh, false, true);
                }
                else if (frame == f2 + 1) {
                    TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                    TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                    TV dp = p.Cross(t).Normalized() * w;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(3*j) = p - dp + t * w * j;
                        cutting_tri_mesh->particles.X(3*j+1) = p + t * w * j;
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp + t * w * j;
                    }
                    mcut->Cut(*cutting_tri_mesh, true, true);
                }
                else if (frame == f3) {
                    for (int i = 0; i < mcut->cutting_particle_material_space.m; ++i) {
                        mcut->cutting_particle_material_space(i) = mcut->volume->particles.X(i);
                    }
                    
                    f = f4 - f3;
                    init_theta = -pi * 3 / 4;
                    dtheta_cut = 2 * pi / f / n;
                    theta = init_theta;
                    
                    init_phi = asin(-0.3 / 0.6);
                    final_phi = asin(0.3 / 0.6);
                    dphi_cut = (final_phi - init_phi) / f / n;
                    phi = init_phi;
                    
                    pt = (final_phi - init_phi) / 2 / pi;
                    
                    TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                    TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                    TV dp = p.Cross(t).Normalized() * w;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(3*j) = p - dp - t * w * (n - j);
                        cutting_tri_mesh->particles.X(3*j+1) = p - t * w * (n - j);
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp - t * w * (n - j);
                    }
                    mcut->Cut(*cutting_tri_mesh, false, true);
                }
                else if (frame > f3 && frame <= f4) {
                    for (int j = 0; j < n+1; ++j) {
                        TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                        TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt);
                        TV dp = p.Cross(t).Normalized() * w;
                        cutting_tri_mesh->particles.X(3*j) = p - dp;
                        cutting_tri_mesh->particles.X(3*j+1) = p;
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp;
                        if (j != n) {
                            theta += dtheta_cut;
                            phi += dphi_cut;
                        }
                    }
                    mcut->Cut(*cutting_tri_mesh, false, true);
                }
                else if (frame == f4 + 1) {
                    TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                    TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                    TV dp = p.Cross(t).Normalized() * w;
                    for (int j = 0; j < n+1; ++j) {
                        cutting_tri_mesh->particles.X(3*j) = p - dp + t * w * j;
                        cutting_tri_mesh->particles.X(3*j+1) = p + t * w * j;
                        cutting_tri_mesh->particles.X(3*j+2) = p + dp + t * w * j;
                    }
                    mcut->Cut(*cutting_tri_mesh, true, true);
                }
                
                VS::start_timer();
                PhysBAM::MATRIX<double,2> rotation(cos(dtheta),sin(dtheta),-sin(dtheta),cos(dtheta));
                for (int r = 0; r < ratio; r++) {
                    int i = 0;
                    for (HASHTABLE_ITERATOR<int> it(mcut->diri_nodes); it.Valid(); it.Next()) {
                        int fixed_node = it.Key();
                        T oldx=mcut->deformable_object->X(fixed_node)(0);
                        T oldy=mcut->deformable_object->X(fixed_node)(1);
                        T oldz=mcut->deformable_object->X(fixed_node)(2);
                        PhysBAM::VECTOR<double,2> old(oldy,oldz);
                        PhysBAM::VECTOR<double,2> newp=rotation*old;
                        mcut->my_constrained->Set_Value(3*i,3*fixed_node);
                        mcut->my_constrained->Set_Value(3*i+1,3*fixed_node+1);
                        mcut->my_constrained->Set_Value(3*i+2,3*fixed_node+2);
                        mcut->my_constrained_locations->Set_Value(3*i,oldx);
                        mcut->my_constrained_locations->Set_Value(3*i+1,newp(0));
                        mcut->my_constrained_locations->Set_Value(3*i+2,newp(1));
                        i++;
                    }
                    mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
                    mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                }
                VS::stop_timer();
                printf("frame %d time:    %f\n", frame, VS::get_time());
                for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                    for (int k = 0; k<3; k++){
                        mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                        mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
                
                generateAndSaveRefinedVolume(refined_volume, frame, outputDir, "cutting_volume");
                WriteToPovRay(refined_volume, outputDir, frame);
            }
            break;
        }
        default:
            break;
    }
}
