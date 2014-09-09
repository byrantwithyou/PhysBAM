//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <sstream>
#include "CONSISTENT_INTERSECTIONS.h"
using namespace PhysBAM;
using namespace std;
typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<T,3> T3;
typedef VECTOR<int,TV::m> TV_INT;
typedef VECTOR<int,2> I2;
typedef VECTOR<int,3> I3;
typedef VECTOR<int,4> I4;

string output_directory="output";
TV shiftV;
template<class TV>
void Dump(const SOLID_BODY_COLLECTION<TV>& solid_body_collection,DEBUG_PARTICLES<TV>& debug_particles,int frame,const char* frame_title)
{
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
    FILE_UTILITIES::Write_To_Text_File(output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),frame_title);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    debug_particles.Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    solid_body_collection.Write(STREAM_TYPE((RW())),output_directory,frame,0,true,true,true,true,false);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    int caseN=1;
    parse_args.Add("-case",&caseN,"number","case number to run");
    parse_args.Add("-o",&output_directory,"dir","output directory");
    parse_args.Parse();

    SOLID_BODY_COLLECTION<TV> solid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    DEBUG_PARTICLES<TV> debug_particles;
    SOLIDS_STANDARD_TESTS<TV> tests(STREAM_TYPE((RW())),getenv("PHYSBAM_DATA_DIRECTORY"),solid_body_collection);

    TRIANGULATED_AREA<T>* triangulated_area=0;
    SEGMENTED_CURVE_2D<T>* curve=0;
    T shift=0;
    int frames=0;
    
    switch (caseN) {
        case 1:
            shift=4.4;
            frames=480;
            triangulated_area=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),1),1);
            curve=SEGMENTED_CURVE_2D<T>::Create();
            curve->mesh.elements=triangulated_area->Get_Boundary_Object().mesh.elements;
            curve->particles.Append(triangulated_area->Get_Boundary_Object().particles);
            curve->particles.X+=TV(-2.2,0);
            shiftV=TV(shift/frames,0);
            break;
        case 2:
            frames=60;
            shift=0.7;
            triangulated_area=TRIANGULATED_AREA<T>::Create();
            curve=SEGMENTED_CURVE_2D<T>::Create();
            triangulated_area->particles.Add_Elements(4);
            triangulated_area->particles.X(0)=TV(0,0.5);
            triangulated_area->particles.X(1)=TV(0,0);
            triangulated_area->particles.X(2)=TV(0.5,0);
            triangulated_area->particles.X(3)=TV(0.5,0.5);
            triangulated_area->Update_Number_Nodes();
            triangulated_area->mesh.elements.Append(I3(0,1,2));
            triangulated_area->mesh.elements.Append(I3(0,2,3));
            curve->particles.Add_Elements(2);
            curve->particles.X(0)=TV(-0.1,0.3);
            curve->particles.X(1)=TV(-0.1,-0.001);
            curve->Update_Number_Nodes();
            curve->mesh.elements.Append(I2(0,1));
            break;
        case 3://tri_seg1
            frames=180;
            shift=2;
            triangulated_area=TRIANGULATED_AREA<T>::Create();
            curve=SEGMENTED_CURVE_2D<T>::Create();
            triangulated_area->particles.Add_Elements(3);
            triangulated_area->particles.X(0)=TV(0,0.5);
            triangulated_area->particles.X(1)=TV(-1/sqrt(3),-0.25);
            triangulated_area->particles.X(2)=TV(1/sqrt(3),-0.25);
            triangulated_area->Update_Number_Nodes();
            triangulated_area->mesh.elements.Append(I3(0,1,2));
            curve->particles.Add_Elements(2);
            curve->particles.X(0)=TV(-1,0.125);
            curve->particles.X(1)=TV(-1,-0.375);
            curve->Update_Number_Nodes();
            curve->mesh.elements.Append(I2(0,1));
            shiftV=TV(shift/frames,0);
            break;
        case 4://tri_seg2
            frames=180;
            shift=2;
            triangulated_area=TRIANGULATED_AREA<T>::Create();
            curve=SEGMENTED_CURVE_2D<T>::Create();
            triangulated_area->particles.Add_Elements(3);
            triangulated_area->particles.X(0)=TV(0,0.5);
            triangulated_area->particles.X(1)=TV(-1/sqrt(3),-0.25);
            triangulated_area->particles.X(2)=TV(1/sqrt(3),-0.25);
            triangulated_area->Update_Number_Nodes();
            triangulated_area->mesh.elements.Append(I3(0,1,2));
            curve->particles.Add_Elements(2);
            curve->particles.X(0)=TV(-1,0.495);
            curve->particles.X(1)=TV(-1,-0.249);
            curve->Update_Number_Nodes();
            curve->mesh.elements.Append(I2(0,1));
            shiftV=TV(shift/frames,0);
            break;
        case 5://tri_seg3
            frames=180;
            shift=1;
            triangulated_area=TRIANGULATED_AREA<T>::Create();
            curve=SEGMENTED_CURVE_2D<T>::Create();
            triangulated_area->particles.Add_Elements(3);
            triangulated_area->particles.X(0)=TV(0,0.5);
            triangulated_area->particles.X(1)=TV(-1/sqrt(3),-0.25);
            triangulated_area->particles.X(2)=TV(1/sqrt(3),-0.25);
            triangulated_area->Update_Number_Nodes();
            triangulated_area->mesh.elements.Append(I3(0,1,2));
            curve->particles.Add_Elements(3);
            curve->particles.X(0)=TV(-1,-0.4);
            curve->particles.X(1)=TV(0,-0.4);
            curve->particles.X(2)=TV(1,-0.4);
            curve->Update_Number_Nodes();
            curve->mesh.elements.Append(I2(0,1));
            curve->mesh.elements.Append(I2(1,2));
            shiftV=TV(0,shift/frames);
            break;
        default:
            break;
    }
    
    TRIANGULATED_AREA<T>& ta=tests.Copy_And_Add_Structure(*triangulated_area);
    SEGMENTED_CURVE_2D<T>& sc=tests.Copy_And_Add_Structure(*curve);

    sc.Update_Number_Nodes();
    ta.Update_Number_Nodes();
    ARRAY<int> cp(sc.mesh.elements.Flattened());
    cp.Prune_Duplicates();


    particles.mass.Fill(1);
    //particles.X.Subset(cp)+=TV(3,0);

    for(int i=0;i<frames;i++){
        LOG::printf("START %i\n",i);
        CONSISTENT_INTERSECTIONS<TV> ci(ta,sc);
        ci.Set_Tol();
        ci.Compute();

        for(HASHTABLE<I2>::ITERATOR it(ci.hash_vv);it.Valid();it.Next()){
            TV A=particles.X(it.Key().x),B=particles.X(it.Key().y);
            Add_Debug_Particle(A,T3(1,0,0));
            Add_Debug_Particle(B,T3(1,0,0));
            Add_Debug_Object(VECTOR<TV,2>(A,B),T3(1,0,0));}

        LOG::printf("ve %P\n",ci.hash_ve);
        for(HASHTABLE<I3,T>::ITERATOR it(ci.hash_ve);it.Valid();it.Next()){
            TV A=particles.X(it.Key().x),B=particles.X(it.Key().y),C=particles.X(it.Key().z),P=B+(C-B)*it.Data();
            Add_Debug_Particle(A,T3(0,1,0));
            Add_Debug_Particle(P,T3(0,1,0));
            Add_Debug_Object(VECTOR<TV,2>(A,P),T3(0,1,0));}

        for(HASHTABLE<I3,T>::ITERATOR it(ci.hash_ev);it.Valid();it.Next()){
            TV A=particles.X(it.Key().z),B=particles.X(it.Key().x),C=particles.X(it.Key().y),P=B+(C-B)*it.Data();
            Add_Debug_Particle(A,T3(1,1,0));
            Add_Debug_Particle(P,T3(1,1,0));
            Add_Debug_Object(VECTOR<TV,2>(A,P),T3(1,1,0));}

        for(HASHTABLE<I4,TV>::ITERATOR it(ci.hash_ee);it.Valid();it.Next()){
            TV A=particles.X(it.Key()(0)),B=particles.X(it.Key()(1)),C=particles.X(it.Key()(2)),D=particles.X(it.Key()(3));
            TV P=A+(B-A)*it.Data().x,Q=C+(D-C)*it.Data().y;
            Add_Debug_Particle(P,T3(1,0,1));
            Add_Debug_Particle(Q,T3(1,0,1));
            Add_Debug_Object(VECTOR<TV,2>(P,Q),T3(1,0,1));}

        LOG::printf("fv %P\n",ci.hash_fv);
        for(HASHTABLE<I4,T3>::ITERATOR it(ci.hash_fv);it.Valid();it.Next()){
            TV A=particles.X(it.Key()(3));
            Add_Debug_Particle(A,T3(1,1,1));}

        Dump(solid_body_collection,debug_particles,i,"dump");
        particles.X.Subset(cp)+=shiftV;
    }


    return 0;
}

