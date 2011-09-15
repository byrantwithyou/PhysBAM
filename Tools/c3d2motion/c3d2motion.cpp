//#####################################################################
// Copyright 2005, Andrew Selle, Eftychios Sifakis
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Motion/MOTION_SEQUENCE.h>
#include <sstream>
#include <string>
#include "C3D_FILE.h"
#include <math.h>

using namespace PhysBAM;

template<class T> void
Convert(const std::string& input_filename,const std::string& output_filename)
{
    C3D_FILE c3d_file;c3d_file.Initialize(input_filename.c_str());
    int marker_count=c3d_file.Number_Channels();
    int frame_count=c3d_file.Number_Frames();
    MOTION_SEQUENCE<T,VECTOR<T, 3> > sequence;
    sequence.Initialize(marker_count,GRID_1D<T>(c3d_file.Number_Frames(),0,(c3d_file.Number_Frames()-1)/(T)c3d_file.Frame_Rate()));
    c3d_file.Get_Labels(sequence.names);
    ARRAY<C3D_POINT> points;
    for(int i=1;i<=marker_count;i++) for(int frame=1;frame<=frame_count;frame++){
        C3D_POINT& point=c3d_file.frames(frame,i);
        sequence.trajectories(i)(frame)=(T)1e-3*VECTOR<T,3>(point.x,point.y,point.z);
        sequence.valid(i)(frame)=(bool)point.camera_visible;}


    std::cout<<"Markers "<<marker_count<<" Frames "<<frame_count<<std::endl;
    std::cout<<"Marker names ";
    for(int i=1;i<=marker_count;i++) std::cout<<" "<<sequence.names(i);
    std::cout<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(output_filename,sequence);
}

template<class T> void
Get_Joint_Angles(const std::string& input_filename,const std::string& output_filename, const char* joint_filename)
{
    C3D_FILE c3d_file;c3d_file.Initialize(input_filename.c_str());
    
    int marker_count=c3d_file.Number_Channels();
    int frame_count=c3d_file.Number_Frames();

    MOTION_SEQUENCE<T,VECTOR<T, 3> > sequence;
    sequence.Initialize(marker_count,GRID_1D<T>(c3d_file.Number_Frames(),0,(c3d_file.Number_Frames()-1)/(T)c3d_file.Frame_Rate()));
    c3d_file.Get_Descriptions(sequence.names);
    ARRAY<C3D_POINT> points;
    for(int i=1;i<=marker_count;i++) for(int frame=1;frame<=frame_count;frame++){
        C3D_POINT& point=c3d_file.frames(frame,i);
        sequence.trajectories(i)(frame)=(T)1e-3*VECTOR<T,3>(point.x,point.y,point.z);
        sequence.valid(i)(frame)=(bool)point.camera_visible;}
    Clear_Spaces(marker_count, sequence);
    std::cout<<"Markers "<<marker_count<<" Frames "<<frame_count<<std::endl;
    std::cout<<"Marker names ";
    for(int i=1;i<=marker_count;i++) std::cout<<" "<<sequence.names(i);
    std::cout<<std::endl;

    FILE *fp=fopen(joint_filename,"rt");
    if(!fp) return;
    int joint_count=0;
    char string[100];
    while (!feof(fp)){
        fgets(string, 100, fp);
        joint_count++;}
    fclose(fp);
    fp=fopen(joint_filename, "rt");
    MOTION_SEQUENCE<T,QUATERNION<T> > joint_angles;
    joint_angles.Initialize(joint_count,GRID_1D<T>(c3d_file.Number_Frames(),0,(c3d_file.Number_Frames()-1)/(T)c3d_file.Frame_Rate()));
    int markerids[8];
    char joint_name[100];
    char joints[8][100];
    for (int i=1;i<=joint_count;i++){
        fgets(string,100,fp);
        sscanf(string,"%s %s %s %s %s %s %s %s %s",joint_name,joints[0],joints[1],joints[2],joints[3],joints[4],joints[5],joints[6],joints[7]);
        joint_angles.names(i)=std::string(joint_name);
        std::cout<<"Joint name is "<< joint_angles.names(i)<<"\n";
        for (int id=0;id<8;id++) markerids[id]=Find_Name(std::string(joints[id]),marker_count,sequence);
        for (int j=1;j<=frame_count;j++){
            if (markerids[0]==-1||markerids[2]==-1||markerids[4]==-1||markerids[6]==-1){
                joint_angles.trajectories(i)(j)=QUATERNION<T>(0,0,0,0);
                joint_angles.valid(i)(j)=true;}
            else{
                if (joint_angles.names(i)=="joint_r_glenohumeral"||joint_angles.names(i)=="joint_r_glenohumeral_left"){
                        VECTOR<T,3> shoulder=Vector(markerids,0,2,j,sequence);
                        VECTOR<T,3> temp=shoulder;
                        shoulder=VECTOR<T,3>(temp.y,-temp.x,-temp.z);
                        VECTOR<T,3> elbow=Vector(markerids,4,6,j,sequence);
                        temp=elbow;
                        elbow=VECTOR<T,3>(temp.y,-temp.x,-temp.z);
                        int original=markerids[2];
                        if (joint_angles.names(i)=="joint_r_glenohumeral") markerids[2]=Find_Name(std::string("rarm"),marker_count,sequence);
                        else markerids[2]=Find_Name(std::string("larm"),marker_count,sequence);
                        VECTOR<T,3> arm=Vector(markerids,0,2,j,sequence);
                        temp=arm;
                        arm=VECTOR<T,3>(temp.y,-temp.x,-temp.z);
                        VECTOR<T,3> projection=arm-VECTOR<T,3>::Dot_Product(arm,elbow)/elbow.Magnitude_Squared()*elbow;
                        QUATERNION<T> rotation=QUATERNION<T>::Rotation_Quaternion(shoulder,elbow);
                        QUATERNION<T> arm_rotation=QUATERNION<T>::Rotation_Quaternion(projection,rotation.v);
                        T angle;
                        if (joint_angles.names(i)=="joint_r_glenohumeral_left") angle=-arm_rotation.Angle();
                        else angle=pi-arm_rotation.Angle();
                        if(VECTOR<T,3>::Dot_Product(elbow,arm_rotation.v)<0) angle = -angle;
                        arm_rotation=QUATERNION<T>(angle,VECTOR<T,3>(0,0,1));
                        markerids[2]=original;
                        shoulder=arm_rotation.Inverse().Rotate(shoulder);
                        elbow=arm_rotation.Inverse().Rotate(elbow);
                        rotation=QUATERNION<T>::Rotation_Quaternion(shoulder,elbow);
                        rotation=QUATERNION<T>(rotation.Angle()-pi/2, rotation.v);
                        joint_angles.trajectories(i)(j)=rotation*arm_rotation;
                }
                else{
                    if (joint_angles.names(i)=="joint_r_hip"||joint_angles.names(i)=="joint_r_hip_left"){
                        VECTOR<T,3> init_angle=Vector(markerids,0,2,j,sequence);
                        VECTOR<T,3> temp=init_angle;
                        init_angle=VECTOR<T,3>(temp.y,temp.z,-temp.x);
                        VECTOR<T,3> final_angle=Vector(markerids,4,6,j,sequence);
                        temp=final_angle;
                        final_angle=VECTOR<T,3>(temp.y,temp.z,-temp.x);
                        joint_angles.trajectories(i)(j)=QUATERNION<T>::Rotation_Quaternion(init_angle,final_angle);  
                        if (joint_angles.names(i)=="joint_r_hip_left") joint_angles.trajectories(i)(j).v.y=-joint_angles.trajectories(i)(j).v.y;} 
                    else{
                        if (joint_angles.names(i)=="joint_l1"){
                            VECTOR<T,3> init_angle = Vector(markerids,0,2,j,sequence);
                            VECTOR<T,3> temp = init_angle;
                            init_angle=VECTOR<T,3>(-temp.y,-temp.x,temp.z);
                            VECTOR<T,3> final_angle = Vector(markerids,4,6,j,sequence);
                            temp = final_angle;
                            final_angle=VECTOR<T,3>(-temp.y,-temp.x,temp.z);
                            joint_angles.trajectories(i)(j)=QUATERNION<T>::Rotation_Quaternion(init_angle,final_angle);}
                        else joint_angles.trajectories(i)(j)=Line_Up(markerids,VECTOR<T,3>(0,0,1),VECTOR<T,3>(1,0,0),j,sequence);}}
                if (joint_angles.names(i)!="joint_r_glenohumeral"&&joint_angles.names(i)!="joint_r_glenohumeral_left")
                    joint_angles.trajectories(i)(j)=joint_angles.trajectories(i)(j).Inverse();
                joint_angles.valid(i)(j)=true;}}}
    fclose(fp);
    FILE_UTILITIES::Write_To_File<T>(output_filename,joint_angles);
}

template <class T> int Find_Name(std::string to_find, int marker_count, MOTION_SEQUENCE<T,VECTOR<T,3> > sequence)
{
    for (int i=1;i<=marker_count;i++) if (sequence.names(i)==to_find) return i;
    return -1;
}

template <class T> void Clear_Spaces(int marker_count, MOTION_SEQUENCE<T,VECTOR<T,3> > &sequence)
{
    for (int i=1;i<=marker_count;i++){
        std::string no_spaces=sequence.names(i);
        for (unsigned int j=0;j<no_spaces.size();){
            if (no_spaces[j]==' '){
                no_spaces=no_spaces.erase(j,1);
                j--;}
            j++;}
        sequence.names(i)=no_spaces;
    }
}

template <class T> QUATERNION<T> Line_Up(int markerids[], VECTOR<T,3> first, VECTOR<T,3> second, int frame, MOTION_SEQUENCE<T,VECTOR<T,3> > sequence)
{
    VECTOR<T,3> parent=Vector(markerids,0,2,frame,sequence);
    VECTOR<T,3> child=Vector(markerids,4,6,frame,sequence);
    QUATERNION<T> line_up_first=QUATERNION<T>::Rotation_Quaternion(child,first);
    line_up_first.Normalize();
    parent=line_up_first.Rotate(parent);
    child=line_up_first.Rotate(child);
    QUATERNION<T> angle=QUATERNION<T>::Rotation_Quaternion(parent,child);
    QUATERNION<T> line_up_second=QUATERNION<T>::Rotation_Quaternion(angle.v,second);
    parent = line_up_second.Rotate(parent);
    child = line_up_second.Rotate(child);
    return QUATERNION<T>::Rotation_Quaternion(parent, child);
}

template <class T> VECTOR<T,3> Vector(int markerids[], int first_index, int second_index, int frame, MOTION_SEQUENCE<T,VECTOR<T,3> > sequence)
{
    VECTOR<T,3> point1=sequence.Frame_X(markerids[first_index],frame);
    if (markerids[first_index+1]!=-1){
        point1+=sequence.Frame_X(markerids[first_index+1],frame);
        point1/=2;}
    VECTOR<T,3> point2=sequence.Frame_X(markerids[second_index],frame);
    if (markerids[second_index+1]!=-1){
        point2+=sequence.Frame_X(markerids[second_index+1],frame);
        point2/=2;}
    return point2-point1;
}

int main(int argc,char *argv[])
{
    if (argc == 3)
        Convert<double>(argv[1],argv[2]);
    else
    {
        if (argc == 4)
            Get_Joint_Angles<double>(argv[1], argv[2], argv[3]);
        else {std::cerr<<"Usage is: "<<argv[0]<<" <input c3d> <output motion> (<input linkages>)"<<std::endl;return 1;}
    }
    return 0;
}
