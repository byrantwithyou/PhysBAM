//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>

using namespace PhysBAM;

typedef float T;

enum modType {SPLIT, MERGE, MOVE};
struct SplitMerge
{
    modType m;
    int pID[2];
};

template<class TV>
void Narain_To_PhysBAM(PARSE_ARGS& parse_args)
{
    typedef VECTOR<int,TV::m> TV_INT;
    static const int d=TV::m;
    std::string input_directory;
    std::string input_name;
    std::string output_directory;
    parse_args.Add("-i",&input_directory,"dir","input simulation directory");
    parse_args.Add("-n",&input_name,"name","input simulation name");
    parse_args.Add("-o",&output_directory,"dir","output simulation directory");
    parse_args.Parse();

    GRID<TV> grid;
    DEFORMABLE_BODY_COLLECTION<TV> deformable_body_collection(0,0);
    DEBUG_PARTICLES<TV> debug_particles;
    ARRAY_VIEW<MATRIX<T,3> > particle_A;
    deformable_body_collection.particles.Add_Array(ATTRIBUTE_ID(100),&particle_A);
    ARRAY<T,FACE_INDEX<d> > face_velocities;

    for(int frame=0;/*TODO*/;frame++)
    {
        std::string simfile = LOG::sprintf("%s/%s-%04i.sand",input_directory,input_name,frame);
        std::fstream file(simfile.c_str(), std::ios::in|std::ios::binary);
        PHYSBAM_ASSERT(file);
        int32_t frame_=0,zb=0,np=0,nb=0,nsm=0;
        file.read((char*)&frame_, 4);
        PHYSBAM_ASSERT(frame==frame_);
        T dt=0,simtime=0,dx=0,rho_max=0,r=0,mu=0;
        TV gravity;
        TV_INT size;
        char buff[100];
        file.read((char*)&dt, 4);
        file.read((char*)&simtime, 4);
        file.read((char*)&gravity, d*4);
        file.read((char*)&size, d*4);
        file.read((char*)&dx, 4);
        file.read((char*)&rho_max, 4);
        file.read((char*)&r, 4);
        file.read((char*)&mu, 4);
        file.read((char*)&zb, 4);
        grid.Initialize(size,RANGE<TV>(TV(),dx*TV(size)),true);
        face_velocities.Resize(grid);
        file.read((char*)&np, 4);
        deformable_body_collection.particles.Resize(np);
        for(int i = 0; i < np; i++)
        {
            file.read((char*)&deformable_body_collection.particles.X(i), d*4);
            file.read((char*)&deformable_body_collection.particles.V(i), d*4);
            file.read((char*)&deformable_body_collection.particles.mass(i), 4);
            file.read((char*)&particle_A(i), d*d*4);
        }
        file.read((char*)&nb, 4);
        for(int i = 0; i < nb; i++)
        {
            file.read(buff, 8);
            if(buff[0]==0)
                file.read(buff, 4+4*d*(d-1));
            int32_t len;
            char src[256];
            file.read((char*)&len, 4);
            file.read(src, len+8);
        }
        file.read((char*)&nsm, 4);
        ARRAY<SplitMerge> smmap(nsm);
        for(int i = 0; i < nsm; i++)
            file.read((char*)&smmap(i),sizeof(SplitMerge));

        if(d==2) file.read(buff, 4*7);

        std::string f=LOG::sprintf("%d",frame);
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),output_directory+"/common/grid",grid);
        if(!system(LOG::sprintf("rm -f %s/%d/mpm_particles.gz ;  ln -s ./deformable_object_particles.gz %s/%d/mpm_particles.gz",output_directory.c_str(),frame,output_directory.c_str(),frame).c_str())){}
        deformable_body_collection.Write(STREAM_TYPE((T)0),output_directory,output_directory,frame,-1,frame==0,false);
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),LOG::sprintf("%s/%d/mac_velocities",output_directory.c_str(),frame),face_velocities);
        for(int i=0;i<deformable_body_collection.particles.X.m;i++){
            Add_Debug_Particle(deformable_body_collection.particles.X(i),VECTOR<T,3>(1,0,1));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,deformable_body_collection.particles.V(i));}
        debug_particles.Write_Debug_Particles(STREAM_TYPE((T)0),output_directory,frame);
    }
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);
    
    if(use_3d)
        Narain_To_PhysBAM<VECTOR<T,3> >(parse_args);
    else
        Narain_To_PhysBAM<VECTOR<T,2> >(parse_args);

    return 0;
}



/*

  common/last_frame
  common/grid.gz
  common/first_frame
  common/deformable_object_structures.gz (??)

  ###/debug_particles.gz
  ###/frame_title
  ###/face_velocities.gz
  ###/deformable_object_particles.gz
  ###/mpm_particles.gz -> deformable_object_particles.gz


*/



