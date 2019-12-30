//#####################################################################
// Copyright 2015, Andre Pradhana and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>

using namespace PhysBAM;

template<class T,int d> void
Convert(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,d> TV;
    int start_at=0,end_at=INT_MAX;
    std::string input_dir,output_dir;
    bool surface_only=false;
    ARRAY<int> list;
    parse_args.Add("-start",&start_at,"start_at","Start conversion from this frame");
    parse_args.Add("-end",&end_at,"end_at","End conversion at this frame");
    parse_args.Add("-o",&output_dir,"dir","Output directory");
    parse_args.Add("-i",&input_dir,"dir","Input directory");
    parse_args.Add("-s",&list,"list","list of structures to convert");
    parse_args.Add("-surface",&surface_only,"dump surface in place of volumetric objects");
    parse_args.Parse();

    PHYSBAM_ASSERT(input_dir!="");
    PHYSBAM_ASSERT(output_dir!="");

    if(end_at==INT_MAX)
    {
        std::string filename=LOG::sprintf("%s/common/last_frame",input_dir);
        PHYSBAM_ASSERT(File_Exists(filename));
        Read_From_Text_File(filename,end_at);
    }

    Create_Directory(output_dir);

    DEFORMABLE_BODY_COLLECTION<TV> deformable_body_collection(0,0);
    VIEWER_DIR viewer_dir(input_dir);
    while(viewer_dir.Find_Next_Directory(0,false) && viewer_dir.frame_stack(0)<=end_at)
    {
        LOG::printf("Frame %P\n",viewer_dir.frame_stack);
        deformable_body_collection.Read(viewer_dir,false);
        if(viewer_dir.frame_stack(0)==start_at)
        {
            if(!list.m) list=IDENTITY_ARRAY<>(deformable_body_collection.structures.m);
            for(auto i:list)
                Create_Directory(LOG::sprintf("%s/%d",output_dir,i));
        }

        for(auto i:list)
        {
            LOG::printf("Object %i\n",i);
            std::string file=LOG::sprintf("%s/%d/object_%04d.obj",output_dir,i,viewer_dir.frame_stack(0));
            STRUCTURE<TV>* structure=deformable_body_collection.structures(i);
            assert(structure);
            TRIANGULATED_SURFACE<T> * ts=0, copy;
            if(auto* obj=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)) ts=obj;
            else if(auto* obj=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure))
            {
                obj->Initialize_Triangulated_Surface();
                ts=obj->triangulated_surface;
            }
            PHYSBAM_ASSERT(ts);
            copy.Compact_Copy(*ts);
            copy.Write_Obj(file);
        }
    }
}

int main(int argc,char *argv[])
{
    bool type_double=true;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Parse(true);

    if(type_double)
        Convert<double,3>(parse_args);
    else
        Convert<float,3>(parse_args);

    return 0;
}
