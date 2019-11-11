//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_LIST
//#####################################################################
#ifndef __STRUCTURE_LIST__
#define __STRUCTURE_LIST__

#include <Core/Data_Structures/DYNAMIC_LIST.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV,class ID>
class STRUCTURE_LIST:public DYNAMIC_LIST<STRUCTURE<TV>,ID>
{
    typedef typename TV::SCALAR T;
    typedef DYNAMIC_LIST<STRUCTURE<TV>,ID> BASE;
public:
    using BASE::Size;using BASE::Needs_Write;using BASE::Active_Element;using BASE::Read;using BASE::Write;
    using BASE::Element_Index;using BASE::Set_Active_Element;using BASE::Is_Active;
    mutable ARRAY<std::string> names;

    STRUCTURE_LIST()
        :BASE()
    {}

    void Read(const VIEWER_DIR& viewer_dir,const std::string& list_name)
    {
        ARRAY<ID> needs_init;
        BASE::Read(viewer_dir,list_name,needs_init);
        Read_From_File(viewer_dir.output_directory+"/common/"+list_name+"_key",names);
        for(int i=0;i<needs_init.m;i++){
            ID id=needs_init(i),index=Element_Index(id);
            Set_Active_Element(index,STRUCTURE<TV>::Create_From_Name(names(id)));
            Read_From_File(LOG::sprintf("%s/common/%s_%d.%s",
                    viewer_dir.output_directory,list_name,id,
                    Active_Element(index)->Extension()),*Active_Element(index));}

        for(ID id=0;id<Size();id++)
            if(Is_Active(id)){
                STRUCTURE<TV>& structure=*Active_Element(Element_Index(id));
                std::string filename=LOG::sprintf("%s/%s_%d.%s",viewer_dir.current_directory,list_name,id,structure.Extension());
                if(File_Exists(filename)) {
                    Read_From_File(filename,structure);
                    structure.update_every_frame=true;}
                else structure.update_every_frame=false;}
    }

    void Write(STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir, const std::string& list_name) const
    {
        BASE::Write(stream_type,viewer_dir,list_name);
        names.Resize(Value(Size()));
        ARRAY<ID>& needs_write=Needs_Write();
        for(int i=0;i<needs_write.Size();i++){
            ID id=needs_write(i);
            int index=Element_Index(id);
            assert(index>=0);
            names(id)=Active_Element(index)->Name();
            Write_To_File(stream_type,LOG::sprintf("%s/common/%s_%d.%s",
                    viewer_dir.output_directory,list_name,id,
                    Active_Element(index)->Extension()),*Active_Element(index));}
        if(viewer_dir.First_Frame() || needs_write.Size())
            Write_To_File(stream_type,viewer_dir.output_directory+"/common/"+list_name+"_key",names);
        needs_write.Remove_All();

        for(ID id=0;id<Size();id++)
            if(Is_Active(id)){
                STRUCTURE<TV>& structure=*Active_Element(Element_Index(id));
                if(structure.update_every_frame){
                    std::string filename=LOG::sprintf("%s/%s_%d.%s",viewer_dir.current_directory,list_name,id,structure.Extension());
                    Write_To_File(stream_type,filename,structure);}}
    }
};
}
#endif
