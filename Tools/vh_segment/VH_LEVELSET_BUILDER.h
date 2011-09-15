//#####################################################################
// Copyright 2005-2006, Andy Selle, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VH_LEVELSET_BUILDER__
#define __VH_LEVELSET_BUILDER__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <string>

class Fl_Hold_Browser;

namespace PhysBAM{

template<class T>
class VH_LEVELSET_BUILDER
{
public:
    std::string data_directory;
    int box_padding_voxels;
    T box_padding;
    bool include_all_tissues;
    RANGE<VECTOR<int,3> > original_box_size;
    ARRAYS<VECTOR<unsigned short,2> > image;
    ARRAYS<VECTOR<unsigned short,2> > image_next;
    GRID<VECTOR<T,3> > grid;
    ARRAYS<VECTOR<T,3> > phi;
    LEVELSET_3D<GRID<VECTOR<T,3> > > levelset;
    VECTOR<int,3> image_resolution;
    VECTOR<T,3> image_size;
    ARRAYS<VECTOR<RANGE<VECTOR<int,3> >,1> > bounding_boxes;
    ARRAYS<VECTOR<std::string,1> > labels;
    
//#####################################################################
    VH_LEVELSET_BUILDER(std::string data_directory_in="../../Private_Data/VH_Raw/");
    void Create_Levelset_Interp(const ARRAY<int>& tissues,RANGE<VECTOR<int,3> >& box,VECTOR<int,3>& granularity,const std::string filename);
    void Create_Levelset(const ARRAY<int>& tissues,RANGE<VECTOR<int,3> >& box,VECTOR<int,3>& granularity,const std::string filename);
    void Dump_Images();
    void Read_Slice(int slice);
    void Read_Slice(int slice,ARRAYS<VECTOR<unsigned short,2> >& image);
    void Convert_To_Image(ARRAYS<VECTOR<VECTOR<T,3> ,2> >& image);
    void Initialize_Grid(RANGE<VECTOR<int,3> >& image_box,VECTOR<int,3>& grid_size);
    VECTOR<int,3> Get_Grid_Size(RANGE<VECTOR<int,3> >& image_box);
    void Make_Sub_Data_Set(std::string output_directory,RANGE<VECTOR<int,3> > box);
    void All_Tissues(Fl_Hold_Browser* tissues, HASHTABLE<int,int>& tissues_hash);
    void All_Tissues_Array(ARRAY<int>& tissues);
    RANGE<VECTOR<int,3> > Get_Bounding_Box(const ARRAY<int>& tissues);
    int Tissue_Index(std::string tissue_label);
private:
    void Compute_Bounding_Boxes();
    void Close_Boundaries();
//#####################################################################
};
}
#endif
