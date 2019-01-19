//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Vectors/VECTOR.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Computations/VORTICITY_UNIFORM.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <limits>
#include <string>
using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

struct FRAME_DATA
{
    GRID<TV> grid;
    ARRAY<T,FACE_INDEX<TV::m> > velocity;
    ARRAY<T,TV_INT> vorticity;
};

void Read_Output_Files(FRAME_DATA& fd,const std::string& input,int frame)
{
    Read_From_File(LOG::sprintf("%s/common/grid",input.c_str()),fd.grid);
    Read_From_File(LOG::sprintf("%s/%d/mac_velocities",input.c_str(),frame),fd.velocity);
}

VECTOR<T,2> Compute_Vorticity(FRAME_DATA& fd)
{
    VECTOR<T,2> r(std::numeric_limits<T>::max(),std::numeric_limits<T>::min());
    RANGE<TV_INT> domain_indices(fd.grid.Domain_Indices());
    domain_indices.Change_Size(-TV_INT::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<TV> lookup(fd.velocity);
    fd.vorticity.Resize(fd.grid.Domain_Indices());
    for(CELL_ITERATOR<TV> it(fd.grid,domain_indices);it.Valid();it.Next()){
        VECTOR<int,2> index=it.Cell_Index();
        T vort_mag=VORTICITY_UNIFORM<TV>::Vorticity(fd.grid,lookup,index).Magnitude();
        fd.vorticity(index)=vort_mag;
        r(0)=min(r(0),vort_mag);
        r(1)=max(r(1),vort_mag);}
    return r;
}

void Dump(const FRAME_DATA& fd,const std::string& output,T vmin,T vmax)
{
    OPENGL_COLOR_RAMP<T>* color_map=OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1);
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > image(fd.vorticity.domain);
    for(RANGE_ITERATOR<2> it(image.domain);it.Valid();it.Next()){
        T v=fd.vorticity(it.index);
        auto c=color_map->Lookup((v-vmin)/(vmax-vmin));
        image(it.index)=VECTOR<T,3>(c.rgba[0],c.rgba[1],c.rgba[2]);}
    PNG_FILE<T>::Write(output+".png",image);
    delete color_map;
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    int frame=0;
    T vmin=std::numeric_limits<T>::min(),vmax=std::numeric_limits<T>::max();
    std::string input="input",output="vort";
    parse_args.Extra(&input,"input","Input directory");
    parse_args.Add("-o",&output,"output","Output image");
    parse_args.Add("-frame",&frame,"frame","Frame");
    parse_args.Add("-vmin",&vmin,"vorticity","Minimum vorticity");
    parse_args.Add("-vmax",&vmax,"vorticity","Maximum vorticity");
    parse_args.Parse();

    FRAME_DATA fd;
    Read_Output_Files(fd,input,frame);
    auto bound=Compute_Vorticity(fd);
    LOG::printf("Min vort: %P Max vort: %P\n",bound(0),bound(1));
    Dump(fd,output,vmin,vmax);
    return 0;
}

