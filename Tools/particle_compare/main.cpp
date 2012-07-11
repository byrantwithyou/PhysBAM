//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function My_Sort
//#####################################################################
template<class T,int d>
struct VECTOR_COMPARATOR{
    bool operator()(const VECTOR<T,d>& v1, const VECTOR<T,d>& v2){
        return v1.Magnitude_Squared()<v2.Magnitude_Squared();}
};
template<class T,int d> void
My_Sort(ARRAY_VIEW<VECTOR<T,d> >& attr)
{
    VECTOR_COMPARATOR<T,d> vec_comp;
    Sort(attr,vec_comp);
}
template<class T> void
My_Sort(ARRAY_VIEW<T>& attr)
{
    Sort(attr);
}
//#####################################################################
// Function Compare_Attribute
//#####################################################################
template<class TV,class T2> bool
Compare_Attribute(ARRAY_VIEW<T2>& attr1,ARRAY_VIEW<T2>& attr2)
{
    if(attr1.m!=attr2.m){LOG::cout<<"WARNING: Number of particles differ"<<std::endl;return false;}
    My_Sort(attr1);My_Sort(attr2);
    for(int i=0;i<attr1.m;i++) if(attr1(i)!=attr2(i)){LOG::cout<<"WARNING: Particles differ starting at "<<i<<" of "<<attr1.m<<std::endl;return false;}
    return true;
}
//#####################################################################
// Function Compare_Attributes
//#####################################################################
template<class TV> bool
Compare_Attributes(PARTICLE_LEVELSET_PARTICLES<TV>& pls_particles_1,PARTICLE_LEVELSET_PARTICLES<TV>& pls_particles_2)
{
    typedef typename TV::SCALAR T;
    
    bool success=true;
    success&=Compare_Attribute<TV,TV>(pls_particles_1.X,pls_particles_2.X);
    if(!success) LOG::cout<<"Failed comparing position "<<std::endl;
    success&=Compare_Attribute<TV,T>(pls_particles_1.radius,pls_particles_2.radius);
    if(!success) LOG::cout<<"Failed comparing radius "<<std::endl;
    success&=Compare_Attribute<TV,T>(pls_particles_1.age,pls_particles_2.age);
    if(!success) LOG::cout<<"Failed comparing age "<<std::endl;
    success&=Compare_Attribute<TV,short unsigned int>(pls_particles_1.quantized_collision_distance,pls_particles_2.quantized_collision_distance);
    if(!success) LOG::cout<<"Failed comparing quantized collision distance "<<std::endl;
    return success;
}
//#####################################################################
// Function Compare_Attributes
//#####################################################################
template<class TV> bool
Compare_Attributes(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& pls_particles_1,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& pls_particles_2)
{
    bool success=true;
    success&=Compare_Attributes((PARTICLE_LEVELSET_PARTICLES<TV>)pls_particles_1,(PARTICLE_LEVELSET_PARTICLES<TV>)pls_particles_2);
    if(!success) return false;
    success&=Compare_Attribute<TV,TV>(pls_particles_1.V,pls_particles_2.V);
    if(!success) LOG::cout<<"Failed comparing velocity "<<std::endl;
    return success;
}
//#####################################################################
// Function Compare_Attributes
//#####################################################################
template<class TV,class T_PARTICLES> bool
Compare_Attributes(GRID<TV>& grid,ARRAY<T_PARTICLES*,VECTOR<int,TV::dimension> >& pls_particles_1,ARRAY<T_PARTICLES*,VECTOR<int,TV::dimension> >& pls_particles_2)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    bool success=true;
    for(typename GRID<TV>::NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        PARTICLE_LEVELSET_PARTICLES<TV> *particles_1=pls_particles_1(iterator.Node_Index()),*particles_2=pls_particles_2(iterator.Node_Index());
        if(!particles_1 && !particles_2) continue;
        if(!particles_1){LOG::cout<<"WARNING: 1 at "<<iterator.Node_Index()<<" is MISSING"<<std::endl;success=false;continue;}
        if(!particles_2){LOG::cout<<"WARNING: 2 at "<<iterator.Node_Index()<<" is MISSING"<<std::endl;success=false;continue;}
        success&=Compare_Attributes(*particles_1,*particles_2);
        if(!success) LOG::cout<<"Failed at index "<<iterator.Node_Index()<<std::endl;break;}
    return success;
}
//#####################################################################
// Function Compare_PLS_Particles
//#####################################################################
template<class TV,class T_PARTICLES> bool
Compare_PLS_Particles(std::string& input_directory_1,std::string& input_directory_2,const std::string& filename,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;

    ARRAY<T_PARTICLES*,TV_INT> pls_particles_1,pls_particles_2;
    FILE_UTILITIES::Read_From_File<T>(input_directory_1+filename,pls_particles_1);FILE_UTILITIES::Read_From_File<T>(input_directory_2+filename,pls_particles_2);
    bool success=Compare_Attributes(grid,pls_particles_1,pls_particles_2);
    if(!success) LOG::cout<<"Failed at frame "<<frame<<std::endl;
    return success;
}
//#####################################################################
// Function Compare_Particles
//#####################################################################
template<class TV> bool
Compare_Particles(std::string& input_directory_1,std::string& input_directory_2,GRID<TV>& grid,int frame,int type,bool verbose=true)
{
    typedef typename TV::SCALAR T;
    
    if(type==-1){
        bool success=true;
        for(int i=0;i<4;i++){
            bool local_success=Compare_Particles(input_directory_1,input_directory_2,grid,frame,i);
            if(local_success) LOG::cout<<" ....Passed"<<std::endl;
            else{std::cout<<"FAILED"<<std::endl;success=false;}}
        return success;}
    std::string filename="";if(verbose) LOG::cout<<"Comparing ";
    if(type==1){filename="positive_particles";if(verbose) LOG::cout<<"Positive Particles";}
    else if(type==2){filename="negative_particles";if(verbose) LOG::cout<<"Negative Particles";}
    else if(type==3){filename="removed_positive_particles";if(verbose) LOG::cout<<"Positive Removed Particles";}
    else if(type==4){filename="removed_negative_particles";if(verbose) LOG::cout<<"Negative Removed Particles";}

    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory_1+"/common/last_frame",last_frame);bool success=true;for(int i=0;i<last_frame;i++) success&=Compare_Particles(input_directory_1,input_directory_2,grid,i,type,false);return success;}

    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    bool success=false;
    if(type==1||type==2) success=Compare_PLS_Particles<TV,PARTICLE_LEVELSET_PARTICLES<TV> >(input_directory_1,input_directory_2,"/"+f+"/"+filename,grid,frame);
    else if(type==3||type==4) success=Compare_PLS_Particles<TV,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(input_directory_1,input_directory_2,"/"+f+"/"+filename,grid,frame);
    return success;
}
//#####################################################################
// Function Compare_Levelset
//#####################################################################
template<class TV> bool
Compare_Levelsets(std::string& input_directory_1,std::string& input_directory_2,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    typedef typename LEVELSET_POLICY<GRID<TV> >::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef VECTOR<int,TV::dimension> TV_INT;

    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory_1+"/common/last_frame",last_frame);bool success=true;for(int i=0;i<last_frame;i++){success&=Compare_Levelsets(input_directory_1,input_directory_2,grid,i);if(!success){LOG::cout<<"Failed at Frame "<<i<<std::endl;break;}}return success;}

    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    bool success=true;ARRAY<T,TV_INT> phi1,phi2;
    T_FAST_LEVELSET l1(grid,phi1);FILE_UTILITIES::Read_From_File<T>(input_directory_1+"/"+f+"/levelset",l1);
    T_FAST_LEVELSET l2(grid,phi2);FILE_UTILITIES::Read_From_File<T>(input_directory_2+"/"+f+"/levelset",l2);
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        if(l1.phi(iterator.Cell_Index())!=l2.phi(iterator.Cell_Index())){LOG::cout<<"WARNING: Phi don't match at index:["<<iterator.Cell_Index()<<std::endl;success=false;break;}}
    return success;
}
//#####################################################################
// Function Compare_Velocities
//#####################################################################
template<class TV> bool
Compare_Velocities(std::string& input_directory_1,std::string& input_directory_2,GRID<TV>& grid,int frame)
{
    typedef typename TV::SCALAR T;
    
    if(frame==-1){int last_frame=0;FILE_UTILITIES::Read_From_Text_File(input_directory_1+"/common/last_frame",last_frame);bool success=true;for(int i=0;i<last_frame;i++){success&=Compare_Velocities(input_directory_1,input_directory_2,grid,i);if(!success){LOG::cout<<"Failed at Frame "<<i<<std::endl;break;}}return success;}

    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    bool success=true;
    ARRAY<T,FACE_INDEX<TV::dimension> > u1;FILE_UTILITIES::Read_From_File<T>(input_directory_1+"/"+f+"/mac_velocities",u1);
    ARRAY<T,FACE_INDEX<TV::dimension> > u2;FILE_UTILITIES::Read_From_File<T>(input_directory_2+"/"+f+"/mac_velocities",u2);
    for(typename GRID<TV>::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        if(u1(iterator.Full_Index())!=u2(iterator.Full_Index())){LOG::cout<<"WARNING: Velocities don't match at axis:["<<iterator.Axis()<<"] index:["<<iterator.Face_Index()<<"]"<<std::endl;success=false;break;}}
    return success;
}
//#####################################################################
// Function Write_Output
//#####################################################################
template<class TV> void
Write_Output(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;

    std::string input_directory_1=parse_args.Extra_Arg(0),input_directory_2=parse_args.Extra_Arg(1);
    int frame=parse_args.Get_Integer_Value("-frame");
    int particle_type=parse_args.Get_Integer_Value("-type");
    GRID<TV> grid;FILE_UTILITIES::Read_From_File<T>(input_directory_1+"/common/grid",grid);
    GRID<TV> sanity;FILE_UTILITIES::Read_From_File<T>(input_directory_2+"/common/grid",sanity);assert(grid==sanity);
    bool success=Compare_Velocities(input_directory_1,input_directory_2,grid,frame);
    if(success) LOG::cout<<"Velocities Passed"<<std::endl;
    else{LOG::cout<<"Velocities FAILED"<<std::endl;return;}
    success=Compare_Levelsets(input_directory_1,input_directory_2,grid,frame);
    if(success) LOG::cout<<"Levelset Passed"<<std::endl;
    else{LOG::cout<<"Levelset FAILED"<<std::endl;return;}
    success=Compare_Particles(input_directory_1,input_directory_2,grid,frame,particle_type);
    if(success) LOG::cout<<"Particles Passed"<<std::endl;
    else{LOG::cout<<"Particles FAILED"<<std::endl;return;}
}
//#####################################################################
// Function Find_Dimension
//#####################################################################
template<class T> void
Find_Dimension(PARSE_ARGS& parse_args)
{
    if(parse_args.Is_Value_Set("-3d")){
        Write_Output<VECTOR<T,3> >(parse_args);}
    else if(parse_args.Is_Value_Set("-2d")){
        Write_Output<VECTOR<T,2> >(parse_args);}
    else{
        Write_Output<VECTOR<T,1> >(parse_args);}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Integer_Argument("-frame",-1,"frame output");
    parse_args.Add_Integer_Argument("-type",-1,"particle type");
    parse_args.Add_Option_Argument("-double","input data is in doubles");
    parse_args.Add_Option_Argument("-2d","input data is 2-D");
    parse_args.Add_Option_Argument("-3d","input data is 3-D");
    parse_args.Set_Extra_Arguments(1,"<input_directory_1>");
    parse_args.Set_Extra_Arguments(2,"<input_directory_2>");
    parse_args.Parse();

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(parse_args.Is_Value_Set("-double")) Find_Dimension<double>(parse_args); 
    else Find_Dimension<float>(parse_args);
#else
    Find_Dimension<float>(parse_args);
#endif
}
//#####################################################################
