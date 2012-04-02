//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_CUTTING.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

using namespace PhysBAM;

// * Handle interface
// * Handle other boundary conditions
// * Handle 
typedef float RW;

std::string output_directory="output";
template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

template<class T,int d>
void Dump_Frame(const ARRAY<T,FACE_INDEX<d> >& u,const char* title)
{
    typedef VECTOR<T,d> TV;
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<TV>().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

template<class TV>
void Flush_Frame()
{
    Dump_Frame(ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),"flush");
}

template<class T,class TV,class T_OBJECT,class CM>
void Dump_Frame(const VECTOR_ND<T>& v,const CELL_MAPPING<TV>& index_map_p,CM* (index_map_u[TV::m]),const T_OBJECT& object,const char* title)
{
    static const int d=TV::m;
    INTERVAL<int> index_range_u[d] = {INTERVAL<int>(0,index_map_u[0]->next_index)};
    for(int i=1;i<d;i++) index_range_u[i]=INTERVAL<int>(index_range_u[i-1].max_corner,index_range_u[i-1].max_corner+index_map_u[i]->next_index);
    INTERVAL<int> index_range_p(index_range_u[d-1].max_corner,index_range_u[d-1].max_corner+index_map_p.next_index);
    INTERVAL<int> index_range_q[d] = {INTERVAL<int>(index_range_p.max_corner,index_range_p.max_corner+object.mesh.elements.m)};
    for(int i=1;i<d;i++) index_range_q[i]=INTERVAL<int>(index_range_q[i-1].max_corner,index_range_q[i-1].max_corner+object.mesh.elements.m);
    const GRID<TV>& grid=index_map_p.grid;

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(1e-11,10,true,true,true);

    char buff[100];
    for(int dir=0;dir<2;dir++)
    {
        sprintf(buff, title, dir?'-':'+');
        ARRAY<T,FACE_INDEX<TV::m> > error(grid);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            int index=index_map_u[it.Axis()]->Get_Index_Fixed(it.index,dir);
            if(index>=0)
                error(it.Full_Index())=v(index+index_range_u[it.Axis()].min_corner);}

        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
            int index=index_map_p.Get_Index_Fixed(it.index,dir);
            if(index>=0)
                Add_Debug_Particle(grid.Center(it.index),color_map(v(index_range_p.min_corner+index)));}

        for(int i=0;i<object.mesh.elements.m;i++){
            Add_Debug_Particle(object.Get_Element(i).Center(),VECTOR<T,3>(0,1,1));
            TV N;
            for(int j=0;j<d;j++) N(j)=v(index_range_q[j].min_corner+i);
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,N);}

        Dump_Frame(error,buff);}
}

template<class TV>
void Integration_Test(int argc,char* argv[])
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    Get_Debug_Particles<TV>();
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-viscosity",1,"viscosity");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-s",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Parse(argc,argv);
    T m=parse_args.Get_Double_Value("-m");
    T s=parse_args.Get_Double_Value("-s");
    T kg=parse_args.Get_Double_Value("-kg");
    T mu=parse_args.Get_Double_Value("-viscosity")*kg/(s*(d==3?m:1));

    TV_INT counts=TV_INT()+6;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*m,true);
    Global_Grid(&grid);
    GRID<TV> coarse_grid(grid.counts/2,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    SPHERE<TV> sphere(TV()+(T).43*m,(T).31*m);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=sphere.Signed_Distance(it.Location());

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    for(int i=0;i<object.mesh.elements.m;i++){
        Add_Debug_Particle(object.Get_Element(i).X(0),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,object.Get_Element(i).X(1)-object.Get_Element(i).X(0));}
    Flush_Frame<TV>();

    printf("\n");
    for(int i=0;i<d;i++) printf("%c [%i %i) ", "uvw"[i], index_range_u[i].min_corner, index_range_u[i].max_corner);
    printf("p [%i %i) ", index_range_p.min_corner, index_range_p.max_corner);
    for(int i=0;i<d;i++) printf("%cq [%i %i) ", "uvw"[i], index_range_q[i].min_corner, index_range_q[i].max_corner);
    printf("\n");

    ARRAY<T,TV_INT> f_body[d][2];
    VECTOR_ND<T> f_interface[d];
    for(int i=0;i<d;i++) f_interface[i].Resize(object.mesh.elements.m);
    for(int i=0;i<d;i++) for(int s=0;s<2;s++) f_body[i][s].Resize(grid.Domain_Indices());

    // TODO: fill in f_interface and f_body 

    VECTOR_ND<T> units(finish);
    VECTOR_ND<T> null[d],null_p(finish),z[d],z_p(finish);
    for(int i=0;i<d;i++){
        null[i].Resize(finish);
        z[i].Resize(finish);
        for(int j=index_range_u[i].min_corner;j<index_range_u[i].max_corner;j++)
            null[i](j)=units(j)=m/s;
        matrix.Times(null[i],z[i]);}

    for(int i=index_range_p.min_corner;i<index_range_p.max_corner;i++)
        null_p(i)=units(i)=kg/(s*s*(d==3?m:1));
    for(int i=0;i<d;i++)
        for(int j=index_range_q[i].min_corner;j<index_range_q[i].max_corner;j++){
            units(j)=kg/(s*s*(d==3?m:1));
            null_p(j)=-object.Get_Element(j-index_range_q[i].min_corner).Normal()(i)*units(j);}

    matrix.Times(null_p,z_p);

    for(int i=0;i<d;i++)
        LOG::cout<<z[i].Magnitude()<<std::endl;
    LOG::cout<<z_p.Magnitude()<<std::endl;

    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"finish setup");
    for(int i=0;i<d;i++){
        char buff[100];
        sprintf(buff, "null %c %%c", "uvw"[i]);
        Dump_Frame(null[i],index_map_p,index_map_u,object,buff);
        sprintf(buff, "error %c %%c", "uvw"[i]);
        Dump_Frame(z[i],index_map_p,index_map_u,object,buff);}

    Dump_Frame(null_p,index_map_p,index_map_u,object,"null p %c");
    Dump_Frame(z_p,index_map_p,index_map_u,object,"error p %c");

    OCTAVE_OUTPUT<T>("M.txt").Write("M",matrix);
//    LOG::cout<<units<<std::endl;
    OCTAVE_OUTPUT<T>("unit.txt").Write("u",units);
}

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    Integration_Test<VECTOR<double,3> >(argc,argv);

    return 0;
}
