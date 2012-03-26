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
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
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
    for(int dir=0;dir<TV::m;dir++)
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

    GRID<TV> grid(TV_INT()+4,RANGE<TV>(TV(),TV()+1)*m,true);
    Global_Grid(&grid);
    GRID<TV> coarse_grid(grid.counts/2,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    SPHERE<TV> sphere(TV()+(T).5*m,(T).3*m);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=sphere.Signed_Distance(it.Location());

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    BASIS_STENCIL_UNIFORM<TV> p_stencil(grid.dX),*u_stencil[d],*udx_stencil[d][d];
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Dice_Stencil();
    for(int i=0;i<d;i++){
        u_stencil[i]=new BASIS_STENCIL_UNIFORM<TV>(grid.dX);
        u_stencil[i]->Set_Face(i);
        u_stencil[i]->Set_Multilinear_Stencil();
        u_stencil[i]->Dice_Stencil();
        for(int j=0;j<d;j++){
            udx_stencil[i][j]=new BASIS_STENCIL_UNIFORM<TV>(*u_stencil[i]);
            udx_stencil[i][j]->Differentiate(j);
            udx_stencil[i][j]->Dice_Stencil();}}

    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT object;
    MARCHING_CUBES<TV>::Create_Surface(object,coarse_grid,phi);
    BASIS_STENCIL_BOUNDARY_UNIFORM<TV> q_stencil(object);

    for(int i=0;i<object.mesh.elements.m;i++){
        Add_Debug_Particle(object.Get_Element(i).X(0),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,object.Get_Element(i).X(1)-object.Get_Element(i).X(0));}
    Flush_Frame<TV>();

    CELL_MAPPING<TV> index_map_p(grid),*index_map_u[d];
    index_map_p.periodic.Fill(true);
    for(int i=0;i<d;i++){
        index_map_u[i]=new CELL_MAPPING<TV>(grid);
        index_map_u[i]->periodic.Fill(true);}

    SYSTEM_MATRIX_HELPER<T> helper;
    RANGE<TV_INT> boundary_conditions;
    boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);
    boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);

    INTERVAL<int> block_uu[d][d],block_p[d],block_q[d];

    // Diagonal blocks
    for(int i=0;i<d;i++){
        int start=helper.start;
        for(int j=0;j<d;j++){
            BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,*udx_stencil[i][j],*udx_stencil[i][j],*index_map_u[i],*index_map_u[i],phi).Compute_Matrix(helper);
            helper.Scale(mu*(1+(i==j)));}
        block_uu[i][i]=helper.Get_Block();
        block_uu[i][i].min_corner=start;
        helper.New_Block();}

    // Off-diagonal viscosity blocks
    int start_off_diag=helper.start;
    for(int i=0;i<d;i++)
        for(int j=i+1;j<d;j++){
            BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,*udx_stencil[i][j],*udx_stencil[j][i],*index_map_u[i],*index_map_u[j],phi).Compute_Matrix(helper);
            helper.Scale(mu);
            block_uu[i][j]=helper.Get_Block();}

    // Pressure blocks
    for(int i=0;i<d;i++){
        BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,*udx_stencil[i][i],p_stencil,*index_map_u[i],index_map_p,phi).Compute_Matrix(helper);
        block_p[i]=helper.Get_Block();}

    // Traction blocks
    for(int i=0;i<d;i++){
        BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>(grid,*u_stencil[i],q_stencil,*index_map_u[i]).Compute_Matrix(helper);
        block_q[i]=helper.Get_Block();}

    INTERVAL<int> index_range_u[d] = {INTERVAL<int>(0,index_map_u[0]->next_index)};
    for(int i=1;i<d;i++) index_range_u[i]=INTERVAL<int>(index_range_u[i-1].max_corner,index_range_u[i-1].max_corner+index_map_u[i]->next_index);
    INTERVAL<int> index_range_p(index_range_u[d-1].max_corner,index_range_u[d-1].max_corner+index_map_p.next_index);
    INTERVAL<int> index_range_q[d] = {INTERVAL<int>(index_range_p.max_corner,index_range_p.max_corner+object.mesh.elements.m)};
    for(int i=1;i<d;i++) index_range_q[i]=INTERVAL<int>(index_range_q[i-1].max_corner,index_range_q[i-1].max_corner+object.mesh.elements.m);
    int finish=index_range_q[d-1].max_corner;

    printf("\n");
    for(int i=0;i<d;i++) printf("%c [%i %i) ", "uvw"[i], index_range_u[i].min_corner, index_range_u[i].max_corner);
    printf("p [%i %i) ", index_range_p.min_corner, index_range_p.max_corner);
    for(int i=0;i<d;i++) printf("%cq [%i %i) ", "uvw"[i], index_range_q[i].min_corner, index_range_q[i].max_corner);
    printf("\n");

    for(int i=0;i<d;i++)
        for(int j=0;j<d;j++)
            helper.Shift(index_range_u[i].min_corner,index_range_u[j].min_corner,block_uu[i][j]);
    for(int i=0;i<d;i++)
        helper.Shift(index_range_u[i].min_corner,index_range_p.min_corner,block_p[i]);
    for(int i=0;i<d;i++)
        helper.Shift(index_range_u[i].min_corner,index_range_q[i].min_corner,block_q[i]);

    helper.Add_Transpose(INTERVAL<int>(start_off_diag,helper.data.m));

    SPARSE_MATRIX_FLAT_MXN<T> matrix;
    helper.Set_Matrix(finish,finish,matrix,1e-14);

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
            null_p(j)=object.Get_Element(j-index_range_q[i].min_corner).Normal()(i)*units(j);}

    matrix.Times(null_p,z_p);

    LOG::cout<<z_p<<std::endl;

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
};

int main(int argc,char* argv[])
{
    Integration_Test<VECTOR<double,3> >(argc,argv);

    return 0;
}
