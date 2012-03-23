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
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

using namespace PhysBAM;

typedef double T;
static const int d=2;
typedef VECTOR<T,d> TV;
typedef VECTOR<int,d> TV_INT;
typedef float RW;

// * Handle interface
// * Handle other boundary conditions
// * Handle 

std::string output_directory="output";
DEBUG_PARTICLES<TV> debug_particles;
// template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);
// template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
GRID<TV> *global_grid;

void Dump_Frame(const ARRAY<T,FACE_INDEX<d> >& u,const char* title)
{
    typedef VECTOR<T,d> TV;
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    debug_particles.Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

void Flush_Frame()
{
    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(*global_grid),"flush");
}

void Dump_Frame(const VECTOR_ND<T>& v,const CELL_MAPPING<TV>& index_map_u,const CELL_MAPPING<TV>& index_map_v,const CELL_MAPPING<TV>& index_map_p,const SEGMENTED_CURVE_2D<T>& curve,const char* title)
{
    int start_v=index_map_u.next_index;
    int start_p=index_map_v.next_index+start_v;
    int start_uq=index_map_p.next_index+start_p;
    int start_vq=curve.mesh.elements.m+start_uq;
//    int total=curve.mesh.elements.m+start_vq;
    const GRID<TV>& grid=index_map_p.grid;

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(1e-11,10,true,true,true);

    char buff[100];
    const CELL_MAPPING<TV>* cell_mappings[2]={&index_map_u, &index_map_v};
    for(int dir=0;dir<TV::m;dir++)
    {
        sprintf(buff, title, dir?'-':'+');
        ARRAY<T,FACE_INDEX<TV::m> > error(grid);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            int index=cell_mappings[it.Axis()]->Get_Index_Fixed(it.index,dir);
            if(index>=0)
                error(it.Full_Index())=v(index+(it.Axis()?start_v:0));}

        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
            int index=index_map_p.Get_Index_Fixed(it.index,dir);
            if(index>=0)
                Add_Debug_Particle(grid.Center(it.index),color_map(v(start_p+index)));}

        for(int i=0;i<curve.mesh.elements.m;i++){
            Add_Debug_Particle(curve.Get_Element(i).Center(),VECTOR<T,3>(0,1,1));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(v(start_uq+i),v(start_vq+i)));}

        Dump_Frame(error,buff);}    
}

void Integration_Test(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-viscosity",1,"viscosity");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-s",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Parse(argc,argv);
    T m=parse_args.Get_Double_Value("-m");
    T s=parse_args.Get_Double_Value("-s");
    T kg=parse_args.Get_Double_Value("-kg");
    T mu=parse_args.Get_Double_Value("-viscosity")*kg/s;
    (void)m;

    GRID<TV> grid(TV_INT()+8,RANGE<TV>(TV(),TV()+1)*m,true);
    global_grid=&grid;
    GRID<TV> coarse_grid(grid.counts,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    SPHERE<TV> sphere(TV()+(T).5*m,(T).3*m);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=sphere.Signed_Distance(it.Location());

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    BASIS_STENCIL_UNIFORM<TV> p_stencil(grid.dX),u_stencil(grid.dX),v_stencil(grid.dX);
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Dice_Stencil();
    u_stencil.Set_Face(0);
    u_stencil.Set_Multilinear_Stencil();
    u_stencil.Dice_Stencil();
    v_stencil.Set_Face(1);
    v_stencil.Set_Multilinear_Stencil();
    v_stencil.Dice_Stencil();
    BASIS_STENCIL_UNIFORM<TV> udx_stencil(u_stencil),udy_stencil(u_stencil),vdx_stencil(v_stencil),vdy_stencil(v_stencil);
    udx_stencil.Differentiate(0);
    udx_stencil.Dice_Stencil();
    udy_stencil.Differentiate(1);
    udy_stencil.Dice_Stencil();
    vdx_stencil.Differentiate(0);
    vdx_stencil.Dice_Stencil();
    vdy_stencil.Differentiate(1);
    vdy_stencil.Dice_Stencil();

    SEGMENTED_CURVE_2D<T> curve;
    MARCHING_CUBES<TV>::Create_Surface(curve,coarse_grid,phi);
    BASIS_STENCIL_BOUNDARY_UNIFORM<TV> q_stencil(curve);

    for(int i=0;i<curve.mesh.elements.m;i++){
        Add_Debug_Particle(curve.Get_Element(i).X(0),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,curve.Get_Element(i).X(1)-curve.Get_Element(i).X(0));}
    Flush_Frame();

    CELL_MAPPING<TV> index_map_p(grid),index_map_u(grid),index_map_v(grid);
    index_map_p.periodic.Fill(true);
    index_map_u.periodic.Fill(true);
    index_map_v.periodic.Fill(true);

    SYSTEM_MATRIX_HELPER<T> helper;
    RANGE<TV_INT> boundary_conditions;
    boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);
    boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udx_stencil,udx_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"udx udx");
    helper.Scale(2*mu);
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udy_stencil,udy_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"udy udy");
    helper.Scale(mu);
    helper.start=0;
    INTERVAL<int> block_uu=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdx_stencil,vdx_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"vdx vdx");
    helper.Scale(2*mu);
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdy_stencil,vdy_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"vdy vdy");
    helper.Scale(mu);
    helper.start=block_uu.max_corner;
    INTERVAL<int> block_vv=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udy_stencil,vdx_stencil,index_map_u,index_map_v,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"udy vdx");
    helper.Scale(mu);
    INTERVAL<int> block_uv=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udx_stencil,p_stencil,index_map_u,index_map_p,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"udx p");
    INTERVAL<int> block_up=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdy_stencil,p_stencil,index_map_v,index_map_p,phi).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"vdy p");
    INTERVAL<int> block_vp=helper.Get_Block();

    BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>(grid,u_stencil,q_stencil,index_map_u).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"u q");
    INTERVAL<int> block_uq=helper.Get_Block();

    BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>(grid,v_stencil,q_stencil,index_map_v).Compute_Matrix(helper);
//    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"v q");
    INTERVAL<int> block_vq=helper.Get_Block();

    int start_v=index_map_u.next_index;
    int start_p=index_map_v.next_index+start_v;
    int start_uq=index_map_p.next_index+start_p;
    int start_vq=curve.mesh.elements.m+start_uq;
    int total=curve.mesh.elements.m+start_vq;

    printf("u %i v %i p %i uq %i vq %i end %i\n", 0, start_v, start_p, start_uq, start_vq, total);

    printf("ind %i %i\n", index_map_u.Get_Index_Fixed(TV_INT(1,4),0), index_map_u.Get_Index_Fixed(TV_INT(7,4),0));
    printf("ind %i %i\n", index_map_u.Get_Index_Fixed(TV_INT(1,4),1), index_map_u.Get_Index_Fixed(TV_INT(7,4),1));

    helper.Shift(start_v,start_v,block_vv);
    helper.Shift(0,start_v,block_uv);
    helper.Shift(0,start_p,block_up);
    helper.Shift(start_v,start_p,block_vp);
    helper.Shift(0,start_uq,block_uq);
    helper.Shift(start_v,start_vq,block_vq);

    helper.Add_Transpose(INTERVAL<int>(block_uv.min_corner,block_vq.max_corner));

    SPARSE_MATRIX_FLAT_MXN<T> matrix;
    helper.Set_Matrix(total,total,matrix,1e-14);

    VECTOR_ND<T> units(total);
    VECTOR_ND<T> null_u(total),null_v(total),null_p(total),z_u(total),z_v(total),z_p(total);
    for(int i=0;i<start_v;i++){null_u(i)=1;units(i)=m/s;}
    for(int i=start_v;i<start_p;i++){null_v(i)=1;units(i)=m/s;}
    for(int i=start_p;i<start_uq;i++){null_p(i)=1;units(i)=kg/(s*s);}
    for(int i=start_uq;i<start_vq;i++){null_p(i)=curve.Get_Element(i-start_uq).Normal().x;units(i)=kg/(s*s);}
    for(int i=start_vq;i<total;i++){null_p(i)=curve.Get_Element(i-start_vq).Normal().y;units(i)=kg/(s*s);}
    null_u*=units;
    null_v*=units;
    null_p*=units;

    matrix.Times(null_u,z_u);
    matrix.Times(null_v,z_v);
    matrix.Times(null_p,z_p);
    z_u*=units;z_u/=kg*m*m/s/s/s;
    z_v*=units;z_v/=kg*m*m/s/s/s;
    z_p*=units;z_p/=kg*m*m/s/s/s;

    LOG::cout<<z_u<<std::endl;
    LOG::cout<<z_v<<std::endl;
    LOG::cout<<z_p<<std::endl;

    for(int i=0;i<curve.mesh.elements.m;i++) Add_Debug_Particle(curve.particles.X(curve.mesh.elements(i).x),VECTOR<T,3>(1,0,0));
    for(int i=0;i<curve.mesh.elements.m;i++) Add_Debug_Particle(curve.particles.X(curve.mesh.elements(i).x)*.7+curve.particles.X(curve.mesh.elements(i).y)*.3,VECTOR<T,3>(0,1,0));
    Dump_Frame(ARRAY<T,FACE_INDEX<d> >(grid),"setup debug particles");
    Dump_Frame(null_u,index_map_u,index_map_v,index_map_p,curve,"null u %c");
    Dump_Frame(null_v,index_map_u,index_map_v,index_map_p,curve,"null v %c");
    Dump_Frame(null_p,index_map_u,index_map_v,index_map_p,curve,"null p %c");
    Dump_Frame(z_u,index_map_u,index_map_v,index_map_p,curve,"error u %c");
    Dump_Frame(z_v,index_map_u,index_map_v,index_map_p,curve,"error v %c");
    Dump_Frame(z_p,index_map_u,index_map_v,index_map_p,curve,"error p %c");

    OCTAVE_OUTPUT<T>("M.txt").Write("M",matrix);
};

int main(int argc,char* argv[])
{
    Integration_Test(argc,argv);

    return 0;
}
