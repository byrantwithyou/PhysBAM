//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Finite_Elements/TRIPLE_JUNCTION_CORRECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>

using namespace PhysBAM;

typedef double T;
//typedef VECTOR<T,d> TV;
//typedef VECTOR<int,d> TV_INT;
std::string output_directory="output";
typedef float RW;
VECTOR<T,3> colors[3]={VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,1,0),VECTOR<T,3>(0,0,1)};

template<class T>
T Weight_Function(T A) {return sqrt(abs(A))*sign(A);}

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class T, int d>
void Dump_Frame(const ARRAY<T,FACE_INDEX<d> >& u,const char* title)
{
    typedef VECTOR<T,d> TV;
    static int frame=0;
    if(frame==0){
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",*Global_Grid<TV>());}

    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<TV>().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

template<class T,class TV>
void Flush_Frame(const char* title)
{
    Dump_Frame(ARRAY<T,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),title);
}

template<class T,class TV_INT>
void Dump_Interface(ARRAY_VIEW<T,TV_INT> p,const ARRAY<VECTOR<TV_INT,3> >& stencils,const VECTOR<T,3>& col)
{
    T default_phi=(T)1.12349871352389e-10;
    typedef VECTOR<T,TV_INT::m> TV;
    const GRID<TV>& grid=*Global_Grid<TV>();
    for(int i=0;i<stencils.m;i++){
        VECTOR<T,3> phi(p.Subset(stencils(i)));
        if(phi.Contains(default_phi)) continue;
        for(int j=0;j<3;j++){
            int k=(j+1)&1;
            int m=3-k-j;
            if((phi(j)>0)==(phi(k)>0)) continue;
            if((phi(j)>0)==(phi(m)>0)) continue;
            TV Xj=grid.Node(stencils(i)(j)),Xk=grid.Node(stencils(i)(k)),Xm=grid.Node(stencils(i)(m));
            TV X=Xj+phi(j)/(phi(j)-phi(k))*(Xk-Xj);
            TV Y=Xj+phi(j)/(phi(j)-phi(m))*(Xm-Xj);
            Add_Debug_Object(VECTOR<TV,2>(X,Y),col);}}
}

template<class T,class TV_INT>
void Dump_Interface(ARRAY_VIEW<T,TV_INT> p,const ARRAY<VECTOR<TV_INT,4> >& stencils,const VECTOR<T,3>& col)
{
    typedef VECTOR<T,TV_INT::m> TV;
    const GRID<TV>& grid=*Global_Grid<TV>();
    for(int i=0;i<stencils.m;i++){
        VECTOR<T,4> phi(p.Subset(stencils(i)));
        for(int o=0;o<4;o++){
            int k=(o+1)&3;
            int m=(o+2)&3;
            int n=(o+3)&3;
            if((phi(o)>0)==(phi(k)>0)) continue;
            if((phi(o)>0)==(phi(m)>0)) continue;
            if((phi(o)>0)==(phi(n)>0)) continue;
            TV Xo=grid.Node(stencils(i)(o)),Xk=grid.Node(stencils(i)(k));
            TV Xm=grid.Node(stencils(i)(m)),Xn=grid.Node(stencils(i)(n));
            TV X=Xo+phi(o)/(phi(o)-phi(k))*(Xk-Xo);
            TV Y=Xo+phi(o)/(phi(o)-phi(m))*(Xm-Xo);
            TV Z=Xo+phi(o)/(phi(o)-phi(n))*(Xn-Xo);
            if((phi(o)>0)^(TETRAHEDRON<T>::Signed_Volume(Xo,Xk,Xm,Xn)>0))
                exchange(X,Y);
            Add_Debug_Object(VECTOR<TV,3>(Y,X,Z),col,col);}
        int k=0;
        for(int m=1;m<4;m++)
            for(int n=1;n<4;n++){
                int o=6-k-m-n;
                if(m==n || o<=n) continue;
                if((phi(k)>0)==(phi(n)>0)) continue;
                if((phi(k)>0)==(phi(o)>0)) continue;
                if((phi(m)>0)==(phi(n)>0)) continue;
                if((phi(m)>0)==(phi(o)>0)) continue;
                TV Xo=grid.Node(stencils(i)(o)),Xk=grid.Node(stencils(i)(k));
                TV Xm=grid.Node(stencils(i)(m)),Xn=grid.Node(stencils(i)(n));
                TV W=Xk+phi(k)/(phi(k)-phi(n))*(Xn-Xk);
                TV X=Xk+phi(k)/(phi(k)-phi(o))*(Xo-Xk);
                TV Y=Xm+phi(m)/(phi(m)-phi(n))*(Xn-Xm);
                TV Z=Xm+phi(m)/(phi(m)-phi(o))*(Xo-Xm);
                if((phi(o)>0)^(TETRAHEDRON<T>::Signed_Volume(Xo,Xk,Xm,Xn)>0))
                    exchange(X,Y);
                    Add_Debug_Object(VECTOR<TV,3>(X,W,Y),col,col);
                    Add_Debug_Object(VECTOR<TV,3>(X,Y,Z),col,col);}}
}

template<class T>
void Evolve(VECTOR<VECTOR<T,3>,3> phi)
{
    for(int t=0;t<20;t++){
        T A=Bad_Fraction(phi);
        LOG::cout<<A<<std::endl;
        phi-=(T).9*Weight_Function(A)+VECTOR<T,3>();}
}

template<class TV,class T,class TV_INT>
void Evolve_Step(GRID<TV>& grid,ARRAY_VIEW<T,TV_INT> q[3],const ARRAY_VIEW<T,TV_INT> p[3],const ARRAY<VECTOR<TV_INT,TV_INT::m+1> >& stencils)
{
    for(int i=0;i<3;i++) q[i].Fill(0);
    for(int i=0;i<stencils.m;i++){
        VECTOR<VECTOR<T,TV::m+1>,3> phi;
        for(int c=0;c<3;c++)
            phi(c)=VECTOR<T,TV::m+1>(p[c].Subset(stencils(i)));
        T A=Bad_Fraction(phi);
        for(int c=0;c<3;c++)
            q[c].Subset(stencils(i))+=A;}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<3;c++)
            q[c](it.index)=p[c](it.index)-(T).25*Weight_Function(q[c](it.index)*grid.dX.Product());
}

template<class TV,class TV_INT>
void Initialize(const GRID<TV>& grid,ARRAY<VECTOR<TV_INT,3> >& stencils,int resolution)
{
    GRID<TV> grid_ta(TV_INT()+resolution,RANGE<TV>::Unit_Box()*resolution,false);
    TRIANGULATED_AREA<T> ta;
    ta.Initialize_Square_Mesh_And_Particles(grid_ta,true);

    stencils.Resize(ta.mesh.elements.m);
    for(int i=0;i<ta.mesh.elements.m;i++)
        for(int v=0;v<TV::m+1;v++)
            stencils(i)(v)=TV_INT(rint(ta.particles.X(ta.mesh.elements(i)(v))));
}

template<class TV,class TV_INT>
void Initialize(const GRID<TV>& grid,ARRAY<VECTOR<TV_INT,4> >& stencils,ARRAY<T,TV_INT> p[3],int resolution)
{
    GRID<TV> grid_tv(TV_INT()+resolution,RANGE<TV>::Unit_Box()*resolution,false);
    TETRAHEDRALIZED_VOLUME<T> tv;
    tv.Initialize_Cube_Mesh_And_Particles(grid_tv);

    stencils.Resize(tv.mesh.elements.m);
    for(int i=0;i<tv.mesh.elements.m;i++)
        for(int v=0;v<TV::m+1;v++)
            stencils(i)(v)=TV_INT(rint(tv.particles.X(tv.mesh.elements(i)(v))));

    PLANE<T> planes[3];
    planes[0]=PLANE<T>(TV(.31,.31,0).Normalized(),TV(.51,.51,.5));
    planes[1]=PLANE<T>(-TV(.31,.11,.1).Normalized(),TV(.51,.41,.5));
    planes[2]=PLANE<T>(-TV(.11,.41,0).Normalized(),TV(.41,.61,.5));

    for(int c=0;c<3;c++){
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next())
            p[c](it.index)=planes[c].Signed_Distance(it.Location());}
}

ROTATION<VECTOR<T,2> > rot(ROTATION<VECTOR<T,2> >::From_Angle(-2*(T)pi/3));
ROTATION<VECTOR<T,2> > rot2(ROTATION<VECTOR<T,2> >::From_Angle((T)pi/10));

template <class TV> typename TV::SCALAR
Levelset(TV X)
{
    typedef typename TV::SCALAR T;
    T angle=atan2(X.y,X.x);
    if(angle>=-2*(T)pi/3 && angle<=-(T)pi/2) return X.Magnitude();
    return max(-X.y,rot.Rotate(X).y);
}

template <class TV> typename TV::SCALAR
Levelset(TV X,int i)
{
    typedef typename TV::SCALAR T;
    const T R=0.35;
    T r=X.Magnitude();
    switch(i){
        case 0:{
            if(X.y>=0){
                if(r<=R) return max(-X.y,r-R);
                else return r-R;}
            else{
                if(X.x<=-R) return (X-TV(-R,0)).Magnitude();
                if(X.x>=R) return (X-TV(R,0)).Magnitude();
                return -X.y;}
        } break;
        case 1:{
            return R-r;
        } break;
        default: PHYSBAM_FATAL_ERROR();
    }
}


template<class TV>
void Compute(PARSE_ARGS& parse_args)
{
    int resolution=8;
    parse_args.Add("-resolution",&resolution,"res","grid resolution");
    parse_args.Parse();

    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    Get_Debug_Particles<TV>();
    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    GRID<TV> grid(TV_INT()+resolution,RANGE<TV>::Centered_Box());
    Global_Grid(&grid);

    ARRAY<VECTOR<TV_INT,TV::m+1> > stencils;
    Initialize(grid,stencils,resolution);

    ARRAY<ARRAY<T,TV_INT> > color_phi(3);
    for(int i=0;i<3;i++)
        color_phi(i).Resize(grid.Node_Indices(3));

    if(1) // triple
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next()){
        TV X=rot2.Rotate(it.Location())+(T).0001;
        color_phi(0)(it.index)=Levelset(X);
        color_phi(1)(it.index)=Levelset(rot.Rotate(X+.01));
        color_phi(2)(it.index)=Levelset(rot.Rotate(rot.Rotate(X)));}
    else // circle
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next()){
        TV X=it.Location();
        color_phi(0)(it.index)=Levelset(X,0);
        color_phi(1)(it.index)=Levelset(-X,0);
        color_phi(2)(it.index)=Levelset(X,1);}

    for(int i=0;i<3;i++){
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next()){
            T p=color_phi(i)(it.index);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        Flush_Frame<T,TV>("Initial level set");}

    ARRAY<ARRAY<ARRAY<T,TV_INT> > > pairwise_phi;
    TRIPLE_JUNCTION_CORRECTION<TV> tjc(grid,color_phi,3);
    tjc.Compute_Pairwise_Level_Set_Data();

    typedef typename MARCHING_CUBES_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    HASHTABLE<TV_INT,CELL_ELEMENTS> index_to_cell_data;
    tjc.Cut_Interface(index_to_cell_data);

    VECTOR<T,3> color_map[4]={VECTOR<T,3>(0,0.7,0),VECTOR<T,3>(0.8,0.8,0),VECTOR<T,3>(0,0.4,1),VECTOR<T,3>(0.8,0.2,0)};
    T sep=(T).08;
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(index_to_cell_data);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.y>=0){
                if(V.color_pair.y>=0) Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.y]);
                if("Alexey was here") Add_Debug_Object(V.face.X+V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.x]);}
            else if(V.color_pair.x>=0) Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.x]);}}
    Flush_Frame<T,TV>("interfaces");
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(index_to_cell_data);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<BOUNDARY_ELEMENT>& boundary_elements=cell_elements.boundary;
        for(int i=0;i<boundary_elements.m;i++){
            const BOUNDARY_ELEMENT& V=boundary_elements(i);
            Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color]);}}
    Flush_Frame<T,TV>("boundarys");

    LOG::Finish_Logging();
}

int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3d");
    parse_args.Parse(true);

//    if(use_3d) Compute<VECTOR<T,3> >(parse_args);
//    else 
    Compute<VECTOR<T,2> >(parse_args);
}

template void Dump_Interface<double,VECTOR<int,2> >(ARRAY_VIEW<double,VECTOR<int,2> >,ARRAY<VECTOR<VECTOR<int,2>,3>,int> const&,VECTOR<double,3> const&);
template void Dump_Interface<double,VECTOR<int,3> >(ARRAY_VIEW<double,VECTOR<int,3> >,ARRAY<VECTOR<VECTOR<int,3>,4>,int> const&,VECTOR<double,3> const&);
template void Dump_Interface<float,VECTOR<int,2> >(ARRAY_VIEW<float,VECTOR<int,2> >,ARRAY<VECTOR<VECTOR<int,2>,3>,int> const&,VECTOR<float,3> const&);
template void Dump_Interface<float,VECTOR<int,3> >(ARRAY_VIEW<float,VECTOR<int,3> >,ARRAY<VECTOR<VECTOR<int,3>,4>,int> const&,VECTOR<float,3> const&);
template void Flush_Frame<double,VECTOR<double,2> >(char const*);
template void Flush_Frame<double,VECTOR<double,3> >(char const*);
template void Flush_Frame<float,VECTOR<float,2> >(char const*);
template void Flush_Frame<float,VECTOR<float,3> >(char const*);
template GRID<VECTOR<double,2> >* Global_Grid<VECTOR<double,2> >(GRID<VECTOR<double,2> >*);
template GRID<VECTOR<float,2> >* Global_Grid<VECTOR<float,2> >(GRID<VECTOR<float,2> >*);
