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
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
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
    typedef VECTOR<T,TV_INT::m> TV;
    const GRID<TV>& grid=*Global_Grid<TV>();
    for(int i=0;i<stencils.m;i++){
        VECTOR<T,3> phi(p.Subset(stencils(i)));
        if(phi.Contains((T)pi/1024)) continue;
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

template<class T,int dp1>
T Bad_Fraction(const VECTOR<VECTOR<T,dp1>,3>& phi)
{
    int pos_mask=0,neg_mask=0;
    for(int c=0;c<3;c++)
        for(int j=0;j<dp1;j++){
            if(phi(c)(j)>0) pos_mask|=1<<c;
            else if(phi(c)(j)<0) neg_mask|=1<<c;}

    if(pos_mask==7){
        if(neg_mask==0) return 1;}
    else if(neg_mask==7){
        if(pos_mask==0) return -1;}
    else return 0;

    for(int c=0;c<3;c++)
        for(int j=0;j<dp1;j++)
            if(T pj=phi(c)(j))
                for(int k=j+1;k<dp1;k++)
                    if(T pk=phi(c)(k))
                        if((pj>0)!=(pk>0)){
                            VECTOR<VECTOR<T,dp1>,3> phi0=phi,phi1=phi;
                            T th=pj/(pj-pk);
                            for(int e=0;e<3;e++){
                                T ej=phi(e)(j),ek=phi(e)(k),el=ej+th*(ek-ej);
                                phi0(e)(j)=ej;
                                phi0(e)(k)=el;
                                phi1(e)(j)=el;
                                phi1(e)(k)=ek;}
                            phi0(c)(k)=phi1(c)(j)=0;
                            return th*Bad_Fraction(phi0)+(1-th)*Bad_Fraction(phi1);}

    PHYSBAM_FATAL_ERROR();
}

template<class T>
struct PAIRWISE_LEVEL_SET_DATA
{
    T phi;
    int valid_flags;
    VECTOR<short,2> trust;

    PAIRWISE_LEVEL_SET_DATA()
        :phi(FLT_MAX),valid_flags(0),trust(-1,-1)
    {}
};

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

template<class T,class TV,class TV_INT>
void Compute_Pairwise_Level_Set_Data(const GRID<TV>& grid,const ARRAY<ARRAY<T,TV_INT> >& phi,int ghost,const ARRAY<VECTOR<TV_INT,TV::m+1> >& stencils,ARRAY<ARRAY<ARRAY<T,TV_INT> > >& pairwise_phi)
{
    T trust_buffer=grid.dX.Max(),valid_width=ghost*grid.dX.Max(),extent=4*valid_width;
    int extrap_width=2*ghost+3;
    ARRAY<PAIRWISE_LEVEL_SET_DATA<T>,TV_INT> pairwise_data(grid.Node_Indices(ghost));

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA<T>& data=pairwise_data(it.index);
        VECTOR<T,3> min(FLT_MAX,FLT_MAX,FLT_MAX);
        for(int i=0;i<phi.m;i++){
            T p=phi(i)(it.index);
            if(p<min(0)){min=min.Remove_Index(2).Insert(p,0);data.trust=data.trust.Remove_Index(1).Insert(i,0);}
            else if(p<min(1)){min=min.Remove_Index(2).Insert(p,1);data.trust(1)=i;}
            else if(p<min(2)) min(2)=p;
            if(p<=extent) data.valid_flags|=1<<i;}
        min-=(min(0)+min(1))/2; // In case of boundary conditions
        if(min(1)>valid_width){
            data=PAIRWISE_LEVEL_SET_DATA<T>();
            continue;}
        if(min(1)>min(2)-trust_buffer){
            data.trust=VECTOR<short,2>(-1,-1);
            continue;}
        data.phi=min.x;
        if(data.trust.y<data.trust.x){
            exchange(data.trust.x,data.trust.y);
            data.phi=min.y;}}

    pairwise_phi.Resize(phi.m);
    for(int i=0;i<phi.m;i++){
        pairwise_phi(i).Resize(phi.m);
        for(int j=i+1;j<phi.m;j++){
            pairwise_phi(i)(j).Resize(grid.Node_Indices(extrap_width),false);
            pairwise_phi(i)(j).Fill((T)pi/1024);}}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA<T>& data=pairwise_data(it.index);
        if(data.trust.x>=0){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(data.trust.y==1,data.trust.Sum()==2,data.trust.x==1)/(1+(data.phi<0)));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(data.phi));
            pairwise_phi(data.trust.x)(data.trust.y)(it.index)=data.phi;}}

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),stencils,VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),stencils,VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),stencils,VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("Initial pairwise level sets");

    for(int i=0;i<phi.m;i++)
        for(int j=i+1;j<phi.m;j++){
            struct LOCAL_MASK:public EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::MASK
            {
                const ARRAY<PAIRWISE_LEVEL_SET_DATA<T>,TV_INT>& pairwise_data;
                VECTOR<short,2> current_pair;
                LOCAL_MASK(const ARRAY<PAIRWISE_LEVEL_SET_DATA<T>,TV_INT>& pairwise_data,const VECTOR<short,2>& current_pair)
                    :pairwise_data(pairwise_data),current_pair(current_pair)
                {}
                bool Inside(const TV_INT& index) PHYSBAM_OVERRIDE
                {return pairwise_data(index).trust==current_pair;}
            } inside_mask(pairwise_data,VECTOR<short,2>(i,j));
            EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Node(grid,inside_mask,extrap_width,pairwise_phi(i)(j),3,extrap_width);}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        T p=pairwise_phi(0)(1)(it.index);
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    Flush_Frame<T,TV>("level set 01");

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),stencils,VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),stencils,VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),stencils,VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("Extrapolated pairwise level sets");

    for(int t=0;t<20;t++){
        HASHTABLE<TRIPLE<int,int,TV_INT>,T> total_size;
        for(int i=0;i<stencils.m;i++){
            int mask=~0;
            for(int j=0;j<stencils(i).m;j++)
                mask&=pairwise_data(stencils(i)(j)).valid_flags;
            for(int a=0;a<phi.m;a++)
                if(mask&(1<<a))
                    for(int b=a+1;b<phi.m;b++)
                        if(mask&(1<<b))
                            for(int c=b+1;c<phi.m;c++)
                                if(mask&(1<<c)){
                                    VECTOR<T,TV::m+1> ab(pairwise_phi(a)(b).Subset(stencils(i)));
                                    VECTOR<T,TV::m+1> bc(pairwise_phi(b)(c).Subset(stencils(i)));
                                    VECTOR<T,TV::m+1> ca(-pairwise_phi(a)(c).Subset(stencils(i)));
                                    T A=Bad_Fraction(VECTOR<VECTOR<T,TV::m+1>,3>(ab,bc,ca));
                                    for(int j=0;j<stencils(i).m;j++){
                                        total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,b,stencils(i)(j)))+=A;
                                        total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(b,c,stencils(i)(j)))+=A;
                                        total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,c,stencils(i)(j)))-=A;}}}
        for(typename HASHTABLE<TRIPLE<int,int,TV_INT>,T>::ITERATOR it(total_size);it.Valid();it.Next())
            pairwise_phi(it.Key().x)(it.Key().y)(it.Key().z)-=(T).25*Weight_Function(it.Data()*grid.dX.Product());}
}

template<class TV,class TV_INT>
void Initialize(const GRID<TV>& grid,ARRAY<VECTOR<TV_INT,3> >& stencils,ARRAY<T,TV_INT> p[3],int resolution)
{
    GRID<TV> grid_ta(TV_INT()+resolution,RANGE<TV>::Unit_Box()*resolution,false);
    TRIANGULATED_AREA<T> ta;
    ta.Initialize_Square_Mesh_And_Particles(grid_ta,true);

    stencils.Resize(ta.mesh.elements.m);
    for(int i=0;i<ta.mesh.elements.m;i++)
        for(int v=0;v<TV::m+1;v++)
            stencils(i)(v)=TV_INT(rint(ta.particles.X(ta.mesh.elements(i)(v))));

    LINE_2D<T> lines[3];
    lines[0]=LINE_2D<T>(TV(.31,.31).Normalized(),TV(.51,.51));
    lines[1]=LINE_2D<T>(-TV(.31,.11).Normalized(),TV(.51,.41));
    lines[2]=LINE_2D<T>(-TV(.11,.41).Normalized(),TV(.41,.61));

    for(int c=0;c<3;c++){
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next())
            p[c](it.index)=lines[c].Signed_Distance(it.Location());}
    
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

template <class TV> typename TV::SCALAR
Levelset(TV X)
{
    typedef typename TV::SCALAR T;
    T angle=atan2(X.y,X.x);
    if(angle>=-2*(T)pi/3 && angle<=-(T)pi/2) return X.Magnitude();
    return max(-X.y,rot.Rotate(X).y);
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

    ARRAY<T,TV_INT> p[3],q[3];
    for(int c=0;c<3;c++){
        p[c].Resize(grid.Node_Indices(3));
        q[c].Resize(grid.Node_Indices(3));}
    ARRAY_VIEW<T,TV_INT> pv[3]={ARRAY_VIEW<T,TV_INT>(p[0]),ARRAY_VIEW<T,TV_INT>(p[1]),ARRAY_VIEW<T,TV_INT>(p[2])};

    ARRAY<VECTOR<TV_INT,TV::m+1> > stencils;
    Initialize(grid,stencils,p,resolution);

    for(int c=0;c<3;c++)
        Dump_Interface(pv[c],stencils,colors[c]);
    Flush_Frame<T,TV>("Initial");

    ARRAY<ARRAY<T,TV_INT> > color_phi(3);
    for(int i=0;i<3;i++)
        color_phi(i).Resize(grid.Node_Indices(3));

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next()){
        TV X=it.Location()+(T).0001;
        color_phi(0)(it.index)=Levelset(X);
        color_phi(1)(it.index)=Levelset(rot.Rotate(X));
        color_phi(2)(it.index)=Levelset(rot.Rotate(rot.Rotate(X)));}

    for(int i=0;i<3;i++){
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,3);it.Valid();it.Next()){
            T p=color_phi(i)(it.index);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
        Flush_Frame<T,TV>("Initial level set");}

    ARRAY<ARRAY<ARRAY<T,TV_INT> > > pairwise_phi;
    Compute_Pairwise_Level_Set_Data(grid,color_phi,3,stencils,pairwise_phi);

    // for(int i=0;i<500;i++){
    //     Evolve_Step(grid,qv,pv,stencils);

    //     for(int c=0;c<3;c++)
    //         pv[c].Exchange(qv[c]);

    //     for(int c=0;c<3;c++)
    //         Dump_Interface(pv[c],stencils,colors[c]);
    //     Flush_Frame<T,TV>("Phew");}
    // Flush_Frame<T,TV>("Phew");

    LOG::Finish_Logging();
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3d");
    parse_args.Parse(true);

//    if(use_3d) Compute<VECTOR<T,3> >(parse_args);
//    else 
    Compute<VECTOR<T,2> >(parse_args);
}

