//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
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
void Dump_Interface(const ARRAY<T,TV_INT>& p,const ARRAY<VECTOR<TV_INT,3> >& stencils,const VECTOR<T,3>& col)
{
    typedef VECTOR<T,TV_INT::m> TV;
    const GRID<TV>& grid=*Global_Grid<TV>();
    for(int i=0;i<stencils.m;i++){
        VECTOR<T,3> phi(p.Subset(stencils(i)));
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
void Dump_Interface(const ARRAY<T,TV_INT>& p,const ARRAY<VECTOR<TV_INT,4> >& stencils,const VECTOR<T,3>& col)
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
VECTOR<T,2> Meet_Phi(const VECTOR<VECTOR<T,4>,2>& phi)
{
    typedef VECTOR<T,2> TV;
    QUADRATIC<T> quad(-phi(0)(0)*phi(1)(3)+phi(0)(3)*phi(1)(0)+phi(0)(1)*phi(1)(2)-phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(3)-phi(0)(2)*phi(1)(1)+phi(0)(0)*phi(1)(1)-phi(0)(3)*phi(1)(2),-2*phi(0)(0)*phi(1)(1)+phi(0)(0)*phi(1)(3)-phi(0)(1)*phi(1)(2)+2*phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(1)-phi(0)(3)*phi(1)(0),-phi(0)(1)*phi(1)(0)+phi(0)(0)*phi(1)(1));
    quad.Compute_Roots();
    if(quad.roots==0){
        quad.roots=1;
        quad.root1=-quad.b/(2*quad.a);}
    T ya[2]={quad.root1,quad.root2},ph[2]={0,FLT_MAX};
    TV X[2];
    for(int i=0;i<quad.roots;i++){
        X[i]=TV((-phi(1)(0)-ya[i]*phi(1)(2)+ya[i]*phi(1)(0))/(phi(1)(1)-phi(1)(0)+ya[i]*phi(1)(3)-ya[i]*phi(1)(2)-ya[i]*phi(1)(1)+ya[i]*phi(1)(0)),ya[i]);
        ph[i]=phi(0)(0)+(phi(0)(1)-phi(0)(0))*X[i].x+(phi(0)(2)-phi(0)(0))*X[i].y+(phi(0)(3)-phi(0)(2)-phi(0)(1)+phi(0)(0))*X[i].x*X[i].y;}

    return X[abs(ph[1])<abs(ph[0])];
}

template<class T>
VECTOR<T,3> Meet_Phi(const VECTOR<VECTOR<T,8>,2>& phi){PHYSBAM_FATAL_ERROR();}

template<class T,int d>
VECTOR<T,(d==4?2:3)> Zero_Phi(const VECTOR<VECTOR<T,d>,3>& phi,VECTOR<T,3>& p)
{
    typedef VECTOR<T,(d==4?2:3)> TV;
    p=VECTOR<T,3>()+FLT_MAX;
    TV bestX,X;
    T pp=0;
    X=Meet_Phi(VECTOR<VECTOR<T,d>,2>(phi.x-phi.y,phi.x-phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,pp);bestX=X;}
    X=Meet_Phi(VECTOR<VECTOR<T,d>,2>(phi.x+phi.y,phi.x-phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,-pp,pp);bestX=X;}
    X=Meet_Phi(VECTOR<VECTOR<T,d>,2>(phi.x-phi.y,phi.x+phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,-pp);bestX=X;}
    X=Meet_Phi(VECTOR<VECTOR<T,d>,2>(-phi.x+phi.y,-phi.x+phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(-pp,pp,pp);bestX=X;}
    return bestX;
}

template<class TV,class T,class TV_INT>
void Evolve_Step(GRID<TV>& grid,ARRAY<T,TV_INT> q[3],const ARRAY<T,TV_INT> p[3],const ARRAY<VECTOR<TV_INT,TV_INT::m+1> >& stencils)
{
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    for(int i=0;i<3;i++) q[i].Fill(0);
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        VECTOR<T,3> pp;
        VECTOR<VECTOR<T,(1<<TV::m)>,3> phi;
        for(int i=0;i<3;i++) phi(i)=VECTOR<T,(1<<TV::m)>(p[i].Subset(bits+it.index));
        // int num_change=0;
        // for(int i=0;i<3;i++) num_change+=(phi(i).Min()<0 && phi(i).Max()>0);
        // if(num_change<2) continue;
        if(!RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<VECTOR<T,(1<<TV::m)>,2>(phi.x,phi.y))) && !RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<VECTOR<T,(1<<TV::m)>,2>(phi.x,phi.z))) && !RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<VECTOR<T,(1<<TV::m)>,2>(phi.y,phi.z)))) continue;

        TV X=Zero_Phi(phi,pp);
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(pp(0)));
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        LOG::cout<<pp<<std::endl;
        for(int i=0;i<3;i++) q[i].Subset(bits+it.index)+=pp(i);}

    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<3;c++)
            q[c](it.index)=p[c](it.index)-(T).25*q[c](it.index);
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
        for(NODE_ITERATOR<TV> it(grid,3);it.Valid();it.Next())
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
        for(NODE_ITERATOR<TV> it(grid,3);it.Valid();it.Next())
            p[c](it.index)=planes[c].Signed_Distance(it.Location());}
}

template<class TV>
void Compute(PARSE_ARGS& parse_args)
{
    int resolution=8;
    parse_args.Add("-resolution",&resolution,"res","grid resolution");
    parse_args.Parse();

    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    GRID<TV> grid(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
    Global_Grid(&grid);

    ARRAY<T,TV_INT> p[3],q[3];
    for(int c=0;c<3;c++){
        p[c].Resize(grid.Node_Indices(3));
        q[c].Resize(grid.Node_Indices(3));}

    ARRAY<VECTOR<TV_INT,TV::m+1> > stencils;
    Initialize(grid,stencils,p,resolution);

    for(int c=0;c<3;c++)
        Dump_Interface(p[c],stencils,colors[c]);
    Flush_Frame<T,TV>("Initial");

    for(int i=0;i<50;i++){
        Evolve_Step(grid,q,p,stencils);

        for(int c=0;c<3;c++)
            p[c].Exchange(q[c]);

        for(int c=0;c<3;c++)
            Dump_Interface(p[c],stencils,colors[c]);
        Flush_Frame<T,TV>("Phew");}
    Flush_Frame<T,TV>("Phew");

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

