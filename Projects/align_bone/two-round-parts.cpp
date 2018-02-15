//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef double RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,TV::m> TV_INT;

void Compute_Samples(ARRAY<int>& sample_points,ARRAY_VIEW<const TV> X,const KD_TREE<TV>& kd_tree,T dist)
{
    ARRAY<int> todo(IDENTITY_ARRAY<>(X.m));
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(123);
    random.Random_Shuffle(todo);
    ARRAY<bool> done(X.m);
    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    for(int i=0;i<todo.m;i++)
    {
        int p=todo(i);
        if(done(p)) continue;
        sample_points.Append(p);
            
        kd_tree.Find_Points_Within_Radius(X(p),sqr(dist),points_found,
            distance_squared_of_points_found,X);
        done.Subset(points_found).Fill(true);
//            Add_Debug_Particle(X,TV(0,1,0));
//            Add_Debug_Object(VECTOR<TV,2>(X,X+normals(p)),TV(1,0,0));
            
        points_found.Remove_All();
        distance_squared_of_points_found.Remove_All();
    }
    printf("samples %i\n",sample_points.m);
}

void Compute_Normals(ARRAY<TV>& normals,ARRAY_VIEW<const TV> X,const KD_TREE<TV>& kd_tree,T dist)
{
    normals.Resize(X.m);

    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    for(int i=0;i<X.m;i++)
    {
        kd_tree.Find_Points_Within_Radius(X(i),sqr(dist),points_found,
            distance_squared_of_points_found,X);

        TV A=X.Subset(points_found).Average();
        SYMMETRIC_MATRIX<T,TV::m> C;
        for(int j=0;j<points_found.m;j++)
            C+=Outer_Product(X(points_found(j))-A);
        C/=points_found.m;

        DIAGONAL_MATRIX<T,TV::m> eigenvalues;
        MATRIX<T,TV::m> eigenvectors;
        C.Solve_Eigenproblem(eigenvalues,eigenvectors);
        int k=eigenvalues.x.Arg_Min();
        TV n=eigenvectors.Column(k);
        normals(i)=n;
        points_found.Remove_All();
        distance_squared_of_points_found.Remove_All();
    }
}



void Compute_Curvature(ARRAY<TV2>& curvatures,ARRAY<TV>& eig0,ARRAY<TV>& eig1,ARRAY_VIEW<const TV> X,const KD_TREE<TV>& kd_tree,const ARRAY<TV>& normals,T dist)
{
    curvatures.Resize(X.m);
    eig0.Resize(X.m);
    eig1.Resize(X.m);

    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    for(int i=0;i<X.m;i++)
    {
        kd_tree.Find_Points_Within_Radius(X(i),sqr(dist),points_found,
            distance_squared_of_points_found,X);
        TV n=normals(i);
        TV u=n.Unit_Orthogonal_Vector(),v=n.Cross(u);
        MATRIX<T,3,2> R(u,v);
            
        MATRIX<T,6> M;
        VECTOR<T,6> b;
        for(int j=0;j<points_found.m;j++)
        {
            TV Z=X(points_found(j))-X(i);
            TV2 Y=R.Transpose_Times(Z);
            VECTOR<T,6> c(sqr(Y.x),Y.x*Y.y,sqr(Y.y),Y.x,Y.y,1);
            M+=Outer_Product(c,c);
            b+=c*n.Dot(Z);
        }
        try
        {
            VECTOR<T,6> s=M.Householder_QR_Solve(b);
            SYMMETRIC_MATRIX<T,2> K(s(0),s(1),s(2));
            DIAGONAL_MATRIX<T,2> eig2;
            MATRIX<T,2> evec2;
            K.Solve_Eigenproblem(eig2,evec2);
            eig2=-eig2;
            curvatures(i)=eig2.x.Sorted();
            int k0=eig2.x.Arg_Min(),k1=1-k0;
            TV K0=R*evec2.Column(k0);
            TV K1=R*evec2.Column(k1);
            eig0(i)=K0;
            eig1(i)=K1;
            points_found.Remove_All();
            distance_squared_of_points_found.Remove_All();
        }
        catch(...)
        {
        }
    }
}

void Orient_Normals(ARRAY<TV>& normals,ARRAY_VIEW<const TV> X,const KD_TREE<TV>& kd_tree,T dist)
{
    ARRAY<TRIPLE<T,int,int> > edges;
    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    for(int i=0;i<X.m;i++)
    {
        kd_tree.Find_Points_Within_Radius(X(i),sqr(dist),points_found,
            distance_squared_of_points_found,X);
        for(int j=0;j<points_found.m;j++){
            int p=points_found(j);
            edges.Append({abs(normals(i).Dot(normals(p))),i,p});}
        points_found.Remove_All();
        distance_squared_of_points_found.Remove_All();
    }
    edges.Sort([](const auto& a,const auto& b){return a.x>b.x;});
    UNION_FIND<> uf(X.m);
    ARRAY<ARRAY<int> > adj(X.m);
    for(int i=0;i<edges.m;i++){
        if(uf.Find(edges(i).y)==uf.Find(edges(i).z)) continue;
        uf.Union(edges(i).y,edges(i).z);
        adj(edges(i).y).Append(edges(i).z);
        adj(edges(i).z).Append(edges(i).y);}
    ARRAY<bool> done(X.m);
    ARRAY<int> todo;
    for(int i=0;i<adj.m;i++)
        if(!done(i))
        {
            todo.Append(i);
            done(i)=true;
            while(todo.m)
            {
                int a=todo.Pop();
                for(int j=0;j<adj(a).m;j++)
                {
                    int b=adj(a)(j);
                    if(done(b)) continue;
                    if(normals(a).Dot(normals(b))<0)
                        normals(b)=-normals(b);
                    done(b)=true;
                    todo.Append(b);
                }
            }
        }
}

int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    T min_dist = 8;
    T seed_dist = 1;
    T cent_dist = 16;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-r",&min_dist,"radius","radius of influence");
    parse_args.Add("-s",&seed_dist,"radius","radius of separation for samples");
    parse_args.Add("-c",&cent_dist,"radius","radius of points for centroid check");
    parse_args.Parse();
    
    GRID<TV> grid(TV_INT()+1,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    const char* scans[]={
        "scans/M1-00819_FemaleRight.obj",
        "scans/M1-02865_FemaleRight.obj",
        "scans/M1-03013_MaleLeft.obj",
        "scans/M10M-7-2_MaleLeft.obj",
        "scans/M70-2985_FemaleLeft.obj",
        "scans/M70-4443_MaleRight.obj"
    };

    for(int s=0;s<6;s++){

    TRIANGULATED_SURFACE<T> tsa;
    tsa.Read_Obj(scans[s]);
    TRIANGULATED_SURFACE<T> &ts=*tsa.Create_Compact_Copy();
    ts.Fill_Holes(true);
    ts.mesh.Make_Orientations_Consistent();
    ts.Update_Vertex_Normals();

    KD_TREE<TV> kd_tree;
    kd_tree.Create_KD_Tree(ts.particles.X);

    ARRAY<TV> normals(ts.particles.X.m);
    
    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    ARRAY<int> sample_points;
    ARRAY<TV2> curvatures;
    ARRAY<TV> eig0,eig1;
    Compute_Samples(sample_points,ts.particles.X,kd_tree,seed_dist);
    printf("SAMPLES DONE\n");
    Compute_Normals(normals,ts.particles.X,kd_tree,min_dist);
    printf("NORMALS DONE\n");
    Orient_Normals(normals,ts.particles.X,kd_tree,min_dist);
    printf("ORIENT DONE\n");
    Compute_Curvature(curvatures,eig0,eig1,ts.particles.X,kd_tree,normals,min_dist);
    printf("CURVATURE DONE\n");

    TV2 cn=curvatures(0),cx=curvatures(0);
    for(int i=0;i<curvatures.m;i++)
    {
        cn=cn.Componentwise_Min(curvatures(i));
        cx=cx.Componentwise_Max(curvatures(i));
    }
    LOG::printf("%P  %P\n",cn,cx);
    
    Dump_Surface(ts,TV(0,.5,0));
    Flush_Frame<TV>("C");
    
    INTERPOLATED_COLOR_MAP<T> cm;
    cm.Initialize_Colors(-.1,.1,false,true,false);
    for(int i=0;i<curvatures.m;i++)
    {
        TV c=cm(curvatures(i).x);
        Add_Debug_Particle(ts.particles.X(i),c);
    }
    Flush_Frame<TV>("C");
    for(int i=0;i<curvatures.m;i++)
    {
        TV c=cm(curvatures(i).y);
        Add_Debug_Particle(ts.particles.X(i),c);
    }
    Flush_Frame<TV>("C");
    cm.Initialize_Colors(-1,1,false,true,false);
    for(int i=0;i<curvatures.m;i++)
    {
        int im=abs(curvatures(i).x)>abs(curvatures(i).y);
        TV c=cm(curvatures(i)(im)/curvatures(i)(1-im));
        Add_Debug_Particle(ts.particles.X(i),c);
    }
    Flush_Frame<TV>("C");
    }    
    return 0;
}

