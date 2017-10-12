//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef double RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+1,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    T min_dist = 1.5;

    TRIANGULATED_SURFACE<T> tsa;
    tsa.Read_Obj("scans/M10M-7-2_MaleLeft.obj");
    TRIANGULATED_SURFACE<T> &ts=*tsa.Create_Compact_Copy();
    ts.mesh.Make_Orientations_Consistent();
    ts.Update_Vertex_Normals();
    
    KD_TREE<TV> kd_tree;
    kd_tree.Create_KD_Tree(ts.particles.X);

    ARRAY<TV> normals(ts.particles.X.m);
    
    ARRAY<int> points_found;
    ARRAY<T> distance_squared_of_points_found;
    for(int i=0;i<ts.particles.X.m;i++)
    {
        TV X=ts.particles.X(i);
        kd_tree.Find_Points_Within_Radius(X,sqr(min_dist),points_found,
            distance_squared_of_points_found,ts.particles.X);

        TV A=ts.particles.X.Subset(points_found).Average();
        SYMMETRIC_MATRIX<T,TV::m> C;
        for(int j=0;j<points_found.m;j++)
            C+=Outer_Product(ts.particles.X(points_found(j))-A);
        C/=points_found.m;

        DIAGONAL_MATRIX<T,TV::m> eigenvalues;
        MATRIX<T,TV::m> eigenvectors;
        C.Solve_Eigenproblem(eigenvalues,eigenvectors);
        int k=eigenvalues.x.Arg_Min();
        TV n=eigenvectors.Column(k);
        TV N=(*ts.vertex_normals)(i);
        if(N.Dot(n)<0) n=-n;
        normals(i)=n;
        // Add_Debug_Particle(X,TV(0,1,0));
        // Debug_Particle_Set_Attribute<TV,TV>("V",n);

        points_found.Remove_All();
        distance_squared_of_points_found.Remove_All();
    }

    for(int i=0;i<ts.particles.X.m;i++)
    {
        TV X=ts.particles.X(i);
        TV N=normals(i),u=N.Unit_Orthogonal_Vector(),v=N.Cross(u);
        MATRIX<T,3,2> R(u,v);
        kd_tree.Find_Points_Within_Radius(X,sqr(min_dist),points_found,
            distance_squared_of_points_found,ts.particles.X);

        MATRIX<T,6> C;
        VECTOR<T,6> b;
        for(int j=0;j<points_found.m;j++)
        {
            TV Z=ts.particles.X(points_found(j))-X;
            TV2 Y=R.Transpose_Times(Z);
            VECTOR<T,6> c(sqr(Y.x),Y.x*Y.y,sqr(Y.y),Y.x,Y.y,1);
            C+=Outer_Product(c,c);
            b+=c*N.Dot(Z);
        }
        VECTOR<T,6> s=C.Inverse_Times(b);
        SYMMETRIC_MATRIX<T,2> K(s(0),s(1),s(2));
        DIAGONAL_MATRIX<T,2> eigenvalues;
        MATRIX<T,2> eigenvectors;
        K.Solve_Eigenproblem(eigenvalues,eigenvectors);
        int k0=eigenvalues.x.Arg_Min(),k1=1-k0;
        T c0=eigenvalues.x(k0),c1=eigenvalues.x(k1);
        TV K0=R*eigenvectors.Column(k0);
        TV K1=R*eigenvectors.Column(k1);
        // Add_Debug_Object(VECTOR<TV,2>(X,X+K0*.2),TV(0,1,0));
        // Add_Debug_Object(VECTOR<TV,2>(X,X+K1*.2),TV(0,0,1));
        Add_Debug_Particle(X,TV(c0>0,c1>0,1));
        
        // SYMMETRIC_MATRIX<T,TV::m> P=(T)1-Outer_Product(N);
        // TV A=P*ts.particles.X.Subset(points_found).Average();
        // SYMMETRIC_MATRIX<T,TV::m> C;
        // for(int j=0;j<points_found.m;j++)
        //     C+=Outer_Product(normals(points_found(j))-N);
        // C/=points_found.m;
        // SYMMETRIC_MATRIX<T,2> K=R.Transpose_Times(C*R).Symmetric_Part();
        
        // DIAGONAL_MATRIX<T,TV::m> eigenvalues;
        // MATRIX<T,TV::m> eigenvectors;
        // K.Solve_Eigenproblem(eigenvalues,eigenvectors);
        // int k0=eigenvalues.x.Arg_Min(),k1=1-k0;
        // TV K0=R*eigenvectors.Column(k0);
        // TV K1=R*eigenvectors.Column(k1);
        
        
        // MATRIX<T,2> K=C.Inverse()*D;
        // LOG::printf("%P\n",K);
        



        
#if 0
        DIAGONAL_MATRIX<T,2> eigenvalues;
        MATRIX<T,2> eigenvectors;
            
        Add_Debug_Particle(X,TV(0,1,0));
//        Debug_Particle_Set_Attribute<TV,TV>("V",n);
#endif
        points_found.Remove_All();
        distance_squared_of_points_found.Remove_All();
    }
    
    // for(int i=0;i<ts.mesh.elements.m;i++)
    //     Add_Debug_Object(ts.Get_Element(i).X,TV(1,0,0));


    

    

    Flush_Frame<TV>("A");
    Flush_Frame<TV>("B");
    
    return 0;
}

