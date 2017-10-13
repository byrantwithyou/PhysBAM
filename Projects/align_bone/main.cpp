//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
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
    T min_dist = 8;
    T seed_dist = 1;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-r",&min_dist,"radius","radius of influence");
    parse_args.Add("-s",&seed_dist,"radius","radius of separation for samples");
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

    for(int s=0;s<sizeof(scans)/sizeof(*scans);s++)
    {
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(123);
        TRIANGULATED_SURFACE<T> tsa;
        tsa.Read_Obj(scans[s]);
        TRIANGULATED_SURFACE<T> &ts=*tsa.Create_Compact_Copy();
        ts.mesh.Make_Orientations_Consistent();
        ts.Update_Vertex_Normals();
    
        KD_TREE<TV> kd_tree;
        kd_tree.Create_KD_Tree(ts.particles.X);

        ARRAY<TV> normals(ts.particles.X.m);
        ARRAY<TV2> curvatures(ts.particles.X.m);
        ARRAY<TV> eig0(ts.particles.X.m);
        ARRAY<TV> eig1(ts.particles.X.m);
    
        ARRAY<int> points_found;
        ARRAY<T> distance_squared_of_points_found;
        ARRAY<int> sample_points,todo(IDENTITY_ARRAY<>(ts.particles.X.m));
        ARRAY<bool> done(ts.particles.X.m);
        random.Random_Shuffle(todo);
        for(int i=0;i<todo.m;i++)
        {
            int p=todo(i);
            if(done(p)) continue;
            sample_points.Append(p);
            
            TV X=ts.particles.X(p);
            kd_tree.Find_Points_Within_Radius(X,sqr(seed_dist),points_found,
                distance_squared_of_points_found,ts.particles.X);
            done.Subset(points_found).Fill(true);
            Add_Debug_Particle(X,TV(0,1,0));
            Add_Debug_Object(VECTOR<TV,2>(X,X+normals(p)),TV(1,0,0));
            
            points_found.Remove_All();
            distance_squared_of_points_found.Remove_All();
        }
        printf("samples %i\n",sample_points.m);
        
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

            TV u=n.Unit_Orthogonal_Vector(),v=n.Cross(u);
            MATRIX<T,3,2> R(u,v);
            
            MATRIX<T,6> M;
            VECTOR<T,6> b;
            for(int j=0;j<points_found.m;j++)
            {
                TV Z=ts.particles.X(points_found(j))-X;
                TV2 Y=R.Transpose_Times(Z);
                VECTOR<T,6> c(sqr(Y.x),Y.x*Y.y,sqr(Y.y),Y.x,Y.y,1);
                M+=Outer_Product(c,c);
                b+=c*n.Dot(Z);
            }
            VECTOR<T,6> s=M.Inverse_Times(b);
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

        TV colors[16]={
            TV(1,0,0),
            TV(1,.5,0),
            TV(1,1,0),
            TV(0,1,0),
            TV(0,1,1),
            TV(0,0,1),
            TV(1,0,1),
            TV(.5,0,0),
            TV(.5,.5,0),
            TV(0,.5,0),
            TV(0,.5,.5),
            TV(0,0,.5),
            TV(.5,0,.5),
            TV(.5,.5,.5),
            TV(1,1,1),
            TV(.2,.2,.2)
        };

        auto compute_score =[](TV2 K)
            {
                if(K(0)<0) return (T)-1;
                if(K(1)<.02) return (T)-1;
//                if(K(1)>3*K(0)) return (T)-1;
                return K(0);

                if(K(1)>=0) return (T)-1;
                return -K.Dot(TV2(1,0));
            };

        ARRAY<T> scores;
        for(int i=0;i<ts.particles.X.m;i++)
        {
            TV2 K=curvatures(i);
            T s=compute_score(K);
            if(s>=0) scores.Append(s);
        }
        PHYSBAM_ASSERT(scores.m);
        scores.Sort();
        T thresh_score=scores(scores.m*.8);

        for(int i=0;i<ts.particles.X.m;i++)
        {
            TV2 K=curvatures(i);
            TV X=ts.particles.X(i);
            int i0=(K(0)>-.1)+(K(0)>-.05)+(K(0)>.05);
            int i1=(K(1)>.05)+(K(1)>.1);
            int j=i0+i1*4;
            continue;
            T s=compute_score(K);
            if(s>thresh_score) Add_Debug_Particle(ts.particles.X(i),TV(0,0,0));
            else if(s>0) Add_Debug_Particle(ts.particles.X(i),TV(1,1,1));
            else Add_Debug_Particle(ts.particles.X(i),colors[j]);
        }
        
        if(1)
        for(int i=0;i<ts.mesh.elements.m;i++)
            Add_Debug_Object(ts.Get_Element(i).X,TV(1,0,0));

        Flush_Frame<TV>(scans[s]);
    }
    

    

    Flush_Frame<TV>("B");
    
    return 0;
}

