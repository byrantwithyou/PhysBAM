//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Math_Tools/cyclic_shift.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

HASHTABLE<TV_INT,int> h;
ARRAY<TV_INT> pts(8+12);

int Safe_Map(const TV& pt)
{
    TV_INT ipt(rint(pt));
    PHYSBAM_ASSERT((TV(ipt)-pt).Magnitude()<1e-4);
    return h.Get(ipt);
}

int Six_Vol(TV_INT A,TV_INT B,TV_INT C)
{
    return (B.Remove_Index(2)-C.Remove_Index(2)).Cross(A.Remove_Index(2)-C.Remove_Index(2)).x*(A.z+B.z+C.z);
}

int Six_Vol(TV_INT a)
{
    return Six_Vol(pts(a.x),pts(a.y),pts(a.z));
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+2,RANGE<TV>::Unit_Box()*2,false);
    ARRAY<T,TV_INT> phi(TV_INT()+2);

    MARCHING_CUBES_CASE<3>::Case_Table();
    for(int i=0;i<8;i++)
    {
        TV_INT j(i&1,i/2&1,i/4&1);
        pts(i)=2*j;
    }
    for(int i=0;i<8;i++)
        for(int j=i+1;j<8;j++)
        {
            int k=i^j;
            if(k&(k-1)) continue;
            int e=MARCHING_CUBES_CASE<3>::edge_lookup[i][j];
            if(e<0) continue;
            e+=8;
            pts(e)=(pts(i)+pts(j))/2;
            LOG::printf("%P %P %P %P %P %P\n",i,j,e,pts(e),pts(i),pts(j));
        }
    for(int i=0;i<pts.m;i++)
        h.Set(pts(i),i);

    for(int c=0;c<256;c++)
    {
        LOG::printf("CASE: %i\n",c);
        for(int i=0;i<8;i++)
        {
            TV_INT j(i&1,i/2&1,i/4&1);
            phi(j)=(c&(1<<i))?-1:1;
        }
        TETRAHEDRALIZED_VOLUME<T> tv;
        MARCHING_CUBES<TV>::Create_Interior(tv,grid,phi,true);
        ARRAY<TV_INT> tris;
        T vol0=0,vol1=0;
        for(int i=0;i<tv.mesh.elements.m;i++)
        {
            vol0+=tv.Get_Element(i).Signed_Volume();
            VECTOR<TV,4> v(tv.particles.X.Subset(tv.mesh.elements(i)));
            VECTOR<int,4> iv;
            for(int j=0;j<4;j++)
                iv(j)=Safe_Map(v(j));
            tris.Append(TV_INT(iv(0),iv(1),iv(2)));
            tris.Append(TV_INT(iv(2),iv(1),iv(3)));
            tris.Append(TV_INT(iv(0),iv(2),iv(3)));
            tris.Append(TV_INT(iv(1),iv(0),iv(3)));
        }

        {
            int vol2=0;
            for(int i=0;i<tris.m;i++)
                vol2+=Six_Vol(tris(i));
            LOG::printf("A %P %P\n",vol0*6,vol2);
        }

        HASHTABLE<TV_INT> keep;
        for(int i=0;i<tris.m;i++)
        {
            TV_INT t=tris(i);
            while(t(0)!=t.Min()) cyclic_shift(t);
            TV_INT u=t;
            exchange(u(1),u(2));
            if(keep.Contains(u))
                keep.Delete(u);
            else
                keep.Set(t);
        }
        tris.Remove_All();
        for(auto it:keep)
        {
            int sv=Six_Vol(it);
            if(sv)
                tris.Append(it);
            vol1+=sv;
        }
        LOG::printf("%P %P\n",vol0*6,vol1);
        LOG::printf("%P\n",tris);
    }

    return 0;
}

