#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<T,3> TV;

    int seed=time(0),size=4;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-s",&seed,"value","Seed");
    parse_args.Add("-n",&size,"value","Domain size");
    parse_args.Parse();

#if 0
    ARRAY<VECTOR<int,4> > m;
    m.Append(VECTOR<int,4>(0,1,2,3));
    m.Append(VECTOR<int,4>(0,1,4,2));
    ARRAY<T> p0(5);
    printf("seed %i\n", seed);
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(seed);
    random.Fill_Uniform(p0,-5,5);
    for(int i=0;i<p0.m;i++) p0(i)=::PhysBAM::rint(p0(i));
    ARRAY<VECTOR<int,4> > c=m,sp;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    LOG::cout<<m<<"  a  "<<c<<"    "<<p0<<std::endl;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,c,sp,p0,weights);
    LOG::cout<<m<<"  b  "<<c<<"    "<<sp<<"    "<<weights<<std::endl;

    HASHTABLE<VECTOR<int,3>,int> hash;
    for(int i=0;i<c.m;i++)
        for(int j=0;j<4;j++)
            hash.Get_Or_Insert(c(i).Remove_Index(j).Sorted())+=m(i)(3)*10+1;
    for(HASHTABLE<VECTOR<int,3>,int>::ITERATOR it(hash);it.Valid();it.Next())
        if(it.Data()==52)
            LOG::cout<<"SHARE "<<it.Key()<<std::endl;
    LOG::cout<<hash<<std::endl;

    return 0;

#else

    GRID<TV> grid(TV_INT()+(size+1),RANGE<TV>::Unit_Box());

    TETRAHEDRALIZED_VOLUME<T>& tv=*TETRAHEDRALIZED_VOLUME<T>::Create();
    tv.Initialize_Cube_Mesh_And_Particles(grid);
    tv.mesh.Initialize_Boundary_Nodes();

    tv.mesh.Initialize_Boundary_Mesh();
    int mm=tv.mesh.boundary_mesh->elements.m;

    ARRAY<T> phi0(tv.particles.number),phi1(tv.particles.number);
    RANDOM_NUMBERS<T> random;
    printf("seed %i\n", seed);
    random.Set_Seed(seed);
    random.Fill_Uniform(phi0,-5,5);
    random.Fill_Uniform(phi1,-5,5);
    for(int i=0;i<phi0.m;i++) phi0(i)=::PhysBAM::rint(phi0(i));
    for(int i=0;i<phi1.m;i++) phi1(i)=::PhysBAM::rint(phi1(i));
    phi0.Subset(*tv.mesh.boundary_nodes).Fill(1);
    phi1.Subset(*tv.mesh.boundary_nodes).Fill(1);

    ARRAY<VECTOR<int,4> > m=tv.mesh.elements,sp,tmp0,tmp1;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    ARRAY<bool> side;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,tv.mesh.elements,tmp0,tmp1,sp,side,phi0,weights);
    m=tmp0;
    tv.mesh.elements=tmp1;
    tv.particles.Add_Elements(weights.m);
    for(int i=0;i<weights.m;i++)
        tv.particles.X(i+phi0.m)=tv.particles.X(weights(i).x.x)*(1-weights(i).y)+tv.particles.X(weights(i).x.y)*weights(i).y;
    phi0.Resize(tv.particles.X.m);

    tv.Update_Number_Nodes();
    tv.mesh.Initialize_Boundary_Nodes();

    tv.mesh.Set_Number_Nodes(tv.mesh.elements.Flattened().Max()+1);
    tv.mesh.Initialize_Boundary_Mesh();
    LOG::cout<<"mesh sizes "<<tv.mesh.boundary_mesh->elements.m<<"  "<<mm<<std::endl;

    for(int i=0;i<tv.mesh.elements.m;i++){
        if(tv.Signed_Volume(i)<-1e-15)
            LOG::cout<<"NEGATIVE VOLUME  "<<tv.Signed_Volume(i)<<"   "<<tv.mesh.elements(i)<<std::endl;}

    for(int i=0;i<tv.mesh.elements.m;i++){
        VECTOR<T,4> a(phi0.Subset(tv.mesh.elements(i)));
        if(a.Min()<0 && a.Max()>0)
            LOG::cout<<"CROSSINGS "<<a<<std::endl;}

    delete &tv;
#endif
    return 0;
}
