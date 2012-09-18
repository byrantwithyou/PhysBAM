#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include "LEVELSET_MESH_CUTTING_3D.h"

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
    ARRAY<LEVELSET_MESH_CUTTING_3D::TET> cut_mesh;

    // ARRAY<VECTOR<int,4> > m;
    // m.Append(VECTOR<int,4>(0,1,2,3));
    // ARRAY<T> p0,p1;
    // p0.Append(0.5);
    // p0.Append(0.2);
    // p0.Append(-0.1);
    // p0.Append(-0.6);
    // p1.Append(0.6);
    // p1.Append(-0.6);
    // p1.Append(-0.4);
    // p1.Append(0.5);

    // LEVELSET_MESH_CUTTING_3D::Subdivide(m,p0,p1,cut_mesh);

    // TETRAHEDRON_MESH tm;
    // for(int i=0;i<cut_mesh.m;i++){
    //     tm.elements.Append(cut_mesh(i).indices);
    //     LOG::cout<<cut_mesh(i).parent<<"  "<<cut_mesh(i).indices<<"  "<<cut_mesh(i).weights<<std::endl;}

    // tm.elements=m;
    // tm.Set_Number_Nodes(tm.elements.Flattened().Max()+1);
    // tm.Initialize_Boundary_Mesh();

    // LOG::cout<<tm.boundary_mesh->elements<<std::endl;

    // return 0;



    GRID<TV> grid(TV_INT()+(size+1),RANGE<TV>::Unit_Box());

    TETRAHEDRALIZED_VOLUME<T>& tv=*TETRAHEDRALIZED_VOLUME<T>::Create();
    tv.Initialize_Cube_Mesh_And_Particles(grid);

    // int k=0;
    // for(int i=0;i<tv.mesh.elements.m;i++)
    //     if(tv.mesh.elements(i).Contains(13))
    //         tv.mesh.elements(k++)=tv.mesh.elements(i);
    // tv.mesh.elements.Resize(k);

//    LOG::cout<<tv.mesh.elements<<std::endl;
//    tv.mesh.elements.Remove_All();
//    tv.mesh.elements.Append(VECTOR<int,4>(30,31,28,19));


    tv.mesh.Initialize_Boundary_Mesh();
    tv.mesh.Initialize_Boundary_Nodes();

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

//    LOG::cout<<tv.mesh.elements<<std::endl;

    LEVELSET_MESH_CUTTING_3D::Subdivide(tv.mesh.elements,phi0,phi1,cut_mesh);

    // for(int i=0;i<cut_mesh.m;i++){
    //     LOG::cout<<cut_mesh(i).parent<<"  "<<cut_mesh(i).indices<<"  "<<cut_mesh(i).weights<<std::endl;}

    TETRAHEDRON_MESH tm;
    for(int i=0;i<cut_mesh.m;i++)
        tm.elements.Append(cut_mesh(i).indices);

    tm.Set_Number_Nodes(tm.elements.Flattened().Max()+1);
    tm.Initialize_Boundary_Mesh();
    LOG::cout<<"mesh sizes "<<tm.boundary_mesh->elements.m<<"  "<<tv.mesh.boundary_mesh->elements.m<<std::endl;
    // LOG::cout<<tm.elements<<std::endl;
    LOG::cout<<tm.boundary_mesh->elements<<std::endl;
    LOG::cout<<tv.mesh.boundary_mesh->elements<<std::endl;

    for(int i=0;i<cut_mesh.m;i++){
        VECTOR<TV,4> pts;
        for(int j=0;j<4;j++)
            for(int k=0;k<4;k++)
                pts(j)+=tv.particles.X(cut_mesh(i).parent(k))*cut_mesh(i).weights(j)(k);
        if(TETRAHEDRON<T>::Signed_Size(pts)<-1e-15)
            LOG::cout<<"NEGATIVE VOLUME  "<<TETRAHEDRON<T>::Signed_Size(pts)<<"   "<<cut_mesh(i).parent<<"   "<<cut_mesh(i).indices<<std::endl;}

    for(int i=0;i<cut_mesh.m;i++){
        VECTOR<T,4> a(phi0.Subset(cut_mesh(i).indices)),b(phi1.Subset(cut_mesh(i).indices));
        if(a.Min()<0 && a.Max()>0 && b.Min()<0 && b.Max()>0)
            LOG::cout<<"CROSSINGS "<<a<<"  "<<b<<std::endl;}

    delete &tv;
    return 0;
}
