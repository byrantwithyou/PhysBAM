//#####################################################################
// Copyright 2003,Neil Molino,Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRINT
//#####################################################################
//
//#####################################################################
#ifndef __PRINT__
#define __PRINT__

#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

template <class T>
class PRINT
{
public:

static void Print_Tetrahedron(const TETRAHEDRON_MESH& mesh,int t)
{std::cout << "tetrahedron " << t << "=(" << mesh.tetrahedrons(1,t) << "," << mesh.tetrahedrons(2,t) << "," << mesh.tetrahedrons(3,t) << "," << mesh.tetrahedrons(4,t) << ")" << std::endl;
}

static void Print_Tetrahedron_Positions(const TETRAHEDRALIZED_VOLUME<T>& tet_vol,int t)
{
    std::cout << "xi=" << tet_vol.particles.X(tet_vol.tetrahedron_mesh.tetrahedrons(1,t)) << "  ";
    std::cout << "xj=" << tet_vol.particles.X(tet_vol.tetrahedron_mesh.tetrahedrons(2,t)) << "  ";
    std::cout << "xk=" << tet_vol.particles.X(tet_vol.tetrahedron_mesh.tetrahedrons(3,t)) << "  ";
    std::cout << "xl=" << tet_vol.particles.X(tet_vol.tetrahedron_mesh.tetrahedrons(4,t)) << std::endl;;
}

static void Print_Tetrahedron_Velocities(const TETRAHEDRALIZED_VOLUME<T>& tet_vol,int t)
{
    std::cout << "vi=" << tet_vol.particles.V(tet_vol.tetrahedron_mesh.tetrahedrons(1,t)) << "  ";
    std::cout << "vj=" << tet_vol.particles.V(tet_vol.tetrahedron_mesh.tetrahedrons(2,t)) << "  ";
    std::cout << "vk=" << tet_vol.particles.V(tet_vol.tetrahedron_mesh.tetrahedrons(3,t)) << "  ";
    std::cout << "vl=" << tet_vol.particles.V(tet_vol.tetrahedron_mesh.tetrahedrons(4,t)) << std::endl;;
}

static void Print_Triangle(const TRIANGLE_MESH & mesh,int t)
{std::cout << "triangle " << t << "=(" << mesh.triangles(1,t) << "," << mesh.triangles(2,t) << "," << mesh.triangles(3,t) <<  ")" << std::endl;
}

static void Print_Parent_Particles(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv, int embedded_particle)
{std::cout << "parent_particles of embedded_particle " << embedded_particle << "=(" << etv.parent_particles(1,embedded_particle) << "," << etv.parent_particles(2,embedded_particle) << ")" << std::endl;}

static void Print_Embedded_Triangles_In_Tetrahedron(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv, int tetrahedron)
{
    int tri1,tri2,tri3,tri4;
    etv.Embedded_Triangles_In_Tetrahedron(tetrahedron,tri1,tri2,tri3,tri4);
    std::cout << "triangles in tetrahedron " << tetrahedron << " are {" << tri1 << "," << tri2 << "," << tri3 << "," << tri4 << "}" << std::endl;
    if(tri1) {std::cout << "   "; Print_Triangle(etv.embedded_surface.triangle_mesh,tri1);}    
    if(tri2) {std::cout << "   "; Print_Triangle(etv.embedded_surface.triangle_mesh,tri2);}    
    if(tri3) {std::cout << "   "; Print_Triangle(etv.embedded_surface.triangle_mesh,tri3);}   
    if(tri4) {std::cout << "   "; Print_Triangle(etv.embedded_surface.triangle_mesh,tri4);}    
}

static void Print_Node_In_Tetrahedron_Is_Material(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,int t)
{std::cout << "node in tetrahedron is material" << t << "=(" << etv.node_in_tetrahedron_is_material(1,t) << "," << etv.node_in_tetrahedron_is_material(2,t) << 
                                     "," << etv.node_in_tetrahedron_is_material(3,t) << "," << etv.node_in_tetrahedron_is_material(4,t) << ")" << std::endl;
}

static void Print_Tetrahedron_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,int tetrahedron)
{
    int i,j,k,l;tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(tetrahedron,i,j,k,l);
    std::cout << "tetrahedron " << tetrahedron << " has signed volume " << TETRAHEDRON<T>::Signed_Volume(tetrahedralized_volumeparticles.X(i),tetrahedralized_volumeparticles.X(j),
                                                                                                         tetrahedralized_volumeparticles.X(k),tetrahedralized_volumeparticles.X(l)) << std::endl;
}

static void Print_Tetrahedrons_Corresponding_To_Tetrahehdron_In_Reference(int tetrahedron_in_reference_configuration)
{
    std::cout << "following tetrahedra correspond to " << tetrahedron_in_reference_configuration << " in the reference configuration" << std::endl;
    for(int t=0;t<mesh.tetrahedrons.m;t++){
        if(Corresponding_Tetrahedron_In_Reference_ETV(t) == tetrahedron_in_reference_configuration){
            Print_Tetrahedron(mesh,t);
            Print_Tetrahedron_Volume(t);
        }
    }
    std::cout << std::endl;
}

static void Print_Marked(const int center,ARRAY<int> &marked)
{
    std::cout << "center=" << center << std::endl;
    for(int p=0;p<(*reference_mesh.neighbor_nodes)(center).m;p++){
        int node=(*reference_mesh.neighbor_nodes)(center)(p);
        std::cout << "(node,marked(node))=(" << node << "," << marked(node) << ")" << std::endl;}
}

static bool Verify_Embedded_Particle_On_Segment(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    bool segment_mesh_defined=etv.tetrahedralized_volume.tetrahedron_mesh.segment_mesh!=0;
    if(!segment_mesh_defined)etv.tetrahedralized_volume.tetrahedron_mesh.Initialize_Segment_Mesh();
    for(int s=0;s<etv.tetrahedralized_volume.tetrahedron_mesh.segment_mesh->segments.m;s++){
        int a,b;etv.tetrahedralized_volume.tetrahedron_mesh.segment_mesh->segments.Get(s,a,b);
        int emb_node=etv.Embedded_Particle_On_Segment(a,b);
        if(emb_node){
            if(!etv.Are_Parents(a,b,emb_node)){std::cout << "problem in Emb Partilce On segment\n" << std::endl;assert(false);exit(1);}
        }
    }
    if(!segment_mesh_defined){delete etv.tetrahedralized_volume.tetrahedron_mesh.segment_mesh;etv.tetrahedralized_volume.tetrahedron_mesh.segment_mesh=0;}
    return true;
}

static bool Verify_Embedded_Triangle_Is_In_Tetrahedron(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,int emb_tri,int tetrahedron)
{
    int a,b,c;etv.embedded_surface.triangle_mesh.triangles.Get(emb_tri,a,b,c);
    int ppa1,ppa2,ppb1,ppb2,ppc1,ppc2;etv.parent_particles.Get(a,ppa1,ppa2);etv.parent_particles.Get(b,ppb1,ppb2);etv.parent_particles.Get(c,ppc1,ppc2);
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppa1,tetrahedron))return false;
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppa2,tetrahedron))return false;
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppb1,tetrahedron))return false;
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppb2,tetrahedron))return false;
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppc1,tetrahedron))return false;
    if(!etv.tetrahedralized_volume.tetrahedron_mesh.Node_In_Tetrahedron(ppc2,tetrahedron))return false;
    return true;
}

static bool Verify_Tetrahedron_Containing_Triangle(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    for(int t=0;t<etv.embedded_surface.triangle_mesh.triangles.m;t++){ 
        int tetrahedron=etv.Tetrahedron_Containing_Triangle(t);
        if(!Verify_Embedded_Triangle_Is_In_Tetrahedron(etv,t,tetrahedron)){std::cout << "problem with tet containing tri" << std::endl;assert(false);exit(1);}
    }
    return true;
}

static bool Verify_Embedded_Triangles_In_Tetrahedron(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    if(etv.embedded_sub_elements_in_parent_element_index->m != etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){assert(false);return false;}
    int total_number_of_emb_tri=0;
    for(int tet=0;tet<etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;tet++){
        int emb_tri1,emb_tri2,emb_tri3,emb_tri4;total_number_of_emb_tri+=etv.Embedded_Triangles_In_Tetrahedron(tet,emb_tri1,emb_tri2,emb_tri3,emb_tri4);
        if(emb_tri1 && !Verify_Embedded_Triangle_Is_In_Tetrahedron(etv,emb_tri1,tet)){std::cout << "problem with emb tri in tet" << std::endl;assert(false);exit(1);}
        if(emb_tri2 && !Verify_Embedded_Triangle_Is_In_Tetrahedron(etv,emb_tri2,tet)){std::cout << "problem with emb tri in tet" << std::endl;assert(false);exit(1);}
        if(emb_tri3 && !Verify_Embedded_Triangle_Is_In_Tetrahedron(etv,emb_tri3,tet)){std::cout << "problem with emb tri in tet" << std::endl;assert(false);exit(1);}
        if(emb_tri4 && !Verify_Embedded_Triangle_Is_In_Tetrahedron(etv,emb_tri4,tet)){std::cout << "problem with emb tri in tet" << std::endl;assert(false);exit(1);}
    }
    if(total_number_of_emb_tri != etv.embedded_surface.triangle_mesh.triangles.m){assert(false);return false;}
    return true;
}

static bool Verify_Children_Of_This_Node(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv, int tet_node)
{
    for(int c=0;c<etv.Number_Of_Children(tet_node);c++){
        int child=etv.Child(tet_node,c);
        if(!etv.Is_Parent(tet_node,child)){assert(false);return false;}
    }
    return true;
}

static bool Verify_Child_Structure(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    if(!etv.embedded_children_index)return true;
    if(etv.embedded_children_index->m != etv.tetrahedralized_volume.tetrahedron_mesh.number_nodes){assert(false);return false;}
    if(etv.tetrahedralized_volume.particles.Size() != etv.tetrahedralized_volume.tetrahedron_mesh.number_nodes){assert(false);return false;}
    int total_number_of_children=0;
    for(int tet_node=0;tet_node<etv.tetrahedralized_volume.particles.Size();tet_node++){
        if(!Verify_Children_Of_This_Node(etv,tet_node)){assert(false);return false;}
        total_number_of_children+=etv.Number_Of_Children(tet_node);
    }
    if(total_number_of_children != 2*etv.embedded_particles.Size()){assert(false);return false;}
    return true;
}

static bool Verify_Array_Sizes_In_ETV(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    if(etv.tetrahedralized_volume.tetrahedron_mesh.number_nodes != etv.tetrahedralized_volume.particles.Size()){assert(false);return false;}
    if(etv.embedded_particles.Size() != etv.embedded_triangle_mesh.number_nodes){assert(false);return false;}
    if(etv.parent_particles.m != etv.embedded_particles.Size()){assert(false);return false;}
    if(etv.interpolation_fraction.m != etv.embedded_particles.Size()){assert(false);return false;}
    if(etv.node_in_tetrahedron_is_material.m != etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){assert(false);return false;}
    if(etv.embedded_sub_elements_in_parent_element_index->m != etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){assert(false);return false;}
    if(etv.embedded_children_index && etv.embedded_children_index->m != etv.tetrahedralized_volume.particles.Size()){assert(false);return false;}
    if(etv.embedded_sub_elements_in_parent_element->m != etv.number_of_embedded_sub_elements_in_parent_element->m){assert(false);return false;}
    return true;
}

static bool Emb_Node_In_Triangle(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T> &etv,const int emb_node,const int emb_tri)
{
    assert(emb_node != 0);
    return etv.embedded_surface.triangle_mesh.Node_In_Triangle(emb_node,emb_tri);
}

static bool Emb_Node_In_Quad(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T> &etv,const int emb_node,const int emb_tri1,const int emb_tri2)
{
    assert(emb_node != 0);
    return etv.embedded_surface.triangle_mesh.Node_In_Triangle(emb_node,emb_tri1) || etv.embedded_surface.triangle_mesh.Node_In_Triangle(emb_node,emb_tri2);
}

static bool Verify_ETV(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv)
{
    if(!Verify_Embedded_Particle_On_Segment(etv)){assert(false);exit(1);return false;}
    if(!Verify_Tetrahedron_Containing_Triangle(etv)){assert(false);exit(1);return false;}
    if(!Verify_Embedded_Triangles_In_Tetrahedron(etv)){assert(false);exit(1);return false;}
    if(!Verify_Child_Structure(etv)){assert(false);exit(1);return false;}
    if(!Verify_Array_Sizes_In_ETV(etv)){assert(false);exit(1);return false;}
    return true;
}

static bool Verify_FTV(FRACTURE_TETRAHEDRALIZED_VOLUME<T>& ftv)
{
    if(!Verify_ETV(ftv.embedded_tetrahedralized_volume)){assert(false);return false;}
    if(!ftv.fracture_bias_stress_scaling.m == ftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){assert(false);return false;}
    return true;
}


static bool Verify_Four_Emb_Triangle_Case(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>&etv,int i,int j,int k,int l,int isolated_emb_tri1,int isolated_emb_tri2,int emb_tri3,int emb_tri4)
{
    assert(isolated_emb_tri1&&isolated_emb_tri2&&emb_tri3&&emb_tri4);
    int ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l),
        jk=etv.Embedded_Particle_On_Segment(j,k),jl=etv.Embedded_Particle_On_Segment(j,l),
        ij=etv.Embedded_Particle_On_Segment(i,j),kl=etv.Embedded_Particle_On_Segment(k,l);
    if(!(ij&&ik&&il&&jk&&jl&&kl)){assert(false);exit(1);return false;}
    if(!Emb_Node_In_Quad(etv,ik,emb_tri3,emb_tri4) || !Emb_Node_In_Quad(etv,il,emb_tri3,emb_tri4) || 
       !Emb_Node_In_Quad(etv,jk,emb_tri3,emb_tri4) || !Emb_Node_In_Quad(etv,jl,emb_tri3,emb_tri4)){assert(false);exit(1);return false;}
    if(Emb_Node_In_Triangle(etv,ij,isolated_emb_tri1)){
        assert(Emb_Node_In_Triangle(etv,kl,isolated_emb_tri2));
        if((i != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri1)) && (j != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri1))){assert(false);exit(1);return false;} 
        if((k != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri2)) && (l != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri2))){assert(false);exit(1);return false;} 
    }
    else{
        assert(Emb_Node_In_Triangle(etv,ij,isolated_emb_tri2));
        assert(Emb_Node_In_Triangle(etv,kl,isolated_emb_tri1));
        if((i != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri2)) && (j != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri2))){assert(false);exit(1);return false;} 
        if((k != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri1)) && (l != etv.Node_Separated_By_Embedded_Triangle(isolated_emb_tri1))){assert(false);exit(1);return false;} 
    }
    return true;
}

static bool Verify_Orientation_Index_Tet(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,int tetrahedron,int orientation_index)
{
    if(!(1<=orientation_index || orientation_index<=24)){assert(false);exit(1);return false;}
    int oi,oj,ok,ol,i,j,k,l;
    etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(tetrahedron,oi,oj,ok,ol);
    permute_four_ints(oi,oj,ok,ol,i,j,k,l,orientation_index);
    int emb_tri1,emb_tri2,emb_tri3,emb_tri4;
    etv.Embedded_Triangles_In_Tetrahedron(tetrahedron,emb_tri1,emb_tri2,emb_tri3,emb_tri4);
    if(!emb_tri1){
        return true;
    }
    else if(emb_tri1 && !emb_tri2){
        assert(!emb_tri3 && !emb_tri4);
        assert(!emb_tri2 && !emb_tri3 && !emb_tri4);
        int ij=etv.Embedded_Particle_On_Segment(i,j),ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l);
        if(!(ij&&ik&&il)){assert(false);exit(1);return false;}
        if(!Emb_Node_In_Triangle(etv,ij,emb_tri1) || !Emb_Node_In_Triangle(etv,ik,emb_tri1) || !Emb_Node_In_Triangle(etv,il,emb_tri1)){assert(false);exit(1);return false;}
        if(!(emb_tri1 == etv.embedded_triangle_mesh.Triangle(ij,ik,il))){assert(false);exit(1);return false;}
        return true;
    }
    else if(emb_tri1 && emb_tri2 && !emb_tri3){
        assert(!emb_tri4);
        if(etv.Cut_By_Quad(tetrahedron)){
            int ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l),
                jk=etv.Embedded_Particle_On_Segment(j,k),jl=etv.Embedded_Particle_On_Segment(j,l);
            if(!(ik&&il&&jk&&jl)){assert(false);exit(1);return false;}
            if(!Emb_Node_In_Quad(etv,ik,emb_tri1,emb_tri2) || !Emb_Node_In_Quad(etv,il,emb_tri1,emb_tri2) || 
               !Emb_Node_In_Quad(etv,jk,emb_tri1,emb_tri2) || !Emb_Node_In_Quad(etv,jl,emb_tri1,emb_tri2)){assert(false);exit(1);return false;}
            return true;
        }
        else {
            int ij=etv.Embedded_Particle_On_Segment(i,j),ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l);
            int jk=etv.Embedded_Particle_On_Segment(j,k),jl=etv.Embedded_Particle_On_Segment(j,l);
            if(!(ij&&ik&&il)){assert(false);exit(1);return false;}
            if(!Emb_Node_In_Triangle(etv,ij,emb_tri1) || !Emb_Node_In_Triangle(etv,ik,emb_tri1) || !Emb_Node_In_Triangle(etv,il,emb_tri1)){assert(false);exit(1);return false;}
            if(!(ij&&jk&&jl)){assert(false);exit(1);return false;}
            if(!Emb_Node_In_Triangle(etv,ij,emb_tri2) || !Emb_Node_In_Triangle(etv,jk,emb_tri2) || !Emb_Node_In_Triangle(etv,jl,emb_tri2)){assert(false);exit(1);return false;}
            return true;
        }
    }
    else if(emb_tri1 && emb_tri2 && emb_tri3 && !emb_tri4){
        int isolated_node1=etv.Node_Separated_By_Embedded_Triangle(emb_tri1);
        int isolated_node2=etv.Node_Separated_By_Embedded_Triangle(emb_tri2);
        int isolated_node3=etv.Node_Separated_By_Embedded_Triangle(emb_tri3);
        if(isolated_node1 && isolated_node2 && isolated_node3){ // 3 subtets + [oct+subtet]
            int ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l),
                jk=etv.Embedded_Particle_On_Segment(j,k),jl=etv.Embedded_Particle_On_Segment(j,l),
                ij=etv.Embedded_Particle_On_Segment(i,j),kl=etv.Embedded_Particle_On_Segment(k,l);
            if(!(ij&&ik&&il&&jk&&jl&&kl)){assert(false);exit(1);return false;}
            assert(i == isolated_node1);
            assert(j == isolated_node2);
            assert(k == isolated_node3 || l == isolated_node3);
            if(!Emb_Node_In_Triangle(etv,ij,emb_tri1) || !Emb_Node_In_Triangle(etv,ik,emb_tri1) || !Emb_Node_In_Triangle(etv,il,emb_tri1)){assert(false);exit(1);return false;}
            if(!Emb_Node_In_Triangle(etv,ij,emb_tri2) || !Emb_Node_In_Triangle(etv,jk,emb_tri2) || !Emb_Node_In_Triangle(etv,jl,emb_tri2)){assert(false);exit(1);return false;}
            if(k == isolated_node3){
                if(!Emb_Node_In_Triangle(etv,ik,emb_tri3) || !Emb_Node_In_Triangle(etv,jk,emb_tri3) || !Emb_Node_In_Triangle(etv,kl,emb_tri3)){assert(false);exit(1);return false;}
            }
            else{
                assert(l == isolated_node3);
                if(!Emb_Node_In_Triangle(etv,il,emb_tri3) || !Emb_Node_In_Triangle(etv,jl,emb_tri3) || !Emb_Node_In_Triangle(etv,kl,emb_tri3)){assert(false);exit(1);return false;}
            }
        }
        else{
            assert(etv.Cut_By_Quad(tetrahedron)); // wedge + [oct/2 + subtet] + subtet
            int ik=etv.Embedded_Particle_On_Segment(i,k),il=etv.Embedded_Particle_On_Segment(i,l),
                jk=etv.Embedded_Particle_On_Segment(j,k),jl=etv.Embedded_Particle_On_Segment(j,l),
                ij=etv.Embedded_Particle_On_Segment(i,j);
            if(!(ij&&ik&&il&&jk&&jl)){assert(false);exit(1);return false;}

            assert(!isolated_node1 && !isolated_node2 && isolated_node3 || isolated_node1 && !isolated_node2 && !isolated_node3);
            if(isolated_node3){
                if(!Emb_Node_In_Triangle(etv,ij,emb_tri3) || !Emb_Node_In_Triangle(etv,ik,emb_tri3) || !Emb_Node_In_Triangle(etv,il,emb_tri3)){assert(false);exit(1);return false;}
                if(!Emb_Node_In_Quad(etv,ik,emb_tri1,emb_tri2) || !Emb_Node_In_Quad(etv,il,emb_tri1,emb_tri2) || 
                   !Emb_Node_In_Quad(etv,jk,emb_tri1,emb_tri2) || !Emb_Node_In_Quad(etv,jl,emb_tri1,emb_tri2)){assert(false);exit(1);return false;}
            }
            else{
                assert(isolated_node1);
                if(!Emb_Node_In_Triangle(etv,ij,emb_tri1) || !Emb_Node_In_Triangle(etv,ik,emb_tri1) || !Emb_Node_In_Triangle(etv,il,emb_tri1)){assert(false);exit(1);return false;}
                if(!Emb_Node_In_Quad(etv,ik,emb_tri2,emb_tri3) || !Emb_Node_In_Quad(etv,il,emb_tri2,emb_tri3) || 
                   !Emb_Node_In_Quad(etv,jk,emb_tri2,emb_tri3) || !Emb_Node_In_Quad(etv,jl,emb_tri2,emb_tri3)){assert(false);exit(1);return false;}
            }
        }
    }
    else if(emb_tri1 && emb_tri2 && emb_tri3 && emb_tri4){
        assert(etv.Cut_By_Quad(tetrahedron));
        int isolated_node1=etv.Node_Separated_By_Embedded_Triangle(emb_tri1);
        int isolated_node2=etv.Node_Separated_By_Embedded_Triangle(emb_tri2);
        int isolated_node3=etv.Node_Separated_By_Embedded_Triangle(emb_tri3);
        int isolated_node4=etv.Node_Separated_By_Embedded_Triangle(emb_tri4);

        // isolated_node triangles are first two in call
        if(isolated_node1 && isolated_node2 && !isolated_node3 && !isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri1,emb_tri2,emb_tri3,emb_tri4)){assert(false);exit(1);return false;}
        }
        else if (isolated_node1 && !isolated_node2 && isolated_node3 && !isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri1,emb_tri3,emb_tri2,emb_tri4)){assert(false);exit(1);return false;}
        }
        else if (isolated_node1 && !isolated_node2 && !isolated_node3 && isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri1,emb_tri4,emb_tri2,emb_tri3)){assert(false);exit(1);return false;}
        }
        else if (!isolated_node1 && isolated_node2 && isolated_node3 && !isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri2,emb_tri3,emb_tri1,emb_tri4)){assert(false);exit(1);return false;}
        }
        else if (!isolated_node1 && isolated_node2 && !isolated_node3 && isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri2,emb_tri4,emb_tri1,emb_tri3)){assert(false);exit(1);return false;}
        }
        else if (!isolated_node1 && !isolated_node2 && isolated_node3 && isolated_node4){
            if(!Verify_Four_Emb_Triangle_Case(etv,i,j,k,l,emb_tri3,emb_tri4,emb_tri1,emb_tri2)){assert(false);exit(1);return false;}
        }
        return true;
    }
    return true;
}

static bool Verify_Orientation_Index(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,ARRAY<int>& orientation_index)
{
    if(orientation_index.m != etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){assert(false);return false;}
    for(int t=0;t<orientation_index.m;t++){
        if(!Verify_Orientation_Index_Tet(etv,t,orientation_index(t))){assert(false);return false;}
    }
    return true;
}


static bool Verify_Node_In_Tetrahedron_Is_Material(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,
                                                   EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& retv,
                                                   VIRTUAL_NODE_ALGORITHM<T>& vna)
{
    ARRAY<bool>to_check(4,1,retv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    for(int t=0;t<etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++){
        int ref_t=vna.corresponding_tetrahedron_in_reference(t);
        for(int a=0;a<4;a++){
            int node=etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons(a,t);
            if(etv.Node_In_Tetrahedron_Is_Material(node,t)){
                if(to_check(a,ref_t)){std::cout << "problem here" << "a=" << a << "  t=" << t << std::endl;assert(false);exit(1);return false;}
                else to_check(a,ref_t)=true;
            }
        }
    }
    for(int t=0;t<to_check.m;t++) for(int a=0;a<4;a++) if(!to_check(a,t)){std::cout << "problem here" << "a=" << a << "  t=" << t << std::endl;assert(false);exit(1);return false;}



    return true;
}









};
}
#endif
