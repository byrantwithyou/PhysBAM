//
//  mesh_cutting.h
//  
//
//  Created by Yuting Wang on 5/23/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#ifndef _mesh_cutting_h
#define _mesh_cutting_h

#include <set>
#include <cassert>
#include <string>

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "GEOMETRY.h"

template<class T> class DEFORMABLE_OBJECT_3D;
template<class T> class BACKWARD_EULER_TIME_STEPPING_3D;
template<class T> class FIXED_COROTATED_ELASTICITY_3D;
template<class T> class HYPERELASTICITY_CONSTITUTIVE_MODEL_3D;
template<class T> class FEM_HYPERELASTICITY_3D;

namespace PhysBAM{
template<class T>
class MESH_CUTTING
{ 
public:
    enum WORKAROUND1 { NumNodesPerTriangle = 3 };
    enum WORKAROUND2 { NumNodesEdgesTriangle = 3 };
    enum WORKAROUND3 { NumFacesPerTet = 4 };
    enum WORKAROUND4 { NumNodesPerTet = 4 };
    enum WORKAROUND5 { NumEdgesPerTet = 6 };
    enum WORKAROUND6 { NumMaterialsPerFace = 6 };
    enum WORKAROUND8 { NumSplitTurnOns = 36 };
    enum WORKAROUND9 { MaxPieces = 24 };
    enum DRAWING_MODE{ normal,picking };
    typedef std::set<int> set_int;
    typedef VECTOR<int, NumNodesPerTriangle> TV_INT;
    typedef VECTOR<T, 3> TV;
    typedef VECTOR<T, NumNodesPerTet> CENTER;
    typedef VECTOR<int, NumNodesPerTet> TET;
    typedef VECTOR<int,2> EDGE;
    typedef VECTOR<T,2> T2;
    typedef VECTOR<T,2> T3;
    typedef VECTOR<T,4> T4;
    typedef VECTOR<int,2> I2;
    typedef VECTOR<int,3> I3;
    typedef VECTOR<int,4> I4;
    typedef VECTOR<int,5> I5;
    typedef VECTOR<int,7> I7;
    typedef HASHTABLE<I4,T4> H;
    
    struct SMOOTH_NORMAL{
        TV normal;
        int n;
        SMOOTH_NORMAL() {n = 0;}
        SMOOTH_NORMAL(TV normal_input) {normal = normal_input; n = 1;}
        void add(const TV& nor) {normal+=nor; ++n;}
        TV value()
        {
            if(n == 0) {
                std::cout << "SMOOTH NORMAL not initialize!\n";
                exit(1);
            }
            else
                return normal/n;
        }
    };
    
    struct TET_CUTTING {
        ARRAY<int> material_ids;
        VECTOR<bool, MaxPieces> has_material;
        VECTOR<bool, NumSplitTurnOns> turned_on;//turned on or not
        VECTOR<bool, 24> no_merge;//merge or not
        
        struct POTENTIAL_CENTER{
            CENTER sum;
            bool set;
            int n;
            
            POTENTIAL_CENTER():set(0), n(0){}
            
            void Add(const CENTER& c){
                if(!set){
                    if (!n)
                        sum=c;
                    else
                        sum+=c;
                    n++;
                }
            }
            
            CENTER Value() const
            {
                if(!n||!set) return sum;
                return sum/n;
            }
        };

        POTENTIAL_CENTER tet_center;
        VECTOR<POTENTIAL_CENTER, NumFacesPerTet> face_centers;
        VECTOR<POTENTIAL_CENTER, NumEdgesPerTet> edge_centers;
        
        TET_CUTTING();
        TET_CUTTING Generate_Sub_Tet(const ARRAY<int>& material_ids_input);
        ARRAY<int> Find_CC(const int& material_id, VECTOR<bool, MaxPieces>& picked);
    };
    
    struct PARENT{
        int id;
        CENTER weight;
        
        PARENT(){}
        PARENT(const int id_input, const CENTER& weight_input):id(id_input), weight(weight_input){}
    };
    
    struct CUT_FACE{
        TV_INT face;
        int sim_tet_id;
        int new_sim_tet_id;
        
        CUT_FACE(){}
        CUT_FACE(const TV_INT& face_input, const int sim_tet_id_input, const int new_sim_tet_id_input):face(face_input), sim_tet_id(sim_tet_id_input), new_sim_tet_id(new_sim_tet_id_input){}
    };
    
    typedef VECTOR<CUT_FACE, NumMaterialsPerFace> TV_TV_INT;
public:
    //meshes
    TETRAHEDRALIZED_VOLUME<T>* sim_volume;
    
    TETRAHEDRALIZED_VOLUME<T>* volume;//cutting volume
    ARRAY<PARENT> weights_in_sim;//decides the particles in cutting volume
    ARRAY<int> ctet2stet;//cutting tet to sim tet mapping
    ARRAY<TET_CUTTING> tet_cuttings;

    //levelset as cutter
    ARRAY<ARRAY<T> > levelsets;
    
    //for drawing
    bool interactive;
    ARRAY<bool> is_blue;
    ARRAY<ARRAY<int> > tet_cc, node_cc;
    ARRAY<int> element2cc;
    ARRAY<int> node2cc;
    ARRAY<int> dragging_elements;
    ARRAY<TV> dragging_weights;
    ARRAY<TV> dragging_targets;
    ARRAY<VECTOR<int,2> > faces_to_draw;
    ARRAY<VECTOR<int,2> > interfaces_to_draw;
    HASHTABLE<VECTOR<int, 3> > cuttingFaces;
    
    //for simulation
    //elasticity
    ARRAY<TV> cutting_particle_material_space;
    DEFORMABLE_OBJECT_3D<T>* deformable_object;
    GEOMETRY::TETRAHEDRON_MESH* our_mesh;
    BACKWARD_EULER_TIME_STEPPING_3D<T>* be;
    HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>* le; 
    FEM_HYPERELASTICITY_3D<T>* fem;
    HASHTABLE<int> diri_nodes;
    ALGEBRA::VECTOR<T>* nodal_volumes;
    ALGEBRA::VECTOR<int>* my_constrained;
    ALGEBRA::VECTOR<T>* my_constrained_locations;
    T timestep; int ratio;
    ARRAY<ALGEBRA::MATRIX_3X3<T> > undeformed_config_copy;
    void Initialize_Elasticity();
    void Reinitialize_Elasticity();
    
    MESH_CUTTING();
    MESH_CUTTING(TETRAHEDRALIZED_VOLUME<T>* volume_input, T timestep_input, int ratio_input, bool interactive_input);
    ~MESH_CUTTING();
    
    void Initialize_Cutting_Volume();
    CENTER Weight_In_Sim_Tet(const CENTER& weight_in_element, const TET& element, const int sim_tet_id);
    CENTER Weight_In_Sim_Tet(const CENTER& weight_in_element, const TET& element, const int sim_tet_id, const ARRAY<TET>& original_mesh);
    void Refine_Cutting_Volume();
    void Update_Cutting_Particles();
    
    void Split_By_Levelset(const int& tet_id, ARRAY<int>& lsi);
    void Cut_By_Levelset(ARRAY<int>& lsi);
    void Become_Sphere();
    void Keep_CC(int cc_id);
    
    void Cut(TRIANGULATED_SURFACE<T>& cutting_surface, bool refine = true);
    int Compute_Intersection(const T& x, const T& y);
    void Update_For_Draw();
    void Draw_For_Picking();
    void Translate_CC(int picked_cc_id, TV translation);
    void Draw(bool drawing_cutting_mesh, const DRAWING_MODE& mode = normal) const;
    void Write_To_File(const std::string& writing_directory, int frame) const;
    void Write_Boundary_Mesh_To_File(const std::string& writing_directory, int frame) const;
    void Refine_And_Save_To(TETRAHEDRALIZED_VOLUME<T>*& refined_volume);
private:
    ARRAY<VECTOR<int, NumNodesPerTet> > original_elements,original_sim_elements;
    HASHTABLE<VECTOR<int,3> > cutFaces;
    ARRAY<int> original_ctet2stet;
    ARRAY<int> sim_tet_from;
            
    void Subdivide_Cutting_Mesh_Into_Eyeball();
    void Set_Dirichlet_Nodes_For_Eyeball();
    T levelset_eyeball(const TV& p);
            
    void Partial_Refine();
    void Fix_Orientation();
    void Split(const int& tet_id, HASHTABLE<int,H>& tri2inter, ARRAY<int>& sim_node_from, ARRAY<bool>& sim_tet_split, ARRAY<int>& material_nodes_from);
    void Union(const TV_TV_INT& cface,  TV_TV_INT& cface_found, UNION_FIND<int>& node_classes, UNION_FIND<int>& sim_node_classes, const int shift, const bool inverse);
    ARRAY<int> Find_Tet_CC(int m, ARRAY<bool>& picked, ARRAY<ARRAY<int> >& itot, ARRAY<ARRAY<int> >& node_cc);
    ARRAY<int> Find_Cutting_Tet_CC(int m, ARRAY<bool>& picked, ARRAY<ARRAY<int> >& itot, ARRAY<ARRAY<int> >& node_cc);
    TV weight2vec(int tet_id, CENTER c) const;
    TV weight2vec_sim(int tet_id, CENTER c) const;
    void Faces_For_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& color_vec, ARRAY<TV>& normal_vec, ARRAY<TV>& refine_edges, const DRAWING_MODE& mode) const;
    void Interfaces_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& refine_edges) const;
    void Interfaces_For_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& refine_edges, ARRAY<TV>& color_vec, const DRAWING_MODE& mode) const;
    void Colors_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& color_vec, const DRAWING_MODE& mode) const;
    ARRAY<TV> Faces_For_OpenGL(const VECTOR<int,2>& f) const;
    ARRAY<TV> Interfaces_For_OpenGL(const VECTOR<int,2>& f) const;
    inline TV& Position(const int element_id, const int node_id) const;
    inline int Index(const int element_id, const int node_id) const;
    int Sorted_Id(const I3& sorted_tri, const I3& tri, int material_id);

};

template<class TV>
void Fit_In_Box(ARRAY_VIEW<TV>& X,const RANGE<TV>& range)
{
    RANGE<TV> bb=RANGE<TV>::Bounding_Box(X);
    typename TV::SCALAR scale=(range.Edge_Lengths()/bb.Edge_Lengths()).Min();
    TV shift=range.Center()-bb.Center()*scale;
    X=X*scale+shift;
}

template<typename T>
VECTOR<T,3> weight_in_tet(const VECTOR<T,3>& p, const VECTOR<T,3>& v1, const VECTOR<T,3>& v2, const VECTOR<T,3>& v3, const VECTOR<T,3>& v4)
{
    MATRIX<T,3> ls(v1-v4, v2-v4, v3-v4);
    if (ls.Determinant() != 0) {
        return ls.Inverse() * (p-v4);
    }
    else {
        std::cout << "bad elemeent from weight_in_tet" << std::endl;
        exit(1);
    }
}

template<typename T>
void weight_in_triangle(const T& x, const T& y, const VECTOR<T,3>& v1, const VECTOR<T,3>& v2, const VECTOR<T,3>& v3, VECTOR<T,2>& w)
{
    VECTOR<T,3> v13 = v1 - v3;
    VECTOR<T,3> v23 = v2 - v3;
    MATRIX<T,2> ls(VECTOR<T,2>(v13[0], v13[1]), VECTOR<T,2>(v23[0], v23[1]));
    if (ls.Determinant() != 0) {
        w = ls.Inverse() * VECTOR<T,2>(x-v3[0], y-v3[1]);
    }
    else {
        w = VECTOR<T,2>(-1,-1);
    }
}

//ray shot from (x,y) on the screen intersects triangle (v1, v2, v3)
template<typename T>
bool intersects(const T& x, const T& y, const VECTOR<T,3>& v1, const VECTOR<T,3>& v2, const VECTOR<T,3>& v3, T& z)
{
    VECTOR<T,2> w;
    weight_in_triangle(x, y, v1, v2, v3, w);
    bool b = (w(0)>=0) && (w(1)>=0) && ((w(0)+w(1))<=1);
    if(b) z = v1[2]*w(0) + v2[2]*w(1) + v3[2]*(1-w(0)-w(1));
    return b;
}
}

#endif
