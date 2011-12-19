//#####################################################################
// Copyright 2004-2007, Unnur Gretarsdottir, Geoffrey Irving, Mike Rodgers, Craig Schroeder, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERIC_EXAMPLE
//#####################################################################
#ifndef __GENERIC_EXAMPLE__
#define __GENERIC_EXAMPLE__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "../Dynamic/DYNAMIC_IMPLICIT_SURFACE.h"
#include "../MESHING_EXAMPLE.h"
#include "../Sphere/SPHERE_IMPLICIT_SURFACE.h"
namespace PhysBAM{

template<class T>
class TETRAHEDRALIZED_VOLUME_FIELD
{
    typedef VECTOR<T,3> TV;
public:
    TETRAHEDRON_MESH mesh;
    PARTICLES<TV> particles;
    TETRAHEDRALIZED_VOLUME<T> field_volume; 
    ARRAY<T> field;
    T default_value;
private:
    T thickness;
    T max_depth;
public:

    TETRAHEDRALIZED_VOLUME_FIELD(const STREAM_TYPE stream_type,const std::string& filename)
        :field_volume(mesh,particles),default_value(0)
    {
        FILE_UTILITIES::Read_From_File(stream_type,filename,field_volume,field);
        field_volume.Update_Tetrahedron_List();
        field_volume.Initialize_Hierarchy();
        field_volume.Initialize_Triangulated_Surface();
        field_volume.triangulated_surface->Update_Triangle_List();
        field_volume.triangulated_surface->Initialize_Hierarchy();
        mesh.Initialize_Adjacent_Elements();
        thickness=field_volume.Minimum_Altitude()/100;
        max_depth=2*field_volume.triangulated_surface->Maximum_Edge_Length();
    }

    T operator()(const TV& X) const
    {PHYSBAM_ASSERT(field_volume.tetrahedron_list && field_volume.triangulated_surface);
    // first search for tetrahedron containing X
    ARRAY<int> intersection_list;field_volume.hierarchy->Intersection_List(X,intersection_list,thickness);
    for(int i=1;i<=intersection_list.m;i++){int t=intersection_list(i);
        if((*field_volume.tetrahedron_list)(t).Inside(X,thickness)){
            const VECTOR<int,4>& nodes=field_volume.mesh.elements(t);
            TV w=(*field_volume.tetrahedron_list)(t).First_Three_Barycentric_Coordinates(X);
            return w[1]*field(nodes[1])+w[2]*field(nodes[2])+w[3]*field(nodes[3])+(1-w.Sum())*field(nodes[4]);}}
    // if that fails, project X onto the surface
    TRIANGULATED_SURFACE<T>& surface=*field_volume.triangulated_surface;
    PHYSBAM_ASSERT(surface.triangle_list);
    int t;T distance;TV surface_X=surface.Surface(X,max_depth,thickness,&t,&distance);
    if(distance>=max_depth) return default_value;
    const VECTOR<int,3>& nodes=surface.mesh.elements(t);
    TV w=(*surface.triangle_list)(t).Barycentric_Coordinates(surface_X);
    return w[1]*field(nodes[1])+w[2]*field(nodes[2])+w[3]*field(nodes[3]);}
};

template<class T>
class TETRAHEDRALIZED_VOLUME_SURFACE_FIELD
{
    typedef VECTOR<T,3> TV;
public:
    TETRAHEDRON_MESH mesh;
    PARTICLES<TV> particles;
    TETRAHEDRALIZED_VOLUME<T> field_volume; 
    ARRAY<T> field;
    T default_value;
private:
    T thickness;
    T max_depth_mult;
    T max_depth;
public:

    TETRAHEDRALIZED_VOLUME_SURFACE_FIELD(const STREAM_TYPE stream_type,const std::string& filename)
        :field_volume(mesh,particles),default_value(0)
    {
        FILE_UTILITIES::Read_From_File(stream_type,filename,field_volume,field);
        field_volume.Update_Tetrahedron_List();
        //field_volume.Initialize_Hierarchy();
        field_volume.Initialize_Triangulated_Surface();
        field_volume.triangulated_surface->Update_Triangle_List();
        field_volume.triangulated_surface->Initialize_Hierarchy();
        mesh.Initialize_Adjacent_Elements();
        thickness=field_volume.Minimum_Altitude()/100;
        max_depth_mult=0.5;
        max_depth=max_depth_mult*field_volume.triangulated_surface->Maximum_Edge_Length();
    }

    void Set_Max_Depth_Mult(T mult=0.5)
    {max_depth_mult=mult;max_depth=max_depth_mult*field_volume.triangulated_surface->Maximum_Edge_Length();}

    T operator()(const TV& X) const
    {PHYSBAM_ASSERT(field_volume.tetrahedron_list && field_volume.triangulated_surface);
    // project X onto the surface
    TRIANGULATED_SURFACE<T>& surface=*field_volume.triangulated_surface;
    PHYSBAM_ASSERT(surface.triangle_list);
    int t;T distance;TV surface_X=surface.Surface(X,max_depth,thickness,&t,&distance);
    if(distance>=max_depth) return default_value;
    const VECTOR<int,3>& nodes=surface.mesh.elements(t);
    TV w=(*surface.triangle_list)(t).Barycentric_Coordinates(surface_X);
    return w[1]*field(nodes[1])+w[2]*field(nodes[2])+w[3]*field(nodes[3]);}
};

template<class T> class GENERIC_EXAMPLE;

template<class T>
class GENERIC_EXAMPLE_EXTRA_REFINEMENT_CRITERIA_HELPER:public EXTRA_REFINEMENT_CRITERIA_HELPER
{
public:
    GENERIC_EXAMPLE<T>& example;

    GENERIC_EXAMPLE_EXTRA_REFINEMENT_CRITERIA_HELPER(GENERIC_EXAMPLE<T>& _example):
        example(_example){};
    int operator()(int index)
    {
        return example.Extra_Refinement_Criteria(index);
    }
};

template<class T>
class GENERIC_EXAMPLE:public MESHING_EXAMPLE<T>
{
    typedef VECTOR<T,3> TV;
public:
    typedef MESHING_EXAMPLE<T> BASE;
    using BASE::tetrahedral_meshing;using BASE::implicit_surface_filename;
    using BASE::use_octree_implicit_surface;using BASE::bcc_lattice_cell_size;using BASE::use_adaptive_refinement;
    using BASE::max_subdivision_levels;using BASE::use_optimization;using BASE::use_dynamics;
    using BASE::number_of_initial_optimization_steps;using BASE::number_of_final_optimization_steps;
    using BASE::force_attraction_coefficient;using BASE::velocity_attraction_coefficient;
    using BASE::allow_tangential_velocity_slip;using BASE::use_finite_element_forces;using BASE::youngs_modulus;
    using BASE::poissons_ratio;using BASE::rayleigh_coefficient;using BASE::edge_spring_stiffness;
    using BASE::edge_spring_overdamping_fraction;using BASE::altitude_spring_stiffness;
    using BASE::altitude_spring_overdamping_fraction;using BASE::time_step;using BASE::number_of_force_steps;
    using BASE::number_of_velocity_steps;using BASE::use_multiple_levelset_implicit_surface;
    using BASE::use_global_quality_criteria_for_early_exit;using BASE::replace_green_refinement_with_embedded_t_junctions;
    using BASE::allow_coarsening_to_non_graded_mesh;using BASE::use_aggressive_tet_pruning;

    std::string data_directory;
    std::string expression;
    RANGE<TV> user_box;

    TETRAHEDRALIZED_VOLUME_FIELD<T>* variable_max_edge_length_tet_field;
    TETRAHEDRALIZED_VOLUME_SURFACE_FIELD<T>* variable_max_edge_length_tetsurf_field;
    bool subdivide_tets_near_boundary_only;
    T torus_outer_radius,torus_inner_radius;
    T snake_radius,snake_length;
    ARRAY<IMPLICIT_OBJECT<TV>*> levelsets;

    GENERIC_EXAMPLE(const STREAM_TYPE stream_type,const std::string& parameter_file)
        :BASE(stream_type),data_directory("../../Public_Data"),variable_max_edge_length_tet_field(0),variable_max_edge_length_tetsurf_field(0)
    {
        PARAMETER_LIST parameter_list;parameter_list.Begin_Parse(parameter_file);

        // levelset
        implicit_surface_filename=parameter_list.Get_Parameter("levelset","<unknown>");
        use_octree_implicit_surface=parameter_list.Get_Parameter("use_octree_implicit_surface",false);
        use_multiple_levelset_implicit_surface=parameter_list.Get_Parameter("use_multiple_levelset_implicit_surface",false);
        if(implicit_surface_filename=="<unknown>"){std::cerr<<"No levelset specified."<<std::endl;exit(1);}

        // parameters
        use_optimization=parameter_list.Get_Parameter("use_optimization",false);
        use_dynamics=parameter_list.Get_Parameter("use_dynamics",true);
        bcc_lattice_cell_size=parameter_list.Get_Parameter("bcc_lattice_cell_size",(T)0);
        use_adaptive_refinement=parameter_list.Get_Parameter("use_adaptive_refinement",true);
        max_subdivision_levels=parameter_list.Get_Parameter("max_subdivision_levels",7);
        T attraction_coefficient=parameter_list.Get_Parameter("attraction_coefficient",(T).1);
        force_attraction_coefficient=parameter_list.Get_Parameter("force_attraction_coefficient",attraction_coefficient);
        velocity_attraction_coefficient=parameter_list.Get_Parameter("velocity_attraction_coefficient",attraction_coefficient);
        tetrahedral_meshing.dynamic_ether_viscosity=parameter_list.Get_Parameter("dynamic_ether_viscosity",(T).5);
        tetrahedral_meshing.use_constant_mass=parameter_list.Get_Parameter("use_constant_mass",true);
        tetrahedral_meshing.boundary_mass_multiplicative_factor=parameter_list.Get_Parameter("boundary_mass_multiplicative_factor",(T)1);
        time_step=parameter_list.Get_Parameter("time_step",(T).25);
        number_of_force_steps=parameter_list.Get_Parameter("number_of_force_steps",0);
        number_of_velocity_steps=parameter_list.Get_Parameter("number_of_velocity_steps",360);
        number_of_initial_optimization_steps=parameter_list.Get_Parameter("number_of_initial_optimization_steps",0);
        number_of_final_optimization_steps=parameter_list.Get_Parameter("number_of_final_optimization_steps",0);
        allow_tangential_velocity_slip=parameter_list.Get_Parameter("allow_tangential_velocity_slip",false);
        tetrahedral_meshing.curvature_subdivision_threshold=parameter_list.Get_Parameter("curvature_subdivision_threshold",(T).7);
        tetrahedral_meshing.curvature_subdivision_threshold*=parameter_list.Get_Parameter("curvature_subdivision_threshold_multiplier",(T)1);
        tetrahedral_meshing.interpolation_error_subdivision_threshold=parameter_list.Get_Parameter("interpolation_error_subdivision_threshold",(T).08);
        tetrahedral_meshing.interpolation_error_subdivision_threshold*=parameter_list.Get_Parameter("interpolation_error_subdivision_threshold_multiplier",(T)1);
        tetrahedral_meshing.maximum_boundary_edge_length=parameter_list.Get_Parameter("maximum_boundary_edge_length",(T)FLT_MAX);
        edge_spring_stiffness=parameter_list.Get_Parameter("edge_spring_stiffness",(T)1e-4);
        altitude_spring_stiffness=parameter_list.Get_Parameter("altitude_spring_stiffness",(T)1e-4);
        T spring_stiffness_multiplier=parameter_list.Get_Parameter("spring_stiffness_multiplier",(T)1);
        use_global_quality_criteria_for_early_exit=parameter_list.Get_Parameter("use_global_quality_criteria_for_early_exit",true);
        replace_green_refinement_with_embedded_t_junctions=parameter_list.Get_Parameter("replace_green_refinement_with_embedded_t_junctions",false);
        allow_coarsening_to_non_graded_mesh=parameter_list.Get_Parameter("allow_coarsening_to_non_graded_mesh",false);
        edge_spring_stiffness*=spring_stiffness_multiplier;altitude_spring_stiffness*=spring_stiffness_multiplier;
        use_aggressive_tet_pruning=parameter_list.Get_Parameter("use_aggressive_tet_pruning",false);
        torus_inner_radius=parameter_list.Get_Parameter("torus_inner_radius",(T)1);
        torus_outer_radius=parameter_list.Get_Parameter("torus_outer_radius",(T)2);
        snake_radius=parameter_list.Get_Parameter("snake_radius",(T).5);
        snake_length=parameter_list.Get_Parameter("snake_length",(T)10);
        tetrahedral_meshing.symmetric_initial_grid=parameter_list.Get_Parameter("symmetric_initial_grid",false);

        expression=parameter_list.Get_Parameter("expression","0");
        user_box.min_corner.x=parameter_list.Get_Parameter("box_min_x",-(T)1);
        user_box.min_corner.y=parameter_list.Get_Parameter("box_min_y",-(T)1);
        user_box.min_corner.z=parameter_list.Get_Parameter("box_min_z",-(T)1);
        user_box.max_corner.x=parameter_list.Get_Parameter("box_max_x",(T)1);
        user_box.max_corner.y=parameter_list.Get_Parameter("box_max_y",(T)1);
        user_box.max_corner.z=parameter_list.Get_Parameter("box_max_z",(T)1);

        // extra refinement criteria
        tetrahedral_meshing.extra_refinement_criteria=new GENERIC_EXAMPLE_EXTRA_REFINEMENT_CRITERIA_HELPER<T>(*this);
        std::string variable_max_edge_length_file=parameter_list.Get_Parameter("variable_max_edge_length_file","");
        if(!variable_max_edge_length_file.empty()){
            subdivide_tets_near_boundary_only=parameter_list.Get_Parameter("variable_max_edge_length_boundary_only",false);
            if(subdivide_tets_near_boundary_only){
                variable_max_edge_length_tetsurf_field=new TETRAHEDRALIZED_VOLUME_SURFACE_FIELD<T>(stream_type,variable_max_edge_length_file);
                variable_max_edge_length_tetsurf_field->field*=parameter_list.Get_Parameter("variable_max_edge_length_scale",(T)1);
                variable_max_edge_length_tetsurf_field->field+=parameter_list.Get_Parameter("variable_max_edge_length_offset",(T)0);
                variable_max_edge_length_tetsurf_field->default_value=(T)FLT_MAX;
                LOG::cout<<"variable max edge length field = ["<<variable_max_edge_length_tetsurf_field->field.Min()
                         <<", "<<variable_max_edge_length_tetsurf_field->field.Max()
                         <<"], subdivide tetrahedra near boundary only: "<<subdivide_tets_near_boundary_only<<std::endl;}
            else{
                variable_max_edge_length_tet_field=new TETRAHEDRALIZED_VOLUME_FIELD<T>(stream_type,variable_max_edge_length_file);
                variable_max_edge_length_tet_field->field*=parameter_list.Get_Parameter("variable_max_edge_length_scale",(T)1);
                variable_max_edge_length_tet_field->field+=parameter_list.Get_Parameter("variable_max_edge_length_offset",(T)0);
                variable_max_edge_length_tet_field->default_value=(T)FLT_MAX;
                LOG::cout<<"variable max edge length field = ["<<variable_max_edge_length_tet_field->field.Min()
                         <<", "<<variable_max_edge_length_tet_field->field.Max()
                         <<"], subdivide tetrahedra near boundary only: "<<subdivide_tets_near_boundary_only<<std::endl;}}
        parameter_list.End_Parse();
    }

    ~GENERIC_EXAMPLE()
    {levelsets.Delete_Pointers_And_Clean_Memory();delete variable_max_edge_length_tet_field;delete variable_max_edge_length_tetsurf_field;}

    void Initialize_Implicit_Surface() PHYSBAM_OVERRIDE
    {if(implicit_surface_filename=="<sphere>"){
        T cell_size=bcc_lattice_cell_size/(1<<max_subdivision_levels);
        tetrahedral_meshing.Initialize(new SPHERE_IMPLICIT_SURFACE<T>(cell_size));
        return;}
    if(implicit_surface_filename=="<torus>"){
        ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >& torus=*new ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >(TORUS<T>(TV(),TV(0,0,1),torus_inner_radius,torus_outer_radius));
        torus.cell_size=bcc_lattice_cell_size/(1<<max_subdivision_levels);
        tetrahedral_meshing.Initialize(&torus);
        return;}
    if(implicit_surface_filename=="<snake>"){
        T half_cylinder_length=(T).5*snake_length-snake_radius;
        ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* sphere1=new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(SPHERE<TV>(TV(-half_cylinder_length,0,0),snake_radius));
        ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* sphere2=new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(SPHERE<TV>(TV(half_cylinder_length,0,0),snake_radius));
        ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >* cylinder=new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(CYLINDER<T>(TV(-half_cylinder_length,0,0),TV(half_cylinder_length,0,0),snake_radius));
        sphere1->cell_size=sphere2->cell_size=cylinder->cell_size=bcc_lattice_cell_size/(1<<max_subdivision_levels);
        levelsets.Append(sphere1);
        levelsets.Append(sphere2);
        levelsets.Append(cylinder);
        tetrahedral_meshing.Initialize(new MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>(&levelsets));
        return;}
    if(implicit_surface_filename=="<dynamic>"){
        DYNAMIC_IMPLICIT_SURFACE<T>* dis=new DYNAMIC_IMPLICIT_SURFACE<T>(user_box);
        dis->Set_Expression(expression);
        tetrahedral_meshing.Initialize(dis);
        return;}
    if(!FILE_UTILITIES::File_Exists(implicit_surface_filename)){
        std::string new_filename=data_directory+"/"+implicit_surface_filename;
        if(!FILE_UTILITIES::File_Exists(new_filename)){std::cerr<<"Can't find "<<implicit_surface_filename<<std::endl;exit(1);}
        implicit_surface_filename=new_filename;}
    LOG::Time("reading levelset");
    LOG::cout<<"levelset file = "<<implicit_surface_filename<<std::endl;
    MESHING_EXAMPLE<T>::Initialize_Implicit_Surface();
    LOG::Stop_Time();}

    int Test(int) const
    {return 0;}

    int Extra_Refinement_Criteria(const int tet_index) const
    {TETRAHEDRALIZED_VOLUME<T> &target_volume=tetrahedral_meshing.solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
     if(!target_volume.tetrahedron_list)target_volume.Update_Tetrahedron_List();
     TV centroid=target_volume.Centroid(tet_index);
     if(variable_max_edge_length_tetsurf_field){
         T tet_min_edge_length=(*target_volume.tetrahedron_list)(tet_index).Minimum_Edge_Length();
         T target_length=(*variable_max_edge_length_tetsurf_field)(centroid);
         return tet_min_edge_length > target_length;
     }else if(variable_max_edge_length_tet_field){
         T tet_min_edge_length=(*target_volume.tetrahedron_list)(tet_index).Minimum_Edge_Length();
         T target_length=(*variable_max_edge_length_tet_field)(centroid);
         return tet_min_edge_length > target_length;}
     return 0;}
//#####################################################################
};
}
#endif

