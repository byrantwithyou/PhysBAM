#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_IMPLICIT_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string>
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"
#include "../../Public_Library/Level_Sets/OCTREE_LEVELSET.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"
#include <Collisions_And_Interactions/FLUID_COLLISION_BODY_LIST_DYADIC.h>
#include <Level_Sets/OCTREE_REMOVED_PARTICLE_PROCESSING.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>

using namespace PhysBAM;

typedef float T;
typedef float RW;

bool Coarsen_Cell(const OCTREE_CELL<T>& cell,const ARRAY<T>& phi,const T band_width_times_small_dx)
{
    T dx=(T)1.01*band_width_times_small_dx+(T)root_three*cell.DX().x;
    for(int i=0;i<8;i++) if(fabs(phi(cell.Node(i))) > dx) return true;
    return false;
}

bool Coarsen_Tree(OCTREE_CELL<T>* cell,const ARRAY<T>& phi,const T refinement_distance)
{
    if(!cell->Has_Children()){return Coarsen_Cell(*cell,phi,refinement_distance);}
    bool can_coarsen=true;for(int i=0;i<8;i++) can_coarsen&=Coarsen_Tree(cell->Child(i),phi,refinement_distance);
    if(can_coarsen){cell->Delete_Children();return Coarsen_Cell(*cell,phi,refinement_distance);}
    else return false;
}

void Fast_March(OCTREE_GRID<T>& octree_grid,ARRAY<T>& octree_phi,PARSE_ARGS& parse_args)
{
    OCTREE_LEVELSET<T> union_octree_levelset(octree_grid,octree_phi);
    if(!parse_args.Get_Option_Value("-no_fmm")){
        int fmm_band=parse_args.Get_Integer_Value("-band");
        std::cout<<"Fast_Marching_Method(band="<<fmm_band<<")..."<<std::flush;
        int fast_marching_method_timer=TIMER::Singleton()->Register_Timer();
        union_octree_levelset.Fast_Marching_Method(0,fmm_band*octree_grid.Minimum_Edge_Length());
        std::cout<<"done"<<std::endl;
        std::cout<<"time for fast_marching_method:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(fast_marching_method_timer)<<std::endl;
    }
}

void Fast_March_And_Write(OCTREE_GRID<T>& octree_grid,ARRAY<T>& octree_phi,const std::string& filename,PARSE_ARGS& parse_args)
{
    OCTREE_LEVELSET<T> union_octree_levelset(octree_grid,octree_phi);
    if(!parse_args.Get_Option_Value("-no_fmm")){
        int fmm_band=parse_args.Get_Integer_Value("-band");
        std::cout<<"Fast_Marching_Method(band="<<fmm_band<<")..."<<std::flush;
        int fast_marching_method_timer=TIMER::Singleton()->Register_Timer();
        union_octree_levelset.Fast_Marching_Method(0,fmm_band*octree_grid.Minimum_Edge_Length());
        std::cout<<"done"<<std::endl;
        std::cout<<"time for fast_marching_method:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(fast_marching_method_timer)<<std::endl;
    }

    std::cout<<"Writing to "<<filename<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(filename,octree_phi);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-s", .25, "scale", "scale particle radius");
    parse_args.Add_Integer_Argument("-d", 2, "depth", "octree maximum depth");
    parse_args.Add_Integer_Argument("-band", 0, "fmm_band", "fast marching bandwidth");
    parse_args.Add_Integer_Argument("-phi_clamp",0); // clamp input phi instead of letting it be a full signed distance
    parse_args.Add_Option_Argument("-phi_clamp_test"); // clamp input phi instead of letting it be a full signed distance
    parse_args.Add_Option_Argument("-no_fmm");
    parse_args.Add_Option_Argument("-no_merge");
    parse_args.Add_Option_Argument("-coarsen");
    parse_args.Add_Option_Argument("-use_velocity_scaling");
    parse_args.Add_Option_Argument("-preserve_volume");
    parse_args.Add_Double_Argument("-dt", (double)1/60, "dt", "dt (time step)");
    parse_args.Add_Option_Argument("-write_merged");
    parse_args.Add_Option_Argument("-write_particle_levelset");
    parse_args.Add_Option_Argument("-write_refined_levelset");
    parse_args.Add_Option_Argument("-write_raw_particles");
    parse_args.Add_Option_Argument("-write_union");
    parse_args.Add_Double_Argument("-blend_test",3);
    parse_args.Add_String_Argument("-o", "", "directory", "output directory");
    parse_args.Add_Double_Argument("-b", .8, "blending", "blending parameter");
    parse_args.Add_Double_Argument("-p", 1, "power", "particle power");
    parse_args.Add_String_Argument("-deformable_objects","");
    parse_args.Add_String_Argument("-rigid_bodies","");
    parse_args.Add_Vector_3D_Argument("-subdomain_minimum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_maximum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_2_minimum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_2_maximum_corner",VECTOR<double,3>());
    parse_args.Set_Extra_Arguments(2, "<input directory> <frame number>");
    parse_args.Parse(argc,argv);
    if (parse_args.Num_Extra_Args() != 2) 
    {
        parse_args.Print_Usage();
        return 1;
    }

    bool merge=!parse_args.Get_Option_Value("-no_merge");
    std::string input_directory(parse_args.Extra_Arg(1));
    std::string output_directory=input_directory;
    if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");

    int frame=atoi(parse_args.Extra_Arg(2).c_str());
    std::string particle_filename=STRING_UTILITIES::string_sprintf("%s/removed_negative_particles.%d",input_directory.c_str(),frame);
    std::string octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/octree_levelset.%d",input_directory.c_str(),frame);
    std::string octree_grid_filename=STRING_UTILITIES::string_sprintf("%s/octree_grid.%d",input_directory.c_str(),frame);
    std::string merged_octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/merged_octree_levelset.%d",output_directory.c_str(),frame);
    std::string merged_octree_grid_filename=STRING_UTILITIES::string_sprintf("%s/merged_octree_grid.%d",output_directory.c_str(),frame);
    std::string particles_only_octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/particles_octree_levelset.%d",output_directory.c_str(),frame);
    std::string refined_octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/refined_octree_levelset.%d",output_directory.c_str(),frame);
    std::string raw_particles_octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/raw_particles_octree_levelset.%d",output_directory.c_str(),frame);
    std::string union_octree_levelset_filename=STRING_UTILITIES::string_sprintf("%s/union_octree_levelset.%d",output_directory.c_str(),frame);

    std::cout<<"Reading particles"<<std::endl;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >*> particles_array;
    FILE_UTILITIES::Read_From_File<RW>(particle_filename,particles_array);

    if(parse_args.Is_Value_Set("-subdomain_minimum_corner") && parse_args.Is_Value_Set("-subdomain_maximum_corner")){
        BOX_3D<T> subdomain(VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_minimum_corner")),VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_maximum_corner")));
        BOX_3D<T> subdomain2=subdomain;
        if(parse_args.Is_Value_Set("-subdomain_2_minimum_corner") && parse_args.Is_Value_Set("-subdomain_2_maximum_corner")){
            subdomain2=BOX_3D<T>(VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_2_minimum_corner")),VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_2_maximum_corner")));}
        std::cout << "Pruning particles outside " << subdomain << " and " << subdomain2 << std::endl;
        int count=0;
        for(int i=0;i<particles_array.m;i++) if(particles_array(i)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles=particles_array(i);
            for(int p=particles->number;p>=1;p--){if(subdomain.Lazy_Outside(particles->X(p)) && subdomain2.Lazy_Outside(particles->X(p))){count++;particles->Delete_Particle(p);}}
            if(!particles->number){delete particles;particles_array(i)=0;}}
        std::cout << "Deleted " << count << " particles" << std::endl;
    }

    std::cout<<"Reading octree grid"<<std::flush;
    OCTREE_GRID<T> octree_grid;
    FILE_UTILITIES::Read_From_File<RW>(octree_grid_filename,octree_grid);
    std::cout<<" (original maximum_depth=" << octree_grid.maximum_depth << std::endl;

    std::cout<<"Reading levelset"<<std::endl;
    ARRAY<T> water_phi;
    FILE_UTILITIES::Read_From_File<RW>(octree_levelset_filename,water_phi);

    T minimum_surface_roughness=1e-5;
    RIGID_BODY_LIST_3D<T> rigid_body_list;
    DEFORMABLE_OBJECT_LIST_3D<T> deformable_object_list;
    FLUID_COLLISION_BODY_LIST_DYADIC<T,OCTREE_GRID<T> > collision_body_list(octree_grid);
    if(parse_args.Is_Value_Set("-rigid_bodies")){
        rigid_body_list.Read<T>(STRING_UTILITIES::string_sprintf("%s/",input_directory.c_str()),frame,true,false);
        ARRAY<int> id_list;STRING_UTILITIES::Parse_Integer_List(parse_args.Get_String_Value("-rigid_bodies"),id_list);
        for(int i=0;i<id_list.m;i++){int id=id_list(i);
            std::cout << "Adding rigid body " << id << " to collision body list" << std::endl;
            rigid_body_list(id)->surface_roughness=minimum_surface_roughness;
            collision_body_list.Add_Body(rigid_body_list(id));}
    }
    if(parse_args.Is_Value_Set("-deformable_objects")){
        std::string prefix=input_directory+"/";
        deformable_object_list.Read_Static_Variables<T>(prefix);
        deformable_object_list.Read_Dynamic_Variables<T>(prefix,frame);
        ARRAY<int> id_list;STRING_UTILITIES::Parse_Integer_List(parse_args.Get_String_Value("-deformable_objects"),id_list);
        for(int i=0;i<id_list.m;i++){int id=id_list(i);
            std::cout << "Adding deformable object " << id << " to collision body list" << std::endl;
            deformable_object_list(id).Initialize_Collision_Geometry();
            deformable_object_list(id).collisions.triangulated_surface->Update_Bounding_Box();
            deformable_object_list(id).collisions.triangulated_surface->Initialize_Triangle_Hierarchy();
            deformable_object_list(id).collisions.triangulated_surface->Update_Triangle_List();
            deformable_object_list(id).collisions.roughness=minimum_surface_roughness;
            collision_body_list.Add_Body(&deformable_object_list(id).collisions);}
    }

    OCTREE_REMOVED_PARTICLE_PROCESSING<T> removed_particle_processing(octree_grid,water_phi,particles_array);
    removed_particle_processing.blending_parameter=(T)parse_args.Get_Double_Value("-b");
    removed_particle_processing.scale=(T)parse_args.Get_Double_Value("-s");
    removed_particle_processing.octree_maximum_depth=parse_args.Get_Integer_Value("-d");
    removed_particle_processing.particle_power=(T)parse_args.Get_Double_Value("-p");
    removed_particle_processing.use_velocity_scaling=parse_args.Get_Option_Value("-use_velocity_scaling");
    removed_particle_processing.preserve_volume=parse_args.Get_Option_Value("-preserve_volume");
    removed_particle_processing.dt=(T)parse_args.Get_Double_Value("-dt");

    if(collision_body_list.collision_bodies.m) removed_particle_processing.Set_Collision_Aware(&collision_body_list);

    removed_particle_processing.Refine_And_Create_Particle_Phi();

    if(parse_args.Get_Option_Value("-write_refined_levelset")){
        std::cout<<"Writing refined water levelset"<<std::endl;
        FILE_UTILITIES::Write_To_File<T>(refined_octree_levelset_filename,water_phi);}

    if(parse_args.Get_Option_Value("-write_raw_particles")){
        std::cout<<"Writing raw particles field"<<std::endl;
        FILE_UTILITIES::Write_To_File<T>(raw_particles_octree_levelset_filename,removed_particle_processing.particle_octree_phi);}

    std::cout<<"Writing octree grid file"<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(merged_octree_grid_filename,octree_grid);

    if(parse_args.Is_Value_Set("-phi_clamp")){
        T phi_clamp=octree_grid.Minimum_Edge_Length()*parse_args.Get_Integer_Value("-phi_clamp");
        std::cout << "clamping phi: " << phi_clamp << std::endl;
        for(int i=0;i<octree_grid.number_of_nodes;i++) water_phi(i)=clamp(water_phi(i),-phi_clamp,phi_clamp);
    }
    else if(parse_args.Is_Value_Set("-phi_clamp_test")){
        T phi_clamp=octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*(T)removed_particle_processing.particle_power;
        std::cout << "clamping phi: " << phi_clamp << std::endl;
        for(int i=0;i<octree_grid.number_of_nodes;i++) water_phi(i)=clamp(water_phi(i),-phi_clamp,phi_clamp);
    }

    ARRAY<T>& particle_octree_phi=removed_particle_processing.particle_octree_phi;

    if(merge){
        std::cout<<"Merging levelset and particle levelset"<<std::endl;
        ARRAY<T> octree_phi;
        removed_particle_processing.Merge_Phi(octree_phi);

        OCTREE_LEVELSET<T> merged_octree_levelset(octree_grid,octree_phi);
        if(!parse_args.Get_Option_Value("-no_fmm")){
            int fmm_band=parse_args.Get_Integer_Value("-band");
            std::cout<<"merged_octree_levelset.Fast_Marching_Method(band="<<fmm_band<<")..."<<std::flush;
            int fast_marching_method_timer=TIMER::Singleton()->Register_Timer();
            merged_octree_levelset.Fast_Marching_Method(0,fmm_band*octree_grid.Minimum_Edge_Length());
            std::cout<<"done"<<std::endl;
            std::cout<<"time for fast_marching_method:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(fast_marching_method_timer)<<std::endl;
        }

        // coarsening
        if(parse_args.Get_Option_Value("-coarsen")){
            std::cout << "Coarsening" << std::endl;

            // TODO: Should we base refinement distance on new min cell size or old one???
            T levelset_refinement_bandwidth=6;
            T refinement_distance=(T).5*levelset_refinement_bandwidth*octree_grid.Minimum_Edge_Length();

            int old_number_of_cells=octree_grid.number_of_cells,old_number_of_nodes=octree_grid.number_of_nodes,old_number_of_faces=octree_grid.number_of_faces;

            for(int i=0;i<octree_grid.uniform_grid.m;i++)for(int j=0;j<octree_grid.uniform_grid.n;j++)for(int ij=0;ij<octree_grid.uniform_grid.mn;ij++)
                if(octree_grid.cells(i,j,ij)->Has_Children()) Coarsen_Tree(octree_grid.cells(i,j,ij),octree_phi,refinement_distance);

            ARRAY<int> cell_mapping_array,node_mapping_array,face_mapping_array;
            octree_grid.Compact_Array_Indices(&cell_mapping_array,&node_mapping_array,&face_mapping_array);

            ARRAY<T> temp(max(old_number_of_cells,old_number_of_nodes,old_number_of_faces),false);  
            ARRAY<T>::Compact_Array_Using_Compaction_Array(octree_phi,node_mapping_array,&temp);
            octree_phi.Resize(octree_grid.number_of_nodes);
            ARRAY<T>::Compact_Array_Using_Compaction_Array(particle_octree_phi,node_mapping_array,&temp);
            particle_octree_phi.Resize(octree_grid.number_of_nodes);
            ARRAY<T>::Compact_Array_Using_Compaction_Array(water_phi,node_mapping_array,&temp);
            water_phi.Resize(octree_grid.number_of_nodes);

            octree_grid.Tree_Topology_Changed();
        }

        std::cout<<"Writing octree levelset file"<<std::endl;
        FILE_UTILITIES::Write_To_File<T>(merged_octree_levelset_filename,merged_octree_levelset.phi);
    }

    if(parse_args.Get_Option_Value("-write_particle_levelset")){
        std::cout<<"Generating particle-only levelset"<<std::endl;
        ARRAY<T> octree_phi(particle_octree_phi);
        octree_phi+=octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
        Fast_March_And_Write(octree_grid,octree_phi,particles_only_octree_levelset_filename,parse_args);}

    // This is the one used for the thin shells fluid coupling paper's funnel example!
    if(parse_args.Get_Option_Value("-write_union")){
        std::cout<<"Generating union levelset"<<std::endl;
   
        // Important to first fast march the particles phi before taking min with water phi because away from its interface the particle phi
        // might have a positive value close to zero. (Without FMM it used to create fringe lines/stairstepping in resulting unioned levelset)
        ARRAY<T> octree_phi(particle_octree_phi);
        octree_phi+=octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
        Fast_March(octree_grid,octree_phi,parse_args);
        for(int i=0;i<octree_grid.number_of_nodes;i++) octree_phi(i)=min(octree_phi(i),water_phi(i));
        Fast_March_And_Write(octree_grid,octree_phi,union_octree_levelset_filename,parse_args);
    }

    if(parse_args.Is_Value_Set("-blend_test")){
        std::cout<<"BLEND TEST"<<std::endl;

        ARRAY<T> octree_phi(octree_grid.number_of_nodes,false);

#if 0
        // special test of plain averaging
        for(int i=0;i<octree_grid.number_of_nodes;i++){
            T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
            octree_phi(i)=0.5*(water_phi(i)+particles_only_phi);
        }

        Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/plain_average_octree_levelset.%d",output_directory.c_str(),frame),parse_args);

        // minmag_1
        for(int i=0;i<octree_grid.number_of_nodes;i++){
            T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
            octree_phi(i)=minmag(particles_only_phi,water_phi(i));
        }

        Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/minmag_1_octree_levelset.%d",output_directory.c_str(),frame),parse_args);

        // minmag_2
        for(int i=0;i<octree_grid.number_of_nodes;i++){
            T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
            T union_phi=min(particles_only_phi,water_phi(i));
            octree_phi(i)=minmag(union_phi,water_phi(i));
        }

        Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/minmag_2_octree_levelset.%d",output_directory.c_str(),frame),parse_args);

        // minmag_3
        for(int i=0;i<octree_grid.number_of_nodes;i++){
            T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
            T union_phi=min(particles_only_phi,water_phi(i));
            octree_phi(i)=union_phi;
        }

        Fast_March(octree_grid,octree_phi,parse_args);

        for(int i=0;i<octree_grid.number_of_nodes;i++){
            octree_phi(i)=minmag(octree_phi(i),water_phi(i));
        }

        Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/minmag_3_octree_levelset.%d",output_directory.c_str(),frame),parse_args);

#endif
        // special blendings
//        for(int cell=0;cell<5;cell++){
        for(int cell=5;cell<=5;cell++){
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)+particle_octree_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
     
        // special blendings
        {int cell=8;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)+particle_octree_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
        {int cell=10;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)+particle_octree_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
        {int cell=15;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)+particle_octree_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }

#if 0
        // special blendings
        for(int cell=0;cell<5;cell++){
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/water_blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
        
        // special blendings
        {int cell=8;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/water_blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
        {int cell=10;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/water_blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
        {int cell=15;
            for(int i=0;i<octree_grid.number_of_nodes;i++){
                T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;
                T alpha=clamp(water_phi(i)/(cell*octree_grid.Minimum_edge_Length()),(T)0,(T)1);
                octree_phi(i)=(1-alpha)*(water_phi(i)) + alpha*particles_only_phi;}
            Fast_March_And_Write(octree_grid,octree_phi,STRING_UTILITIES::string_sprintf("%s/water_blend_%d_octree_levelset.%d",output_directory.c_str(),cell,frame),parse_args);
        }
#endif
    }

    particles_array.Delete_Pointers_And_Clean_Memory();
    return 0;
}

