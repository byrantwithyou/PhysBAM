#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_IMPLICIT_SURFACE.h>
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"
#include "PARTICLE_BLENDER.h"
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_SURFACE_MAKER.h>

using namespace PhysBAM;

typedef float T;
typedef float RW;
typedef ARRAYS_3D<PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR_3D<T> >*> T_PARTICLES_ARRAY_3D;
typedef ARRAYS_3D<PARTICLE_LEVELSET_PARTICLES<T,VECTOR_3D<T> >*> T_NEGATIVE_PARTICLES_ARRAY_3D;

void Get_Particle_Bounding_Boxes(const T_PARTICLES_ARRAY_3D& particles_array, ARRAYS_3D<ARRAY<BOX_3D<T> > >& particles_bounding_box,double scale,PARTICLE_BLENDER<T>& particle_blender,int& number_of_particles)
{
    for(int i=particles_array.m_start;i<=particles_array.m_end;i++) for(int j=particles_array.n_start;j<=particles_array.n_end;j++) for(int ij=particles_array.mn_start;ij<=particles_array.mn_end;ij++) {
        PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR_3D<T> >* particles = particles_array(i,j,ij);
        if(particles) for(int p=1;p<=particles->number;p++) {
            T radius=scale*particles->radius(p);
            particles_bounding_box(i,j,ij).Append(particle_blender.Get_Bounding_Box(3*radius,radius,particles->X(p),particles->V(p)));
            number_of_particles++;
        }
    }
    
}

void Refine_Octree_Cell_To_Box(OCTREE_GRID<T>& octree_grid, OCTREE_CELL<T>* cell, const BOX_3D<T>& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells)
{
    if(!cell->Intersection(box)) return;
    if(cell->Depth_Of_This_Cell()>=octree_grid.maximum_depth){intersecting_cells.Append(cell);return;}
    if(!cell->Has_Children()) cell->Create_Children(octree_grid.number_of_cells,0,octree_grid.number_of_nodes,0,octree_grid.number_of_faces,0,&octree_grid);
    for(int i=0;i<8;i++) Refine_Octree_Cell_To_Box(octree_grid,cell->Child(i),box,intersecting_cells);
}

void Pad_Grid(GRID_3D<T>& grid)
{
    if(grid.number_of_cells_x%2!=0){grid.m+=1;grid.xmax+=grid.dx;grid.number_of_cells_x+=1;std::cout<<"Padded uniform grid m"<<std::endl;}
    if(grid.number_of_cells_y%2!=0){grid.n+=1;grid.ymax+=grid.dy;grid.number_of_cells_y+=1;std::cout<<"Padded uniform grid m"<<std::endl;}
    if(grid.number_of_cells_z%2!=0){grid.mn+=1;grid.zmax+=grid.dz;grid.number_of_cells_z+=1;std::cout<<"Padded uniform grid m"<<std::endl;}
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-s", .25, "scale", "scale particle radius");
    parse_args.Add_Integer_Argument("-d", 2, "depth", "octree maximum depth");
    parse_args.Add_String_Argument("-o", "output", "directory", "output directory");
    parse_args.Add_Double_Argument("-b", .8, "blending", "blending parameter");
    parse_args.Add_Double_Argument("-p", 1, "power", "particle power");
    parse_args.Add_Double_Argument("-e", 1, "object_expansion", "number of grid cells to expand the objects");
    parse_args.Set_Extra_Arguments(2, "<input directory> <frame number>");
    parse_args.Parse(argc,argv);
    if (parse_args.Num_Extra_Args() != 2) 
    {
        parse_args.Print_Usage();
        return 1;
    }
    double scale=parse_args.Get_Double_Value("-s");
    double blending_parameter=parse_args.Get_Double_Value("-b");
    int octree_maximum_depth=parse_args.Get_Integer_Value("-d");
    double particle_power=parse_args.Get_Double_Value("-p");
    double object_expansion=parse_args.Get_Double_Value("-e");
    std::string output_directory=parse_args.Get_String_Value("-o");
    std::string input_directory(parse_args.Extra_Arg(1));
    int frame=atoi(parse_args.Extra_Arg(2).c_str());
    std::string particle_filename=STRING_UTILITIES::string_sprintf("%s/removed_negative_particles.%d",input_directory.c_str(),frame);
    std::string negative_particle_filename=STRING_UTILITIES::string_sprintf("%s/negative_particles.%d",input_directory.c_str(),frame);
    std::string levelset_filename=STRING_UTILITIES::string_sprintf("%s/levelset.%d",input_directory.c_str(),frame);
    std::string deformable_object_filename=STRING_UTILITIES::string_sprintf("%s/deformable_object_1_embedded_tetrahedralized_volume_boundary_surface.%d",input_directory.c_str(),frame);

    std::cout<<"PARAMETERS:\n\tscale:"<<scale<<"\tblending:"<<blending_parameter<<"\toctree depth:"<<octree_maximum_depth<<"\tobject expantion:"<<object_expansion<<std::endl;
    PARTICLE_BLENDER<T> particle_blender(blending_parameter);

    std::cout<<"Reading particles"<<std::endl;
    T_PARTICLES_ARRAY_3D particles_array;
    FILE_UTILITIES::Read_From_File<RW>(particle_filename,particles_array);

    std::cout<<"Reading negative particles"<<std::endl;
    T_NEGATIVE_PARTICLES_ARRAY_3D negative_particles_array;
    FILE_UTILITIES::Read_From_File<RW>(negative_particle_filename,negative_particles_array);

    std::cout<<"Reading levelset"<<std::endl;
    ARRAYS<VECTOR<T,3> >* phi=new ARRAYS<VECTOR<T,3> >;GRID_3D<T>* grid=new GRID_3D<T>;
    RENDERING_IMPLICIT_SURFACE<T>* rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(*grid,*phi);
    LEVELSET_IMPLICIT_SURFACE<T> &implicit_surface = *(LEVELSET_IMPLICIT_SURFACE<T>*)rendering_implicit_surface->implicit_surface;
    FILE_UTILITIES::Read_From_File<RW>(levelset_filename,implicit_surface);

    ARRAY<TRIANGULATED_SURFACE<T>*> triangulated_surfaces;

    std::cout<<"Reading deformable object"<<std::endl;
    DEFORMABLE_OBJECT_LIST_3D<T> deformable_object_list;
    deformable_object_list.Read_Static_Variables<RW>(input_directory+"/",frame);
    deformable_object_list.Read_Dynamic_Variables<RW>(input_directory+"/",frame);
//    assert(deformable_object_list.deformable_objects.m==1);
    std::cout<<"-- Added "<<deformable_object_list.deformable_objects.m<<" objects"<<std::endl;
    for(int i=1;i<=deformable_object_list.deformable_objects.m;i++){
        DEFORMABLE_OBJECT_3D<T>& deformable_object=deformable_object_list(i);
        deformable_object.embedded_tetrahedralized_volume->tetrahedralized_volume.tetrahedron_mesh.Initialize_Neighbor_Nodes();
        FILE_UTILITIES::Read_From_File<RW>(deformable_object_filename,deformable_object.embedded_tetrahedralized_volume_boundary_surface->boundary_surface);
        triangulated_surfaces.Append(&deformable_object.embedded_tetrahedralized_volume_boundary_surface->boundary_surface);}

    std::cout<<"Reading rigid objects"<<std::endl;
    RIGID_BODY_LIST_3D<T> rigid_body_list;
    rigid_body_list.Read<RW>(input_directory+"/",frame);
    std::cout<<"-- Added "<<rigid_body_list.rigid_bodies.m<<" objects"<<std::endl;
    for(int i=6;i<=rigid_body_list.rigid_bodies.m;i++){
        RIGID_BODY_3D<T>& rigid_body=*rigid_body_list.rigid_bodies(i);//*rigid_body_list(i);
        for(int i=1;i<=rigid_body.triangulated_surface->particles.number;i++) rigid_body.triangulated_surface->particles.X(i)=rigid_body.World_Space_Point(rigid_body.triangulated_surface->particles.X(i));
        rigid_body.triangulated_surface->Update_Bounding_Box();rigid_body.triangulated_surface->Refresh_Auxiliary_Structures();
        triangulated_surfaces.Append(rigid_body.triangulated_surface);}
    
    LEVELSET_3D<T>& levelset=implicit_surface.levelset;
    ARRAYS<VECTOR<T,3> > phi_objects(*grid,3);
    ARRAYS<VECTOR<T,3> >::copy(5*grid->max_dx_dy_dz,phi_objects);
    LEVELSET_MAKER_DYADIC<T> levelset_maker;
    for(int object=1;object<=triangulated_surfaces.m;object++){
        TRIANGULATED_SURFACE<T>& boundary_surface=*triangulated_surfaces(object);
        if(boundary_surface.triangle_mesh.triangles.m==0)continue;
        boundary_surface.Refresh_Auxiliary_Structures();boundary_surface.Update_Bounding_Box();
        boundary_surface.Initialize_Triangle_Hierarchy();
        BOX_3D<T> bounding_box=*boundary_surface.bounding_box;
        if(!bounding_box.Intersection(grid->Domain()))continue;
        bounding_box.Scale_About_Center((T)1.2);
        VECTOR_3D<int> min_corner,max_corner;grid->Clamped_Index_End_Minus_One(bounding_box.Minimum_Corner(),min_corner);
        grid->Clamped_Index_End_Minus_One(bounding_box.Maximum_Corner(),max_corner);max_corner+=VECTOR_3D<int>(1,1,1);
        bounding_box.xmin=grid->x(min_corner.x);bounding_box.ymin=grid->y(min_corner.y);bounding_box.zmin=grid->z(min_corner.z);
        bounding_box.xmax=grid->x(max_corner.x);bounding_box.ymax=grid->y(max_corner.y);bounding_box.zmax=grid->z(max_corner.z);
        GRID_3D<T> rasterization_grid(max_corner.x-min_corner.x,max_corner.y-min_corner.y,max_corner.z-min_corner.z,bounding_box);
        ARRAYS<VECTOR<T,3> > phi_object(rasterization_grid);
        if(levelset_maker.Compute_Level_Set(boundary_surface,rasterization_grid,phi_object)){
            for(int ii=1;ii<=rasterization_grid.m;ii++)for(int jj=1;jj<=rasterization_grid.n;jj++)for(int iijj=1;iijj<=rasterization_grid.mn;iijj++){
                int i=ii+min_corner.x-1,j=jj+min_corner.y-1,ij=iijj+min_corner.z-1;
                if(phi_object(ii,jj,iijj)<phi_objects(i,j,ij)){
                    phi_objects(i,j,ij)=phi_object(ii,jj,iijj);}}}}

    {RENDERING_IMPLICIT_SURFACE<T>* rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(*grid,phi_objects);
    LEVELSET_IMPLICIT_SURFACE<T> &implicit_surface = *(LEVELSET_IMPLICIT_SURFACE<T>*)rendering_implicit_surface->implicit_surface;
    std::cout<<"Writing object levelset file"<<std::endl;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/object_levelset.%d",output_directory.c_str(),frame).c_str(),implicit_surface);}

    phi_objects-=object_expansion*grid->max_dx_dy_dz;
    for(int i=1;i<=grid->m;i++)for(int j=1;j<=grid->n;j++)for(int ij=1;ij<=grid->mn;ij++)
        levelset.phi(i,j,ij)=min(levelset.phi(i,j,ij),phi_objects(i,j,ij));

    std::cout<<"Adding negative particles to removed negative particles"<<std::endl;
    for(int i=negative_particles_array.m_start;i<=negative_particles_array.m_end;i++) for(int j=negative_particles_array.n_start;j<=negative_particles_array.n_end;j++) for(int ij=negative_particles_array.mn_start;ij<=negative_particles_array.mn_end;ij++) {
        PARTICLE_LEVELSET_PARTICLES<T,VECTOR_3D<T> >* particles = negative_particles_array(i,j,ij);
        if(particles) for(int p=1;p<=particles->number;p++){
            if(implicit_surface.levelset.Phi(particles->X(p))>0) {
                if(!particles_array(i,j,ij)) particles_array(i,j,ij)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR_3D<T> >();
                PARTICLES::Move_Particle(*particles,*particles_array(i,j,ij),p);}}}

    std::cout<<"Creating octree grid"<<std::endl;
    Pad_Grid(implicit_surface.levelset.grid);
    OCTREE_GRID<T> octree_grid;
    octree_grid.Initialize(implicit_surface.levelset.grid,octree_maximum_depth,3,true,false);
    std::cout<<"octree grid number of cells:"<<octree_grid.number_of_cells<<std::endl;

    // compute particle bounding boxes
    ARRAYS<VECTOR<ARRAY<BOX_3D<T> >,3> > particles_bounding_box(particles_array.m_start,particles_array.m_end,particles_array.n_start,particles_array.n_end,particles_array.mn_start,particles_array.mn_end);
    int number_of_particles=0;
    Get_Particle_Bounding_Boxes(particles_array,particles_bounding_box,scale,particle_blender,number_of_particles);
    ARRAYS<VECTOR<ARRAY<ARRAY<OCTREE_CELL<T>*> >,3> > particles_intersecting_cells(particles_array.m_start,particles_array.m_end,particles_array.n_start,particles_array.n_end,particles_array.mn_start,particles_array.mn_end);

    std::cout<<"Refining octree to particles"<<std::endl;
    printf("Number of cells=%d\n",octree_grid.number_of_cells);
    PROGRESS_INDICATOR refinement_progress(number_of_particles);
    int refinement_timer=TIMER::Singleton()->Register_Timer();
    std::cout <<"Number of particles:" << number_of_particles<<std::endl;
    for(int i=particles_array.m_start;i<=particles_array.m_end;i++) for(int j=particles_array.n_start;j<=particles_array.n_end;j++) for(int ij=particles_array.mn_start;ij<=particles_array.mn_end;ij++) {
        PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR_3D<T> >* particles = particles_array(i,j,ij);
        if(particles){
            for(int p=1;p<=particles->number;p++){
                BOX_3D<T> particle_bounding_box=particles_bounding_box(i,j,ij)(p);
                VECTOR_3D<T> particle_box_size(particle_bounding_box.Size());
                int cells_x=(int)ceil(particle_box_size.x/2/octree_grid.uniform_grid.dx);
                int cells_y=(int)ceil(particle_box_size.y/2/octree_grid.uniform_grid.dy);
                int cells_z=(int)ceil(particle_box_size.z/2/octree_grid.uniform_grid.dz);
                int m_start=max(octree_grid.cells.m_start,i-cells_x);int m_end=min(octree_grid.cells.m_end,i+cells_x);
                int n_start=max(octree_grid.cells.n_start,j-cells_y);int n_end=min(octree_grid.cells.n_end,j+cells_y);
                int mn_start=max(octree_grid.cells.mn_start,ij-cells_z);int mn_end=min(octree_grid.cells.mn_end,ij+cells_z);
                particles_intersecting_cells(i,j,ij).Append(ARRAY<OCTREE_CELL<T>* >());
                for(int m=m_start;m<=m_end;m++) for(int n=n_start;n<=n_end;n++) for(int mn=mn_start;mn<=mn_end;mn++) 
                    Refine_Octree_Cell_To_Box(octree_grid,octree_grid.cells(m,n,mn),particle_bounding_box,particles_intersecting_cells(i,j,ij)(p));
                refinement_progress.Progress();}}
    }
    std::cout<<"time for refinement:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(refinement_timer)<<std::endl;
    printf("Number of cells=%d\n",octree_grid.number_of_cells);

    std::cout<<"Initializing particle phi on octree"<<std::endl;
    ARRAY<T> particle_octree_phi(octree_grid.number_of_nodes,true);
    PROGRESS_INDICATOR particle_phi_progress(number_of_particles);
    int particle_field_timer=TIMER::Singleton()->Register_Timer();
    ARRAY<VECTOR_3D<T> >& node_locations=octree_grid.Node_Locations();
    OPERATION_HASH node_operations(octree_grid.number_of_nodes);
    for(int i=particles_array.m_start;i<=particles_array.m_end;i++) for(int j=particles_array.n_start;j<=particles_array.n_end;j++) for(int ij=particles_array.mn_start;ij<=particles_array.mn_end;ij++) {
        PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR_3D<T> >* particles = particles_array(i,j,ij);
        if(particles){
            for(int p=1;p<=particles->number;p++){
                BOX_3D<T> particle_bounding_box=particles_bounding_box(i,j,ij)(p);
                T radius=scale*particles->radius(p);
                T velocity_magnitude_squared=particles->V(p).Magnitude_Squared();
                ARRAY<OCTREE_CELL<T>*> intersecting_cells=particles_intersecting_cells(i,j,ij)(p);
                MATRIX_4X4<T> rotation; if(velocity_magnitude_squared!=0) rotation=MATRIX_4X4<T>::Rotation_Matrix(particles->V(p),VECTOR_3D<T>(1,0,0));
                for(int c=1;c<=intersecting_cells.m;c++) for(int k=0;k<8;k++){
                    int node_index=intersecting_cells(c)->Node(k); if(node_operations.Is_Marked_Current(node_index)) continue;
                    T distance=(velocity_magnitude_squared==0)?particle_blender.Get_Distance(3*radius,particles->X(p),node_locations(node_index))
                        :particle_blender.Get_Distance(3*radius,radius,particles->X(p),particles->V(p),node_locations(node_index),rotation);
                    particle_octree_phi(node_index)-=octree_grid.minimum_cell_size*particle_blender.C(distance)*particle_power;
                    node_operations.Mark(node_index);}
                particle_phi_progress.Progress();node_operations.Next_Operation();
            }
        }
    }
    std::cout<<"time for particle field:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(particle_field_timer)<<std::endl;

    std::cout<<"Merging levelset and particle levelset"<<std::endl;
    ARRAY<T> octree_phi(octree_grid.number_of_nodes,false);
    for(int i=1;i<=octree_grid.number_of_nodes;i++) octree_phi(i)=implicit_surface.levelset.Phi(node_locations(i));//+particle_octree_phi(i);

    std::cout<<"merged_octree_levelset.Fast_Marching_Method()..."<<std::flush;
    int fast_marching_method_timer=TIMER::Singleton()->Register_Timer();
    OCTREE_LEVELSET<T> merged_octree_levelset(octree_grid,octree_phi);
    merged_octree_levelset.Fast_Marching_Method();
    std::cout<<"done"<<std::endl;
    std::cout<<"time for fast_marching_method:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(fast_marching_method_timer)<<std::endl;

    std::cout<<"Writing octree levelset file"<<std::endl;
    std::ofstream octree_levelset_file(STRING_UTILITIES::string_sprintf("%s/octree_levelset.%d",output_directory.c_str(),frame).c_str());
    merged_octree_levelset.phi.Write<T>(octree_levelset_file);
    octree_levelset_file.close();

    std::cout<<"Writing octree grid file"<<std::endl;
    std::ofstream octree_grid_file(STRING_UTILITIES::string_sprintf("%s/octree_grid.%d",output_directory.c_str(),frame).c_str());
    octree_grid.Write<T>(octree_grid_file);
    octree_grid_file.close();

    delete rendering_implicit_surface;delete phi;delete grid;
    return 0;
}

