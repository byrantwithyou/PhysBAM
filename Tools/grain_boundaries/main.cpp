#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TETRAHEDRA.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>

using namespace PhysBAM;

int main(int argc,const char *argv[])
{
    typedef float T;
    typedef float RW;

    RW rw=RW();
    STREAM_TYPE stream_type(rw);
    typedef VECTOR<T,3> TV;
    
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-tri","","tri input filename");
    parse_args.Add_String_Argument("-tet","","tet input filename");
    parse_args.Add_String_Argument("-seed_centers","","txt file of seed region centers");
    parse_args.Add_Integer_Argument("-refine_levels",0,"number of refinement levels");
    parse_args.Add_Integer_Argument("-seeds",10,"number grain boundary regions");
    parse_args.Add_Integer_Argument("-max_coarse_bcc_res",10,"resolution of largest bounding box dimension");
    parse_args.Add_Integer_Argument("-num_smoothing_steps",0,"number of smooting steps");

    parse_args.Parse(argc,argv);

    int max_refinement_levels=0;
    int max_dim_resolution=20;
    int num_seeds=500;
    int num_interface_surface_smooting_steps=1;
    std::string input_geom_filename,seed_center_filename;
    bool using_existing_seed_centers=false;
    bool tri_input_geom=false;
    bool tet_input_geom=false;

    if(parse_args.Is_Value_Set("-refine_levels")) max_refinement_levels=parse_args.Get_Integer_Value("-refine_levels");
    if(parse_args.Is_Value_Set("-max_coarse_bcc_res")) max_dim_resolution=parse_args.Get_Integer_Value("-max_coarse_bcc_res");
    if(parse_args.Is_Value_Set("-seeds")) num_seeds=parse_args.Get_Integer_Value("-seeds");
    if(parse_args.Is_Value_Set("-num_smoothing_steps")) num_interface_surface_smooting_steps=parse_args.Get_Integer_Value("-num_smoothing_steps");
    if(parse_args.Is_Value_Set("-tet")){input_geom_filename=parse_args.Get_String_Value("-tet");tet_input_geom=true;tri_input_geom=false;}
    if(parse_args.Is_Value_Set("-tri")){input_geom_filename=parse_args.Get_String_Value("-tri");tri_input_geom=true;tet_input_geom=false;}
    if(parse_args.Is_Value_Set("-seed_centers")){
        seed_center_filename=parse_args.Get_String_Value("-seed_centers");
        using_existing_seed_centers=true;
        max_refinement_levels=0;}

    std::cout<<"Input geom file name = "<<input_geom_filename<<std::endl;
    if(using_existing_seed_centers)
        std::cout<<"Input grain centers file = "<<seed_center_filename<<std::endl;
    std::cout<<"Number of refinement levels = "<<max_refinement_levels<<std::endl;
    std::cout<<"Max dim resolution in unrefined mesh= "<<max_dim_resolution<<std::endl;
    std::cout<<"Number of grain boudnary regions = "<<num_seeds<<std::endl;
    std::cout<<"Number of smoothing steps = "<<num_interface_surface_smooting_steps<<std::endl;

    if(!tri_input_geom && !tet_input_geom) PHYSBAM_FATAL_ERROR();

    DEFORMABLE_PARTICLES<TV> fracturing_geom_particles;
    TETRAHEDRON_MESH fracturing_geom_mesh;
    TETRAHEDRALIZED_VOLUME<T> fracturing_geom(fracturing_geom_mesh,fracturing_geom_particles);

    DEFORMABLE_PARTICLES<TV> model_geom_particles;
    TETRAHEDRON_MESH model_geom_mesh;
    TETRAHEDRALIZED_VOLUME<T> model_geom(model_geom_mesh,model_geom_particles);

    TRIANGLE_MESH tri_model_geom_mesh;
    TRIANGULATED_SURFACE<T> tri_model_geom(tri_model_geom_mesh,model_geom_particles);

    if(tet_input_geom)
        FILE_UTILITIES::Read_From_File(stream_type,input_geom_filename,model_geom);
    else if(tri_input_geom)
        FILE_UTILITIES::Read_From_File(stream_type,input_geom_filename,tri_model_geom);


    ARRAY<TV> fracture_points;
    if(using_existing_seed_centers){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(seed_center_filename,false);
        int number_of_points;
        (*input)>>number_of_points;
        fracture_points.Resize(number_of_points);
        LOG::cout<<"Number fracture points = "<<number_of_points<<std::endl;
        for (int k=1;k<=number_of_points;k++) (*input) >> fracture_points(k);
        delete input;}

    if(tet_input_geom) model_geom.Update_Bounding_Box();
    else if(tri_input_geom) tri_model_geom.Update_Bounding_Box();

    BOX<TV> model_geom_bounding_box;

    if(tet_input_geom)
        model_geom_bounding_box=*(model_geom.bounding_box);
    else if(tri_input_geom)
        model_geom_bounding_box=*(tri_model_geom.bounding_box);

    TV model_geom_dimensions=model_geom_bounding_box.Edge_Lengths();

    int m,n,mn;
    T max_dim=model_geom_dimensions.Max();
    m=(int)(max_dim_resolution*(model_geom_dimensions.x/max_dim));
    n=(int)(max_dim_resolution*(model_geom_dimensions.y/max_dim));
    mn=(int)(max_dim_resolution*(model_geom_dimensions.z/max_dim));
    std::cout<<"Fracture mesh dimensions "<<m<<" , "<<n<<" , "<<mn<<" ."<<std::endl;
    model_geom_bounding_box.Scale_About_Center((T)1.25);
    GRID<TV> grid(m,n,mn,model_geom_bounding_box);
    fracturing_geom.Initialize_Octahedron_Mesh_And_Particles(grid);

    //FILE_UTILITIES::Write_To_File(stream_type,"unrefined_fracture_geometry.tet",fracturing_geom);
    
    RED_GREEN_TETRAHEDRA<T> redgreen(fracturing_geom);
    
    T fracture_point_thickness=(T).01;
    for(int level=0;level<max_refinement_levels;level++){
        fracturing_geom.Update_Tetrahedron_List();
        fracturing_geom.Initialize_Hierarchy();
        ARRAY<bool> refine_tet(fracturing_geom.mesh.elements.m);
        for(int p=0;p<fracture_points.m;p++){
            TV& X=fracture_points(p);
            ARRAY<int> intersection_list;
            (*fracturing_geom.hierarchy).Intersection_List(X,intersection_list,fracture_point_thickness);
            refine_tet.Subset(intersection_list).Fill(true);}
        ARRAY<int> tets_to_refine;
        for(int t=0;t<fracturing_geom.mesh.elements.m;t++) if(refine_tet(t)) tets_to_refine.Append(t);
        redgreen.Refine_Simplex_List(tets_to_refine);}

    //FILE_UTILITIES::Write_To_File(stream_type,"refined_fracture_geometry.tet",fracturing_geom);

    
    assert(num_seeds>0);
    RANDOM_NUMBERS<T> random_number;
    random_number.Set_Seed(0); // ensure repeatability for testing purposes
    ARRAY<int> tet_color(fracturing_geom_mesh.elements.m);
    tet_color.Fill(0);
    ARRAY<TV> grain_centers(num_seeds);
    QUEUE<int> frontier_tets(num_seeds);

    //read seed centers from file
    if(using_existing_seed_centers){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(seed_center_filename,false);
        (*input)>>num_seeds;
        grain_centers.Resize(num_seeds);
        LOG::cout<<"Number seeds from file = "<<num_seeds<<std::endl;
        for (int k=1;k<=num_seeds;k++) (*input) >> grain_centers(k);
        std::cout<<"initializing seed tets..."<<std::endl;
        delete input;

        fracturing_geom.Update_Tetrahedron_List();
        fracturing_geom.Initialize_Hierarchy();
        ARRAY<bool> seen_this_tet_before(fracturing_geom.mesh.elements.m);
        for(int s=0;s<num_seeds;s++){
            int seed_tet;
            ARRAY<int> intersection_list;
            (*fracturing_geom.hierarchy).Intersection_List(grain_centers(s),intersection_list,(T)1e-1);
            for(int i=intersection_list.m;i>=1;i--)
                if(!(*fracturing_geom.tetrahedron_list)(intersection_list(i)).Inside(grain_centers(s)))
                    intersection_list.Remove_Index_Lazy(i);
            if(intersection_list.m<1){
                std::cout<<"Couldn't find a tet for seed "<<s<<" of "<<num_seeds<<" total input seeds"<<std::endl;
                break;}
            seed_tet=intersection_list(1);
            if(!seen_this_tet_before(seed_tet)){
                frontier_tets.Safe_Enqueue(seed_tet);
                tet_color(seed_tet)=s;}
            else
                std::cout<<"Seed "<<s<<" was in the same tet as some other seed"<<std::endl;
            seen_this_tet_before(seed_tet)=true;}}
    else{
        std::cout<<"initializing seed tets..."<<std::endl;
        for(int s=0;s<num_seeds;s++){
            int random_tet;
            do {
                random_tet=random_number.Get_Uniform_Integer(1,fracturing_geom_mesh.elements.m);
            } while(tet_color(random_tet)>0);
            tet_color(random_tet)=s;
            int i,j,k,l;fracturing_geom_mesh.elements(random_tet).Get(i,j,k,l);
            grain_centers(s)=(T).25*(fracturing_geom_particles.X(i)+fracturing_geom_particles.X(j)+fracturing_geom_particles.X(k)+fracturing_geom_particles.X(l));
            frontier_tets.Safe_Enqueue(random_tet);}

        std::cout<<"Writing grain centers to file"<<std::endl;

        std::string grain_output_file="grain_centers.txt";std::ostream* output=FILE_UTILITIES::Safe_Open_Output(grain_output_file,false);
        (*output)<<num_seeds<<std::endl;
        for (int k=1;k<=num_seeds;k++) (*output)<<grain_centers(k)<<std::endl;
        delete output;}

    

    std::cout<<"initializing mesh structures..."<<std::endl;
    fracturing_geom_mesh.Initialize_Incident_Elements();
    fracturing_geom_mesh.Initialize_Adjacent_Elements();
    fracturing_geom_mesh.Initialize_Triangle_Mesh();
    fracturing_geom_mesh.Initialize_Tetrahedron_Faces();

    const ARRAY< ARRAY<int> >& adjacent_elements=*(fracturing_geom_mesh.adjacent_elements);
    const TRIANGLE_MESH& triangle_mesh=*(fracturing_geom_mesh.triangle_mesh);
    const ARRAY<VECTOR<int,4> >& tetrahedron_faces=*(fracturing_geom_mesh.tetrahedron_faces);

    QUEUE<int> interface_tets(0);
    ARRAY<bool> interface_tets_contains(fracturing_geom_mesh.elements.m);
    interface_tets_contains.Fill(false);

    std::cout<<"growing seed regions via BFS..."<<std::endl;
    // at the same time, queue up the tets along the interfaces between regions
    while(!frontier_tets.Empty()){
        int tet=frontier_tets.Dequeue();
        int color=tet_color(tet);
        const ARRAY<int>& adj_tets=adjacent_elements(tet);
        for(int t=0;t<adj_tets.m;t++){
            int adj_tet=adj_tets(t);
            int& adj_color=tet_color(adj_tet);
            if(adj_color==0){
                adj_color=color;
                frontier_tets.Safe_Enqueue(adj_tet);
            }
            else if(adj_color!=color){
                if(!interface_tets_contains(tet)){
                    interface_tets.Safe_Enqueue(tet);
                    interface_tets_contains(tet)=true;
                }
                if(!interface_tets_contains(adj_tet)){
                    interface_tets.Safe_Enqueue(adj_tet);
                    interface_tets_contains(adj_tet)=true;
                }
            }
        }
    }

    // debug check that all tets are colored
    for(int t=0;t<fracturing_geom_mesh.elements.m;t++)
        assert(tet_color(t)>0);



    std::cout<<"attempting to smooth out interface..."<<std::endl;
    while(!interface_tets.Empty()){
        int tet=interface_tets.Dequeue();
        assert(interface_tets_contains(tet));
        interface_tets_contains(tet)=false;
        int& color=tet_color(tet);
        HASHTABLE<int,T> color_area_hash;
        const ARRAY<int>& adj_tets=adjacent_elements(tet);
        // associating each face of the tet with the color of its adjacent tet,
        //  total the face area of each color
        for(int t=0;t<adj_tets.m;t++){
            int adj_tet=adj_tets(t);
            int adj_color=tet_color(adj_tet);
            T face_area=TRIANGLE_3D<T>(fracturing_geom_particles.X.Subset(triangle_mesh.elements(tetrahedron_faces(tet)(t)))).Area();
            assert(face_area>T(0.0));
            if(color_area_hash.Contains(adj_color))
                color_area_hash.Get(adj_color) += face_area;
            else
                color_area_hash.Insert(adj_color,face_area);
        }
        // determine the color that accounts for the most face area
        T max_area=T(0.0);
        int max_area_color=0;
        for(HASHTABLE<int,T>::iterator i_color_area_hash=color_area_hash.begin();i_color_area_hash!=color_area_hash.end();++i_color_area_hash){
            if(i_color_area_hash->state!=ENTRY_ACTIVE)
                continue;
            int area_color=i_color_area_hash->key;
            T area=i_color_area_hash->data;
            if(area>max_area){
                max_area=area;
                max_area_color=area_color;
            }
        }
        assert(max_area_color>0.0);
        if(max_area_color!=color){
            // color the tet the color which accounts for the most face area,
            //  thereby attempting to minimize the area of the interfaces
            //  between regions
            // also enqueue tets newly exposed to the interface
            color=max_area_color;
            for(int t=0;t<adj_tets.m;t++){
                int adj_tet=adj_tets(t);
                int adj_color=tet_color(adj_tet);
                if(adj_color!=color && !interface_tets_contains(adj_tet)){
                    interface_tets.Safe_Enqueue(adj_tet);
                    interface_tets_contains(adj_tet)=true;
                }
            }
        }
    }



    std::cout<<"building interface mesh..."<<std::endl;
    ARRAY<VECTOR<int,3> > interface_tris;
    for(int tet=0;tet<fracturing_geom_mesh.elements.m;tet++){
        const VECTOR<int,4>& tet_vertices=fracturing_geom_mesh.elements(tet);
        int color=tet_color(tet);
        const ARRAY<int>& adj_tets=adjacent_elements(tet);
        for(int t=0;t<adj_tets.m;t++){
            int adj_tet=adj_tets(t);
            int adj_color=tet_color(adj_tet);
            if(adj_color<=color) // "<=" to prevent double-counting tris
                continue;
            // since tet and adj_tet have differing colors, add the tri common
            //  to both, oriented such that its normal points out from tet
            const VECTOR<int,4>& adj_tet_vertices=fracturing_geom_mesh.elements(adj_tet);
            // find the vertex i of tet which is not in the common tri
            int i=1;
            for(;i<=4;++i){
                int vertex_i=tet_vertices(i);
                int j=1;
                for(;j<=4 && vertex_i!=adj_tet_vertices(j);++j);
                if(j>4)
                    break;
            }
            assert(i<=4);
            // pull the other 3 vertices (not including i) from tet
            VECTOR<int,3> tri_vertices;
            for(int j=0;j<4;j++){
                if(j<i)
                    tri_vertices(j)=tet_vertices(j);
                else if(j>i)
                    tri_vertices(j-1)=tet_vertices(j);
            }
            // re-orient if removed vertex is even (that's the rule; something
            //  to do with even and odd permutations, probably)
            if(i%2==0){
                int temp=tri_vertices(1);
                tri_vertices(1)=tri_vertices(2);
                tri_vertices(2)=temp;
            }
            interface_tris.Append(tri_vertices);
        }
    }



    TRIANGLE_MESH interface_mesh(fracturing_geom_particles.number,interface_tris);
    TRIANGULATED_SURFACE<T> interface_surface(interface_mesh,fracturing_geom_particles);

    interface_mesh.Initialize_Neighbor_Nodes();

    //smooth interface surface
    ARRAY<TV> smoothed_vertices(fracturing_geom_particles.number),smoothed_vertices_save(fracturing_geom_particles.number);

    //don't smootht the boundary nodes
    interface_mesh.Initialize_Node_On_Boundary();

    for(int p=0;p<fracturing_geom_particles.number;p++) smoothed_vertices(p)=fracturing_geom_particles.X(p);
    
    for(int s=0;s<num_interface_surface_smooting_steps;s++){
        for(int p=0;p<fracturing_geom_particles.number;p++) smoothed_vertices_save(p)=smoothed_vertices(p);
        for(int p=0;p<fracturing_geom_particles.number;p++){
            ARRAY<int>& neighbor_nodes=(*interface_mesh.neighbor_nodes)(p);
            TV average_neighbor_location=smoothed_vertices_save(p);
            if(!(*interface_mesh.node_on_boundary)(p)){
                for(int n=0;n<neighbor_nodes.m;n++) average_neighbor_location+=smoothed_vertices_save(neighbor_nodes(n));
                average_neighbor_location/=(T)(neighbor_nodes.m+1);}
            smoothed_vertices(p)=average_neighbor_location;}}

    for(int p=0;p<fracturing_geom_particles.number;p++) fracturing_geom_particles.X(p)=smoothed_vertices(p);
    

    FILE_UTILITIES::Write_To_File(stream_type,"interface_surface.tri",interface_surface);

    return 0;
}
