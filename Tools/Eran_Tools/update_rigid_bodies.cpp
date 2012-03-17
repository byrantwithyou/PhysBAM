#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/TETRAHEDRALIZED_VOLUME_LIST.h"

using namespace PhysBAM;

int main(int argc, char *argv[])
{
    std::string filename;
    if(argc > 1) filename=argv[1];
    else filename="rigid_body_key";

    ARRAY<int> rigid_body_id_to_triangulated_surface_id;
    ARRAY<int> rigid_body_id_to_implicit_surface_id;
    ARRAY<int> rigid_body_id_to_tetrahedralized_volume_id;

    FILE_UTILITIES::Read_From_File<float>(filename,rigid_body_id_to_triangulated_surface_id,rigid_body_id_to_implicit_surface_id);
    if(rigid_body_id_to_triangulated_surface_id.m != rigid_body_id_to_implicit_surface_id.m){
        std::cout << "Error: rigid_body_id_to_triangulated_surface_id.m != rigid_body_id_to_implicit_surface_id.m" << std::endl;
        exit(1);
    }

    char version=1;
    rigid_body_id_to_tetrahedralized_volume_id.Resize(rigid_body_id_to_triangulated_surface_id.m);
    std::cout << "Writing " << rigid_body_id_to_triangulated_surface_id.m << " entries" << std::endl;
    FILE_UTILITIES::Write_To_File<float>(filename,version,rigid_body_id_to_triangulated_surface_id,rigid_body_id_to_implicit_surface_id,rigid_body_id_to_tetrahedralized_volume_id);

    TETRAHEDRALIZED_VOLUME_LIST<float> tetrahedralized_volume_list;
    tetrahedralized_volume_list.Write<float>("rigid_body_tetrahedralized_volume_",0);
}
