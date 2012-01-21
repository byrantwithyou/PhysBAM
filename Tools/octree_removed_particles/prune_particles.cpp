#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <fstream>
#include <iostream>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>

using namespace PhysBAM;

typedef float T;
typedef float RW;

ARRAY<int>* Get_Id_Array(PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles)
{
    PARTICLE_ATTRIBUTE<int>* id_attribute=particles->attribute_collection->Get_Id(particles->Attribute_Map());
    if(id_attribute) return &id_attribute->array;
    else {std::cerr << "Can't find ID array" << std::endl;exit(1);}
}
    
int main(int argc, char* argv[])
{
    PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >::Attribute_Map().Store(PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >::Attribute_Map().id);

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","");
    parse_args.Add_Vector_3D_Argument("-subdomain_minimum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_maximum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_2_minimum_corner",VECTOR<double,3>());
    parse_args.Add_Vector_3D_Argument("-subdomain_2_maximum_corner",VECTOR<double,3>());

    parse_args.Add_Vector_3D_Argument("-cylinder_bottom",VECTOR<double,3>(1,.7,2));
    parse_args.Add_Vector_3D_Argument("-cylinder_top",VECTOR<double,3>(1,4,2));
    parse_args.Add_Double_Argument("-radius",.35);
    parse_args.Add_Vector_3D_Argument("-cylinder_2_bottom",VECTOR<double,3>(1,.7,2));
    parse_args.Add_Vector_3D_Argument("-cylinder_2_top",VECTOR<double,3>(1,4,2));
    parse_args.Add_Double_Argument("-radius_2",.35);
    parse_args.Add_Option_Argument("-use_cylinder");
    parse_args.Add_Option_Argument("-use_cylinder_2");

    parse_args.Add_Option_Argument("-delete_outside");
    parse_args.Add_Option_Argument("-delete_inside");
    parse_args.Add_Option_Argument("-get_ids_outside");
    parse_args.Add_Option_Argument("-get_ids_inside");
    parse_args.Add_String_Argument("-ids_to_keep","");
    parse_args.Add_String_Argument("-ids_to_delete","");
    parse_args.Set_Extra_Arguments(1, "<filename>");
    parse_args.Parse(argc,argv);
    if (parse_args.Num_Extra_Args() != 1) 
    {
        parse_args.Print_Usage();
        return 1;
    }

    std::string filename=parse_args.Extra_Arg(1);
    std::string output_filename=filename+"_pruned";
    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");

    bool get_ids_outside=parse_args.Get_Option_Value("-get_ids_outside");
    bool get_ids_inside=parse_args.Get_Option_Value("-get_ids_inside");
    bool get_ids=get_ids_outside||get_ids_inside;

    if(!get_ids) std::cout<<"Reading particles from " << filename <<std::endl;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >*> particles_array;
    FILE_UTILITIES::Read_From_File<RW>(filename,particles_array);

    ARRAY<int> ids;ids.Preallocate(500);

    bool ids_to_keep=parse_args.Is_Value_Set("-ids_to_keep");
    bool ids_to_delete=parse_args.Is_Value_Set("-ids_to_delete");
    bool delete_outside=!parse_args.Get_Option_Value("-delete_inside");

    if(ids_to_keep||ids_to_delete){
        std::string ids_filename;
        if(ids_to_keep) ids_filename=parse_args.Get_String_Value("-ids_to_keep");
        else ids_filename=parse_args.Get_String_Value("-ids_to_delete");
        std::ifstream input_stream(ids_filename.c_str());
        int id,max_id=0;while(input_stream >> id){ids.Append(id);max_id=max(max_id,id);}
        ARRAY<bool> id_in_list(max_id);for(int i=0;i<ids.m;i++) id_in_list(ids(i))=true;

        int count=0;
        for(int i=0;i<particles_array.m;i++) if(particles_array(i)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles=particles_array(i);
            ARRAY<int>* id_array=Get_Id_Array(particles);
            for(int p=particles->number;p>=1;p--){
                int id=(*id_array)(p);
                bool in_list=(id<=id_in_list.m && id_in_list(id));
                if((ids_to_keep && !in_list) || (ids_to_delete && in_list)){count++;particles->Delete_Particle(p);}
            }
            if(!particles->number){delete particles;particles_array(i)=0;}}
        std::cout << "Deleted " << count << " particles" << std::endl;
    }
    else if(parse_args.Is_Value_Set("-subdomain_minimum_corner") && parse_args.Is_Value_Set("-subdomain_maximum_corner")){
        BOX_3D<T> subdomain(VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_minimum_corner")),VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_maximum_corner")));
        BOX_3D<T> subdomain2=subdomain;
        if(parse_args.Is_Value_Set("-subdomain_2_minimum_corner") && parse_args.Is_Value_Set("-subdomain_2_maximum_corner")){
            subdomain2=BOX_3D<T>(VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_2_minimum_corner")),VECTOR<T,3>(parse_args.Get_Vector_3D_Value("-subdomain_2_maximum_corner")));}
        if(!get_ids) std::cout << "Pruning particles outside " << subdomain << " and " << subdomain2 << std::endl;
        int count=0;
        for(int i=0;i<particles_array.m;i++) if(particles_array(i)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles=particles_array(i);
            ARRAY<int>* id_array=0;if(get_ids) id_array=Get_Id_Array(particles);
            if(get_ids_outside){
                for(int p=particles->number;p>=1;p--)
                    if(subdomain.Lazy_Outside(particles->X(p)) && subdomain2.Lazy_Outside(particles->X(p))) ids.Append((*id_array)(p));
            }
            else if(get_ids_inside){
                for(int p=particles->number;p>=1;p--)
                    if(!(subdomain.Lazy_Outside(particles->X(p)) && subdomain2.Lazy_Outside(particles->X(p)))) ids.Append((*id_array)(p));
            }
            else if(delete_outside){
                for(int p=particles->number;p>=1;p--){
                    if(subdomain.Lazy_Outside(particles->X(p)) && subdomain2.Lazy_Outside(particles->X(p))){
                        count++;particles->Delete_Particle(p);}}
            }
            else{
                for(int p=particles->number;p>=1;p--){
                    if(!(subdomain.Lazy_Outside(particles->X(p)) && subdomain2.Lazy_Outside(particles->X(p)))){
                        count++;particles->Delete_Particle(p);}}
            }
            if(!particles->number){delete particles;particles_array(i)=0;}}
        if(!get_ids) std::cout << "Deleted " << count << " particles" << std::endl;
    }
    else if(parse_args.Get_Option_Value("-use_cylinder")){
        VECTOR<T,3> bottom(parse_args.Get_Vector_3D_Value("-cylinder_bottom"));
        VECTOR<T,3> top(parse_args.Get_Vector_3D_Value("-cylinder_top"));
        T radius=parse_args.Get_Double_Value("-radius");
        CYLINDER<T> cylinder(bottom,top,radius);
        CYLINDER<T> cylinder2=cylinder;
        if(parse_args.Is_Value_Set("-use_cylinder_2")){
            VECTOR<T,3> bottom2(parse_args.Get_Vector_3D_Value("-cylinder_2_bottom"));
            VECTOR<T,3> top2(parse_args.Get_Vector_3D_Value("-cylinder_2_top"));
            T radius2=parse_args.Get_Double_Value("-radius_2");
            cylinder2=CYLINDER<T>(bottom2,top2,radius2);
            std::cerr << "and cylinder " << bottom2 << ", " << top2 << ", r="<< radius2 << std::endl;}
        std::cerr << "cylinder " << bottom << ", " << top << ", r="<< radius << std::endl;
        int count=0;
        for(int i=0;i<particles_array.m;i++) if(particles_array(i)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles=particles_array(i);
            ARRAY<int>* id_array=0;if(get_ids) id_array=Get_Id_Array(particles);
            if(get_ids_outside){
                for(int p=particles->number;p>=1;p--)
                    if(cylinder.Lazy_Outside(particles->X(p)) && cylinder2.Lazy_Outside(particles->X(p))) ids.Append((*id_array)(p));
            }
            else if(get_ids_inside){
                for(int p=particles->number;p>=1;p--)
                    if(!(cylinder.Lazy_Outside(particles->X(p)) && cylinder2.Lazy_Outside(particles->X(p)))) ids.Append((*id_array)(p));
            }
            else if(delete_outside){
                for(int p=particles->number;p>=1;p--){
                    if(cylinder.Lazy_Outside(particles->X(p)) && cylinder2.Lazy_Outside(particles->X(p))){
                        count++;particles->Delete_Particle(p);}}
            }
            else{
                for(int p=particles->number;p>=1;p--){
                    if(!(cylinder.Lazy_Outside(particles->X(p)) && cylinder2.Lazy_Outside(particles->X(p)))){
                        count++;particles->Delete_Particle(p);}}
            }
            if(!particles->number){delete particles;particles_array(i)=0;}}
        if(!get_ids) std::cout << "Deleted " << count << " particles" << std::endl;
    }
    else if(get_ids){
        for(int i=0;i<particles_array.m;i++) if(particles_array(i)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<T,VECTOR<T,3> >* particles=particles_array(i);
            ARRAY<int>* id_array=0;if(get_ids) id_array=Get_Id_Array(particles);
            for(int p=particles->number;p>=1;p--){ids.Append((*id_array)(p));}}
    }

    if(!get_ids){
        std::cout << "Writing particles to " << output_filename << std::endl;
        FILE_UTILITIES::Write_To_File<RW>(output_filename,particles_array);}
    else{
        for(int i=0;i<ids.m;i++) std::cout << ids(i) << std::endl;}
    
    return 0;
}

