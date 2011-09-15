#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

std::string input_filename, output_filename;
double threshold_fraction;

template<class T,class RW>
void Seal_Holes(PARSE_ARGS& parse_args)
{
    bool merge_coincident_vertices=!parse_args.Get_Option_Value("-no_vertex_merging");
    bool fill_holes=!parse_args.Get_Option_Value("-no_hole_filling");

    TRIANGULATED_SURFACE<T>* surface=0;
    FILE_UTILITIES::Create_From_File<RW>(input_filename,surface);
    surface->Remove_Degenerate_Triangles(1e-6);
    surface->Close_Surface(merge_coincident_vertices,threshold_fraction,fill_holes,true);

    LOG::cout<<"Writing to "<<output_filename<<std::endl;
    FILE_UTILITIES::Write_To_File<RW>(output_filename,*surface);
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Double_Argument("-threshold",.1,"threshold","fraction of minimium segment length for threshold");
    parse_args.Add_Option_Argument("-no_vertex_merging");
    parse_args.Add_Option_Argument("-no_hole_filling");
    parse_args.Add_String_Argument("-o","","output filename");
    parse_args.Set_Extra_Arguments(1,"<filename>");

    parse_args.Parse(argc,argv);
    
    if(parse_args.Num_Extra_Args()<1) parse_args.Print_Usage(true);
    else input_filename=parse_args.Extra_Arg(1);

    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    else output_filename=FILE_UTILITIES::Get_Basename(input_filename)+"_sealed."+FILE_UTILITIES::Get_File_Extension(input_filename);

    threshold_fraction=parse_args.Get_Double_Value("-threshold");

    if(!parse_args.Get_Option_Value("-double")) Seal_Holes<float,float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Seal_Holes<double,float>(parse_args);
#else
    else{LOG::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
