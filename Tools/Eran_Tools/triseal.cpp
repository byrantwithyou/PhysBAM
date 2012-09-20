#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

template<class T,class RW>
void Seal_Holes(PARSE_ARGS& parse_args)
{
    std::string input_filename,output_filename;
    T threshold_fraction=.1;
    bool merge_coincident_vertices=true,fill_holes=true;

    parse_args.Add("-threshold",&threshold_fraction,"threshold","fraction of minimium segment length for threshold");
    parse_args.Add_Not("-no_vertex_merging",&merge_coincident_vertices,"merge coincident vertices");
    parse_args.Add_Not("-no_hole_filling",&fill_holes,"do not fill holes");
    parse_args.Add("-o",&output_filename,"file","output filename");
    parse_args.Set_Extra_Arguments(1,"<filename>");
    parse_args.Parse();

    if(parse_args.Num_Extra_Args()<1) parse_args.Print_Usage(true);
    else input_filename=parse_args.Extra_Arg(0);

    if(!output_filename.size()) output_filename=FILE_UTILITIES::Get_Basename(input_filename)+"_sealed."+FILE_UTILITIES::Get_File_Extension(input_filename);

    TRIANGULATED_SURFACE<T>* surface=0;
    FILE_UTILITIES::Create_From_File<RW>(input_filename,surface);
    surface->Remove_Degenerate_Triangles(1e-6);
    surface->Close_Surface(merge_coincident_vertices,threshold_fraction,fill_holes,true);

    LOG::cout<<"Writing to "<<output_filename<<std::endl;
    FILE_UTILITIES::Write_To_File<RW>(output_filename,*surface);
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    bool type_double=false;
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Parse(true);
    
    if(!type_double) Seal_Holes<float,float>(parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Seal_Holes<double,float>(parse_args);
#else
    else{LOG::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
