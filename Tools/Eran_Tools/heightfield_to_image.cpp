#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/EXR_FILE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T>
void Convert(PARSE_ARGS &parse_args)
{
    std::string input_filename=parse_args.Extra_Arg(1);
    std::string output_filename;

    if (parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    else output_filename=FILE_UTILITIES::Get_Basename(input_filename)+".exr";

    ARRAYS<VECTOR<T,2> > heightfield;

    FILE_UTILITIES::Read_From_File<T>(input_filename,heightfield);

    if (FILE_UTILITIES::Get_File_Extension(output_filename)=="exr"){
        ARRAYS<VECTOR<VECTOR<float,3> ,2> > image(1,heightfield.m,1,heightfield.n);
        for(int i=0;i<image.m;i++) for(int j=0;j<image.n;j++){
            image(i,j).x=(float)heightfield(heightfield.m_start+i-1,heightfield.n_start+j-1);
        }
        EXR_FILE::Write_Row_Column_Image(output_filename,image);
    }
    else{
        std::cerr << "Unsupported image file extension: " << FILE_UTILITIES::Get_File_Extension(output_filename) << std::endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_String_Argument("-o", "");
    parse_args.Set_Extra_Arguments(-1, "<filename>");

    parse_args.Parse(argc, argv);

    if (parse_args.Num_Extra_Args() < 1) return 1;

    if (parse_args.Get_Option_Value("-double")) Convert<double>(parse_args);
    else Convert<float>(parse_args);
}
