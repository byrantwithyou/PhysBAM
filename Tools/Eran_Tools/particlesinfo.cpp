#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include <Particles/SOLIDS_PARTICLES.h>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T,class T_PARTICLES>
void print_info(const char *filename, int verbose_level)
{
    ifstream input_file(filename, ios::binary);
    T_PARTICLES particles;
    particles.template Read<T>(input_file);

    cout << "number = " << particles.number << std::endl;
    cout << "NUMBER = " << particles.number << endl;
    cout << "number of attributes = " << particles.number_of_attributes << endl;
    if (verbose_level >= 1)
    {
        for (int i = 1; i <= particles.array_size; i++)
        {
            cout << i << ": ";
            cout << "X=" << particles.X(i) << " ";
            if (particles.store_velocity) cout << "V=" << particles.V(i) << " ";
            if (particles.store_acceleration) cout << "A=" << particles.A(i) << " ";
            if (verbose_level >= 2)
            {
                if (particles.store_mass) cout << "mass=" << particles.mass(i) << " ";
                if (particles.store_radius) cout << "radius=" << particles.radius(i) << " ";
                if (particles.store_temperature) cout << "temperature=" << particles.temperature(i) << " ";
                if (particles.store_density) cout << "density=" << particles.density(i) << " ";
                if (particles.store_age) cout << "age=" << particles.age(i) << " ";
                if (particles.store_id) cout << "id=" << particles.id(i) << " ";
                if (particles.store_material_coordinates) cout << "material_coorinates=" << particles.Xm(i) << " ";
#ifdef USE_NAMED_ATTRIBUTES
                if (particles.store_named_attributes && verbose_level >= 3) {
                    cout << "named attributes:" << endl;
                    particles.named_attributes.Print(std::cout,i);
                }
#endif
            }
            cout << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    int dim = 3;
    bool type_double = false,opt_1d=false,opt_2d=false,opt_3d=false,opt_v=false,opt_vv=false,opt_vvv=false;
    char filename[256];
    int verbose_level = 0;
    std::string filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-1d",&opt_1d,"1d");
    parse_args.Add("-2d",&opt_2d,"2d");
    parse_args.Add("-3d",&opt_3d,"3d");
    parse_args.Add("-v",&opt_v,"v");
    parse_args.Add("-vv",&opt_vv,"vv");
    parse_args.Add("-vvv",&opt_vvv,"vvv");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    if (opt_1d) dim = 1;
    if (opt_2d) dim = 2;
    if (opt_3d) dim = 3;
    if (opt_v) verbose_level = 1;
    if (opt_vv) verbose_level = 2;
    if (opt_vvv) verbose_level = 3;

    if (!FILE_UTILITIES::File_Exists(filename))
    {
        cerr << "File " << filename << " does not exist" << std::endl;
        return 1;
    }

    cout << "Filename: " << filename
         << " [" << dim << "D, " << ((type_double)?"double":"float") << "]" << endl;

    if (type_double)
    {
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        switch (dim)
        {
            case 1: print_info<double,SOLIDS_PARTICLES<double,VECTOR<double,1> > >(filename, verbose_level); break;
            case 2: print_info<double,SOLIDS_PARTICLES<double,VECTOR<double,2> > >(filename, verbose_level); break;
            default: print_info<double,SOLIDS_PARTICLES<double,VECTOR<double,3> > >(filename, verbose_level); break;
        }
#else
        std::cerr<<"Double support not enabled."<<std::endl;exit(1);
#endif
    }
    else
    {
        switch (dim)
        {
            case 1: print_info<float,SOLIDS_PARTICLES<float,VECTOR<float,1> > >(filename, verbose_level); break;
            case 2: print_info<float,SOLIDS_PARTICLES<float,VECTOR<float,2> > >(filename, verbose_level); break;
            default: print_info<float,SOLIDS_PARTICLES<float,VECTOR<float,3> > >(filename, verbose_level); break;
        }
    }
}
