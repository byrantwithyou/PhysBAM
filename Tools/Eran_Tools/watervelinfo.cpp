#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

int main(int argc, char *argv[])
{
    char filename[256];

    PARSE_ARGS parse_args(argc,argv);

    std::string filename;
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    cout << "Filename: " << filename << endl;

    ifstream vel_file(filename, ios::binary);
    if (!vel_file) return -1;

    ARRAYS<VECTOR<double,3> > u, v, w;
    u.Read<double>(vel_file);
    v.Read<double>(vel_file);
    w.Read<double>(vel_file);

    double max_magnitude = 0;
    for (int i = 1; i <= u.m-1; i++)
        for (int j = 1; j <= v.n-1; j++)
            for (int k = 1; k <= w.mn-1; k++)
            {
                VECTOR<double,3> avg_vel = 0.5*VECTOR<double,3>(u(i,j,k)+u(i+1,j,k),
                                                  v(i,j,k)+v(i,j+1,k),
                                                  w(i,j,k)+w(i,j,k+1));
                cout << i << " " << j << " " << k << ": " << avg_vel << endl;
                max_magnitude = PhysBAM::max(max_magnitude, avg_vel.Magnitude());
            }

    cout << endl;
    cout << "Maximum magnitude: " << max_magnitude << endl;
}
