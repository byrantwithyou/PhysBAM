#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

int main(int argc, char *argv[])
{
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse(argc, argv);

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

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
