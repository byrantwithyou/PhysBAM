#include <Tools/Matrices/MATRIX_4X4.h>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    POLYGONAL_SURFACE surface;

    surface.Read(std::cin);

    MATRIX<double,4> xform = MATRIX<double,4>::Rotation_Matrix_Z_Axis(0.1);
//    MATRIX<double,4> xform = MATRIX<double,4>::Scale_Matrix(VECTOR<double,3>(1,2,3));
    for(int i=0;i<surface.vertices.m;i++) surface.vertices(i)=xform*surface.vertices(i);

    surface.Write(std::cout);
}
