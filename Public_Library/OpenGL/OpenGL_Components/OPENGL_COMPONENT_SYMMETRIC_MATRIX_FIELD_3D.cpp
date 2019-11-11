//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Reinitialize
//##################################################################### 
template<class T> void OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>::
Reinitialize()
{
    if(!draw) return;
    valid=false;
    std::string tmp_filename=viewer_dir.current_directory+"/"+field_filename;
    if(File_Exists(tmp_filename))
        Read_From_File(tmp_filename,field);
    else return;
    opengl_symmetric_matrix_field.Update();
    valid=true;
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<float>;
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<double>;
}
