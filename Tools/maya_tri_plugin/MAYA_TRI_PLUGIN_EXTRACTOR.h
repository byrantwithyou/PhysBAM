//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _MAYA_TRI_PLUGIN_EXTRACTOR_
class MDagPath;
class MFnDagNode;
class ofstream;

#include <maya/MDagPath.h>
#include <maya/MPxFileTranslator.h>
#include <maya/MStatus.h>

template<class T>
class MAYA_TRI_PLUGIN_EXTRACTOR
{
private:
    PhysBAM::SOLIDS_PARTICLES<T,PhysBAM::VECTOR_3D<T> > particles;
    PhysBAM::TRIANGLE_MESH triangle_mesh;
    PhysBAM::TRIANGULATED_SURFACE<T> triangulated_surface;
    PhysBAM::ARRAY<PhysBAM::VECTOR_3D<T> > colors;
    int base;

public:
    MAYA_TRI_PLUGIN_EXTRACTOR()
        :triangulated_surface(triangle_mesh,particles),base(1)
    {}

//#####################################################################
    MStatus Extract_All();
    MStatus Extract_Selection();
    MStatus Extract_Mesh(const MDagPath dag_path);
    MStatus Write(std::string filename_base);
//#####################################################################

};
#endif