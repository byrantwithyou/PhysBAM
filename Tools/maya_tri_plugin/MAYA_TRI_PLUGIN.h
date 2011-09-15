//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _MAYA_TRI_PLUGIN_
class MDagPath;
class MFnDagNode;
class ofstream;

#include <maya/MDagPath.h>
#include <maya/MPxFileTranslator.h>
#include <maya/MStatus.h>

template<class T>
class MAYA_TRI_PLUGIN:public MPxFileTranslator {
public:
    static void* creator();

    virtual MStatus writer(const MFileObject& file,const MString& optionsString,MPxFileTranslator::FileAccessMode mode);
    virtual MStatus reader(const MFileObject& file,const MString& optionsString,MPxFileTranslator::FileAccessMode mode);
    virtual bool haveWriteMethod() const;
    virtual bool haveReadMethod() const;
    virtual bool canBeOpened() const;
    MString defaultExtension() const;

//#####################################################################
};
#endif 
