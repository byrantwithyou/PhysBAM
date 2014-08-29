#!/bin/bash

function all_cpp
{
    pushd $1 >& /dev/null
    find -name '*.cpp' | sed 's@^\./@@'
    popd >& /dev/null
}

function all_h
{
    pushd $1 >& /dev/null
    find -name '*.h' | sed 's@^\./@@'
    popd >& /dev/null
}

echo 'D:\ucla\PhysBAM\Public_Library;%(AdditionalIncludeDirectories);C:\local\boost_1_56_0;D:\ucla\Summer 2014\downloads\zlib128-dll\include' > incs.txt

echo 'D:/ucla/Summer 2014/downloads/zlib128-dll/lib;%(AdditionalLibraryDirectories)' > lib-dirs.txt

echo '..\..\..\..\Summer 2014\downloads\zlib128-dll\lib\zdll.lib' > libs.txt

rm -f projs.txt

for p in Compressible Deformables Dynamics Fluids Geometry Incompressible OpenGL Ray_Tracing Rigids Solids Tools ; do
    mkdir -p ../../Public_Library/$p/msvc13
    all_cpp ../../Public_Library/$p > all-cpp.txt
    all_h ../../Public_Library/$p > all-h.txt
    args="Public_Library/$p $p 1 -- all-cpp.txt all-h.txt incs.txt lib-dirs.txt libs.txt --"
    ./gen_msvc13_project.pl $args < template.vcxproj > ../../Public_Library/$p/msvc13/$p.vcxproj
    ./gen_msvc13_project.pl $args < template.vcxproj.filters > ../../Public_Library/$p/msvc13/$p.vcxproj.filters
    echo $p Public_Library/$p >> projs.txt
done

for p in be_evolution ; do
    mkdir -p ../../Projects/$p/msvc13
    all_cpp ../../Projects/$p > all-cpp.txt
    all_h ../../Projects/$p > all-h.txt
    args="Projects/$p $p 0 -- all-cpp.txt all-h.txt incs.txt lib-dirs.txt libs.txt --"
    ./gen_msvc13_project.pl $args < template.vcxproj > ../../Projects/$p/msvc13/$p.vcxproj
    ./gen_msvc13_project.pl $args < template.vcxproj.filters > ../../Projects/$p/msvc13/$p.vcxproj.filters
    echo $p Projects/$p >> projs.txt
done

mkdir -p ../../msvc13
args=". PhysBAM 0 -- projs.txt --"
./gen_msvc13_project.pl $args < template.sln > ../../msvc13/$p.sln

