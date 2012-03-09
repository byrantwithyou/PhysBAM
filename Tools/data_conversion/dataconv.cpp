/* dataconv
 * -----------------
 * Includes converters to convert double files to float files and vice versa.
 * Adding the particle parameter requires putting the DEFORMABLE_PARTICLES.h file included
 * in this directory into Public_Library/Particles/DEFORMABLE_PARTICLES.h and recompiling.
 * The code is organized by program functionality (dtof, ftod, etc.).
 *
 * The exectutable produced takes the following command line:
 * dataconv [dtof/ftod/pmodf/pmodd] [infile] [outfile]
 */

#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <fstream>
#include <iostream>
#include <string>
#include "../../Public_Library/Geometry/IMPLICIT_SURFACE.h"
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include "../../Public_Library/Geometry/OCTREE_IMPLICIT_SURFACE.h"
#include "../../Public_Library/Grids/OCTREE_MESH.h"
#include "../../Public_Library/Level_Sets/OCTREE_LEVELSET.h"
#include "../../Public_Library/Particles/PARTICLE_3D.h"

using namespace PhysBAM;
using namespace std;

int main(int argc, char* argv[])
{

     if (argc < 4) { 
        cout << "Arguments needed: [dtof/ftod/pmodf/pmodd] [infile] [outfile]" << endl;
        return 0;
    }
    string program(argv[1]);
    string infilename(argv[2]);
    string extension(infilename.substr(infilename.length() -3, infilename.length()-1));

    string outname(argv[3]);
    ofstream os(outname.c_str(),ios::binary);
    if(!os) {cerr<<"Could not open "<<outname<<endl;return -1;}

    //*************************
    // DTOF
    //*************************
    if (program == "dtof") {
        cout << "***** converting " << infilename.c_str() <<" from double to float format..." << endl;

        if (extension == string("tri")) {
            TRIANGULATED_SURFACE<double> *triSurface= 
                TRIANGULATED_SURFACE<double>::Create_From_File(infilename.c_str());        
            triSurface->Write_Float(os);
        
        } else if (extension == string("oct")) {
            OCTREE_IMPLICIT_SURFACE<double> *octSurface =
                OCTREE_IMPLICIT_SURFACE<double>::Create_From_File(infilename.c_str());
            octSurface->Write_Float(os);
    
        } else if (extension == string("phi")) {
            LEVELSET_IMPLICIT_SURFACE<double> *impSurf =
                LEVELSET_IMPLICIT_SURFACE<double>::Create_From_File(infilename.c_str());
            impSurf->Write_Float(os);

        } else if (extension == string("tet")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename.c_str()<<endl;return -1;}        
            
            TETRAHEDRALIZED_VOLUME<double> *tetVol = new TETRAHEDRALIZED_VOLUME<double>(*new TETRAHEDRON_MESH,*new PARTICLE_3D<double>);
            tetVol->Read(input);
            tetVol->Write_Float(os);

        } else if (extension == string("rgd")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename<<endl;return -1;}
            RIGID_BODY_3D<double> rb;
            rb.Read(input);
            rb.Write_Float(os);

        } else {
            cout << "Invalid extension.  Cannot convert." << endl;
            return -1;
        }
        cout << "Successfully converted " << infilename << " to " << outname << endl;

    //*************************
    // FTOD
    //*************************
    } else if (program == "ftod") {
        cout << "***** converting " << infilename <<" from float to double format..." << endl;
        if (extension == string("tri")) {
            TRIANGULATED_SURFACE<float> *triSurface= 
                TRIANGULATED_SURFACE<float>::Create_From_File(infilename.c_str());        
            triSurface->Write_Double(os);
        
        } else if (extension == string("oct")) {
            OCTREE_IMPLICIT_SURFACE<float> *octSurface =
                OCTREE_IMPLICIT_SURFACE<float>::Create_From_File(infilename.c_str());
            octSurface->Write_Double(os);
    
        } else if (extension == string("phi")) {
            LEVELSET_IMPLICIT_SURFACE<float> *impSurf =
                LEVELSET_IMPLICIT_SURFACE<float>::Create_From_File(infilename.c_str());
            impSurf->Write_Double(os);

        } else if (extension == string("tet")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename.c_str()<<endl;return -1;}        

            TETRAHEDRALIZED_VOLUME<float> *tetVol = new TETRAHEDRALIZED_VOLUME<float>(*new TETRAHEDRON_MESH,*new PARTICLE_3D<float>);
            tetVol->Read(input);
            tetVol->Write_Double(os);

        } else if (extension == string("rgd")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename<<endl;return -1;}
            RIGID_BODY_3D<float> rb;
            rb.Read(input); //read float
            rb.Write_Double(os);
        } else {
            cout << "Invalid extension.  Cannot convert." << endl;
            return -1;
        }
        cout << "Successfully converted " << infilename << " to " << outname << endl;

    //*************************
    // pmodd
    //*************************
    } else if (program == "pmodd") { // particle mod
        cout << "***** adding particle parameter to " << infilename << endl;

        if (extension == string("tri")) {
            TRIANGULATED_SURFACE<double> *triSurface = 
                TRIANGULATED_SURFACE<double>::Create_From_File(infilename.c_str());

            triSurface->Write(os);
        } else if (extension == string("tet")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename<<endl;return -1;}        
            TETRAHEDRALIZED_VOLUME<double> *tetVol = new TETRAHEDRALIZED_VOLUME<double>(*new TETRAHEDRON_MESH,*new PARTICLE_3D<double>);
            tetVol->Read(input);
            tetVol->Write(os);
        }

    //*************************
    // pmodf
    //*************************
    } else if (program == "pmodf") { // particle mod
        cout << "***** adding particle parameter to " << infilename << endl;
        
        if (extension == string("tri")) {
            TRIANGULATED_SURFACE<float> *triSurface = 
                TRIANGULATED_SURFACE<float>::Create_From_File(infilename.c_str());
            triSurface->Write(os);

        } else if (extension == string("tet")) {
            ifstream input(infilename.c_str(),ios::binary);
            if(!input) {cerr<<"Could not open "<<infilename.c_str()<<endl;return -1;}        
            TETRAHEDRALIZED_VOLUME<float> *tetVol = new TETRAHEDRALIZED_VOLUME<float>(*new TETRAHEDRON_MESH,*new PARTICLE_3D<float>);
            tetVol->Read(input);
            tetVol->Write(os);
        }
    } else {
        cout << "program not recognized." << endl;
    }
    os.close();
    return 0;
}
