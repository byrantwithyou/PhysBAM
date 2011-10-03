#ifndef __DATA_EXCHANGE_POLYGON_MESH__
#define __DATA_EXCHANGE_POLYGON_MESH__

#include "fixed_vector.h"
namespace data_exchange{

/*
  This class describes the surface connectivity of a surface mesh.  Two helper
  functions are provided for populating this class.  Polygons are listed in
  counterclockwise order.

  Members:

  * polygons
  
  This list stores the indices of the faces as a flat array.  The indices of
  each polygon should be listed in counterclockwise order, but the polygons
  themselves can be listed in any order.  The polygons need not be triangles,
  but they will be converted into triangles before simulation.  Specifying
  arbitrary polygons is potentially useful if additional data is later
  associated with polygons, since then the indexing can be carried over to the
  simulator.

  * polygon_counts

  These numbers indicate the number of points in each polygon.  The sum of the
  integers in this list should equal the length of the polygons array.  If the
  first entry is a 3, then the first three integers in the polygons array
  describe a triangle.  If it is a 4, then the first four integers of the
  polygons array describe a quadralateral, and so on.  The order of the integers
  in this list corresponds to the order in which polygons are listed in the
  polygons array.
 */

struct polygon_mesh
{
    std::vector<int> polygons; // counterclockwise order
    std::vector<int> polygon_counts; // number of vertices per polygon

    polygon_mesh() {}

    template<int d>
    void insert_polygon(const fixed_vector<int,d>& vertices)
    {
        insert_polygon(vertices.data, d);
    }

    void insert_polygon(const int* vertices, int n)
    {
        for(int i=0; i<n; i++)
            polygons.push_back(vertices[i]);
        polygon_counts.push_back(n);
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & polygons & polygon_counts;
    }
};

}

#endif
