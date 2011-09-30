#ifndef __DATA_EXCHANGE_GEOMETRY__
#define __DATA_EXCHANGE_GEOMETRY__

#include "primitives.h"
namespace data_exchange{

struct polygon_mesh
{
    std::vector<int> polygons; // counterclockwise order
    std::vector<int> polygon_counts; // number of vertices per polygon

    polygon_mesh() {}

    template<int d>
    void insert_polygon(fixed_vector<int,d>& vertices)
    {
        insert_polygon(vertices.data, d);
    }

    void insert_polygon(int* vertices, int n)
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
