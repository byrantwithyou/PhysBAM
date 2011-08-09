//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Geometry/Close_Mesh_Object.h>
#include <Jeffrey_Utilities/Geometry/Divide_Cell2.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_ONE_OVER_VOLUME_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Grid/VISIT_IF_SIGN_PREDICATE_GRID_VISITOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Math/Equal_Relative_Tolerance.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <Jeffrey_Utilities/Write_Mesh_Object_To_File.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/MAIN_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"

#include "Lagrangify_Level_Set.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T, int D >
struct EMBEDDED_SURFACE;

template< class T >
struct EMBEDDING_CELL_STATISTICS;

template< class T, int D >
struct DIVIDE_CELL_VISITOR;

} // namespace

template< class T, int D >
int Lagrangify_Level_Set(
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index)
{
    namespace qi = boost::spirit::qi;

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    std::cout << "Resolving level set against grid..." << std::endl;
    BASIC_TIMER timer;

    {
        BASIC_TIMER timer;

        const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

        const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
        const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);

        typedef typename EMBEDDED_SURFACE<T,D>::type EMBEDDED_SURFACE_TYPE;
        boost::scoped_ptr< EMBEDDED_SURFACE_TYPE > p_embedded_surface(EMBEDDED_SURFACE_TYPE::Create());
        EMBEDDED_SURFACE_TYPE& embedded_surface = *p_embedded_surface;
        boost::mutex embedded_surface_mutex;

        EMBEDDING_CELL_STATISTICS<T> embedding_cell_statistics;
        boost::mutex embedding_cell_statistics_mutex;

        std::cout << "Dividing embedding cells..." << std::endl;
        timer.Restart();
        Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT<2>(
            main_params.general.n_thread,
            cell_multi_index_bound,
            Make_Compose_Function(
                SIGN_FUNCTION(),
                phi_of_fine_index,
                fine_multi_index_bound
            ),
            Make_Visit_If_Sign_Predicate_Grid_Visitor(
                Make_Equal_Function(0),
                DIVIDE_CELL_VISITOR<T,D>(
                    main_params,
                    phi_of_fine_index,
                    embedded_surface, embedded_surface_mutex,
                    embedding_cell_statistics, embedding_cell_statistics_mutex
                )
            )
        );

        std::cout << "  # of embedding cells = " << embedding_cell_statistics.n << std::endl;
        std::cout << "  min (-) (normalized) volume measure = " << embedding_cell_statistics.min_negative_volume_measure << std::endl;
        std::cout << "  max (-) (normalized) volume measure = " << embedding_cell_statistics.max_negative_volume_measure << std::endl;
        std::cout << "  min (+) (normalized) volume measure = " << embedding_cell_statistics.min_positive_volume_measure << std::endl;
        std::cout << "  max (+) (normalized) volume measure = " << embedding_cell_statistics.max_positive_volume_measure << std::endl;

        embedded_surface.Update_Number_Nodes();

        std::cout << "  # of particles (before closing) = " << embedded_surface.particles.array_collection->Size()     << std::endl;
        std::cout << "  # of elements  (before closing) = " << embedded_surface.mesh.elements.Size() << std::endl;
        Close_Mesh_Object(embedded_surface);
        std::cout << "  # of particles  (after closing) = " << embedded_surface.particles.array_collection->Size()     << std::endl;
        std::cout << "  # of elements   (after closing) = " << embedded_surface.mesh.elements.Size() << std::endl;

        std::cout << "[Dividing embedding cells...] " << timer.Elapsed() << " s" << std::endl;

        {
            const std::string& filename = main_params.output.embedded_surface_filename;
            if(!filename.empty()) {
                std::cout << "Writing embedded surface to file \"" << filename << "\"...";
                timer.Restart();
                Write_Mesh_Object_To_File(embedded_surface, filename);
                std::cout << timer.Elapsed() << " s" << std::endl;
            }
        }
    }

    std::cout << "[Resolving level set against grid...] "
              << timer.Elapsed() << " s"
              << std::endl;

    return 0;
}

namespace
{

template< class T > struct EMBEDDED_SURFACE<T,2> { typedef SEGMENTED_CURVE_2D<T> type; };
template< class T > struct EMBEDDED_SURFACE<T,3> { typedef TRIANGULATED_SURFACE<T> type; };

template< class T >
struct EMBEDDING_CELL_STATISTICS
{
    int n;
    T min_negative_volume_measure;
    T max_negative_volume_measure;
    T min_positive_volume_measure;
    T max_positive_volume_measure;

    EMBEDDING_CELL_STATISTICS()
        : n(0),
          min_negative_volume_measure(std::numeric_limits<T>::infinity()),
          max_negative_volume_measure(0),
          min_positive_volume_measure(std::numeric_limits<T>::infinity()),
          max_positive_volume_measure(0)
    { }
};

template< class T, int D >
struct DIVIDE_CELL_VISITOR
{
    typedef typename EMBEDDED_SURFACE<T,D>::type EMBEDDED_SURFACE_TYPE;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        DIVIDE_CELL_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( MAIN_PARAMS<T,D> )) const &, main_params ))
        (( typename ARRAY_VIEW<const T> const, phi_of_fine_index ))
        (( typename EMBEDDED_SURFACE_TYPE&, embedded_surface ))
        (( /******/ boost::mutex&, embedded_surface_mutex ))
        (( typename EMBEDDING_CELL_STATISTICS<T>&, embedding_cell_statistics ))
        (( /******/ boost::mutex&, embedding_cell_statistics_mutex ))
    )
public:

    typedef void result_type;

    void operator()(const int cell_linear_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;

        const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

        const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
        const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
        const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;
        const T dv = dx.Product();

        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const VECTOR<T,D> x0 = Multi_Index_X(min_x, max_x, fine_multi_index_bound, 2 * cell_multi_index);

        boost::scoped_ptr< EMBEDDED_SURFACE_TYPE > p_cell_local_embedded_surface(EMBEDDED_SURFACE_TYPE::Create());
        EMBEDDED_SURFACE_TYPE& cell_local_embedded_surface = *p_cell_local_embedded_surface;

        T negative_volume_measure = 0;
        T positive_volume_measure = 0;

        Divide_Cell2(
            dx / 2,
            cell_multi_index,
            Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
            Make_Visitor_Sequence(
                POLYTOPE_VISITOR(x0, cell_local_embedded_surface),
                Make_Integrate_One_Over_Volume_Polytope_Visitor(-1, negative_volume_measure),
                Make_Integrate_One_Over_Volume_Polytope_Visitor(+1, positive_volume_measure)
            ),
            -1 // sign_of_zero
        );

        assert(0 <= negative_volume_measure && negative_volume_measure <= dv);
        assert(0 <= positive_volume_measure && positive_volume_measure <= dv);
        const T volume_measure = negative_volume_measure + positive_volume_measure;
        static_cast<void>(volume_measure);
        assert(Equal_Relative_Tolerance<16>(volume_measure, dv));

        cell_local_embedded_surface.Update_Number_Nodes();
        Close_Mesh_Object(cell_local_embedded_surface);

        {
            const std::string& filename_format =
                main_params.output.cell_local_embedded_surface_filename_format;
            if(!filename_format.empty()) {
                boost::format filename_formatter(filename_format);
                for(int d = 1; d <= D; ++d)
                    filename_formatter % cell_multi_index[d];
                const std::string filename = filename_formatter.str();
                Write_Mesh_Object_To_File(cell_local_embedded_surface, filename);
            }
        }

        {
            boost::lock_guard< boost::mutex > _embedded_surface_scoped_lock(embedded_surface_mutex);
            const int particle_offset = embedded_surface.particles.array_collection->Size();
            embedded_surface.particles.array_collection->Add_Elements(
                cell_local_embedded_surface.particles.array_collection->Size()
            );
            for(int p = 1; p <= cell_local_embedded_surface.particles.array_collection->Size(); ++p)
                embedded_surface.particles.X(particle_offset + p) =
                    cell_local_embedded_surface.particles.X(p);
            for(int e = 1; e <= cell_local_embedded_surface.mesh.elements.Size(); ++e)
                embedded_surface.mesh.elements.Append(
                    particle_offset + cell_local_embedded_surface.mesh.elements(e)
                );
        }

        {
            boost::lock_guard< boost::mutex > _embedding_cell_statistics_scoped_lock(embedding_cell_statistics_mutex);
            ++embedding_cell_statistics.n;
            embedding_cell_statistics.min_negative_volume_measure =
                std::min(embedding_cell_statistics.min_negative_volume_measure, negative_volume_measure / dv);
            embedding_cell_statistics.max_negative_volume_measure =
                std::max(embedding_cell_statistics.max_negative_volume_measure, negative_volume_measure / dv);
            embedding_cell_statistics.min_positive_volume_measure =
                std::min(embedding_cell_statistics.min_positive_volume_measure, positive_volume_measure / dv);
            embedding_cell_statistics.max_positive_volume_measure =
                std::max(embedding_cell_statistics.max_positive_volume_measure, positive_volume_measure / dv);
        }
    }

private:
    struct POLYTOPE_VISITOR;
};

template< class T, int D >
struct DIVIDE_CELL_VISITOR<T,D>::
POLYTOPE_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POLYTOPE_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, x0 ))
        (( typename EMBEDDED_SURFACE_TYPE&, cell_local_embedded_surface ))
    )
public:
    typedef void result_type;
    template< class T_VERTICES_X >
    void operator()(
        const int sign,
        const T_VERTICES_X& vertices_x,
        const SUB_CUBE_POLYTOPE<T,D>& polytope) const
    {
        if(sign != -1 || polytope.alignment != SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED)
            return;
        polytope.Add_To(cell_local_embedded_surface, x0, vertices_x);
    }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Lagrangify_Level_Set<T,D>( \
    const MAIN_PARAMS<T,D>& main_params, \
    RAND_MT19937_UNIFORM_REAL<T>::type& rand, \
    const ARRAY_VIEW<const T> phi_of_fine_index);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
