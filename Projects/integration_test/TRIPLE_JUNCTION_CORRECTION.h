#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
namespace PhysBAM{
template<class TV>
class TRIPLE_JUNCTION_CORRECTION
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MARCHING_CUBES_COLOR<TV>::HASH_CELL_DATA HASH_CELL_DATA;
    typedef typename MARCHING_CUBES_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef VECTOR<T,(1<<TV::m)> PHI;
public:
    enum WORKAROUND {num_stencils=TV::m==2?2:5};
    const GRID<TV>& grid;
    ARRAY<ARRAY<T,TV_INT> >& phi;
    int ghost;
    T default_phi;
    T trust_buffer; // distance from sonic points to trusted region
    T valid_width; // distance from interface we care about
    T extent; // compute pairwise level set far enough so we can exprapolate
    int extrap_width; // how far to extrapolate

    struct PAIRWISE_LEVEL_SET_DATA
    {
        T phi;
        int valid_flags;
        VECTOR<short,2> trust;

        PAIRWISE_LEVEL_SET_DATA();
    };

    ARRAY<PAIRWISE_LEVEL_SET_DATA,TV_INT> pairwise_data;
    ARRAY<ARRAY<ARRAY<T,TV_INT> > > pairwise_phi;
    ARRAY<T,TV_INT> combined_phi;
    ARRAY<int,TV_INT> combined_color;

    TRIPLE_JUNCTION_CORRECTION(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi,int ghost);
    void Compute_Pairwise_Data();
    void Initialize_Pairwise_Level_Set();
    void Fill_Valid_Region_With_Exprapolation();
    void One_Step_Triple_Junction_Correction();
    void Update_Color_Level_Sets();
    void Cut_Interface(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data);
    void Cut_Stencil_With_Phi(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const TV_INT& cell,int s);
    void Cut_Stencil_With_Pairwise_Phi(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const TV_INT& cell,int s);
    void Compute_Pairwise_Level_Set_Data();
    void Evolve_Step(ARRAY<T,TV_INT> q[3],const ARRAY<T,TV_INT> p[3]);
    static TV Zero_Phi(const VECTOR<PHI,3>& phi,VECTOR<T,3>& p);
    static TV Meet_Phi(const VECTOR<PHI,2>& phi);
    void Cut_Interface(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data);
    void Cut_Stencil_With_Phi(CELL_ELEMENTS& ce,const TV_INT& cell,const VECTOR<TV_INT,TV::m+1>& st);
    void Cut_Stencil_With_Pairwise_Phi(CELL_ELEMENTS& ce,const TV_INT& cell,const VECTOR<TV_INT,TV::m+1>& st);
    void Fill_Combined_Level_Set();
    int Fill_Combined_Level_Set_At_Index(const TV_INT& node);
    int Fill_Phi_From_Pairwise_Level_Set_At_Index(const TV_INT& node,int color);
};
}
