#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/draw_linear_cell_complex.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;


struct LCC_items
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef std::tuple<Vertex_attrib> Attributes;
  };
};

using LCC_traits = CGAL::Linear_cell_complex_traits<3, K>;
using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, LCC_traits, LCC_items>;
using Dart_descriptor =  LCC::Dart_descriptor;

template < typename GT,
           typename Vb = CGAL::Triangulation_vertex_base_3<GT> >
class Triangulation_vertex_base_with_id_3
  : public Vb
{
  std::size_t _id;
public:

  typedef typename Vb::Cell_handle                   Cell_handle;
  typedef typename Vb::Point                         Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Triangulation_vertex_base_with_id_3<GT, Vb2>   Other;
  };

  Triangulation_vertex_base_with_id_3()
    : Vb() {}

  Triangulation_vertex_base_with_id_3(const Point & p)
    : Vb(p) {}

  Triangulation_vertex_base_with_id_3(const Point & p, Cell_handle c)
    : Vb(p, c) {}

  Triangulation_vertex_base_with_id_3(Cell_handle c)
    : Vb(c) {}

  const std::size_t& id() const { return _id; }
  std::size_t&       id()       { return _id; }
};


using Vb = Triangulation_vertex_base_with_id_3<K>;
using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
using TDS = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using DT = CGAL::Delaunay_triangulation_3<K, TDS>;


// modified version of CGAL::import_from_triangulation_3
template < class LCC, class Triangulation >
typename LCC::Dart_descriptor fill_lcc
(LCC& alcc, const Triangulation &atr,
 std::map<typename Triangulation::Cell_handle,
          typename LCC::Dart_descriptor >* avol_to_dart=nullptr)
{
  static_assert( LCC::dimension>=3 && LCC::ambient_dimension==3 );

  // Case of empty triangulations.
  if (atr.number_of_vertices() == 0) return LCC::null_descriptor;

  // Check the dimension.
  if (atr.dimension() != 3) return LCC::null_descriptor;
  CGAL_assertion(atr.is_valid());

  typedef typename Triangulation::Vertex_handle    TVertex_handle;
  typedef typename Triangulation::Vertex_iterator  TVertex_iterator;
  typedef typename Triangulation::Cell_iterator    TCell_iterator;
  typedef typename std::map
    < TCell_iterator, typename LCC::Dart_descriptor >::iterator itmap_tcell;

  // Create vertices in the map and associate in a map
  // TVertex_handle and vertices in the map.
  std::map< TVertex_handle, typename LCC::Vertex_attribute_descriptor > TV;
  for (TVertex_iterator itv = atr.vertices_begin();
       itv != atr.vertices_end(); ++itv)
  {
    TV[itv] = alcc.create_vertex_attribute(itv->point()/* , itv->id() */);
  }

  // Create the tetrahedron and create a map to link Cell_iterator
  // and tetrahedron.
  TCell_iterator it;

  std::map<typename Triangulation::Cell_handle, typename LCC::Dart_descriptor> TC;
  std::map<typename Triangulation::Cell_handle, typename LCC::Dart_descriptor>*
    mytc = (avol_to_dart==nullptr?&TC:avol_to_dart);

  itmap_tcell maptcell_it;

  typename LCC::Dart_descriptor res=LCC::null_descriptor, dart=LCC::null_descriptor;
  typename LCC::Dart_descriptor cur=LCC::null_descriptor, neighbor=LCC::null_descriptor;

  for (it = atr.cells_begin(); it != atr.cells_end(); ++it)
  {
    {
      res = alcc.make_tetrahedron(TV[it->vertex(0)],
                                  TV[it->vertex(1)],
                                  TV[it->vertex(2)],
                                  TV[it->vertex(3)]);

      if ( dart==LCC::null_descriptor )
      {
        if ( it->vertex(0) == atr.infinite_vertex() )
          dart = res;
        else if ( it->vertex(1) == atr.infinite_vertex() )
          dart = alcc.next(res);
        else if ( it->vertex(2) == atr.infinite_vertex() )
          dart = alcc.previous(res);
        else if ( it->vertex(3) == atr.infinite_vertex() )
          dart = alcc.previous(alcc.template opposite<2>(res));
      }

      for (unsigned int i = 0; i < 4; ++i)
      {
        switch (i)
        {
        case 0: cur = alcc.template opposite<2>(alcc.next(res)); break;
        case 1: cur = alcc.template opposite<2>(alcc.previous(res)); break;
        case 2: cur = alcc.template opposite<2>(res); break;
        case 3: cur = res; break;
        }

        maptcell_it = mytc->find(it->neighbor(i));
        if (maptcell_it != mytc->end())
        {
          switch (atr.mirror_index(it,i) )
          {
          case 0: neighbor = alcc.template opposite<2>(alcc.next(maptcell_it->second));
            break;
          case 1: neighbor = alcc.template opposite<2>(alcc.previous(maptcell_it->second));
            break;
          case 2: neighbor = alcc.template opposite<2>(maptcell_it->second); break;
          case 3: neighbor = maptcell_it->second; break;
          }
          while (alcc.vertex_attribute(neighbor) !=
                 alcc.vertex_attribute(alcc.other_extremity(cur)) )
            neighbor = alcc.next(neighbor);
          alcc.template topo_sew<3>(cur, alcc.other_orientation(neighbor));
        }
      }
      (*mytc)[it] = res;
    }
  }
  CGAL_assertion(dart!=LCC::null_descriptor);
  return dart;
}

int main()
{
  const int nb_points = 10000;

  std::vector<Point> v;
  v.reserve(nb_points);

  CGAL::Random_points_in_sphere_3<Point> gen (2.0);
  for (int i = 0; i < nb_points; ++i)  v.push_back (*gen++);


  DT dt(v.begin(), v.end());
  std::size_t id=0;
  std::vector<DT::Vertex_handle> vertices(nb_points);
  for (DT::Vertex_handle vh : dt.finite_vertex_handles())
  {
    vh->id()=id++;
    vertices[vh->id()]=vh;
  }


  LCC lcc;
  std::map<DT::Cell_handle, Dart_descriptor > vol_to_dart;

  Dart_descriptor d=fill_lcc(lcc, dt, &vol_to_dart);

  CGAL::draw(lcc);

}
