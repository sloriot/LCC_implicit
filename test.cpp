#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/draw_linear_cell_complex.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3,3, CGAL::Linear_cell_complex_traits<3, K> >;
using Dart_descriptor =  LCC::Dart_descriptor;
using Point = K::Point_3;

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

  Dart_descriptor d=CGAL::import_from_triangulation_3(lcc, dt, &vol_to_dart);

  CGAL::draw(lcc);

}
