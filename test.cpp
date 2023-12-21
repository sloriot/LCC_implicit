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

//~ typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
//~ typedef LCC_3::Dart_descriptor Dart_descriptor;
//~ typedef LCC_3::Point           Point;

using DT = CGAL::Delaunay_triangulation_3<K>;


int main()
{
  const int nb_points = 10000;

  std::vector<Point> v;
  v.reserve(nb_points);

  CGAL::Random_points_in_sphere_3<Point> gen (2.0);
  for (int i = 0; i < nb_points; ++i)  v.push_back (*gen++);


  DT dt(v.begin(), v.end());

  LCC lcc;
  std::map<DT::Cell_handle, Dart_descriptor > vol_to_dart;

  Dart_descriptor d=CGAL::import_from_triangulation_3(lcc, dt, &vol_to_dart);

  CGAL::draw(lcc);

}
