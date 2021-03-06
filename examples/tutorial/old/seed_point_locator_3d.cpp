#include "viennamesh/algorithm/viennagrid.hpp"
#include "viennamesh/algorithm/tetgen.hpp"


int main()
{
  typedef viennagrid::plc_3d_mesh GeometryMeshType;

  typedef viennagrid::result_of::point<GeometryMeshType>::type PointType;
  typedef viennagrid::result_of::vertex_handle<GeometryMeshType>::type GeometryVertexHandle;
  typedef viennagrid::result_of::line_handle<GeometryMeshType>::type GeometryLineHandle;

  // creating the geometry mesh
  viennamesh::result_of::parameter_handle< GeometryMeshType >::type geometry_handle = viennamesh::make_parameter<GeometryMeshType>();
  GeometryMeshType & geometry = geometry_handle();


  double s = 10.0;
  GeometryVertexHandle vtx[8];
  GeometryLineHandle lines[12];

  vtx[0] = viennagrid::make_vertex( geometry, PointType(0, 0, 0) );
  vtx[1] = viennagrid::make_vertex( geometry, PointType(0, s, 0) );
  vtx[2] = viennagrid::make_vertex( geometry, PointType(s, 0, 0) );
  vtx[3] = viennagrid::make_vertex( geometry, PointType(s, s, 0) );

  vtx[4] = viennagrid::make_vertex( geometry, PointType(0, 0, s) );
  vtx[5] = viennagrid::make_vertex( geometry, PointType(0, s, s) );
  vtx[6] = viennagrid::make_vertex( geometry, PointType(s, 0, s) );
  vtx[7] = viennagrid::make_vertex( geometry, PointType(s, s, s) );


  // bottom 4 lines
  lines[0] = viennagrid::make_line( geometry, vtx[0], vtx[1] );
  lines[1] = viennagrid::make_line( geometry, vtx[1], vtx[3] );
  lines[2] = viennagrid::make_line( geometry, vtx[3], vtx[2] );
  lines[3] = viennagrid::make_line( geometry, vtx[2], vtx[0] );

  // top 4 lines
  lines[4] = viennagrid::make_line( geometry, vtx[4], vtx[5] );
  lines[5] = viennagrid::make_line( geometry, vtx[5], vtx[7] );
  lines[6] = viennagrid::make_line( geometry, vtx[7], vtx[6] );
  lines[7] = viennagrid::make_line( geometry, vtx[6], vtx[4] );

  // columns
  lines[8] = viennagrid::make_line( geometry, vtx[0], vtx[4] );
  lines[9] = viennagrid::make_line( geometry, vtx[1], vtx[5] );
  lines[10] = viennagrid::make_line( geometry, vtx[2], vtx[6] );
  lines[11] = viennagrid::make_line( geometry, vtx[3], vtx[7] );



  viennagrid::make_plc( geometry, lines+0, lines+4 );
  viennagrid::make_plc( geometry, lines+4, lines+8 );

  {
    GeometryLineHandle cur_lines[4];
    cur_lines[0] = lines[0];
    cur_lines[1] = lines[9];
    cur_lines[2] = lines[4];
    cur_lines[3] = lines[8];
    viennagrid::make_plc( geometry, cur_lines+0, cur_lines+4 );
  }

  {
    GeometryLineHandle cur_lines[4];
    cur_lines[0] = lines[2];
    cur_lines[1] = lines[11];
    cur_lines[2] = lines[6];
    cur_lines[3] = lines[10];
    viennagrid::make_plc( geometry, cur_lines+0, cur_lines+4 );
  }

  {
    GeometryLineHandle cur_lines[4];
    cur_lines[0] = lines[3];
    cur_lines[1] = lines[10];
    cur_lines[2] = lines[7];
    cur_lines[3] = lines[8];
    viennagrid::make_plc( geometry, cur_lines+0, cur_lines+4 );
  }

  {
    GeometryLineHandle cur_lines[4];
    cur_lines[0] = lines[1];
    cur_lines[1] = lines[11];
    cur_lines[2] = lines[5];
    cur_lines[3] = lines[9];
    viennagrid::make_plc( geometry, cur_lines+0, cur_lines+4 );
  }



  // creating the seed point locator algorithm
  viennamesh::algorithm_handle mesher( new viennamesh::tetgen::algorithm() );
  viennamesh::algorithm_handle extract_seed_points( new viennamesh::extract_seed_points::algorithm() );

  mesher->set_input( "default", geometry_handle );
  extract_seed_points->link_input( "default", mesher, "default" );

  mesher->run();
  extract_seed_points->run();

  typedef viennamesh::result_of::seed_point_container<PointType>::type PointContainerType;
  viennamesh::result_of::parameter_handle<PointContainerType>::type point_container = extract_seed_points->get_output<PointContainerType>( "default" );
  if (point_container)
  {
    std::cout << "Number of extracted seed points: " << point_container().size() << std::endl;
    for (PointContainerType::iterator it = point_container().begin(); it != point_container().end(); ++it)
      std::cout << "  " << it->first << " " << it->second << std::endl;
  }
}
