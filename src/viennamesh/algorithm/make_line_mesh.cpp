/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include "viennamesh/algorithm/make_line_mesh.hpp"

#include "viennagrid/algorithm/norm.hpp"
#include "viennagrid/algorithm/extract_seed_points.hpp"


namespace viennamesh
{
  namespace detail
  {
    template<typename MeshT>
    class vertex_handle_point_sorter_1d
    {
    public:
      typedef typename viennagrid::result_of::const_vertex_handle<MeshT>::type ConstVertexHandleType;

      vertex_handle_point_sorter_1d(MeshT const & mesh_) : mesh(mesh_) {}

      bool operator()( ConstVertexHandleType const & v0, ConstVertexHandleType const & v1 ) const
      {
        return
          viennagrid::point(mesh, viennagrid::dereference_handle(mesh, v0) )[0] <
          viennagrid::point(mesh, viennagrid::dereference_handle(mesh, v1) )[0];
      }

    private:
      MeshT const & mesh;
    };

    template<typename MeshT>
    class point_sorter_1d
    {
    public:
      typedef typename viennagrid::result_of::point<MeshT>::type PointType;

      bool operator()( PointType const & p0, PointType const & p1 ) const
      { return p0 [0] < p1[0]; }

      bool operator()( std::pair<PointType, int> const & sp0, std::pair<PointType, int> const & sp1 ) const
      { return sp0.first[0] < sp1.first[0]; }
    };

    template<typename GeometryT, typename OutputMeshT, typename OutputSegmentationT>
    bool generate_line_mesh( GeometryT const & geometry, OutputMeshT & mesh, OutputSegmentationT & segmentation,
                        double cell_size, bool use_different_segment_ids_for_unknown_segments, double min_geometry_point_distance,
                        seed_point_1d_container const & seed_points_, point_1d_container const & hole_points_ )
    {
      typedef typename viennagrid::result_of::point<GeometryT>::type PointType;

      typedef typename viennagrid::result_of::const_vertex_handle<GeometryT>::type GeometryVertexHandleType;
      typedef typename viennagrid::result_of::vertex_handle<OutputMeshT>::type OutputVertexHandleType;

      typedef typename viennagrid::result_of::const_vertex_range<GeometryT>::type ConstVertexRangeType;
      typedef typename viennagrid::result_of::iterator<ConstVertexRangeType>::type ConstVertexIteratorType;

      seed_point_1d_container seed_points(seed_points_);
      point_1d_container hole_points(hole_points_);

      // copy and sort vertices
      ConstVertexRangeType vertices( geometry );
      std::list<GeometryVertexHandleType> sorted_geometry_points;
      for (ConstVertexIteratorType vit = vertices.begin(); vit != vertices.end(); ++vit)
        sorted_geometry_points.push_back( vit.handle() );
      sorted_geometry_points.sort( vertex_handle_point_sorter_1d<GeometryT>(geometry) );

      // remove points which would lead to too short lines
      {
        typename std::list<GeometryVertexHandleType>::iterator vhit0 = sorted_geometry_points.begin();
        typename std::list<GeometryVertexHandleType>::iterator vhit1 = vhit0; ++vhit1;

        while (vhit1 != sorted_geometry_points.end())
        {
          double length = std::abs(viennagrid::point(geometry, *vhit0)[0] - viennagrid::point(geometry, *vhit1)[0]);

          if (length < min_geometry_point_distance)
          {
            sorted_geometry_points.erase(vhit1++);
          }
          else
          {
            ++vhit0;
            ++vhit1;
          }
        }
      }

      // create points in output mesh
      std::vector<OutputVertexHandleType> sorted_points;
      for (typename std::list<GeometryVertexHandleType>::iterator vhit = sorted_geometry_points.begin(); vhit != sorted_geometry_points.end(); ++vhit)
      {
        OutputVertexHandleType vertex_handle = viennagrid::make_vertex(mesh, viennagrid::point(geometry, *vhit) );
        sorted_points.push_back(vertex_handle);
      }


      // determine unused segment ID
      int default_segment_id = -1;
      for (seed_point_1d_container::const_iterator it = seed_points.begin(); it != seed_points.end(); ++it)
        default_segment_id = std::max(it->second, default_segment_id);
      ++default_segment_id;
      info(10) << "Default segment id: " << default_segment_id << std::endl;


      // sort seed points and hole points
      std::sort( seed_points.begin(), seed_points.end(), point_sorter_1d<GeometryT>() );
      std::sort( hole_points.begin(), hole_points.end(), point_sorter_1d<GeometryT>() );

      seed_point_1d_container::const_iterator spit = seed_points.begin();
      point_1d_container::const_iterator hpit = hole_points.begin();


      // iterate over all intervals in sorted point array
      for (unsigned int i = 1; i < sorted_points.size(); ++i)
      {
        // getting start and end point of the current interval
        double start = viennagrid::point( mesh, sorted_points[i-1])[0];
        double end = viennagrid::point( mesh, sorted_points[i])[0];


        // updating the hole point iterator to fit into the current interval
        while ( (hpit != hole_points.end()) && ((*hpit)[0] < start) )
          ++hpit;

        // if the hole point is inside the current interval -> skip current interval
        if ( (hpit != hole_points.end()) && ((*hpit)[0] < end) )
          continue;

        // updating the seed point iterator to fit into the current interval
        while ( (spit != seed_points.end()) && (spit->first[0] < start) )
          ++spit;


        // determining the segment ID of the current interval
        int output_segment_id;
        if ( (spit != seed_points.end()) && (spit->first[0] < end) )
          output_segment_id = spit->second;
        else
        {
          output_segment_id = default_segment_id;
          if (use_different_segment_ids_for_unknown_segments)
            ++default_segment_id;
        }


        // calculate the interval length and the line count in the current interval
        double length = std::abs(end-start);
        unsigned int new_line_count = 1;
        if (cell_size > 0.0)
          new_line_count = static_cast<unsigned int>(length / cell_size + 0.5);


        // if only one line should be created
        if (new_line_count <= 1)
        {
          viennagrid::make_line( segmentation.get_make_segment(output_segment_id), sorted_points[i-1], sorted_points[i] );
        }
        else
        {
          // create all line segments
          unsigned int lines_in_between = new_line_count - 2;
          double delta = (end-start) / new_line_count;

          OutputVertexHandleType last_vertex = viennagrid::make_vertex( mesh, PointType(start+delta) );
          viennagrid::make_line( segmentation.get_make_segment(output_segment_id), sorted_points[i-1], last_vertex );

          for (unsigned int j = 0; j < lines_in_between; ++j)
          {
            OutputVertexHandleType tmp = viennagrid::make_vertex( mesh, PointType(start+delta*(j+2)) );
            viennagrid::make_line( segmentation.get_make_segment(output_segment_id), last_vertex, tmp );

            last_vertex = tmp;
          }

          viennagrid::make_line( segmentation.get_make_segment(output_segment_id), last_vertex, sorted_points[i] );
        }
      }

      return true;
    }
  }




  make_line_mesh::make_line_mesh() :
    input_mesh(*this, parameter_information("mesh","mesh","The input mesh, segmented and non-segmented brep 1d mesh supported")),
    input_seed_points(*this, parameter_information("seed_points","seed_point_1d_container","A container of seed points with corresonding segment names")),
    input_hole_points(*this, parameter_information("hole_points","point_1d_container","A container of hole points")),
    cell_size(*this, parameter_information("cell_size","double","The desired maximum size of tetrahedrons, all tetrahedrons will be at most this size"), -1.0),
    use_different_segment_for_unknown_segments(*this, parameter_information("use_different_segment_for_unknown_segments","bool","Determines, if different segments should be used for elements with unkonwn segment or if the should all be assigned to one segment"), false),
    extract_segment_seed_points(*this, parameter_information("extract_segment_seed_points","bool","Should seed points be extracted from the input segmentation?"), true),
    relative_min_geometry_point_distance(*this, parameter_information("relative_min_geometry_point_distance","double","Relative value for checking if two close points should be merged"), 1e-10),
    absolute_min_geometry_point_distance(*this, parameter_information("absolute_min_geometry_point_distance","double","Absolute value for checking if two close points should be merged, overwrites relative_min_geometry_point_distance")),
    output_mesh(*this, parameter_information("mesh", "mesh", "The output mesh, segmented line 1d")) {}

  std::string make_line_mesh::name() const { return "ViennaGrid 1D line mesher"; }
  std::string make_line_mesh::id() const { return "make_line_mesh"; }



  template<typename GeometrySegmentationT>
  void make_line_mesh::extract_seed_points( GeometrySegmentationT const & segmentation, point_1d_container const & hole_points,
                            seed_point_1d_container & seed_points )
  {
    typedef typename viennagrid::result_of::point<GeometrySegmentationT>::type PointType;

    int highest_segment_id = -1;
    for (seed_point_1d_container::iterator spit = seed_points.begin(); spit != seed_points.end(); ++spit)
      highest_segment_id = std::max( highest_segment_id, spit->second );
    ++highest_segment_id;

    for (typename GeometrySegmentationT::const_iterator sit = segmentation.begin(); sit != segmentation.end(); ++sit)
    {
      viennagrid::line_1d_mesh viennagrid_mesh;
      viennagrid::line_1d_segmentation viennagrid_segmentation(viennagrid_mesh);

      detail::generate_line_mesh<typename GeometrySegmentationT::segment_handle_type, viennagrid::line_1d_mesh, viennagrid::line_1d_segmentation>( *sit, viennagrid_mesh, viennagrid_segmentation, -1, false, -1, seed_point_1d_container(), hole_points );

      std::vector<PointType> local_seed_points;
      viennagrid::extract_seed_points( viennagrid_mesh, local_seed_points );
      for (unsigned int i = 0; i < local_seed_points.size(); ++i)
      {
        info(5) << "Found seed point: " << local_seed_points[i] << std::endl;
        seed_points.push_back( std::make_pair(local_seed_points[i], highest_segment_id) );
      }
      highest_segment_id++;
    }
  }



  template<typename GeometryT, typename GeometrySegmentationT, typename OutputMeshT, typename OutputSegmentationT>
  bool make_line_mesh::generic_run_impl()
  {
    typedef viennagrid::segmented_mesh<GeometryT, GeometrySegmentationT> SegmentedGeometryType;

    typedef typename viennamesh::result_of::const_parameter_handle<GeometryT>::type GeometryHandleType;
    typedef typename viennamesh::result_of::const_parameter_handle<SegmentedGeometryType>::type SegmentedGeometryHandleType;

    typedef typename viennagrid::segmented_mesh<OutputMeshT, OutputSegmentationT> OutputSegmentedMesh;

    typedef typename viennagrid::result_of::const_vertex_handle<GeometryT>::type GeometryVertexHandleType;

    typedef typename viennagrid::result_of::const_vertex_range<GeometryT>::type ConstVertexRangeType;
    typedef typename viennagrid::result_of::iterator<ConstVertexRangeType>::type ConstVertexIteratorType;


    GeometryT const * geometry = NULL;
    GeometrySegmentationT const * geometry_segmentation = NULL;

    // query input parameters
    SegmentedGeometryHandleType segmented_geometry_handle = input_mesh.get<SegmentedGeometryType>();
    if (segmented_geometry_handle)
    {
      geometry = &(segmented_geometry_handle().mesh);
      geometry_segmentation = &(segmented_geometry_handle().segmentation);
    }
    else
    {
      GeometryHandleType geometry_handle = input_mesh.get<GeometryT>();
      if (!geometry_handle)
        return false;

      geometry = &(geometry_handle());
    }


    // query possible output parameters: mesh and segmented mesh
    output_parameter_proxy<OutputSegmentedMesh> omp(output_mesh);

    info(10) << "Using cell size: " << cell_size() << std::endl;
    info(10) << "Using different segment IDs for unknown segments: " << std::boolalpha << use_different_segment_for_unknown_segments() << std::endl;

    // query seed points input parameter
    seed_point_1d_container seed_points;
    if (input_seed_points.valid() && !input_seed_points().empty())
    {
      info(10) << "Found seed points -> enabling make_segmented_mesh" << std::endl;
      seed_points = input_seed_points();
    }

    // query hole points input parameter
    point_1d_container hole_points;
    if (input_hole_points.valid() && !input_hole_points().empty())
    {
      info(10) << "Found hole points" << std::endl;
      hole_points = input_hole_points();
    }


    if (extract_segment_seed_points() && geometry_segmentation && geometry_segmentation->size() >= 1)
    {
      extract_seed_points( *geometry_segmentation, hole_points, seed_points );
    }



    // decide, if mesh or segmented mesh is used (seed points available?)
    OutputMeshT & mesh = omp().mesh;
    OutputSegmentationT & segmentation = omp().segmentation;

    // copy and sort vertices
    ConstVertexRangeType vertices( *geometry );
    std::list<GeometryVertexHandleType> sorted_geometry_points;
    for (ConstVertexIteratorType vit = vertices.begin(); vit != vertices.end(); ++vit)
      sorted_geometry_points.push_back( vit.handle() );
    sorted_geometry_points.sort( detail::vertex_handle_point_sorter_1d<GeometryT>(*geometry) );


    // query and determine minimal line length
    double min_distance =
      (viennagrid::point(*geometry, sorted_geometry_points.back())[0] -
      viennagrid::point(*geometry, sorted_geometry_points.front())[0]) * relative_min_geometry_point_distance();

    if (absolute_min_geometry_point_distance.valid())
      min_distance = absolute_min_geometry_point_distance();

    detail::generate_line_mesh( *geometry, mesh, segmentation, cell_size(), use_different_segment_for_unknown_segments(), min_distance, seed_points, hole_points );

    return true;
  }


  bool make_line_mesh::run_impl()
  {
    if (generic_run_impl<viennagrid::brep_1d_mesh, viennagrid::brep_1d_segmentation, viennagrid::line_1d_mesh, viennagrid::line_1d_segmentation>())
      return true;

    return false;
  }

}
