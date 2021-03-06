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

#include "viennamesh/algorithm/project_mesh.hpp"

namespace viennamesh
{

  template<typename InputMeshT, typename InputSegmentationT, typename OutputMeshT, typename OutputSegmentationT>
  void project_mesh_impl( InputMeshT const & input_mesh, InputSegmentationT const & input_segmentation,
                          OutputMeshT & output_mesh, OutputSegmentationT & output_segmentation )
  {
    typedef typename viennagrid::result_of::point<InputMeshT>::type InputPointType;
    typedef typename viennagrid::result_of::point<OutputMeshT>::type OutputPointType;

    typedef typename viennagrid::result_of::cell<InputMeshT>::type InputCellType;
    typedef typename viennagrid::result_of::cell_id<InputMeshT>::type InputCellIDType;

    typedef typename viennagrid::result_of::cell_handle<OutputMeshT>::type OutputCellHandleType;

    typedef typename viennagrid::result_of::vertex_id<InputMeshT>::type InputVertexIDType;
    typedef typename viennagrid::result_of::vertex_handle<OutputMeshT>::type OutputVertexHandleType;

    std::map<InputVertexIDType, OutputVertexHandleType> vertex_map;

    typedef typename viennagrid::result_of::const_vertex_range<InputMeshT>::type ConstInputVertexRangeType;
    typedef typename viennagrid::result_of::iterator<ConstInputVertexRangeType>::type ConstInputVertexIteratorType;

    ConstInputVertexRangeType vertices(input_mesh);
    for (ConstInputVertexIteratorType vit = vertices.begin(); vit != vertices.end(); ++vit)
    {
      InputPointType in_point = viennagrid::point(*vit);
      OutputPointType out_point;
      for (unsigned int i = 0; i != out_point.size(); ++i)
        out_point[i] = in_point[i];

      vertex_map[ vit->id() ] = viennagrid::make_vertex( output_mesh, out_point );
    }

    std::map<InputCellIDType, OutputCellHandleType> cell_map;

    typedef typename viennagrid::result_of::segment_handle<InputSegmentationT>::type InputSegmentHandleType;
    typedef typename viennagrid::result_of::segment_handle<OutputSegmentationT>::type OutputSegmentHandleType;

    for (typename InputSegmentationT::const_iterator sit = input_segmentation.begin(); sit != input_segmentation.end(); ++sit)
    {
      OutputSegmentHandleType & output_segment = output_segmentation.get_make_segment( sit->id() );

      typedef typename viennagrid::result_of::const_cell_range<InputSegmentHandleType>::type ConstInputCellRangeType;
      typedef typename viennagrid::result_of::iterator<ConstInputCellRangeType>::type ConstInputCellIteratorType;

      ConstInputCellRangeType cells(*sit);
      for (ConstInputCellIteratorType cit = cells.begin(); cit != cells.end(); ++cit)
      {
        typedef typename viennagrid::result_of::const_vertex_range<InputCellType>::type ConstVertexOnInputCellRangeType;
        typedef typename viennagrid::result_of::iterator<ConstVertexOnInputCellRangeType>::type ConstVertexOnInputCellIteratorType;

        typename std::map<InputCellIDType, OutputCellHandleType>::iterator cmit = cell_map.find( cit->id() );
        if (cmit != cell_map.end())
        {
          viennagrid::add( output_segment, viennagrid::dereference_handle(output_mesh, cmit->second) );
        }
        else
        {
          ConstVertexOnInputCellRangeType vertices_on_cell( *cit );
          std::vector<OutputVertexHandleType> vertex_handles;
          for (ConstVertexOnInputCellIteratorType vit = vertices_on_cell.begin(); vit != vertices_on_cell.end(); ++vit)
            vertex_handles.push_back( vertex_map[vit->id()] );

          cell_map[ cit->id() ] = viennagrid::make_cell( output_segment, vertex_handles.begin(), vertex_handles.end() );
        }
      }
    }
  }



  project_mesh::project_mesh() :
    input_mesh(*this, parameter_information("mesh","mesh","The input mesh, segmented line 3d mesh and segmented triangular 3d mesh supported")),
    target_dimension(*this, parameter_information("target_dimension","int","The target geometric dimension of the resulting mesh")),
    output_mesh(*this, parameter_information("mesh", "mesh", "The output mesh, same type of mesh as input mesh but the target geometric dimension")) {}

  std::string project_mesh::name() const { return "ViennaGrid Project"; }
  std::string project_mesh::id() const { return "project_mesh"; }


  template<typename InputMeshT, typename InputSegmentationT, typename OutputMeshT, typename OutputSegmentationT>
  bool project_mesh::generic_run( int target_dimension )
  {
    if (viennagrid::result_of::geometric_dimension<OutputMeshT>::value != target_dimension)
      return false;

    typedef viennagrid::segmented_mesh<InputMeshT, InputSegmentationT> InputSegmentedMeshType;
    typename viennamesh::result_of::const_parameter_handle<InputSegmentedMeshType>::type imp = input_mesh.get<InputSegmentedMeshType>();
    if (imp)
    {
      typedef viennagrid::segmented_mesh<OutputMeshT, OutputSegmentationT> OutputSegmentedMeshType;

      output_parameter_proxy<OutputSegmentedMeshType> omp(output_mesh);

      project_mesh_impl( imp().mesh, imp().segmentation, omp().mesh, omp().segmentation );
      return true;
    }

    return false;
  }


  bool project_mesh::run_impl()
  {
    if (generic_run<viennagrid::line_3d_mesh, viennagrid::line_3d_segmentation,
                    viennagrid::line_2d_mesh, viennagrid::line_2d_segmentation>( target_dimension() ))
      return true;

    if (generic_run<viennagrid::triangular_3d_mesh, viennagrid::triangular_3d_segmentation,
                    viennagrid::triangular_2d_mesh, viennagrid::triangular_2d_segmentation>( target_dimension() ))
      return true;

    error(1) << "Input Parameter 'default' (type: mesh) is missing or of non-convertable type" << std::endl;
    return false;
  }


}
