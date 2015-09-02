#ifndef VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_WRITER_HPP
#define VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_WRITER_HPP

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

#include <string>
#include <vector>
#include <boost/container/flat_map.hpp>

#include "H5Cpp.h"

#include "viennautils/exception.hpp"

#include "viennagridpp/mesh/mesh.hpp"
#include "viennagridpp/quantity_field.hpp"
#include "viennagridpp/range.hpp"

namespace viennamesh
{

struct tdr_writer_error : virtual viennautils::exception {};

class sentaurus_tdr_writer {
  typedef viennagrid::const_mesh                                          MeshType;
  typedef viennagrid::result_of::element<MeshType>::type                  ElementType;
  typedef viennagrid::result_of::point<MeshType>::type                    PointType;
  typedef viennagrid::result_of::region<MeshType>::type                   RegionType;
  typedef viennagrid::result_of::id<ElementType>::type                    ElementId;
  typedef viennagrid::result_of::id<RegionType>::type                     RegionId;

  typedef viennagrid::result_of::const_vertex_range<MeshType>::type       MeshVertexRange;
  typedef viennagrid::result_of::iterator<MeshVertexRange>::type          MeshVertexIterator;

  typedef viennagrid::result_of::const_vertex_range<ElementType>::type    BoundaryVertexRange;
  typedef viennagrid::result_of::iterator<BoundaryVertexRange>::type      BoundaryVertexIterator;

  typedef viennagrid::result_of::const_region_range<MeshType>::type       RegionRange;
  typedef viennagrid::result_of::iterator<RegionRange>::type              RegionIterator;

  typedef viennagrid::result_of::const_vertex_range<RegionType>::type     RegionVertexRange;
  typedef viennagrid::result_of::iterator<RegionVertexRange>::type        RegionVertexIterator;

  typedef viennagrid::result_of::const_cell_range<RegionType>::type       CellRange;
  typedef viennagrid::result_of::iterator<CellRange>::type                CellIterator;
  
  MeshType const & mesh;
  std::vector<viennagrid::quantity_field> const & quantities;
  
  H5::H5File file;
  H5::Group geometry;
  unsigned int dimension;
  boost::container::flat_map<ElementId, int32_t> vertex_mapping;
  
  void create_geometry();
  void write_identity_transformation();
  void write_vertices();
  void write_regions();
  void write_quantity_fields();
  
  
public:
  sentaurus_tdr_writer(std::string const & filename, MeshType const & mesh, std::vector<viennagrid::quantity_field> const & quantities)
    : mesh(mesh), quantities(quantities), file(filename, H5F_ACC_TRUNC) {}
  
  void write_to_tdr();
};

void write_to_tdr(std::string const & filename, viennagrid::const_mesh const & mesh, std::vector<viennagrid::quantity_field> const & quantities);

} //end of namespace viennamesh

#endif
