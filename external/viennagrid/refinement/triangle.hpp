#ifndef VIENNAGRID_TOPOLOGY_TRIANGLE_HPP
#define VIENNAGRID_TOPOLOGY_TRIANGLE_HPP

/* =======================================================================
   Copyright (c) 2011-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                     ViennaGrid - The Vienna Grid Library
                            -----------------

   Authors:      Karl Rupp                           rupp@iue.tuwien.ac.at
                 Josef Weinbub                    weinbub@iue.tuwien.ac.at
               
   (A list of additional contributors can be found in the PDF manual)

   License:      MIT (X11), see file LICENSE in the base directory
======================================================================= */

#include "viennagrid/forwards.hpp"
#include "viennagrid/topology/vertex.hpp"
#include "viennagrid/topology/line.hpp"
// #include "viennagrid/detail/element_iterators.hpp"
#include "viennagrid/algorithm/norm.hpp"

/** @file refinement/triangle.hpp
    @brief Provides refinement routines for a triangle
*/

namespace viennagrid
{
    
    
  template<typename CellType, typename DomainType, typename VertexHandleContainer>
  void create_refinement_cell(DomainType & domain, VertexHandleContainer vertex_handle_container, unsigned int i0, unsigned int i1, unsigned int i2)
  {
      typedef typename VertexHandleContainer::iterator VertexHandleIteratorType;
      typedef typename std::iterator_traits<VertexHandleIteratorType>::value_type VertexHandleType;
      storage::static_array< VertexHandleType, element_topology::boundary_cells<triangle_tag, vertex_tag>::num > cellvertices;
      
      cellvertices[0] = *viennagrid::advance(vertex_handle_container.begin(), i0);
      cellvertices[1] = *viennagrid::advance(vertex_handle_container.begin(), i1);
      cellvertices[2] = *viennagrid::advance(vertex_handle_container.begin(), i2);
      
      viennagrid::create_element<CellType>( domain, cellvertices.begin(), cellvertices.end() );
  }

  
  /** @brief Specialization of the refinement class for a triangle */
  template <>
  struct element_refinement<triangle_tag>
  {
    
    /** @brief No refinement. Just put same cell into new domain. */
    template <typename CellType, typename DomainTypeOut>
    static void apply0(CellType const & cell_in, DomainTypeOut & segment_out)
    {
//       typedef typename CellType::config_type        ConfigTypeIn;
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::vertex_tag>::type            VertexOnCellRange;
      typedef typename viennagrid::result_of::iterator<VertexOnCellRange>::type         VertexOnCellIterator;            
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::line_tag>::type            EdgeOnCellRange;
      typedef typename viennagrid::result_of::iterator<EdgeOnCellRange>::type           EdgeOnCellIterator;            
      
      typedef typename viennagrid::result_of::handle<DomainTypeOut, viennagrid::vertex_tag>::type             VertexHandleType;
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, vertex_tag>::type                                      VertexTypeOut;
      typedef typename VertexTypeOut::id_type VertexIDTypeOut;

      const unsigned int num_vertices = element_topology::boundary_cells<triangle_tag, vertex_tag>::num;
      storage::static_array<VertexHandleType, num_vertices> vertex_handles;
//       VertexType * vertices[topology::bndcells<triangle_tag, 0>::num];
      
      //
      // Step 1: Get vertices on the new domain
      //
      
      //grab existing vertices:
      VertexOnCellRange vertices_on_cell = viennagrid::elements<viennagrid::vertex_tag>(cell_in);
      VertexOnCellIterator vocit = vertices_on_cell.begin();
      vertex_handles[0] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[1] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[2] = viennagrid::find_by_id( segment_out, vocit->id() ).handle();

      //
      // Step 2: Add new cells to new domain:
      //
      viennagrid::create_element<CellType>( segment_out, vertex_handles.begin(), vertex_handles.end() );
      
      
//       CellType new_cell;
//       storage::static_array<VertexHandleType, element_topology::boundary_cells<tetrahedron_tag, vertex_tag>::num> cellvertices;
//       
//       //0-3-2:
//       cellvertices[0] = vertices[0];
//       cellvertices[1] = vertices[1];
//       cellvertices[2] = vertices[2];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);

    } //apply0()
    
    
    /** @brief Refinement for one edge to be bisected.
     *
     * Orientation of vertices (established by rotating the triangle appropriately)
     *
     *           2
     *         /   \ 
     *        /     \ 
     *       0 - 3 - 1
     */
    template <typename CellType, typename DomainTypeOut>
    static void apply1(CellType const & cell_in, DomainTypeOut & segment_out)
    {
//       typedef typename CellType::config_type        ConfigTypeIn;
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::vertex_tag>::type            VertexOnCellRange;
      typedef typename viennagrid::result_of::iterator<VertexOnCellRange>::type         VertexOnCellIterator;            
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::line_tag>::type            EdgeOnCellRange;
      typedef typename viennagrid::result_of::iterator<EdgeOnCellRange>::type           EdgeOnCellIterator;            
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, viennagrid::vertex_tag>::type             VertexType;
      typedef typename viennagrid::result_of::handle<DomainTypeOut, viennagrid::vertex_tag>::type             VertexHandleType;
      typedef typename viennagrid::result_of::element<DomainTypeOut, viennagrid::line_tag>::type             EdgeType;
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, vertex_tag>::type                                      VertexTypeOut;
      typedef typename VertexTypeOut::id_type VertexIDTypeOut;


      const unsigned int num_vertices = element_topology::boundary_cells<triangle_tag, vertex_tag>::num;
      storage::static_array<VertexHandleType, num_vertices+1> vertex_handles;
      
      
      //
      // Step 1: Get vertices on the new domain
      //
      
      //grab existing vertices:
      VertexOnCellRange vertices_on_cell = viennagrid::elements<viennagrid::vertex_tag>(cell_in);
      VertexOnCellIterator vocit = vertices_on_cell.begin();
      vertex_handles[0] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[1] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[2] = viennagrid::find_by_id( segment_out, vocit->id() ).handle();

      //add vertices from edge
      EdgeOnCellRange edges_on_cell = viennagrid::elements<viennagrid::line_tag>(cell_in);
      std::size_t offset = 0;
      EdgeOnCellIterator eocit = edges_on_cell.begin();
      EdgeType const & e0 = *eocit; ++eocit;
      EdgeType const & e1 = *eocit; ++eocit;
      EdgeType const & e2 = *eocit;
      
      if ( (viennadata::access<refinement_key, bool>(refinement_key())(e0)) == true )
      {
        vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e0) ).handle();
        offset = 0;
      }
      else if ( (viennadata::access<refinement_key, bool>(refinement_key())(e1)) == true )
      {
        vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e1) ).handle();
        offset = 2;
      }
      else if ( (viennadata::access<refinement_key, bool>(refinement_key())(e2)) == true )
      {
        vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e2) ).handle();
        offset = 1;
      }
      else
      {
        assert( (2 < 2) && "Logic error: Triangle does not have an edges for bisection!");
      }
      
      //
      // Step 2: Add new cells to new domain:
      //
      
      create_refinement_cell<CellType>( segment_out, vertex_handles, (offset + 0) % num_vertices, 3, (offset + 2) % num_vertices );
      
//       CellType new_cell;
//       VertexType * cellvertices[topology::bndcells<triangle_tag, 0>::num];
//       
//       //0-3-2:
//       cellvertices[0] = vertices[(offset + 0) % topology::bndcells<triangle_tag, 0>::num];
//       cellvertices[1] = vertices[3];
//       cellvertices[2] = vertices[(offset + 2) % topology::bndcells<triangle_tag, 0>::num];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);

      
      create_refinement_cell<CellType>( segment_out, vertex_handles, 3, (offset + 1) % num_vertices, (offset + 2) % num_vertices );
      //3-1-2:
//       cellvertices[0] = vertices[3];
//       cellvertices[1] = vertices[(offset + 1) % topology::bndcells<triangle_tag, 0>::num];
//       cellvertices[2] = vertices[(offset + 2) % topology::bndcells<triangle_tag, 0>::num];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);
    } //apply1()
    

    /** @brief Refinement for one edge to be bisected.
     *
     * Orientation of vertices:  (established by rotating the triangle appropriately)
     *
     *            2
     *          /   \
     *         /     4
     *        /       \
     *       0 -- 3 -- 1
    */
    template <typename CellType, typename DomainTypeOut>
    static void apply2(CellType const & cell_in, DomainTypeOut & segment_out)
    {
//       typedef typename CellType::config_type        ConfigTypeIn;
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::vertex_tag>::type            VertexOnCellRange;
      typedef typename viennagrid::result_of::iterator<VertexOnCellRange>::type         VertexOnCellIterator;            
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::line_tag>::type            EdgeOnCellRange;
      typedef typename viennagrid::result_of::iterator<EdgeOnCellRange>::type           EdgeOnCellIterator;            
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, viennagrid::vertex_tag>::type             VertexType;
      typedef typename viennagrid::result_of::handle<DomainTypeOut, viennagrid::vertex_tag>::type             VertexHandleType;
      typedef typename viennagrid::result_of::element<DomainTypeOut, viennagrid::line_tag>::type             EdgeType;
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, vertex_tag>::type                                      VertexTypeOut;
      typedef typename VertexTypeOut::id_type VertexIDTypeOut;


      const unsigned int num_vertices = element_topology::boundary_cells<triangle_tag, vertex_tag>::num;
      storage::static_array<VertexHandleType, num_vertices+2> vertex_handles;

      
//       VertexType * vertices[topology::bndcells<triangle_tag, 0>::num + 2];
      
      //
      // Step 1: Get vertices on the new domain
      //
      
      //grab existing vertices:
      VertexOnCellRange vertices_on_cell = viennagrid::elements<viennagrid::vertex_tag>(cell_in);
      VertexOnCellIterator vocit = vertices_on_cell.begin();
      vertex_handles[0] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[1] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[2] = viennagrid::find_by_id( segment_out, vocit->id() ).handle();

      //Find rotation offset such that first two edges are to be refined
      EdgeOnCellRange edges_on_cell = viennagrid::elements<viennagrid::line_tag>(cell_in);
      std::size_t offset = 0;
      
      EdgeOnCellIterator eocit = edges_on_cell.begin();
      EdgeType const & e0 = *eocit; ++eocit;
      EdgeType const & e1 = *eocit; ++eocit;
      EdgeType const & e2 = *eocit;
      
      if ( (viennadata::access<refinement_key, bool>(refinement_key())(e0) == true)
           && (viennadata::access<refinement_key, bool>(refinement_key())(e1) == true) )
      {
          vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e1) ).handle();
          vertex_handles[4] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e0) ).handle();
          
//         vertices[3] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e1)]);
//         vertices[4] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e0)]);
        offset = 2;
      }
      else if ( (viennadata::access<refinement_key, bool>(refinement_key())(e0) == true)
           && (viennadata::access<refinement_key, bool>(refinement_key())(e2) == true) )
      {
          vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e0) ).handle();
          vertex_handles[4] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e2) ).handle();
          
//         vertices[3] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e0)]);
//         vertices[4] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e2)]);
        offset = 0;
      }
      else if ( (viennadata::access<refinement_key, bool>(refinement_key())(e1) == true)
           && (viennadata::access<refinement_key, bool>(refinement_key())(e2) == true) )
      {
          vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e2) ).handle();
          vertex_handles[4] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(e1) ).handle();
          
//         vertices[3] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e2)]);
//         vertices[4] = &(viennagrid::ncells<0>(segment_out.domain())[viennadata::access<refinement_key, std::size_t>(refinement_key())(e1)]);
        offset = 1;
      }
      else
      {
        assert( (2 < 2) && "Logic error: Triangle does not have two edges for bisection!");
      }
      
      //
      // Step 2: Add new cells to new domain:
      //
//       CellType new_cell;
//       VertexType * cellvertices[topology::bndcells<triangle_tag, 0>::num];
      
      //3-1-4:
//       cellvertices[0] = vertices[3];
//       cellvertices[1] = vertices[(offset + 1) % topology::bndcells<triangle_tag, 0>::num];
//       cellvertices[2] = vertices[4];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);

        create_refinement_cell<CellType>( segment_out, vertex_handles, 3, (offset + 1) % num_vertices, 4 );

      //split second-longest edge
      VertexHandleType vh0 = vertex_handles[(offset + 0) % num_vertices];
      VertexHandleType vh1 = vertex_handles[(offset + 1) % num_vertices];
      VertexHandleType vh2 = vertex_handles[(offset + 2) % num_vertices];
      double len_edge1 = viennagrid::norm( viennagrid::point(segment_out, vh1) - viennagrid::point(segment_out, vh0) );
      double len_edge2 = viennagrid::norm( viennagrid::point(segment_out, vh2) - viennagrid::point(segment_out, vh1) );
      
      if (len_edge1 > len_edge2) //split edge [v0, v1] again
      {
//         //0-3-2:
//         cellvertices[0] = vertices[(offset + 0) % topology::bndcells<triangle_tag, 0>::num];
//         cellvertices[1] = vertices[3];
//         cellvertices[2] = vertices[(offset + 2) % topology::bndcells<triangle_tag, 0>::num];
//         new_cell.vertices(cellvertices);
//         segment_out.push_back(new_cell);
        
        create_refinement_cell<CellType>( segment_out, vertex_handles, (offset + 0) % num_vertices, 3, (offset + 2) % num_vertices );

//         //2-3-4:
//         cellvertices[0] = vertices[(offset + 2) % topology::bndcells<triangle_tag, 0>::num];
//         cellvertices[1] = vertices[3];
//         cellvertices[2] = vertices[4];
//         new_cell.vertices(cellvertices);
//         segment_out.push_back(new_cell);
        
        create_refinement_cell<CellType>( segment_out, vertex_handles, (offset + 2) % num_vertices, 3, 4 );
      }
      else //split edge [v1, v2]
      {
//         //0-3-4:
//         cellvertices[0] = vertices[(offset + 0) % topology::bndcells<triangle_tag, 0>::num];
//         cellvertices[1] = vertices[3];
//         cellvertices[2] = vertices[4];
//         new_cell.vertices(cellvertices);
//         segment_out.push_back(new_cell);
        
        create_refinement_cell<CellType>( segment_out, vertex_handles, (offset + 0) % num_vertices, 3, 4 );

//         //0-4-2:
//         cellvertices[0] = vertices[(offset + 0) % topology::bndcells<triangle_tag, 0>::num];
//         cellvertices[1] = vertices[4];
//         cellvertices[2] = vertices[(offset + 2) % topology::bndcells<triangle_tag, 0>::num];
//         new_cell.vertices(cellvertices);
//         segment_out.push_back(new_cell);
        
        create_refinement_cell<CellType>( segment_out, vertex_handles, (offset + 0) % num_vertices, 4, (offset + 2) % num_vertices );
      }
      
      
    } //apply2()
    


    
    /** @brief Refinement of a triangle with three edges to be refined (uniform refinement) */
    template <typename CellType, typename DomainTypeOut>
    static void apply3(CellType const & cell_in, DomainTypeOut & segment_out)
    {
//       typedef typename CellType::config_type        ConfigTypeIn;
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::vertex_tag>::type            VertexOnCellRange;
      typedef typename viennagrid::result_of::iterator<VertexOnCellRange>::type         VertexOnCellIterator;            
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::line_tag>::type            EdgeOnCellRange;
      typedef typename viennagrid::result_of::iterator<EdgeOnCellRange>::type           EdgeOnCellIterator;            
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, viennagrid::vertex_tag>::type             VertexType;
      typedef typename viennagrid::result_of::handle<DomainTypeOut, viennagrid::vertex_tag>::type             VertexHandleType;
      
      typedef typename viennagrid::result_of::element<DomainTypeOut, vertex_tag>::type                                      VertexTypeOut;
      typedef typename VertexTypeOut::id_type VertexIDTypeOut;

      const unsigned int num_vertices = element_topology::boundary_cells<triangle_tag, vertex_tag>::num;
      const unsigned int num_lines = element_topology::boundary_cells<triangle_tag, line_tag>::num;
      
      storage::static_array<VertexHandleType, num_vertices+num_lines> vertex_handles;
      
//       VertexType * vertices[topology::bndcells<triangle_tag, 0>::num
//                             + topology::bndcells<triangle_tag, 1>::num];
      
      //
      // Step 1: Get vertices on the new domain
      //
      
      //grab existing vertices:
      VertexOnCellRange vertices_on_cell = viennagrid::elements<viennagrid::vertex_tag>(cell_in);
      VertexOnCellIterator vocit = vertices_on_cell.begin();
      vertex_handles[0] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[1] = viennagrid::find_by_id( segment_out, vocit->id() ).handle(); ++vocit;
      vertex_handles[2] = viennagrid::find_by_id( segment_out, vocit->id() ).handle();

      //add vertices from edge
      EdgeOnCellRange edges_on_cell = viennagrid::elements<viennagrid::line_tag>(cell_in);
      EdgeOnCellIterator eocit = edges_on_cell.begin();
      vertex_handles[3] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(*eocit) ).handle(); ++eocit;
      vertex_handles[4] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(*eocit) ).handle(); ++eocit;
      vertex_handles[5] = viennagrid::find_by_id( segment_out, viennadata::access<refinement_key, VertexIDTypeOut>(refinement_key())(*eocit) ).handle();
      
      //
      // Step 2: Add new cells to new domain:
      //
//       CellType new_cell;
//       VertexType * cellvertices[topology::bndcells<triangle_tag, 0>::num];
//       
//       //0-3-4:
//       cellvertices[0] = vertices[0];
//       cellvertices[1] = vertices[3];
//       cellvertices[2] = vertices[4];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);
      
      create_refinement_cell<CellType>( segment_out, vertex_handles, 0, 3, 4 );

//       //3-1-5:
//       cellvertices[0] = vertices[3];
//       cellvertices[1] = vertices[1];
//       cellvertices[2] = vertices[5];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);
      
      create_refinement_cell<CellType>( segment_out, vertex_handles, 3, 1, 5 );

//       //4-5-2:
//       cellvertices[0] = vertices[4];
//       cellvertices[1] = vertices[5];
//       cellvertices[2] = vertices[2];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);
      
      create_refinement_cell<CellType>( segment_out, vertex_handles, 4, 5, 2 );

//       //4-3-5:
//       cellvertices[0] = vertices[4];
//       cellvertices[1] = vertices[3];
//       cellvertices[2] = vertices[5];
//       new_cell.vertices(cellvertices);
//       segment_out.push_back(new_cell);
      
      create_refinement_cell<CellType>( segment_out, vertex_handles, 4, 3, 5 );
      
    } //apply3()


    /** @brief Public entry function for the refinement of a triangle.
     * 
     * @param cell_in       The triangle to be refined
     * @param segment_out   The domain or segment the refined triangles are written to
     */
    template <typename CellType, typename DomainTypeOut>
    static void apply(CellType const & cell_in, DomainTypeOut & segment_out)
    {
      typedef typename viennagrid::result_of::const_element_range<CellType, viennagrid::line_tag>::type            EdgeOnCellRange;
      typedef typename viennagrid::result_of::iterator<EdgeOnCellRange>::type                 EdgeOnCellIterator;            
      
      std::size_t edges_to_refine = 0;
      EdgeOnCellRange edges_on_cell = viennagrid::elements<viennagrid::line_tag>(cell_in);
      for (EdgeOnCellIterator eocit = edges_on_cell.begin();
                              eocit != edges_on_cell.end();
                            ++eocit)
      {
        if (viennadata::access<refinement_key, bool>(refinement_key())(*eocit) == true)
          ++edges_to_refine;
      }
      
      switch (edges_to_refine)
      {
        case 0: apply0(cell_in, segment_out); break;
        case 1: apply1(cell_in, segment_out); break;
        case 2: apply2(cell_in, segment_out); break;
        case 3: apply3(cell_in, segment_out); break;
        default: //nothing to do...
                break;
      }
    } //apply()
  };
  
  
  
  
  
  
    
}

#endif

