#ifndef VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_READER_HPP
#define VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_READER_HPP

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
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>

#include "viennagridpp/mesh/element_creation.hpp"
#include "viennagridpp/quantity_field.hpp"
#include "viennagridpp/algorithm/inclusion.hpp"
#include "viennagridpp/algorithm/distance.hpp"
#include "viennagridpp/point_accessor.hpp"

#include "viennameshpp/logger.hpp"

#include "H5Cpp.h"

using std::string;
using namespace H5;

#define mythrow(a) { std::ostringstream s; s << a; throw std::runtime_error(s.str()); }


namespace viennamesh
{
  int read_int(const H5Object &g, const string name);
  double read_double(const H5Object &g, const string name);
  string read_string(const H5Object &g, const string name);
  
  // get the number of values in the given dataset
  hsize_t get_dataset_size(const DataSet& dataset);
  
  // return a vector of all values in the given dataset
  std::vector<double> read_vector(const DataSet& dataset);

  
  template <typename MeshT>
  struct tdr_geometry;
  class tdr_region;
  
  /*
   * This class is responsible for the construction of elements from the
   * tdr data and adding these elements to the mesh/regions.
   */
  class tdr_element {
    tdr_region* region;
    viennagrid::element_tag element_type;
    std::vector<int> vertices;
    
  public:
    /*
     * Constructs the element from raw data read from the tdr file.
     * The parameter begin points to the first value of that element
     * in the element list (the element type) and is modified so that
     * it points to the next element after this method finishes. The
     * paramter end should point to the end of the element list (not the element)
     * and is just used to prevent reading beyond the end of the vector.
     */
    tdr_element(tdr_region* region, int* & begin, int* end);
    
    /*
     * Constructs viennagrid element and adss it to the mesh
     */
    template<typename MeshT, typename RegionT>
    void add_to_mesh(tdr_geometry<MeshT> * geometry, RegionT & reg);
    
    /*
     * Extrudes this element (line -> triangle, triangle -> tetrahedron).
     * This method is applied to all elements in interface regions if
     * "extrude_contacts" is set.
     */
    template <typename MeshT>
    void extrude(tdr_geometry<MeshT> * geometry);
    
  private:
    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1);

    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1, PointT const & p2);
  };
  
  class tdr_region {
    friend class tdr_element;
    typedef std::vector<tdr_element>::iterator element_iter;
    
    int id;
    string name;
    std::vector<tdr_element> elements;
    
    
    /*
     * A region can have one of the following types:
     *  0: "normal" region, contains elements
     *  1: outside interface region, at the edge of 1 other region, contains only lines
     *  2: inside interface region, between 2 other regions, contains only lines
     */
    int type;
    
    /*
     * contains all vertices in this regions. Needed for quantity fields
     */
    std::set<int> vertices;
    
    /*
     * contains the ids of all elements in this region. Needed for quantity fields.
     */
    std::vector<int> element_indices;
    
    /*
     * contains all vertices that were newly created when extruding elements. Needed for quantity fields.
     */
    std::vector<int> extruded_vertices;
    
  public:
    /*
     * Reads the region including all contained elements from the tdr file.
     * 
     * The main reason this is not done in the constructor is because then the
     * region pointer in tdr_element would be invalidated when the region object
     * is moved/copied into the vector where it is stored.
     */
    void read(int regnr, const Group &reg)
    {
      id = regnr;
      name = read_string(reg,"name");
      type = read_int(reg,"type");

      if(type > 2)
        mythrow("Unknown region type " << type << " found in " << reg.getObjName());
      
      //iterate over all elements_* and read the contained elements
      const int n=reg.getNumObjs();
      for (int i=0; i<n; i++)
      {
        std::string objname = reg.getObjnameByIdx(i);
        if (objname.find("elements_") == 0)
        {
          read_elements(reg.openDataSet(objname));
        }

        if (objname.find("part_") == 0)
        {
          read_elements(reg.openGroup(objname).openDataSet("elements"));
        }
      }
    }
    
    /*
     * Creates the viennagrid region and ads it to the mesh.
     * Then adds all child elements to the mesh.
     */
    template <typename MeshT>
    void add_to_mesh(tdr_geometry<MeshT> * geometry)
    {
      typedef typename tdr_geometry<MeshT>::RegionType RegionType;
      
      //TODO: Handle type 2 regions?
      if(type == 2)
      {
        warning(0) << "Type 2 regions (interfaces) are ignored!" << std::endl;
        return;
      }
      
      RegionType reg = geometry->mesh.get_or_create_region(id);
      reg.set_name(name);
      
      for(element_iter it = elements.begin(); it != elements.end(); ++it)
      {
        tdr_element & e = *it;
        e.add_to_mesh(geometry, reg);
      }
    }
    
    /*
     * Redirects to the private overload of this method with the correct target paramter.
     * 
     * location_type 0 = vertices
     * location_type 3 = elements
     */
    void apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, int location_type)
    {
      switch(location_type) {
        case 0:
          apply_quantity(quantity, values, vertices);
          break;
        case 3:
          apply_quantity(quantity, values, element_indices);
          break;
        default:
          mythrow("Invalid location type in quantity field " << quantity.get_name());
      }
    }
    
    /*
     * When using "extrude_contacts" new elements and vertices are created for which
     * no values exist in the tdr file. This method takes the average of the neighbor
     * verties or element to fill in the missing values.
     */
    template <typename MeshT>
    void extrude_quantities(tdr_geometry<MeshT> * geometry)
    {
      typedef typename tdr_geometry<MeshT>::VertexType VertexType;
      typedef typename tdr_geometry<MeshT>::RegionType RegionType;
      typedef typename viennagrid::result_of::neighbor_range<MeshT>::type NeighborRangeType;
      typedef typename viennagrid::result_of::cell_range<RegionType>::type CellRangeType;
      typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;
      
      //only regions of type 1 contain extruded elements
      //TODO: Handle type 2 regions?
      if(type == 1)
      {
	// iterate over all quantity fields
        std::map<string, viennagrid::quantity_field>::iterator qit = geometry->quantities.begin();
        for(; qit != geometry->quantities.end(); ++qit)
        {
          viennagrid::quantity_field & quantity = qit->second;
          int dimension = quantity.topologic_dimension();
          
          if(dimension == 0)
          {
	    //vertex quantity -> iterate over the extruded vertices
            std::vector<int>::const_iterator it = extruded_vertices.begin();
            for(; it != extruded_vertices.end(); ++it)
            {
	      //set neighbor average as quantity value
              VertexType v = geometry->vertices[*it];
              NeighborRangeType neighbor_vertices(geometry->mesh, v, 1, 0);
              quantity.set(v, get_neighbor_average(quantity, neighbor_vertices));
            }
          }
          else if(dimension == geometry->dim)
          {
	    //element quantity -> iterate over the extruded elements
	    //Note: type 1 regions contain only extruded elements
            RegionType reg = geometry->mesh.get_or_create_region(id);
            CellRangeType cells(reg);
            
            CellIteratorType cell = cells.begin();
            for(; cell != cells.end(); ++cell)
            {
	      //set neighbor average as quantity value
              NeighborRangeType neighbor_cells(geometry->mesh, *cell, dimension-1, dimension);
              quantity.set(*cell, get_neighbor_average(quantity, neighbor_cells));
            }
          }
          else
          {
            mythrow("connot extrude quantity field of topological dimension " << dimension);
          }
        }
      }
    }
    
    int get_type() {return type;}
    
  private:
    /*
     * Read all elements contained in the given dataset
     */
    void read_elements(const DataSet &elem)
    {
      hsize_t size = get_dataset_size(elem);

      boost::scoped_array<int> element_data(new int[size]);      
      elem.read( element_data.get(), PredType::NATIVE_INT);

      int* el = element_data.get();
      int* el_end = el + size;
      
      while(el != el_end)
        elements.push_back(tdr_element(this, el, el_end));
    }
    
    /*
     * For every given vertex/element given in "target" it sets the given value in the quantity_field.
     */
    template <typename ContainerT>
    void apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, const ContainerT & target)
    {
      if(values.size() != target.size())
        mythrow("Quantity field " << quantity.get_name() << " in region " << id << " has invalid size");
      
      int i = 0;
      typename ContainerT::const_iterator t = target.begin();
      for(; t != target.end(); ++t)
      {
        quantity.set(*t, values[i++]);
      }
    }
    
    template <typename NeighborRangeType>
    inline double get_neighbor_average(const viennagrid::quantity_field & quantity, NeighborRangeType &neighbors)
    {
      typedef typename viennagrid::result_of::iterator<NeighborRangeType>::type NeighborIteratorType;
      double value = 0;
      int i = 0;
      
      NeighborIteratorType neighbor = neighbors.begin();
      for(; neighbor != neighbors.end(); ++neighbor, ++i)
      {
        value += quantity.get(*neighbor);
      }
      return value / (double)i;
    }
  };
  
  /*
   * This class is the main class of the tdr reader.
   */
  template<typename MeshT>
  struct tdr_geometry
  {
    typedef typename viennagrid::result_of::point<MeshT>::type PointType;
    typedef typename viennagrid::result_of::element<MeshT>::type VertexType;
    typedef typename viennagrid::result_of::element<MeshT>::type ElementType;
    typedef typename viennagrid::result_of::region<MeshT>::type RegionType;
    
    typedef std::vector<tdr_region>::iterator region_iterator;
    
    MeshT const & mesh;
    std::vector<VertexType> vertices;
    std::vector<tdr_region> regions;
    std::map<string, viennagrid::quantity_field> quantities;
    
    unsigned int nvertices;
    int dim,nregions,ndatasets;
    
    bool extrude_contacts;
    double extrude_contacts_scale;
    bool fill_triangle_contacts;
    
    tdr_geometry(MeshT const & mesh) : mesh(mesh) {}

    /*
     * Reads the transformation Matrix+Vector from the tdr file.
     * 
     * At the moment it only throws an error if the read transformation
     * is not equal to the identity transformation.
     */
    void read_transformation(const Group &trans)
    {
      boost::array<double, 9> trans_matrix, identity_matrix;
      boost::array<double, 3> trans_vector, zero_vector;
      
      // initialize identity tranformation
      zero_vector.fill(0);
      identity_matrix.fill(0);
      identity_matrix[0] = identity_matrix[4] = identity_matrix[8] = 1;
      
      // read transformation from file
      const DataSet &A=trans.openDataSet("A");
      const DataSet &b=trans.openDataSet("b");
      A.read( trans_matrix.elems, PredType::NATIVE_DOUBLE);
      b.read( trans_vector.elems, PredType::NATIVE_DOUBLE);
      
      if(trans_matrix != identity_matrix || trans_vector != zero_vector)
        mythrow("transformation not equal to identity");
    }

    void read_vertices(const DataSet &vert)
    {
      hsize_t size = get_dataset_size(vert);
      
      if (nvertices != size)
        mythrow("nvertices not equal vertices.dim");

      CompType mtype2( sizeof(double)*dim );
      
      mtype2.insertMember( "x", 0, PredType::NATIVE_DOUBLE);

      if (dim>1)
        mtype2.insertMember( "y", sizeof(double), PredType::NATIVE_DOUBLE);

      if (dim>2)
        mtype2.insertMember( "z", sizeof(double)*2, PredType::NATIVE_DOUBLE);

      boost::scoped_array<double> s2(new double[size*dim]);
      vert.read( s2.get(), mtype2 );
      
      for (unsigned int i=0; i < size*dim; i += dim)
      {
        PointType point(dim);
        
        point[0] = s2[i];
        if (dim>1)
        point[1] = s2[i+1];
        if (dim>2)
        point[2] = s2[i+2];

        vertices.push_back(viennagrid::make_vertex(mesh, point));
      }
    }
    
    void read_regions(const Group & group)
    {
      nregions = read_int(group, "number of regions");
      regions.reserve(nregions);
      
      for (int i=0; i<nregions; i++)
      {
        const Group &reg=group.openGroup("region_" + boost::lexical_cast<string>(i));
        regions.push_back(tdr_region());
        regions.back().read(i, reg);
      }
    }

    void read_dataset(const Group &dataset)
    {
      string name = read_string(dataset,"name");
      int regnr = read_int(dataset,"region");
      int location_type = read_int(dataset,"location type");
      string unit=read_string(dataset,"unit:name");
      double conversion_factor = read_double(dataset, "conversion factor");
      
      if (location_type != 0 && location_type != 3)
        mythrow("Location type " << location_type << " in dataset " << name << " not supported");

      if (read_int(dataset,"structure type") != 0)
        mythrow("Dataset " << name << " structure type not 0");

      if (read_int(dataset,"value type") != 2)
        mythrow("Dataset " << name << " value type not 2");
      
      tdr_region & region = regions[regnr];
      
      //TODO:: Handle type 2 regions?
      if(region.get_type() == 2)
        return;
      
      viennagrid::quantity_field & quantity = quantities[name];
      
      //check if quantity field is already initialized
      if(quantity.get_name() != name)
      {
	//and initialize it if not
        if(location_type == 0)
          quantity.init(0, 1);
        else
          quantity.init(dim, 1);
        
        quantity.set_name(name);
        quantity.set_unit(unit);
      }
      
      if(quantity.unit() != unit)
        mythrow("Units not the same in all regions in quantity field " << name << " - " << dataset.getObjName());
      
      std::vector<double> values = read_vector(dataset.openDataSet("values"));
      
      if(conversion_factor != 1.0)
      {
	warning(0) << "Values multiplied with conversion factor " << conversion_factor << " in " << dataset.getObjName() << std::endl;
	
	for(unsigned int i = 0; i < values.size(); ++i)
	{
	  values[i] *= conversion_factor;
	}
      }
      
      region.apply_quantity(quantity, values, location_type);
    }

    void read_datasets(const Group &state)
    {
      if (read_int(state,"number of plots") != 0)
        mythrow("'Number of plots' not equal 0");

      if (read_int(state,"number of string streams") != 0)
        mythrow("'Number of string streams' not equal 0");

      if (read_int(state,"number of xy plots") != 0)
        mythrow("'Number of xy plots' not equal 0");

      ndatasets=read_int(state,"number of datasets");
      
      for (int i=0; i<ndatasets; i++)
      {
        read_dataset(state.openGroup("dataset_" + boost::lexical_cast<string>(i)));
      }
    }

    void read_geometry(const Group &geometry)
    {
      int type = read_int(geometry,"type");
      int dum = read_int(geometry,"number of states");
      dim = read_int(geometry,"dimension");
      nvertices = read_int(geometry,"number of vertices");
      
      /*
       * type 0 = boundary files, contains only polygons and lines and no datasets
       * type 1 = mesh files
       */
      if(type != 0 && type != 1)
        mythrow("Unknown geometry type " << type << " in " << geometry.getObjName());

      if (type == 1 && dum != 1)
        mythrow("Number of states not one");
      
      if(type == 0 && extrude_contacts)
      {
	//extrude contacts does not work with polygons, because is_inside is not implemented for polygons in viennagrid
        warning(0) << "extrude_contacts not available with geometry type 0" << std::endl;
        extrude_contacts = false;
      }

      read_transformation(geometry.openGroup("transformation"));
      read_vertices(geometry.openDataSet("vertex"));
      read_regions(geometry);
      
      for(region_iterator it = regions.begin(); it != regions.end(); ++it)
      {
        it->add_to_mesh(this);
      }
      
      //only type 1 geometries contain datasets
      if(type == 1)
      {
        read_datasets(geometry.openGroup("state_0"));
        
        for(region_iterator it = regions.begin(); it != regions.end(); ++it)
        {
          it->extrude_quantities(this);
        }
      }
    }

    void read_file(const string &filename)
    {
      boost::shared_ptr<H5File> file( new H5File(filename.c_str(), H5F_ACC_RDWR) );

      if (file->getNumObjs()!=1)
        mythrow("File has not exactly one collection (number of collections = " << file->getNumObjs());

      const Group &collection = file->openGroup("collection");
      
      if (read_int(collection,"number of geometries") != 1)
        mythrow("Not only one geometry");
      
      if (read_int(collection,"number of plots") != 0)
        viennamesh::warning(1) << "file contains plots, skip them" << std::endl;
      
      read_geometry(collection.openGroup("geometry_0"));
    }


//     template<typename PointT>
//     PointT normal_vector(PointT const & p0, PointT const & p1)
//     {
//       PointT line = p1-p0;
//       std::swap(line[0], line[1]);
//       line[0] = -line[0];
//       return line;
//     }
// 
//     template<typename PointT>
//     PointT normal_vector(PointT const & p0, PointT const & p1, PointT const & p2)
//     {
//       return viennagrid::cross_prod( p1-p0, p2-p0 );
//     }
    
    /*
     * This function increases the topological dimension of the given element.
     * It is used on contact elements because they have a topological dimension
     * which is less than the cell dimension
     */
//     void extrude_element(element_t & element)
//     {
//       PointType center;
//       PointType normal;
//       double size = 0;
//       
//       /*normal = viennagrid::normal_vector(viennagrid::mesh_point_accessor(mesh), element);
//       center = viennagrid::centroid(viennagrid::mesh_point_accessor(mesh), element);
//       
//       viennagrid_element_type extruded_type;
//       if(element.is_line())
//       {
//         extruded_type = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
//         size = viennagrid::distance(viennagrid::vertices(element)[0], viennagrid::vertices(element)[1]) / 2.0;
//       }
//       else if(element.is_triangle())
//       {
//         extruded_type = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
//         
//         size = viennagrid::volume(viennagrid::mesh_point_accessor(mesh), element);
// //         size = std::max(viennagrid::distance(center, viennagrid::vertices(element)[0]), viennagrid::distance(center, viennagrid::vertices(element)[1]));
// //         size = std::max(size, viennagrid::distance(center, viennagrid::vertices(element)[2]));
//       }
//       else
//       {
//         mythrow("extruding elements of type " << viennagrid_element_type_string(element.element_type) << " is not supported");
//       }*/
//       
//       viennagrid_element_type extruded_type = VIENNAGRID_ELEMENT_TYPE_NO_ELEMENT;
//       if (element.element_type == VIENNAGRID_ELEMENT_TYPE_LINE)
//       {
//         PointType p0 = viennagrid::get_point( vertices[element.vertices[0]] );
//         PointType p1 = viennagrid::get_point( vertices[element.vertices[1]] );
// 
//         center = (p0+p1)/2;
//         normal = normal_vector(p0, p1);
// 
//         size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
//         extruded_type = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
//       }
//       else if (element.element_type == VIENNAGRID_ELEMENT_TYPE_TRIANGLE)
//       {
//         PointType p0 = viennagrid::get_point( vertices[element.vertices[0]] );
//         PointType p1 = viennagrid::get_point( vertices[element.vertices[1]] );
//         PointType p2 = viennagrid::get_point( vertices[element.vertices[2]] );
// 
//         center = (p0+p1+p2)/3;
//         normal = normal_vector(p0, p1, p2);
// 
//         size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
//         size = std::max(size, viennagrid::distance(center, p2));
//         
//         extruded_type = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
//       }
//       else
//       {
//         mythrow("extruding elements of type " << viennagrid_element_type_string(element.element_type) << " is not supported");
//       }
// 
//       normal.normalize();
//       normal *= size * extrude_contacts_scale;
// 
//       PointType new_vertex = center + normal;
// 
//       typedef typename viennagrid::result_of::cell_range<MeshT>::type CellRangeType;
//       typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;
// 
//       CellRangeType cells(mesh);
//       CellIteratorType cit = cells.begin();
//       for (; cit != cells.end(); ++cit)
//       {
//         if ( viennagrid::is_inside(*cit, new_vertex) )
//           break;
//       }
// 
//       if (cit != cells.end())
//         new_vertex = center - normal;
// 
//       VertexType v = viennagrid::make_vertex(mesh, new_vertex);
//       
//       vertices.push_back(v);
//       newly_created_vertices[v.id()] = element.vertices;
//       
//       /*ElementType e = viennagrid::make_element( region, viennagrid::element_tag::from_internal(element.element_type),
//                             cell_vertices.begin(), cell_vertices.end() );*/
//       
//       element.vertices.push_back(v.id());
//       element.element_type = extruded_type;
//     }
    
//     void fill_triangle(ElementType const & element, RegionType const & region)
//     {
//       if(!element.is_triangle())
//         mythrow("fill_triangle_contacts works only with triangles");
//       
//       typedef typename viennagrid::result_of::neighbor_range<RegionType>::type NeighborElementRangeType;
//       typedef typename viennagrid::result_of::iterator<NeighborElementRangeType>::type NeighborElementRangeIterator;
//       
//       NeighborElementRangeType neighbor_triangles(region, element, 0, 2);
//       for (NeighborElementRangeIterator ntit = neighbor_triangles.begin(); ntit != neighbor_triangles.end(); ++ntit)
//       {
//         boost::array<VertexType, 3> triangle;
//         
//         // find shared vertex
//         for (int i = 0; i != 3; ++i)
//           for (int j = 0; j != 3; ++j)
//           {
//             if ( viennagrid::vertices(element)[i] == viennagrid::vertices(*ntit)[j] )
//               triangle[0] = viennagrid::vertices(element)[i];
//           }
// 
//         // find newly created vertex in this triangle
//         for (int i = 0; i != 3; ++i)
//         {
//           if (viennagrid::regions(mesh, viennagrid::vertices(element)[i]).size() == 1)
//             triangle[1] = viennagrid::vertices(element)[i];
//         }
// 
//         // find newly created vertex in neighbor triangle
//         for (int i = 0; i != 3; ++i)
//         {
//           if (viennagrid::regions(mesh, viennagrid::vertices(*ntit)[i]).size() == 1)
//             triangle[2] = viennagrid::vertices(*ntit)[i];
//         }
// 
//         std::sort( triangle.begin(), triangle.end() );
//         
//         viennagrid::make_element(region, viennagrid::element_tag::triangle(), triangle.begin(), triangle.end());
//       }
//     }
// 
//     /* sets quantity field values for vertices and elements created when extruding contacts */
//     void extrude_contacts_quantity_fields()
//     {
//       for (std::map<std::string, viennagrid::quantity_field>::iterator it = quantity_fields_.begin();
//                                                                        it != quantity_fields_.end();
//                                                                      ++it)
//       {
//         viennagrid::quantity_field & quantities = it->second;
//         
//         if(quantities.topologic_dimension() > 0)
//         {
//           for(typename std::vector<RegionType>::const_iterator rit = contact_regions.begin(); rit != contact_regions.end(); ++rit)
//           {
//             typedef typename viennagrid::result_of::cell_range<RegionType>::type CellRangeType;
//             typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIterType;
//             
//             CellRangeType cells(*rit);
//             
//             for(CellIterType cit = cells.begin(); cit != cells.end(); ++cit)
//             {
//               quantities.set((*cit), 0.0);
//             }
//           }
//         }
//         else
//         {
//           for (std::map< int, std::vector<int> >::const_iterator nvit = newly_created_vertices.begin();
//                                                                   nvit != newly_created_vertices.end();
//                                                                 ++nvit)
//           {
//             double val = 0.0;
//             for (std::size_t i = 0; i != (*nvit).second.size(); ++i)
//               val += quantities.get( (*nvit).second[i] );
//             val /= (*nvit).second.size();
// 
//             quantities.set( vertices[nvit->first], val );
//           }
//         }
//       }
//     }

    /* returns a vector of all quantity fields read from the file */
    std::vector<viennagrid::quantity_field> quantity_fields()
    {
      std::vector<viennagrid::quantity_field> results;
      for (std::map<std::string, viennagrid::quantity_field>::iterator it = quantities.begin();
                                                                       it != quantities.end();
                                                                     ++it)
      {
        results.push_back( it->second );
      }

      return results;
    }

  };
  
  template<typename MeshT, typename RegionT>
  void tdr_element::add_to_mesh(tdr_geometry<MeshT> * geometry, RegionT & reg)
  {
    typedef typename tdr_geometry<MeshT>::VertexType VertexType;
    typedef typename tdr_geometry<MeshT>::ElementType ElementType;
    
    if(geometry->extrude_contacts && region->type != 0)
      extrude(geometry);
    
    // wo only have the vertex indices, retrieve the actual vertex objects
    std::vector<VertexType> cell_vertices;
    for(std::vector<int>::iterator it = vertices.begin();
                                  it != vertices.end(); ++it)
    {
      cell_vertices.push_back(geometry->vertices[*it]);
    }
  
    ElementType el = viennagrid::make_element( reg, element_type,
                      cell_vertices.begin(), cell_vertices.end() );
    
    region->element_indices.push_back(el.id());
  }
    
  template <typename MeshT>
  void tdr_element::extrude(tdr_geometry<MeshT> * geometry)
  {
    typedef typename tdr_geometry<MeshT>::PointType    PointType;
    typedef typename viennagrid::result_of::cell_range<MeshT>::type CellRangeType;
    typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;
    
    PointType center;
    PointType normal;
    double size = 0;
    
    if (element_type.is_line())
    {
      PointType p0 = viennagrid::get_point( geometry->vertices[vertices[0]] );
      PointType p1 = viennagrid::get_point( geometry->vertices[vertices[1]] );

      center = (p0+p1)/2;
      normal = normal_vector(p0, p1);

      size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
      
      element_type = viennagrid::element_tag::triangle();
    }
    else if (element_type.is_triangle())
    {
      PointType p0 = viennagrid::get_point( geometry->vertices[vertices[0]] );
      PointType p1 = viennagrid::get_point( geometry->vertices[vertices[1]] );
      PointType p2 = viennagrid::get_point( geometry->vertices[vertices[2]] );

      center = (p0+p1+p2)/3;
      normal = normal_vector(p0, p1, p2);

      size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
      size = std::max(size, viennagrid::distance(center, p2));
      
      element_type = viennagrid::element_tag::tetrahedron();
    }
    else
    {
      mythrow("extruding elements of type " << element_type.name() << " is not supported");
    }

    normal.normalize();
    normal *= size * geometry->extrude_contacts_scale;

    PointType new_vertex = center + normal;

    //check if the new vertex is inside any cell of the mesh
    CellRangeType cells(geometry->mesh);
    CellIteratorType cit = cells.begin();
    for (; cit != cells.end(); ++cit)
    {
      if ( viennagrid::is_inside(*cit, new_vertex) )
        break;
    }

    //if the vertex is inside the mesh, extrude in the other direction
    if (cit != cells.end())
      new_vertex = center - normal;

    //create the new vertex and add it to the element
    vertices.push_back(geometry->vertices.size());
    region->extruded_vertices.push_back(geometry->vertices.size());
    geometry->vertices.push_back(viennagrid::make_vertex(geometry->mesh, new_vertex));
  }
    
  template<typename PointT>
  PointT tdr_element::normal_vector(PointT const & p0, PointT const & p1)
  {
    PointT line = p1-p0;
    std::swap(line[0], line[1]);
    line[0] = -line[0];
    return line;
  }

  template<typename PointT>
  PointT tdr_element::normal_vector(PointT const & p0, PointT const & p1, PointT const & p2)
  {
    return viennagrid::cross_prod( p1-p0, p2-p0 );
  }

}

#endif
