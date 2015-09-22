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
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>

#include "viennautils/exception.hpp"

#include "viennagridpp/mesh/element_creation.hpp"
#include "viennagridpp/quantity_field.hpp"
#include "viennagridpp/algorithm/inclusion.hpp"
#include "viennagridpp/algorithm/distance.hpp"
#include "viennagridpp/point_accessor.hpp"

#include "viennameshpp/logger.hpp"

#include "H5Cpp.h"


namespace tdr
{
  struct read_error : virtual viennautils::exception {};
  
  int read_int(const H5::H5Object &g, const std::string name);
  double read_double(const H5::H5Object &g, const std::string name);
  std::string read_string(const H5::H5Object &g, const std::string name);
  
  // get the number of values in the given dataset
  hsize_t get_dataset_size(const H5::DataSet& dataset);
  
  // return a vector of all values in the given dataset
  std::vector<double> read_vector_double(const H5::DataSet& dataset);
  std::vector<int> read_vector_int(const H5::DataSet& dataset);

  
  template <typename MeshT>
  struct geometry;
  class region;
  
  /*
   * This class is responsible for the construction of elements from the
   * tdr data and adding these elements to the mesh/regions.
   */
  class element {
    viennagrid::element_tag element_type_;
    std::vector<int> vertices_;
    
  public:
    /*
     * Constructs the element from raw data read from the tdr file.
     * The parameter begin points to the first value of that element
     * in the element list (the element type) and is modified so that
     * it points to the next element after this method finishes. The
     * paramter end should point to the end of the element list (not the element)
     * and is just used to prevent reading beyond the end of the vector.
     */
    element(region & region, std::vector<int>::iterator & begin, std::vector<int>::iterator const & end);
    
    /*
     * Constructs viennagrid element and adss it to the mesh
     */
    template<typename MeshT>
    void add_to_mesh(geometry<MeshT> & geometry, region & region);
    
    /*
     * Extrudes this element (line -> triangle, triangle -> tetrahedron).
     * This method is applied to all elements in interface regions if
     * "extrude_contacts" is set.
     */
    template <typename MeshT>
    void extrude(geometry<MeshT> & geometry, region & region);
    
  private:
    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1);

    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1, PointT const & p2);
  };
  
  /*
   * A region can have one of the following types:
   *  0: "normal" region, contains elements
   *  1: outside interface region, at the edge of 1 other region, contains only lines
   *  2: inside interface region, between 2 other regions, contains only lines
   */
  enum region_type {
    normal = 0,
    interface = 1,
    inside_interface = 2
  };
  
  class region {
    int id_;
    std::string name_;
    std::vector<element> elements_;
    
    region_type type_;
    
    /*
     * contains all vertices in this regions. Needed for quantity fields
     */
    std::set<int> vertices_;
    
    /*
     * contains the ids of all elements in this region. Needed for quantity fields.
     */
    std::vector<int> element_indices_;
    
    /*
     * contains all vertices that were newly created when extruding elements. Needed for quantity fields.
     */
    std::vector<int> extruded_vertices_;
    
  public:
    /*
     * Reads the region including all contained elements from the tdr file.
     * 
     * The main reason this is not done in the constructor is because then the
     * region pointer in element would be invalidated when the region object
     * is moved/copied into the vector where it is stored.
     */
    region(int regnr, const H5::Group &reg);
    
    /*
     * Creates the viennagrid region and ads it to the mesh.
     * Then adds all child elements to the mesh.
     */
    template <typename MeshT>
    void add_to_mesh(geometry<MeshT> & geometry);
    
    /*
     * Redirects to the private overload of this method with the correct target paramter.
     * 
     * location_type 0 = vertices
     * location_type 3 = elements
     */
    void apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, int location_type);
    
    /*
     * When using "extrude_contacts" new elements and vertices are created for which
     * no values exist in the tdr file. This method takes the average of the neighbor
     * verties or element to fill in the missing values.
     */
    template <typename MeshT>
    void extrude_quantity(geometry<MeshT> & geometry, viennagrid::quantity_field & quantity);
    
    /*
     * If the fill_triangle_contacts parameter is set, this method will connect extruded vertices
     * using triangles.
     */
    template <typename MeshT>
    void fill_triangle_contacts(geometry<MeshT> & geometry);
    
    region_type type() { return type_; }
    int id() { return id_; }
    
    void add_vertex(int id) { vertices_.insert(id); }
    void add_extruded_vertex(int id) { extruded_vertices_.push_back(id); }
    void add_element(int id) { element_indices_.push_back(id); }
    
  private:
    /*
     * Read all elements contained in the given dataset
     */
    void read_elements(const H5::DataSet &elem);
    
    /*
     * For every given vertex/element given in "target" it sets the given value in the quantity_field.
     */
    template <typename ContainerT>
    void apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, const ContainerT & target);
    
    template <typename NeighborRangeType>
    double get_neighbor_average(const viennagrid::quantity_field & quantity, NeighborRangeType &neighbors);
  };
  
  /*
  * type 0 = boundary files, contains only polygons and lines and no datasets
  * type 1 = mesh files
  */
  enum geometry_type {
    boundary = 0,
    mesh = 1
  };
  
  /*
   * This class is the main class of the tdr reader.
   */
  template<typename MeshT>
  class geometry
  {
  public:
    typedef typename viennagrid::result_of::point<MeshT>::type PointType;
    typedef typename viennagrid::result_of::element<MeshT>::type VertexType;
    typedef typename viennagrid::result_of::element<MeshT>::type ElementType;
    typedef typename viennagrid::result_of::region<MeshT>::type RegionType;
    typedef typename viennagrid::result_of::neighbor_range<MeshT>::type NeighborRangeType;
    typedef typename viennagrid::result_of::iterator<NeighborRangeType>::type NeighborIteratorType;
    typedef typename viennagrid::result_of::cell_range<MeshT>::type CellRangeType;
    typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;
    
    typedef std::vector<boost::shared_ptr<region> >::iterator region_iterator;
    
  private:
    MeshT const & mesh_;
    std::vector<VertexType> vertices_;
    std::vector<boost::shared_ptr<region> > regions_;
    std::map<std::string, viennagrid::quantity_field> quantities_;
    
    int dim_;
    
    bool extrude_contacts_;
    double extrude_contacts_scale_;
    bool fill_triangle_contacts_;
    
  public:
    geometry(MeshT const & mesh) : mesh_(mesh), extrude_contacts_(true),
        extrude_contacts_scale_(1.0), fill_triangle_contacts_(false) {}

    void read_file(const std::string &filename);

    /* returns a vector of all quantity fields read from the file */
    std::vector<viennagrid::quantity_field> quantity_fields();
    
    int dimension() { return dim_; }
    
    bool extrude_contacts() { return extrude_contacts_; }
    double extrude_contacts_scale() { return extrude_contacts_scale_; }
    
    void set_extrude_contacts(bool extrude_contacts);
    void set_extrude_contacts_scale(double extrude_contacts_scale);
    void set_fill_triangle_contacts(bool fill_triangle_contacts);
    
    /*
     * Creates a new vertex and returns its index
     */
    int create_vertex(PointType p)
    {
      vertices_.push_back(viennagrid::make_vertex(mesh_, p));
      return vertices_.size() - 1;
    }
    
    VertexType get_vertex(unsigned int vertex_id)
    {
      assert(vertex_id < vertices_.size());
      return vertices_[vertex_id];
    }
    
    RegionType get_or_create_region(int region_id)
    {
      return mesh_.get_or_create_region(region_id);
    }
    
    template <typename ElementT>
    NeighborRangeType get_neighbors(ElementT e, int connection_dim, int neighbor_dim)
    {
      return NeighborRangeType(mesh_, e, connection_dim, neighbor_dim);
    }
    
    CellRangeType get_cells()
    {
      return CellRangeType(mesh_);
    }
    
  private:
    /*
     * Reads the transformation Matrix+Vector from the tdr file.
     * 
     * At the moment it only throws an error if the read transformation
     * is not equal to the identity transformation.
     */
    void read_transformation(const H5::Group &trans);

    void read_vertices(const H5::DataSet &vert, unsigned int nvertices);
    
    void read_regions(const H5::Group & group);

    void read_dataset(const H5::Group &dataset);

    void read_datasets(const H5::Group &state);

    void read_geometry(const H5::Group &geometry);
  };
  
//***************************  tdr::element  ***************************
  
  template<typename MeshT>
  void element::add_to_mesh(geometry<MeshT> & geometry, region & region)
  {
    typedef typename geometry<MeshT>::VertexType VertexType;
    typedef typename geometry<MeshT>::ElementType ElementType;
    
    if(geometry.extrude_contacts() && region.type() != normal)
      extrude(geometry, region);
    
    // wo only have the vertex indices, retrieve the actual vertex objects
    std::vector<VertexType> cell_vertices;
    for(std::vector<int>::iterator it = vertices_.begin();
                                  it != vertices_.end(); ++it)
    {
      cell_vertices.push_back(geometry.get_vertex(*it));
    }
  
    ElementType el = viennagrid::make_element( geometry.get_or_create_region(region.id()), 
        element_type_, cell_vertices.begin(), cell_vertices.end() );
    
    region.add_element(el.id());
  }
    
  template <typename MeshT>
  void element::extrude(geometry< MeshT > & geometry, region & region)
  {
    typedef typename geometry<MeshT>::PointType    PointType;
    typedef typename geometry<MeshT>::CellRangeType CellRangeType;
    typedef typename geometry<MeshT>::CellIteratorType CellIteratorType;
    
    PointType center;
    PointType normal;
    double size = 0;
    
    if (element_type_.is_line())
    {
      PointType p0 = viennagrid::get_point( geometry.get_vertex(vertices_[0]) );
      PointType p1 = viennagrid::get_point( geometry.get_vertex(vertices_[1]) );

      center = (p0+p1)/2;
      normal = normal_vector(p0, p1);

      size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
      
      element_type_ = viennagrid::element_tag::triangle();
    }
    else if (element_type_.is_triangle())
    {
      PointType p0 = viennagrid::get_point( geometry.get_vertex(vertices_[0]) );
      PointType p1 = viennagrid::get_point( geometry.get_vertex(vertices_[1]) );
      PointType p2 = viennagrid::get_point( geometry.get_vertex(vertices_[2]) );

      center = (p0+p1+p2)/3;
      normal = normal_vector(p0, p1, p2);

      size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
      size = std::max(size, viennagrid::distance(center, p2));
      
      element_type_ = viennagrid::element_tag::tetrahedron();
    }
    else
    {
      viennautils::make_exception<read_error>("extruding elements of type " + element_type_.name() + " is not supported");
    }

    normal.normalize();
    normal *= size * geometry.extrude_contacts_scale();

    PointType new_vertex = center + normal;

    //check if the new vertex is inside any cell of the mesh
    CellRangeType cells = geometry.get_cells();
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
    int new_index = geometry.create_vertex(new_vertex);
    vertices_.push_back(new_index);
    region.add_extruded_vertex(new_index);
  }
    
  template<typename PointT>
  PointT element::normal_vector(PointT const & p0, PointT const & p1)
  {
    PointT line = p1-p0;
    std::swap(line[0], line[1]);
    line[0] = -line[0];
    return line;
  }

  template<typename PointT>
  PointT element::normal_vector(PointT const & p0, PointT const & p1, PointT const & p2)
  {
    return viennagrid::cross_prod( p1-p0, p2-p0 );
  }
  
//***************************  tdr::region  ***************************
  
  template <typename MeshT>
  void region::add_to_mesh(geometry<MeshT> & geometry)
  {
    typedef typename geometry<MeshT>::RegionType RegionType;
    
    //TODO: Handle type 2 regions?
    if(type_ == inside_interface)
    {
      viennamesh::warning(1) << "Type 2 regions (interfaces) are ignored!" << std::endl;
      return;
    }
    
    RegionType reg = geometry.get_or_create_region(id_);
    reg.set_name(name_);
    
    for(std::vector<element>::iterator it = elements_.begin(); it != elements_.end(); ++it)
    {
      element & e = *it;
      e.add_to_mesh(geometry, *this);
    }
  }
  
  template <typename MeshT>
  void region::extrude_quantity(geometry<MeshT> & geometry, viennagrid::quantity_field & quantity)
  {
    typedef typename geometry<MeshT>::VertexType VertexType;
    typedef typename geometry<MeshT>::RegionType RegionType;
    typedef typename geometry<MeshT>::NeighborRangeType NeighborRangeType;
    typedef typename viennagrid::result_of::cell_range<RegionType>::type CellRangeType;
    typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;
    
    //only regions of type 1 contain extruded elements
    //TODO: Handle type 2 regions?
    if(type_ == interface)
    {
      int dimension = quantity.topologic_dimension();
      
      if(dimension == 0)
      {
        //vertex quantity -> iterate over the extruded vertices
        std::vector<int>::const_iterator it = extruded_vertices_.begin();
        for(; it != extruded_vertices_.end(); ++it)
        {
          //set neighbor average as quantity value
          VertexType v = geometry.get_vertex(*it);
          NeighborRangeType neighbor_vertices = geometry.get_neighbors(v, 1, 0);
          quantity.set(v, get_neighbor_average(quantity, neighbor_vertices));
        }
      }
      else if(dimension == geometry.dimension())
      {
        //element quantity -> iterate over the extruded elements
        //Note: type 1 regions contain only extruded elements
        RegionType reg = geometry.get_or_create_region(id_);
        CellRangeType cells(reg);
        
        CellIteratorType cell = cells.begin();
        for(; cell != cells.end(); ++cell)
        {
          //set neighbor average as quantity value
          NeighborRangeType neighbor_cells = geometry.get_neighbors(*cell, dimension-1, dimension);
          quantity.set(*cell, get_neighbor_average(quantity, neighbor_cells));
        }
      }
      else
      {
        viennautils::make_exception<read_error>("connot extrude quantity field of topological dimension " + boost::lexical_cast<std::string>(dimension));
      }
    }
  }
  
  template <typename MeshT>
  void region::fill_triangle_contacts(geometry<MeshT> & geometry)
  {
    typedef typename geometry<MeshT>::VertexType VertexType;
    typedef typename geometry<MeshT>::RegionType RegionType;
    typedef typename geometry<MeshT>::NeighborRangeType NeighborRangeType;
    typedef typename geometry<MeshT>::NeighborIteratorType NeighborIteratorType;
    
    //TODO: Handle type 2 regions?
    if(type_ != interface)
      return;
    
    RegionType reg = geometry.get_or_create_region(id_);
    std::vector<VertexType> triangles;
    
    //iterate over all pairs of extruded vertices
    for(std::vector<int>::iterator it1 = extruded_vertices_.begin(); it1 != extruded_vertices_.end(); ++it1)
      for(std::vector<int>::iterator it2 = it1 + 1; it2 != extruded_vertices_.end(); ++it2)
      {
        VertexType v1 = geometry.get_vertex(*it1);
        VertexType v2 = geometry.get_vertex(*it2);
        
        //get the neighbor vertices
        NeighborRangeType n1 = geometry.get_neighbors(v1, 1, 0);
        NeighborRangeType n2 = geometry.get_neighbors(v2, 1, 0);
        
        for(NeighborIteratorType nit1 = n1.begin(); nit1 != n1.end(); ++nit1)
          for(NeighborIteratorType nit2 = n2.begin(); nit2 != n2.end(); ++nit2)
          {
            if(*nit1 == *nit2)
            {
              //viennagrid::make_triangle(reg, v1, v2, *nit1);
              triangles.push_back(v1);
              triangles.push_back(v2);
              triangles.push_back(*nit1);
            }
          }
      }
    
    for(typename std::vector<VertexType>::iterator t = triangles.begin(); t != triangles.end(); t += 3)
    {
      viennagrid::make_element(reg, viennagrid::element_tag::triangle(), t, t+3);
    }
  }
  
  template <typename NeighborRangeType>
  double region::get_neighbor_average(const viennagrid::quantity_field & quantity, NeighborRangeType &neighbors)
  {
    typedef typename viennagrid::result_of::iterator<NeighborRangeType>::type NeighborIteratorType;
    double value = 0;
    int count = 0;
    
    NeighborIteratorType neighbor = neighbors.begin();
    for(; neighbor != neighbors.end(); ++neighbor)
    {
      //only consider neighbors with a valid value
      if(quantity.valid(*neighbor))
      {
        value += quantity.get(*neighbor);
        ++count;
      }
    }
    return value / (double)count;
  }
  
//***************************  tdr::geometry  ***************************
  
  template <typename MeshT>
  void geometry<MeshT>::read_transformation(const H5::Group &trans)
  {
    boost::array<double, 9> trans_matrix, identity_matrix;
    boost::array<double, 3> trans_vector, zero_vector;
    
    // initialize identity tranformation
    zero_vector.fill(0);
    identity_matrix.fill(0);
    identity_matrix[0] = identity_matrix[4] = identity_matrix[8] = 1;
    
    // read transformation from file
    const H5::DataSet &A=trans.openDataSet("A");
    const H5::DataSet &b=trans.openDataSet("b");
    A.read( trans_matrix.elems, H5::PredType::NATIVE_DOUBLE);
    b.read( trans_vector.elems, H5::PredType::NATIVE_DOUBLE);
    
    if(trans_matrix != identity_matrix || trans_vector != zero_vector)
      viennautils::make_exception<read_error>("transformation not equal to identity");
  }

  template <typename MeshT>
  void geometry<MeshT>::read_vertices(const H5::DataSet &vert, unsigned int nvertices)
  {
    hsize_t size = get_dataset_size(vert);
    
    if (nvertices != size)
      viennautils::make_exception<read_error>("nvertices not equal vertices.dim");

    H5::CompType vtype( sizeof(double)*dim_ );
    
    vtype.insertMember( "x", 0, H5::PredType::NATIVE_DOUBLE);

    if (dim_>1)
      vtype.insertMember( "y", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    if (dim_>2)
      vtype.insertMember( "z", sizeof(double)*2, H5::PredType::NATIVE_DOUBLE);

    std::vector<double> data;
    data.reserve(size*dim_);
    vert.read( &data[0], vtype );
    
    for (unsigned int i=0; i < size*dim_; i += dim_)
    {
      PointType point(dim_);
      
      point[0] = data[i];
      if (dim_>1)
      point[1] = data[i+1];
      if (dim_>2)
      point[2] = data[i+2];

      vertices_.push_back(viennagrid::make_vertex(mesh_, point));
    }
  }
  
  template <typename MeshT>
  void geometry<MeshT>::read_regions(const H5::Group & group)
  {
    int nregions = read_int(group, "number of regions");
    regions_.reserve(nregions);
    
    for (int i=0; i<nregions; i++)
    {
      const H5::Group &reg=group.openGroup("region_" + boost::lexical_cast<std::string>(i));
      regions_.push_back(boost::make_shared<region>(i, reg));
    }
  }

  template <typename MeshT>
  void geometry<MeshT>::read_dataset(const H5::Group &dataset)
  {
    std::string name = read_string(dataset,"name");
    int regnr = read_int(dataset,"region");
    int location_type = read_int(dataset,"location type");
    std::string unit=read_string(dataset,"unit:name");
    double conversion_factor = read_double(dataset, "conversion factor");
    
    if (location_type != 0 && location_type != 3)
      viennautils::make_exception<read_error>("Location type " + boost::lexical_cast<std::string>(location_type) + " in dataset " + name + " not supported");

    if (read_int(dataset,"structure type") != 0)
      viennautils::make_exception<read_error>("Dataset " + name + " structure type not 0");

    if (read_int(dataset,"value type") != 2)
      viennautils::make_exception<read_error>("Dataset " + name + " value type not 2");
    
    boost::shared_ptr<region> reg = regions_[regnr];
    
    //TODO:: Handle type 2 regions?
    if(reg->type() == inside_interface)
      return;
    
    viennagrid::quantity_field & quantity = quantities_[name];
    
    //check if quantity field is already initialized
    if(quantity.get_name() != name)
    {
      //and initialize it if not
      if(location_type == 0)
        quantity.init(0, 1);
      else
        quantity.init(dim_, 1);
      
      quantity.set_name(name);
      quantity.set_unit(unit);
    }
    
    if(quantity.unit() != unit)
      viennautils::make_exception<read_error>("Units not the same in all regions in quantity field " + name + " - " + dataset.getObjName());
    
    std::vector<double> values = read_vector_double(dataset.openDataSet("values"));
    
    if(conversion_factor != 1.0)
    {
      viennamesh::warning(1) << "Values multiplied with conversion factor " << conversion_factor << " in " << dataset.getObjName() << std::endl;
      
      for(unsigned int i = 0; i < values.size(); ++i)
      {
        values[i] *= conversion_factor;
      }
    }
    
    reg->apply_quantity(quantity, values, location_type);
  }

  template <typename MeshT>
  void geometry<MeshT>::read_datasets(const H5::Group &state)
  {
    if (read_int(state,"number of plots") != 0)
      viennautils::make_exception<read_error>("'Number of plots' not equal 0");

    if (read_int(state,"number of string streams") != 0)
      viennautils::make_exception<read_error>("'Number of string streams' not equal 0");

    if (read_int(state,"number of xy plots") != 0)
      viennautils::make_exception<read_error>("'Number of xy plots' not equal 0");

    int ndatasets=read_int(state,"number of datasets");
    
    for (int i=0; i<ndatasets; i++)
    {
      read_dataset(state.openGroup("dataset_" + boost::lexical_cast<std::string>(i)));
    }
  }

  template <typename MeshT>
  void geometry<MeshT>::read_geometry(const H5::Group &geometry)
  {
    //int type = read_int(geometry,"type");
    int dum = read_int(geometry,"number of states");
    dim_ = read_int(geometry,"dimension");
    unsigned int nvertices = read_int(geometry,"number of vertices");
    
    geometry_type type;
    switch(read_int(geometry,"type"))
    {
      case 0:
        type = boundary;
        break;
      case 1:
        type = mesh;
        break;
      default:
        viennautils::make_exception<read_error>("Unknown geometry type " + boost::lexical_cast<std::string>(type) + " in " + geometry.getObjName());
    }

    if (type == mesh && dum != 1)
      viennautils::make_exception<read_error>("Number of states not one");
    
    if(type == boundary && extrude_contacts_)
    {
      //extrude contacts does not work with polygons, because is_inside is not implemented for polygons in viennagrid
      viennamesh::warning(1) << "extrude_contacts not available with geometry type 0" << std::endl;
      set_extrude_contacts(false);
    }

    read_transformation(geometry.openGroup("transformation"));
    read_vertices(geometry.openDataSet("vertex"), nvertices);
    read_regions(geometry);
    
    for(region_iterator it = regions_.begin(); it != regions_.end(); ++it)
    {
      (*it)->add_to_mesh(*this);
      
      if(fill_triangle_contacts_)
        (*it)->fill_triangle_contacts(*this);
    }
    
    //only type 1 geometries contain datasets
    if(type == mesh)
    {
      read_datasets(geometry.openGroup("state_0"));
      
      for(region_iterator reg = regions_.begin(); reg != regions_.end(); ++reg)
      {
        for(std::map<std::string, viennagrid::quantity_field>::iterator q = quantities_.begin(); q != quantities_.end(); ++q)
        {
          (*reg)->extrude_quantity(*this, q->second);
        }
      }
    }
  }

  template <typename MeshT>
  void geometry<MeshT>::read_file(const std::string &filename)
  {
    boost::shared_ptr<H5::H5File> file = boost::make_shared<H5::H5File>(filename.c_str(), H5F_ACC_RDWR);

    if (file->getNumObjs()!=1)
      viennautils::make_exception<read_error>("File has not exactly one collection");

    const H5::Group &collection = file->openGroup("collection");
    
    if (read_int(collection,"number of geometries") != 1)
      viennautils::make_exception<read_error>("File has not exactly one geometry");
    
    if (read_int(collection,"number of plots") != 0)
      viennamesh::warning(1) << "file contains plots, skip them" << std::endl;
    
    read_geometry(collection.openGroup("geometry_0"));
  }

  template <typename MeshT>
  std::vector<viennagrid::quantity_field> geometry<MeshT>::quantity_fields()
  {
    std::vector<viennagrid::quantity_field> results;
    for (std::map<std::string, viennagrid::quantity_field>::iterator it = quantities_.begin();
                                                                      it != quantities_.end();
                                                                    ++it)
    {
      results.push_back( it->second );
    }

    return results;
  }
  
  template <typename MeshT>
  void geometry<MeshT>::set_extrude_contacts(bool extrude_contacts)
  {
    extrude_contacts_ = extrude_contacts;
    
    if(!extrude_contacts && fill_triangle_contacts_)
    {
      fill_triangle_contacts_ = false;
      viennamesh::warning(1) << "extrude_contacts is disabled so fill_triangle_contacts has no effect" << std::endl;
    }
  }
  
  template <typename MeshT>
  void geometry<MeshT>::set_extrude_contacts_scale(double extrude_contacts_scale)
  {
    extrude_contacts_scale_ = extrude_contacts_scale;
  }
  
  template <typename MeshT>
  void geometry<MeshT>::set_fill_triangle_contacts(bool fill_triangle_contacts)
  { 
    if(extrude_contacts_)
    {
      fill_triangle_contacts_ = fill_triangle_contacts;
    }
    else
    {
      viennamesh::warning(1) << "extrude_contacts is disabled so fill_triangle_contacts has no effect" << std::endl;
      fill_triangle_contacts = false;
    }
  }
}

#endif
