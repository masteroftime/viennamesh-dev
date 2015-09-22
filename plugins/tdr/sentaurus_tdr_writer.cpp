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

#include "sentaurus_tdr_writer.hpp"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/container/flat_map.hpp>

#include "viennameshpp/logger.hpp"
#include "viennagridpp/range.hpp"

#include "H5Cpp.h"

namespace tdr
{
namespace
{

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


void write_attribute(H5::H5Object & obj, std::string const & name, std::string const & value)
{
  H5::StrType strdatatype(H5::PredType::C_S1, value.size());
  H5::Attribute attr = obj.createAttribute(name, strdatatype, H5::DataSpace(H5S_SCALAR));
  attr.write(strdatatype, value);
}

void write_attribute(H5::H5Object & obj, std::string const & name, int value)
{
  H5::Attribute attr = obj.createAttribute(name, H5::PredType::NATIVE_INT, H5::DataSpace());
  attr.write(H5::PredType::NATIVE_INT, &value);
}

void write_attribute(H5::H5Object & obj, std::string const & name, double value)
{
  H5::Attribute attr = obj.createAttribute(name, H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
  attr.write(H5::PredType::NATIVE_DOUBLE, &value);
}

template <typename T>
H5::DataSet write_dataset(H5::Group & group, std::string const & name, H5::DataType const & type, hsize_t length, std::vector<T> const & data)
{
  H5::DataSpace dataspace(1, &length);
  H5::DataSet dataset = group.createDataSet(name, type, dataspace);
  dataset.write(&data[0], type);
  return dataset;
}

/*
 * creates the "collection" and "geometry_0" with their most impartant parameters
 */
void sentaurus_tdr_writer::create_geometry()
{
  H5::Group collection = file.createGroup("collection");

  write_attribute(collection, "number of geometries", 1);
  write_attribute(collection, "number of plots", 0);

  geometry = collection.createGroup("geometry_0");

  write_attribute(geometry, "type", 1);
  write_attribute(geometry, "dimension", static_cast<int>(viennagrid::geometric_dimension(mesh)));
  write_attribute(geometry, "number of vertices", static_cast<int>(viennagrid::vertices(mesh).size()));
  write_attribute(geometry, "number of regions", static_cast<int>(mesh.region_count()));
  write_attribute(geometry, "number of states", 1);
}

/*
 * At the moment transformations are not supported by the TDR writer so
 * this methods writes the identity transformation to the TDR file
 */
void sentaurus_tdr_writer::write_identity_transformation()
{
  //transformation Group
  //for some reason this is apparently always 3d...
  H5::Group transformation = geometry.createGroup("transformation");
  write_attribute(transformation, "type", 1);

  std::vector<double> identity_matrix(9, 0.0);
  identity_matrix[0] = identity_matrix[4] = identity_matrix[8] = 1.0;
  write_dataset(transformation, "A", H5::PredType::NATIVE_DOUBLE, 9, identity_matrix);

  std::vector<double> zero_vector(3, 0.0);
  write_dataset(transformation, "b", H5::PredType::NATIVE_DOUBLE, 3, zero_vector);
}

/*
 * Writes the coordinates of all vertices to the "vertex" dataset.
 * 
 * Also creates vertex_mapping which links a vertex id from the mesh
 * to the vertex id in the TDR file.
 */
void sentaurus_tdr_writer::write_vertices()
{
  vertex_mapping.reserve(viennagrid::vertices(mesh).size());

  //store the vertex coordinates in a single vector (e.g.: x,y,z,x,y,z,x,...)
  std::vector<double> vertex_coordinates;
  vertex_coordinates.reserve(dimension*viennagrid::vertices(mesh).size());
  MeshVertexRange mesh_vertices(mesh);
  int id_counter = 0;
  for (MeshVertexIterator it = mesh_vertices.begin(); it != mesh_vertices.end(); ++it, ++id_counter)
  {
    vertex_mapping[(*it).id()] = id_counter;

    PointType const & p = viennagrid::get_point(*it);
    for (PointType::size_type i = 0; i < dimension; ++i)
    {
      vertex_coordinates.push_back(p[i]);
    }
  }

  //create the H5 Datatype based on the dimension of the mesh
  H5::CompType vertex_type( dimension*sizeof(double) );
  vertex_type.insertMember( "x", 0, H5::PredType::NATIVE_DOUBLE);
  if (dimension > 1)
  {
    vertex_type.insertMember( "y", sizeof(double), H5::PredType::NATIVE_DOUBLE);
  }
  if (dimension > 2)
  {
    vertex_type.insertMember( "z", 2*sizeof(double), H5::PredType::NATIVE_DOUBLE);
  }

  write_dataset(geometry, "vertex", vertex_type, vertex_coordinates.size()/dimension, vertex_coordinates);
}

/*
 * creates the region groups ("region_x") and writes the elements of each regions to the TDR file
 */
void sentaurus_tdr_writer::write_regions()
{
  //iterate over all regions
  RegionRange regions(mesh);
  int region_counter = 0;
  for (RegionIterator region_it = regions.begin(); region_it != regions.end(); ++region_it, ++region_counter)
  {
    //write name etc.
    H5::Group region_group = geometry.createGroup("region_" + boost::lexical_cast<std::string>(region_counter));
    write_attribute(region_group, "type", 0);
    write_attribute(region_group, "name", (*region_it).get_name());
    write_attribute(region_group, "material", "unknown"); //TODO
    write_attribute(region_group, "number of parts", 1);

    std::vector<int32_t> region_element_data;

    //iterate over all cells in the region
    CellRange cells(*region_it);
    unsigned int num_elements = 0;
    for (CellIterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it, ++num_elements)
    {
      // the first entry for each element is the element type id
      if((*cell_it).tag().is_line())
      {
        region_element_data.push_back(1);
      }
      else if((*cell_it).tag().is_triangle())
      {
        region_element_data.push_back(2);
      }
      else if((*cell_it).tag().is_quadrilateral())
      {
        region_element_data.push_back(3);
      }
      else if((*cell_it).tag().is_polygon())
      {
        region_element_data.push_back(4);
        
        // for polygons, the number of vertices follows after the type id
        uint32_t vertex_count = 0;
        BoundaryVertexRange vertices(*cell_it);
        for(BoundaryVertexIterator it = vertices.begin(); it != vertices.end(); ++it)
          ++vertex_count;
        
        region_element_data.push_back(vertex_count);
      }
      else if((*cell_it).tag().is_tetrahedron())
      {
        region_element_data.push_back(5);
      }
      else
      {
        throw viennautils::make_exception<writer_error>("Unsupported element type " + (*cell_it).tag().name());
      }

      //then write the ids of all vertices that make up this element
      BoundaryVertexRange vertices(*cell_it);
      for (BoundaryVertexIterator vertex_it = vertices.begin(); vertex_it != vertices.end(); ++vertex_it)
      {
	region_element_data.push_back(vertex_mapping[(*vertex_it).id()]);
      }
    }

    H5::DataSet elements = write_dataset(region_group, "elements_0", H5::PredType::NATIVE_INT32, region_element_data.size(), region_element_data);
    write_attribute(elements, "number of elements", static_cast<int>(num_elements));
  }
}

/*
 * Creates the "state_0" group and writes all quantity fields ("dataset_x").
 */
void sentaurus_tdr_writer::write_quantity_fields()
{
  //datasets
  H5::Group state = geometry.createGroup("state_0");
  write_attribute(state, "name", "state_0");
  write_attribute(state, "number of plots", 0);
  write_attribute(state, "number of string streams", 0);
  write_attribute(state, "number of xy plots", 0);
  
  RegionRange regions(mesh);

  //iterate over all quantity fields
  int num_datasets = 0;
  for (unsigned int i = 0; i < quantities.size(); ++i)
  {
    viennagrid::quantity_field const & quantity = quantities[i];
    
    //iterate over all regions
    int region_num = 0;
    for(RegionIterator region_it = regions.begin(); region_it != regions.end(); ++region_it, ++region_num)
    {
      //store all values which have to be written in a vector
      std::vector<double> values;
      bool valid = true;
      int location_type;
      
      if(quantity.topologic_dimension() == 0)
      {
        location_type = 0;
        
        RegionVertexRange vertices(*region_it);
        for(RegionVertexIterator vertex_it = vertices.begin(); vertex_it != vertices.end(); ++vertex_it)
        {
          if (!quantity.valid(*vertex_it))
          {
            valid = false;
            break;
          }
          
          values.push_back(quantity.get(*vertex_it));
        }
        
        if(!valid)
        {
          viennamesh::warning(1) << "Skipping quantity field " << quantity.get_name() << " in region " << region_num
            << " which is not defined for all vertices in that region" << std::endl;
        }
      }
      else if(quantity.topologic_dimension() == viennagrid::cell_dimension(mesh))
      {
        location_type = 3;
        
        CellRange cells(*region_it);
        for(CellIterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it)
        {
          if (!quantity.valid(*cell_it))
          {
            valid = false;
            break;
          }
          
          values.push_back(quantity.get(*cell_it));
        }
        
        if(!valid)
        {
          viennamesh::warning(1) << "Skipping quantity field " << quantity.get_name() << " in region " << region_num
            << " which is not defined for all cells in that region" << std::endl;
        }
      }
      else
      {
        valid = false;
        
        viennamesh::warning(1) << "Skipping quantity field " << quantity.get_name() << " in region " << region_num
            << ". Only quantity fields on vertices and cells are supported!" << std::endl;
      }

      if (valid) //only write quantities that are defined on the entire region
      {
	H5::Group dataset_group = state.createGroup("dataset_" + boost::lexical_cast<std::string>(num_datasets++));
	write_attribute(dataset_group, "number of values", static_cast<int>(values.size()));
	write_attribute(dataset_group, "location type", location_type);
	write_attribute(dataset_group, "structure type", 0);
	write_attribute(dataset_group, "value type", 2);
	write_attribute(dataset_group, "name", quantity.get_name());
	write_attribute(dataset_group, "quantity", quantity.get_name());
	write_attribute(dataset_group, "conversion factor", 1.0);
	write_attribute(dataset_group, "region", region_num);
	write_attribute(dataset_group, "unit:name", quantity.unit());
	write_dataset(dataset_group, "values", H5::PredType::NATIVE_DOUBLE, values.size(), values);
      }
    }
  }
  write_attribute(state, "number of datasets", num_datasets);
}

void sentaurus_tdr_writer::write_to_tdr()
{
  try 
  {
    dimension = viennagrid::geometric_dimension(mesh);
    if (dimension != 2)
    {
      throw viennautils::make_exception<writer_error>("TDR writer currently supports only two dimensional meshes");
    }
    
    create_geometry();
    write_identity_transformation();
    write_vertices();
    write_regions();
    write_quantity_fields();
  }
  catch(H5::Exception const & e)
  {
    throw viennautils::make_exception<writer_error>("caught HDF5 exception in HDF5 function: " + e.getFuncName() + " - with message: " + e.getDetailMsg());
  }
}

} //end of anonymous namespace

void write_to_tdr(const std::string& filename, const viennagrid::const_mesh& mesh, const std::vector< viennagrid::quantity_field >& quantities)
{
  sentaurus_tdr_writer w(filename, mesh, quantities);
  w.write_to_tdr();
}


} //end of namespace tdr
