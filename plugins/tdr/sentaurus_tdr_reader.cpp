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

#include "sentaurus_tdr_reader.hpp"

namespace tdr
{
  int read_int(const H5::H5Object &g, const std::string name)
  {
    int i;
    H5::Attribute a=g.openAttribute(name);
    if (a.getTypeClass()!=H5T_INTEGER)
      viennautils::make_exception<read_error>("Wrong class in atrribute " + name + " (" + g.getObjName() + ")");
    a.read(a.getDataType(),&i);
    return i;
  }

  double read_double(const H5::H5Object &g, const std::string name)
  {
    double i;
    H5::Attribute a=g.openAttribute(name);
    //if (a.getTypeClass()!=H5T_DOUBLE)
    if (a.getTypeClass()!=H5T_FLOAT)
      viennautils::make_exception<read_error>("Wrong class in atrribute " + name + " (" + g.getObjName() + ")");
    a.read(a.getDataType(),&i);
    return i;
  }

  //string read_string(const Group &g, const string name)
  std::string read_string(const H5::H5Object &g, const std::string name)
  {
    H5::Attribute a=g.openAttribute(name);
    if (a.getTypeClass()!=H5T_STRING)
      viennautils::make_exception<read_error>("Wrong class in atrribute " + name + " (" + g.getObjName() + ")");

    std::string result;
    a.read(a.getDataType(), result);

    return result;
  }
      
  template <typename T>
  std::vector<T> read_vector(const H5::DataSet & dataset, const H5::PredType type)
  {
    std::vector<T> result;
    
    result.resize(get_dataset_size(dataset));
    dataset.read( &result[0], type);
    
    return result;
  }
  
  std::vector<double> read_vector_double(const H5::DataSet& dataset)
  {
    return read_vector<double>(dataset, H5::PredType::NATIVE_DOUBLE);
  }
  
  std::vector<int> read_vector_int(const H5::DataSet& dataset)
  {
    return read_vector<int>(dataset, H5::PredType::NATIVE_INT);
  }
  
  hsize_t get_dataset_size(const H5::DataSet& dataset)
  {
    const H5::DataSpace &dataspace = dataset.getSpace();
    int ndims = dataspace.getSimpleExtentNdims();
    
    if(ndims != 1)
      viennautils::make_exception<read_error>("Dataset " + dataset.getObjName() + " has mor than one dimension");
    
    hsize_t size;
    dataspace.getSimpleExtentDims(&size);
    return size;
  }
  
  
  element::element(region & region, std::vector<int>::iterator & begin, std::vector<int>::iterator const & end)
  {
    int vertex_num;
  
    switch (*begin++)
    {
      case 1:
        element_type_ = viennagrid::element_tag::line();
        vertex_num = 2;
        break;
      case 2:
        element_type_ = viennagrid::element_tag::triangle();
        vertex_num = 3;
        break;
      case 3:
        element_type_ = viennagrid::element_tag::quadrilateral();
        vertex_num = 4;
        break;
      case 4:
        element_type_ = viennagrid::element_tag::polygon();
        vertex_num = *begin++;
        break;
      case 5:
        element_type_ = viennagrid::element_tag::tetrahedron();
        vertex_num = 4;
        break;
      default:
        viennautils::make_exception<read_error>("Element type " + boost::lexical_cast<std::string>(begin[-1]) 
            + " in region_" + boost::lexical_cast<std::string>(region.id()) + " not known");
    }
    
    if(begin+vertex_num > end)
      viennautils::make_exception<read_error>("Incomplete element in region_" + boost::lexical_cast<std::string>(region.id()));

    for (int i=0; i<vertex_num; i++)
    {
      region.add_vertex(*begin);
      vertices_.push_back(*begin++);
    }
  }
  
  region::region(int regnr, const H5::Group &reg)
  {
    id_ = regnr;
    name_ = read_string(reg,"name");
    
    switch(read_int(reg, "type"))
    {
      case 0:
        type_ = normal;
        break;
      case 1:
        type_ = interface;
        break;
      case 2:
        type_ = inside_interface;
        break;
      default:
        throw viennautils::make_exception<read_error>("Unknown region type " + boost::lexical_cast<std::string>(type_) + " found in " + reg.getObjName());
    }
    
    //iterate over all elements_* and read the contained elements
    const int n=reg.getNumObjs();
    for (int i=0; i<n; i++)
    {
      std::string objname = reg.getObjnameByIdx(i);
      if (objname.find("elements_") == 0)
      {
        read_elements(reg.openDataSet(objname));
      }

//       removed because it is untested
//       if (objname.find("part_") == 0)
//       {
//         read_elements(reg.openGroup(objname).openDataSet("elements"));
//       }
    }
  }
  
  void region::read_elements(const H5::DataSet &elem)
  {
    std::vector<int> element_data = read_vector_int(elem);

    std::vector<int>::iterator it = element_data.begin();
    
    while(it != element_data.end())
      elements_.push_back(element(*this, it, element_data.end()));
  }
  
  void region::apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, int location_type)
  {
    switch(location_type) {
      case 0:
        apply_quantity(quantity, values, vertices_);
        break;
      case 3:
        apply_quantity(quantity, values, element_indices_);
        break;
      default:
        throw viennautils::make_exception<read_error>("Invalid location type in quantity field " + quantity.get_name());
    }
  }
  
  template <typename ContainerT>
  void region::apply_quantity(viennagrid::quantity_field & quantity, const std::vector<double> & values, const ContainerT & target)
  {
    if(values.size() != target.size())
      viennautils::make_exception<read_error>("Quantity field " + quantity.get_name() + " in region " + boost::lexical_cast<std::string>(id_) + " has invalid size");
    
    int i = 0;
    for(typename ContainerT::const_iterator t = target.begin(); t != target.end(); ++t)
    {
      quantity.set(*t, values[i++]);
    }
  }
}
