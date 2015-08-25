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

namespace viennamesh
{
  int read_int(const H5Object &g, const string name)
  {
    int i;
    Attribute a=g.openAttribute(name);
    if (a.getTypeClass()!=H5T_INTEGER)
      mythrow("Wrong class in atrribute");
    a.read(a.getDataType(),&i);
    return i;
  }

  double read_double(const H5Object &g, const string name)
  {
    double i;
    Attribute a=g.openAttribute(name);
    //if (a.getTypeClass()!=H5T_DOUBLE)
    if (a.getTypeClass()!=H5T_FLOAT)
      mythrow("Wrong class " << typeid(a.getTypeClass()).name() << " in atrribute");
    a.read(a.getDataType(),&i);
    return i;
  }

  //string read_string(const Group &g, const string name)
  string read_string(const H5Object &g, const string name)
  {
    Attribute a=g.openAttribute(name);
    if (a.getTypeClass()!=H5T_STRING)
      mythrow("Wrong class in atrribute");

    string result;
    a.read(a.getDataType(), result);

    return result;
  }
  
  hsize_t get_dataset_size(const DataSet& dataset)
  {
    const DataSpace &dataspace = dataset.getSpace();
    int ndims = dataspace.getSimpleExtentNdims();
    
    if(ndims != 1)
      mythrow("Dataset " << dataset.getObjName() << " has " << ndims << " dimensions");
    
    hsize_t size;
    dataspace.getSimpleExtentDims(&size);
    return size;
  }

std::vector< double > read_vector(const DataSet& dataset)
  {
    std::vector<double> result;
    
    result.resize(get_dataset_size(dataset));
    dataset.read( &result[0], PredType::NATIVE_DOUBLE);
    
    return result;
  }
  
  
  tdr_element::tdr_element(tdr_region* region, int*& begin, int* end)
    : region(region)
  {
    int vertex_num;
  
    switch (*begin++)
    {
      case 1:
        element_type = viennagrid::element_tag::line();
        vertex_num = 2;
        break;
      case 2:
        element_type = viennagrid::element_tag::triangle();
        vertex_num = 3;
        break;
      case 3:
        element_type = viennagrid::element_tag::quadrilateral();
        vertex_num = 4;
        break;
      case 4:
        element_type = viennagrid::element_tag::polygon();
        vertex_num = *begin++;
        break;
      case 5:
        element_type = viennagrid::element_tag::tetrahedron();
        vertex_num = 4;
        break;
      default:
        mythrow("Element type " << begin[-1] << " in region_" << region->id << " not known");
    }
    
    if(begin+vertex_num > end)
      mythrow("Incomplete element in region_" << region->id);

    for (int i=0; i<vertex_num; i++)
    {
      region->vertices.insert(*begin);
      vertices.push_back(*begin++);
    }
    
    /*if(region_type != 0 && extrude_contacts)
      extrude_element(element);
      
      if(region_type != 0 && extrude_contacts && fill_triangle_contacts) 
        fill_triangle(e, region);*/
    
    
  }
}
