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

    char * buf = new char[a.getDataType().getSize()+1];
    a.read(a.getDataType(),buf);
    buf[a.getDataType().getSize()]='\0';

    string result(buf);
    delete[] buf;

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

}
