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
#include <set>
#include <map>
#include <vector>
#include <typeinfo>
#include <cstdlib>
#include <boost/scoped_array.hpp>

using std::string;

#include <iostream>
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#define FALSE   0

#define mythrow(a) { std::cerr << a << std::endl; throw; }


#include "viennagridpp/mesh/element_creation.hpp"
#include "viennagridpp/quantity_field.hpp"

#include "viennagridpp/algorithm/cross_prod.hpp"
#include "viennagridpp/algorithm/centroid.hpp"
#include "viennagridpp/algorithm/inclusion.hpp"
#include "viennagridpp/algorithm/geometry.hpp"
#include "viennagridpp/algorithm/distance.hpp"
#include "viennagridpp/algorithm/norm.hpp"


namespace viennamesh
{

  // Operator function
  //herr_t get_all_groups(hid_t loc_id, const char *name, void *opdata);
  //herr_t file_info(hid_t loc_id, const char *name, void *opdata);

  int read_int(const H5Object &g, const string name);
  double read_double(const H5Object &g, const string name);
  string read_string(const H5Object &g, const string name);
  hsize_t get_dataset_size(const DataSet& dataset);

  struct element_t
  {
    viennagrid_element_type element_type;
    std::vector<int> vertices;
  };
  
  
  template<typename MeshT>
  struct tdr_geometry
  {
    typedef typename viennagrid::result_of::point<MeshT>::type PointType;
    typedef typename viennagrid::result_of::element<MeshT>::type VertexType;
    typedef typename viennagrid::result_of::region<MeshT>::type RegionType;
    
    unsigned int nvertices;
    int dim,nregions,ndatasets;
    double trans_matrix[9],trans_move[3];
    std::map< int, std::vector<int> > newly_created_vertices;
    
    MeshT const & mesh;
    std::vector<VertexType> vertices;
    std::map<std::string, viennagrid::quantity_field> quantity_fields_;
    
    bool extrude_contacts;
    double extrude_contacts_scale;
    
    tdr_geometry(MeshT const & mesh) : mesh(mesh) {}

    void read_transformation(const Group &trans)
    {
      const DataSet &A=trans.openDataSet("A");
      const DataSet &b=trans.openDataSet("b");
      A.read( trans_matrix, PredType::NATIVE_DOUBLE);
      b.read( trans_move, PredType::NATIVE_DOUBLE);
    }

    typedef struct coord2_t {
      double x;
      double y;
      double z;
    } coord2_t;


    void read_vertex(const DataSet &vert)
    {
      hsize_t size = get_dataset_size(vert);
      
      if (nvertices != size)
        mythrow("nvertices not equal vertices.dim");

      CompType mtype2( sizeof(coord2_t) );
      
      mtype2.insertMember( "x", HOFFSET(coord2_t, x), PredType::NATIVE_DOUBLE);

      if (dim>1)
        mtype2.insertMember( "y", HOFFSET(coord2_t, y), PredType::NATIVE_DOUBLE);

      if (dim>2)
        mtype2.insertMember( "z", HOFFSET(coord2_t, z), PredType::NATIVE_DOUBLE);

      boost::scoped_array<coord2_t> s2(new coord2_t[size]);
      vert.read( s2.get(), mtype2 );
      
      for (unsigned int i=0; i<size; i++)
      {
        PointType point(dim);
        
        point[0] = s2[i].x;
        if (dim>1)
        point[1] = s2[i].y;
        if (dim>2)
        point[2] = s2[i].z;

        vertices.push_back(viennagrid::make_vertex(mesh, point));
      }
    }

    void read_elements(RegionType region, const DataSet &elem, bool extrude_elements)
    {
      hsize_t size = get_dataset_size(elem);

      boost::scoped_array<int> el(new int[size]);
      elem.read( el.get(), PredType::NATIVE_INT);

      unsigned int elct=0;
      while (elct < size)
      {
        element_t element;
        int vertex_num;
        
        switch (el[elct++])
        {
          case 1:
            element.element_type = VIENNAGRID_ELEMENT_TYPE_LINE;
            vertex_num = 2;
            break;
          case 2:
            element.element_type = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
            vertex_num = 3;
            break;
          case 3:
            element.element_type = VIENNAGRID_ELEMENT_TYPE_QUADRILATERAL;
            vertex_num = 4;
            break;
          case 5:
            element.element_type = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
            vertex_num = 4;
            break;
          default:
            mythrow("Element type " << el[elct-1] << " in region " << region.get_name() << " not known");
        }

        for (int i=0; i<vertex_num; i++)
          element.vertices.push_back(el[elct++]);
        
        if(extrude_elements)
          extrude_element(element);
        
        std::vector<VertexType> cell_vertices;
        for(std::vector<int>::iterator it = element.vertices.begin();
                                        it != element.vertices.end(); ++it)
        {
          cell_vertices.push_back(vertices[*it]);
        }
        
        viennagrid::make_element( region, viennagrid::element_tag::from_internal(element.element_type),
                            cell_vertices.begin(), cell_vertices.end() );
      }
    }

    void read_region(const int regnr, const Group &reg)
    {
      string name = read_string(reg,"name");
      int type = read_int(reg,"type");
      
      /*
       * It seems that:
       *  0 -> "normal" region
       *  1 -> contact
       *  2 -> interface between 2 regions
       */
      if(type > 2)
        mythrow("Unknown region type " << type << " found in " << reg.getObjName());
      
      RegionType region = mesh.get_or_create_region(regnr);
      region.set_name(name);
      
      const int n=reg.getNumObjs();
      for (int i=0; i<n; i++)
      {
        std::string objname = reg.getObjnameByIdx(i);
        if (objname.find("elements_") == 0)
        {
          const DataSet &ds=reg.openDataSet(objname);
          read_elements(region, ds, type==0 ? false:true);
        }

        if (objname.find("part_") == 0)
        {
          const Group &part=reg.openGroup(objname);
          read_elements(region, part.openDataSet("elements"), type==0 ? false:true);
        }
      }
    }

    /*region_t &find_region(int regnr)
    {
      for (std::map<string,region_t>::iterator R=region.begin(); R!=region.end(); R++)
        if (R->second.regnr==regnr)
          return R->second;
      
      mythrow("Region " << regnr << " not found");
    }*/

    template <typename RangeT>
    void read_values(const DataSet &values, viennagrid::quantity_field & quantities, RangeT const & el)
    {
      typedef typename viennagrid::result_of::iterator<RangeT>::type IteratorType;
      
      hsize_t size = get_dataset_size(values);
      unsigned int nvalues = el.end() - el.begin();
      
      if (nvalues != size)
        mythrow("Dataset " << values.getObjName() << " should have " << nvalues << " values, but has " << size);

      boost::scoped_array<double> v(new double[size]);
      values.read( v.get(), PredType::NATIVE_DOUBLE);
      
      int i = 0;
      for(IteratorType it = el.begin(); it != el.end(); ++it)
      {
        quantities.set(*it, v[i++]);
      }
    }

    void read_dataset(const Group &dataset)
    {
      string name = read_string(dataset,"name");

      //string quantity = read_string(dataset,"quantity");
      int regnr = read_int(dataset,"region");
      //int nvalues=read_int(dataset,"number of values");
      //double conversion_factor = read_double(dataset,"conversion factor");

      if (read_int(dataset,"location type") != 0)
        mythrow("Dataset " << name << " location type not 0");

      if (read_int(dataset,"structure type") != 0)
        mythrow("Dataset " << name << " structure type not 0");

      if (read_int(dataset,"value type") != 2)
        mythrow("Dataset " << name << " value type not 2");

      // In this dataset we have a group
      // tag_group_0???
      // units: take only the unit name
      // Dataset values: the actual values

      /*string unit;
      int n=dataset.getNumObjs(),i;
      for (i=0; i<n; i++)
      {
        if (dataset.getObjnameByIdx(i)=="unit")
        {
          const Group &u=dataset.openGroup("unit");
          unit=read_string(u,"name");
          break;
        }
      }
      if (i==n)
      {
        unit=read_string(dataset,"unit:name");
      }*/

      /*region_t &region=find_region(regnr);
      region.dataset[name].name=name;
      region.dataset[name].quantity=quantity;
      region.dataset[name].nvalues=nvalues;
      region.dataset[name].conversion_factor=conversion_factor;
      region.dataset[name].unit=unit;*/
      
      RegionType region = mesh.get_region(regnr);
      
      viennagrid::quantity_field & quantities = quantity_fields_[name];
      
      if(quantities.get_name() != name)
      {
        quantities.init(0, 1);
        quantities.set_name(name);
      }
      
      typedef typename viennagrid::result_of::const_vertex_range<RegionType>::type VertexRangeType;

      VertexRangeType vertices(region);

      read_values(dataset.openDataSet("values"), quantities, vertices);
    }

    void read_attribs0(const Group &state)
    {
      if (read_int(state,"number of plots") != 0)
        mythrow("Numberofplots not equal 0");

      if (read_int(state,"number of string streams") != 0)
        mythrow("Numberofstringstreams not equal 0");

      if (read_int(state,"number of xy plots") != 0)
        mythrow("Numberofxyplots not equal 0");

      ndatasets=read_int(state,"number of datasets");
      
      for (int i=0; i<ndatasets; i++)
      {
        char a[100];
        sprintf(a,"dataset_%d",i);
        read_dataset(state.openGroup(a));
      }
    }

    void read_geometry(const Group &geometry)
    {
      int typ,dum;
      char name[100];
      typ=read_int(geometry,"type");
      switch (typ)
      {
        case 1: break;
        default : mythrow(__LINE__ << ": Unknown type " << typ << " in geometry");
      }
      dim=read_int(geometry,"dimension");

      nvertices=read_int(geometry,"number of vertices");

      nregions=read_int(geometry,"number of regions");

      dum=read_int(geometry,"number of states");
      if (dum!=1)
        mythrow("Number of states not one");

      const Group &trans=geometry.openGroup("transformation");
      read_transformation(trans);
      const DataSet &vert=geometry.openDataSet("vertex");
      read_vertex(vert);
      
      for (int i=0; i<nregions; i++)
      {
        sprintf(name,"region_%d",i);
        const Group &reg=geometry.openGroup(name);
        read_region(i,reg);
      }
      
      read_attribs0(geometry.openGroup("state_0"));
      extrude_contacts_quantity_fields();
    }

    void read_collection(const Group &collection)
    {
      int i;

      i=read_int(collection,"number of geometries");
      if (i!=1)
        mythrow("Not only one geometry");
      i=read_int(collection,"number of plots");
      if (i!=0)
        fprintf(stderr,"We have plots, skip them\n");

      read_geometry(collection.openGroup("geometry_0"));
    }


    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1)
    {
      PointT line = p1-p0;
      std::swap(line[0], line[1]);
      line[0] = -line[0];
      return line;
    }

    template<typename PointT>
    PointT normal_vector(PointT const & p0, PointT const & p1, PointT const & p2)
    {
      return viennagrid::cross_prod( p1-p0, p2-p0 );
    }
    
    void extrude_element(element_t & element)
    {
      PointType center;
      PointType normal;
      double size = 0;

      viennagrid_element_type contact_tag = VIENNAGRID_ELEMENT_TYPE_NO_ELEMENT;
      if (element.element_type == VIENNAGRID_ELEMENT_TYPE_LINE)
      {
        PointType p0 = viennagrid::get_point( vertices[element.vertices[0]] );
        PointType p1 = viennagrid::get_point( vertices[element.vertices[1]] );

        center = (p0+p1)/2;
        normal = normal_vector(p0, p1);

        size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
        contact_tag = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
      }
      else if (element.element_type == VIENNAGRID_ELEMENT_TYPE_TRIANGLE)
      {
        PointType p0 = viennagrid::get_point( vertices[element.vertices[0]] );
        PointType p1 = viennagrid::get_point( vertices[element.vertices[1]] );
        PointType p2 = viennagrid::get_point( vertices[element.vertices[2]] );

        center = (p0+p1+p2)/3;
        normal = normal_vector(p0, p1, p2);

        size = std::max(viennagrid::distance(center, p0), viennagrid::distance(center, p1));
        size = std::max(size, viennagrid::distance(center, p2));
        
        contact_tag = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
      }
      else
      {
        std::cout << "NOT SUPPORTED" << std::endl;
      }

      normal.normalize();
      normal *= size * extrude_contacts_scale;

      PointType new_vertex = center + normal;

      typedef typename viennagrid::result_of::cell_range<MeshT>::type CellRangeType;
      typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;

      CellRangeType cells(mesh);
      CellIteratorType cit = cells.begin();
      for (; cit != cells.end(); ++cit)
      {
        if ( viennagrid::is_inside(*cit, new_vertex) )
          break;
      }

      if (cit != cells.end())
        new_vertex = center - normal;

      VertexType v = viennagrid::make_vertex(mesh, new_vertex);
      
      vertices.push_back(v);
      newly_created_vertices[v.id()] = element.vertices;
      
      element.vertices.push_back(v.id());
      element.element_type = contact_tag;
    }

    /* sets quantity field values for vertices and elements created when extruding contacts */
    void extrude_contacts_quantity_fields()
    {
      for (std::map<std::string, viennagrid::quantity_field>::iterator it = quantity_fields_.begin();
                                                                       it != quantity_fields_.end();
                                                                     ++it)
      {
        //TODO: check quantity field type -> quantites on vertices or elements?
        viennagrid::quantity_field & quantities = it->second;
        
        for (std::map< int, std::vector<int> >::const_iterator nvit = newly_created_vertices.begin();
                                                                 nvit != newly_created_vertices.end();
                                                               ++nvit)
        {
          double val = 0.0;
          for (std::size_t i = 0; i != (*nvit).second.size(); ++i)
            val += quantities.get( (*nvit).second[i] );
          val /= (*nvit).second.size();

          quantities.set( (*nvit).first, val );
        }
      }
    }

    std::vector<viennagrid::quantity_field> quantity_fields()
    {
      std::vector<viennagrid::quantity_field> results;
      for (std::map<std::string, viennagrid::quantity_field>::iterator it = quantity_fields_.begin();
                                                                       it != quantity_fields_.end();
                                                                     ++it)
      {
        results.push_back( it->second );
      }

      return results;
    }

  };

}

#endif
