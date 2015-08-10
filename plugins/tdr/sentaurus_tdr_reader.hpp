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
    viennagrid_element_type element_tag;
    std::vector<int> vertex_indices;
  };

  struct region_t
  {
    int regnr;
    string name,material;
    int nelements,npointidx;
    std::vector<element_t> elements;
  };


  template<typename VertexT>
  struct region_contacts
  {
    string region_name;
    std::vector<element_t> elements;
  };

  template<typename MeshT>
  struct tdr_geometry
  {
    typedef typename viennagrid::result_of::point<MeshT>::type PointType;
    typedef typename viennagrid::result_of::element<MeshT>::type VertexType;
    typedef typename viennagrid::result_of::region<MeshT>::type RegionType;
    
    unsigned int nvertices;
    int dim,nregions,ndatasets;
    std::map<string,region_t> region;
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

      /*
      * Read two fields x and y from s2 dataset. Fields in the file
      * are found by their names "x_name" and "y_name".
      */
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

    void read_elements(region_t &region, const DataSet &elem)
    {
      hsize_t size = get_dataset_size(elem);

      boost::scoped_array<int> el(new int[size]);
      elem.read( el.get(), PredType::NATIVE_INT);

      unsigned int elct=0;
      while (elct < size)
      {
        element_t element;
        switch (el[elct++])
        {
          case 1:
            element.element_tag = VIENNAGRID_ELEMENT_TYPE_LINE;
            for (int i=0; i<2; i++)
              element.vertex_indices.push_back(el[elct++]);
            break;
          case 2:
            element.element_tag = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
            for (int i=0; i<3; i++)
              element.vertex_indices.push_back(el[elct++]);
            break;
          case 3:
            element.element_tag = VIENNAGRID_ELEMENT_TYPE_QUADRILATERAL;
            for (int i=0; i<4; i++)
              element.vertex_indices.push_back(el[elct++]);
            break;
          case 5:
            element.element_tag = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
            for (int i=0; i<4; i++)
              element.vertex_indices.push_back(el[elct++]);
            break;
          default:mythrow("Element type " << el[elct-1] << " in region " << region.name << " not known");
        }

        region.elements.push_back(element);
      }
    }

    void read_region(const int regnr, const Group &reg)
    {
      string name=read_string(reg,"name");

      string material;
      std::stringstream ss;
      const int typ=read_int(reg,"type");
      switch (typ)
      {
        case 0: // we have a nD region
          material=read_string(reg,"material");
          break;
        case 1: // we have a nD-1 region (contact... drain, source, ...)
          material="Contact";
          // attribute "part 0", no idea for what
          break;
        case 2:
          ss << "Interface" << read_int(reg,"bulk 0") << read_int(reg,"bulk 1");
          material = ss.str();
          break;
        //default : mythrow(__LINE__ << ": Unknown type " << typ << " in region " << read_string(reg,"name"));
        default: std::cerr << __LINE__ << ": Unknown type " << typ << " in region " << read_string(reg,"name") << std::endl;
         std::cerr << "I do not throw, I only skip\n";
          return;
      }

      region[name].regnr=regnr;
      region[name].name=name;
      region[name].material=material;
      const int n=reg.getNumObjs();
      for (int i=0; i<n; i++)
      {
        std::string objname = reg.getObjnameByIdx(i);
        if (objname.find("elements_") == 0)
        {
          const DataSet &ds=reg.openDataSet(objname);
          region[name].nelements=read_int(ds,"number of elements");
          read_elements(region[name],ds);
        }

        if (objname.find("part_") == 0)
        {
          const Group &part=reg.openGroup(objname);
          region[name].nelements=read_int(part,"number of elements");
          read_elements(region[name],part.openDataSet("elements"));
        }
      }

    }

    region_t &find_region(int regnr)
    {
      for (std::map<string,region_t>::iterator R=region.begin(); R!=region.end(); R++)
        if (R->second.regnr==regnr)
          return R->second;
      mythrow("Region " << regnr << " not found");
    }

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
      if (name.find("Stress")!=name.npos)
        return;

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

      string unit;
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
      }

      /*region_t &region=find_region(regnr);
      region.dataset[name].name=name;
      region.dataset[name].quantity=quantity;
      region.dataset[name].nvalues=nvalues;
      region.dataset[name].conversion_factor=conversion_factor;
      region.dataset[name].unit=unit;*/
      
      region_t &reg=find_region(regnr);
      RegionType region = mesh.get_region(reg.name);
      
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

      for (int i=0; i<nregions; i++)
      {
        sprintf(name,"region_%d",i);
        const Group &reg=geometry.openGroup(name);
        read_region(i,reg);
      }

      const Group &trans=geometry.openGroup("transformation");
      read_transformation(trans);
      const DataSet &vert=geometry.openDataSet("vertex");
      read_vertex(vert);
      // What do the Units???
      
      to_viennagrid();
      
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

    bool to_viennagrid()
    {
      viennagrid_element_type cell_type = VIENNAGRID_ELEMENT_TYPE_VERTEX;
      for (std::map<string,region_t>::iterator S=region.begin(); S!=region.end(); S++)
      {
        std::vector<element_t> &elements=S->second.elements;
        for (std::vector<element_t>::iterator E=elements.begin(); E!=elements.end(); E++)
        {
          element_t const & e = *E;
          cell_type = viennagrid_topological_max( cell_type, e.element_tag );
        }
      }


      std::map<string, region_contacts<VertexType> > contact_elements;

      for (std::map<string,region_t>::iterator S=region.begin(); S!=region.end(); S++)
      {
        string region_name = S->second.name;

        std::vector<element_t> &elements=S->second.elements;
        for (std::vector<element_t>::iterator E=elements.begin(); E!=elements.end(); E++)
        {
          element_t const & e = *E;

          std::vector<VertexType> cell_vertices(e.vertex_indices.size());
          for (std::size_t i = 0; i < e.vertex_indices.size(); ++i)
            cell_vertices[i] = vertices[e.vertex_indices[i]];

          if (e.element_tag == cell_type)
          {
            viennagrid::make_element( mesh.get_or_create_region(region_name),
                                      viennagrid::element_tag::from_internal(e.element_tag),
                                      cell_vertices.begin(), cell_vertices.end() );
          }
          else
          {
            contact_elements[region_name].region_name = region_name + "_contact";
            contact_elements[region_name].elements.push_back(e);
          }
        }
      }

      if (extrude_contacts)
      {
        for (typename std::map<string, region_contacts<VertexType> >::iterator rc = contact_elements.begin(); rc != contact_elements.end(); ++rc)
        {
          for (std::size_t i = 0; i < rc->second.elements.size(); ++i)
          {
            element_t const & element = rc->second.elements[i];
//             std::vector<VertexType> handles = rc->second.vertex_handles[i];

            PointType center;
            PointType normal;
            double size = 0;

            viennagrid_element_type contact_tag = VIENNAGRID_ELEMENT_TYPE_NO_ELEMENT;
            if (element.element_tag == VIENNAGRID_ELEMENT_TYPE_LINE)
            {
              PointType p0 = viennagrid::get_point( vertices[element.vertex_indices[0]] );
              PointType p1 = viennagrid::get_point( vertices[element.vertex_indices[1]] );

              center = (p0+p1)/2;
              normal = normal_vector(p0, p1);

              size = std::max(size, viennagrid::distance(center, p0));
              size = std::max(size, viennagrid::distance(center, p1));
              contact_tag = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
            }
            else if (element.element_tag == VIENNAGRID_ELEMENT_TYPE_TRIANGLE)
            {
              PointType p0 = viennagrid::get_point( vertices[element.vertex_indices[0]] );
              PointType p1 = viennagrid::get_point( vertices[element.vertex_indices[1]] );
              PointType p2 = viennagrid::get_point( vertices[element.vertex_indices[2]] );

              center = (p0+p1+p2)/3;
              normal = normal_vector(p0, p1, p2);

              size = std::max(size, viennagrid::distance(center, p0));
              size = std::max(size, viennagrid::distance(center, p1));
              size = std::max(size, viennagrid::distance(center, p2));
              contact_tag = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
            }
            else
            {
              std::cout << "NOT SUPPORTED" << std::endl;
            }

            normal /= viennagrid::norm_2(normal);
            normal *= size * extrude_contacts_scale;

            PointType other_vertex = center + normal;

            typedef typename viennagrid::result_of::cell_range<MeshT>::type CellRangeType;
            typedef typename viennagrid::result_of::iterator<CellRangeType>::type CellIteratorType;

            CellRangeType cells(mesh);
            CellIteratorType cit = cells.begin();
            for (; cit != cells.end(); ++cit)
            {
              if ( viennagrid::is_inside(*cit, other_vertex) )
                break;
            }

            if (cit != cells.end())
              other_vertex = center - normal;

            std::vector<VertexType> cell_vertices;
            std::vector<int> other_vertex_indices;

//             for (auto index : element.vertex_indices)
            for (std::vector<int>::const_iterator it = element.vertex_indices.begin(); it != element.vertex_indices.end(); ++it)
            {
              cell_vertices.push_back( vertices[*it] );
              other_vertex_indices.push_back( vertices[*it].id() );
            }
            cell_vertices.push_back( viennagrid::make_vertex(mesh, other_vertex) );
            newly_created_vertices[ cell_vertices.back().id() ] = other_vertex_indices;

            viennagrid::make_element( mesh.get_or_create_region(rc->second.region_name),
                                      viennagrid::element_tag::from_internal(contact_tag),
                                      cell_vertices.begin(), cell_vertices.end() );
          }
        }
      }


      return true;
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
