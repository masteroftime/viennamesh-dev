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

#include <memory>

#include "tdr_reader.hpp"
#include "viennameshpp/core.hpp"
#include "sentaurus_tdr_reader.hpp"

namespace viennamesh
{
  tdr_reader::tdr_reader() {}
  std::string tdr_reader::name() { return "tdr_reader"; }


  bool tdr_reader::run(viennamesh::algorithm_handle &)
  {
    string_handle filename = get_required_input<string_handle>("filename");

    std::string path = base_path();
    std::string full_filename;

    if (!path.empty() && filename()[0] != '/')
    {
      info(1) << "Using base path: " << path << std::endl;
      full_filename = path + "/" + filename();
    }
    else
      full_filename = filename();

    
    mesh_handle output_mesh = make_data<mesh_handle>();
    tdr_geometry<mesh_handle::CPPResultType> geometry(output_mesh());
    
    if ( get_input<bool>("extrude_contacts").valid() )
      geometry.extrude_contacts = get_input<bool>("extrude_contacts")();
    else
      geometry.extrude_contacts = true;

    if ( get_input<double>("extrude_contacts_scale").valid() )
      geometry.extrude_contacts_scale = get_input<double>("extrude_contacts_scale")();
    else
      geometry.extrude_contacts_scale = 1.0;
    
    if ( get_input<bool>("fill_triangle_contacts").valid() )
      geometry.fill_triangle_contacts = get_input<bool>("fill_triangle_contacts")();
    else
      geometry.fill_triangle_contacts = false;

    try {
      geometry.read_file(full_filename);
    } catch (std::runtime_error &e) {
      error(1) << e.what() << std::endl;
      return false;
    }


    std::vector<viennagrid::quantity_field> quantity_fields = geometry.quantity_fields();
    if (!quantity_fields.empty())
    {
      quantity_field_handle output_quantity_fields = make_data<viennagrid::quantity_field>();
      output_quantity_fields.set(quantity_fields);

      set_output( "quantities", output_quantity_fields );
    }

    set_output("mesh", output_mesh);

    return true;
  }

}
