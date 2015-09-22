#ifndef VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_WRITER_HPP
#define VIENNAMESH_ALGORITHM_IO_SENTAURUS_TDR_WRITER_HPP

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
#include <vector>

#include "viennautils/exception.hpp"

#include "viennagridpp/mesh/mesh.hpp"
#include "viennagridpp/quantity_field.hpp"

namespace tdr
{
  struct writer_error : virtual viennautils::exception {};
  
  void write_to_tdr(std::string const & filename, viennagrid::const_mesh const & mesh, std::vector<viennagrid::quantity_field> const & quantities);
}

#endif
