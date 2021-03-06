#ifndef VIENNAMESH_CORE_EXCEPTIONS_HPP
#define VIENNAMESH_CORE_EXCEPTIONS_HPP

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

#include <exception>
#include <stdexcept>

namespace viennamesh
{
  class input_parameter_not_found_exception : public std::runtime_error
  {
  public:
    input_parameter_not_found_exception(std::string const & message_) : std::runtime_error(message_) {}
    virtual ~input_parameter_not_found_exception() throw() {}
  };

  class output_not_convertable_to_referenced_value_exception : public std::runtime_error
  {
  public:

    output_not_convertable_to_referenced_value_exception(std::string const & message_) : std::runtime_error(message_) {}
    virtual ~output_not_convertable_to_referenced_value_exception() throw() {}
  };

  class metric_not_implemented_or_supported_exception : public std::runtime_error
  {
  public:

    metric_not_implemented_or_supported_exception(std::string const & message_) : std::runtime_error(message_) {}
    virtual ~metric_not_implemented_or_supported_exception() throw() {}
  };

  class interface_check_failed_exception : public std::runtime_error
  {
  public:

    interface_check_failed_exception(std::string const & message_) : std::runtime_error(message_) {}
    virtual ~interface_check_failed_exception() throw() {}
  };

  class create_sizing_function_exception : public std::runtime_error
  {
  public:
    create_sizing_function_exception(std::string const & message_) : std::runtime_error(message_) {}
    virtual ~create_sizing_function_exception() throw() {}
  };
}

#endif
