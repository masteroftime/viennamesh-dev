#ifndef VIENNAMESH_CORE_PARAMETER_HPP
#define VIENNAMESH_CORE_PARAMETER_HPP

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

#include <vector>

#include "viennamesh/forwards.hpp"
#include "viennamesh/utils/logger.hpp"


namespace viennamesh
{

  template<typename ValueT>
  struct type_information
  {
    // default implementation is empty! spezialise this class if you want something special
    static void init() {}
    static std::string name() { return typeid(ValueT).name(); }
  };


  class base_parameter;

  template<typename T>
  class parameter_wrapper;




  template<typename ObjectT>
  class parameter_handle_t : public shared_ptr<ObjectT>
  {
  public:
    parameter_handle_t() {}

    template<typename AnotherObjectT>
    parameter_handle_t( parameter_handle_t<AnotherObjectT> const & rhs ) : shared_ptr<ObjectT>(rhs) {}

    parameter_handle_t( shared_ptr<ObjectT> const & ptr ) : shared_ptr<ObjectT>(ptr) {}

    explicit parameter_handle_t( ObjectT * ptr ) : shared_ptr<ObjectT>(ptr) {}
  };

  template<typename ValueT>
  class parameter_handle_t< parameter_wrapper<ValueT> > : public shared_ptr< parameter_wrapper<ValueT> >
  {
    typedef parameter_wrapper<ValueT> ObjectType;
  public:
    parameter_handle_t() {}
    parameter_handle_t( shared_ptr<ObjectType> const & ptr ) : shared_ptr<ObjectType>(ptr) {}

    explicit parameter_handle_t( ObjectType * ptr ) : shared_ptr<ObjectType>(ptr) {}

    ValueT & operator()() { return (*this)->value(); }
    ValueT const & operator()() const { return (*this)->value(); }
  };

  template<typename ValueT>
  class parameter_handle_t< const parameter_wrapper<ValueT> > : public shared_ptr<const parameter_wrapper<ValueT> >
  {
    typedef parameter_wrapper<ValueT> NonConstObjectType;
    typedef const parameter_wrapper<ValueT> ObjectType;
  public:
    parameter_handle_t() {}

    parameter_handle_t( shared_ptr<ObjectType> const & ptr ) : shared_ptr<ObjectType>(ptr) {}
    parameter_handle_t( shared_ptr<NonConstObjectType> const & ptr ) : shared_ptr<ObjectType>(ptr) {}

    parameter_handle_t & operator=( parameter_handle_t<ObjectType> const & ptr )
    {
      shared_ptr<const parameter_wrapper<ValueT> >::operator=(ptr);
      return *this;
    }
    parameter_handle_t & operator=( parameter_handle_t<NonConstObjectType> const & ptr )
    {
      shared_ptr<const parameter_wrapper<ValueT> >::operator=(ptr);
      return *this;
    }


    explicit parameter_handle_t( ObjectType * ptr ) : shared_ptr<ObjectType>(ptr) {}

    ValueT const & operator()() { return (*this)->value(); }
    ValueT const & operator()() const { return (*this)->value(); }
  };



  typedef parameter_handle_t<base_parameter> parameter_handle;
  typedef parameter_handle_t<const base_parameter> const_parameter_handle;

  namespace result_of
  {
    template<typename ValueT>
    struct parameter_handle
    {
      typedef parameter_handle_t< parameter_wrapper<ValueT> > type;
    };

    template<typename ValueT>
    struct const_parameter_handle
    {
      typedef parameter_handle_t< const parameter_wrapper<ValueT> > type;
    };


    template<typename ValueT>
    struct parameter_handle<const ValueT>
    {
      typedef typename const_parameter_handle<ValueT>::type type;
    };
  }

  typedef result_of::parameter_handle<bool>::type bool_parameter_handle;
  typedef result_of::const_parameter_handle<bool>::type const_bool_parameter_handle;

  typedef result_of::parameter_handle<int>::type int_parameter_handle;
  typedef result_of::const_parameter_handle<int>::type const_int_parameter_handle;

  typedef result_of::parameter_handle<double>::type double_parameter_handle;
  typedef result_of::const_parameter_handle<double>::type const_double_parameter_handle;

  typedef result_of::parameter_handle<std::string>::type string_parameter_handle;
  typedef result_of::const_parameter_handle<std::string>::type const_string_parameter_handle;


  namespace detail
  {
    template<typename ValueT, typename SourceT>
    struct dynamic_handle_cast_impl
    {
      static typename result_of::parameter_handle<ValueT>::type cast( shared_ptr<SourceT> const & ptr )
      {
        return dynamic_pointer_cast< parameter_wrapper<ValueT> >(ptr);
      }
    };

    template<typename ValueT, typename SourceT>
    struct dynamic_handle_cast_impl<const ValueT, SourceT>
    {
      static typename result_of::const_parameter_handle<ValueT>::type cast( shared_ptr<SourceT> const & ptr )
      {
        return dynamic_pointer_cast< const parameter_wrapper<ValueT> >(ptr);
      }
    };
  }


  template<typename ValueT, typename SourceT>
  typename result_of::parameter_handle<ValueT>::type dynamic_handle_cast( shared_ptr<SourceT> const & ptr )
  { return detail::dynamic_handle_cast_impl<ValueT, SourceT>::cast(ptr); }

  template<typename ValueT>
  typename result_of::parameter_handle<ValueT>::type make_parameter()
  { return typename result_of::parameter_handle<ValueT>::type( new parameter_wrapper<ValueT>() ); }

  template<typename ValueT>
  typename result_of::parameter_handle<ValueT>::type make_parameter( ValueT const & value )
  { return typename result_of::parameter_handle<ValueT>::type( new parameter_wrapper<ValueT>(value) ); }

  template<typename ValueT>
  typename result_of::parameter_handle<ValueT>::type make_reference_parameter( ValueT & value )
  { return typename result_of::parameter_handle<ValueT>::type( new parameter_wrapper<ValueT>(&value) ); }

  class type_properties
  {
  public:

    typedef std::map<std::string, std::string> PropertyMapType;


    PropertyMapType & get_map(type_info_wrapper const & ti)
    {
      return type_properties_map[ti];
    }

    std::pair<std::string, bool> get_property( type_info_wrapper const & ti, std::string const & name ) const;

    void set_property( type_info_wrapper const & ti, std::string const & key, std::string const & value )
    {
      get_map(ti)[key] = value;
    }

    template<typename ValueT>
    void set_property( std::string const & key, std::string const & value )
    {
      set_property( typeid( parameter_wrapper<ValueT> ), key, value );
    }

    bool match_properties( type_info_wrapper const & ti, PropertyMapType const & to_match ) const;
    bool match_property( type_info_wrapper const & ti, std::string const & key, std::string const & value ) const;

    void print() const;


    static type_properties & get()
    {
      static type_properties type_properties_;
      return type_properties_;
    }

  private:
    std::map<type_info_wrapper, PropertyMapType> type_properties_map;
  };









  class base_conversion_function
  {
  public:
    base_conversion_function(int function_depth_) : function_depth(function_depth_) {}
    virtual ~base_conversion_function() {}

    virtual bool convert( const_parameter_handle const &, parameter_handle const & ) const = 0;
    virtual parameter_handle get_converted( const_parameter_handle const & ) const = 0;

    virtual type_info_wrapper input_type() const = 0;
    virtual type_info_wrapper output_type() const = 0;

    int function_depth;
  };

  template<typename InputValueT, typename OutputValueT>
  class conversion_function : public base_conversion_function
  {
  public:
    typedef parameter_wrapper<InputValueT> InputParameterType;
    typedef parameter_wrapper<OutputValueT> OutputParameterType;

    conversion_function() : base_conversion_function(1) {}
    conversion_function( function<bool (InputValueT const &, OutputValueT &)> const & f) :
        base_conversion_function(1), convert_function(f) {}

    virtual bool convert( const_parameter_handle const & input, parameter_handle const & output ) const
    {
#ifdef DEBUG
      return convert_function( dynamic_cast<InputParameterType const &>(*input).value(), dynamic_cast<OutputParameterType &>(*output).value() );
#else
      return convert_function( static_cast<InputParameterType const &>(*input).value(), static_cast<OutputParameterType &>(*output).value() );
#endif
    }

    virtual parameter_handle get_converted( const_parameter_handle const & input ) const
    {
      parameter_handle_t<OutputParameterType> result( new OutputParameterType() );
      if (convert(input, result))
        return result;
      return parameter_handle_t<OutputParameterType>();
    }

    virtual type_info_wrapper input_type() const
    { return type_info_wrapper::make<InputParameterType>(); }

    virtual type_info_wrapper output_type() const
    { return type_info_wrapper::make<OutputParameterType>(); }


    function<bool (InputValueT const &, OutputValueT &)> convert_function;
  };

  class dual_conversion_function : public base_conversion_function
  {
  public:
    dual_conversion_function( parameter_handle_t<base_conversion_function> const & first_, parameter_handle_t<base_conversion_function> const & second_ ) :
        base_conversion_function(first_->function_depth + second_->function_depth), first(first_), second(second_) {}

    virtual bool convert( const_parameter_handle const & input, parameter_handle const & output ) const
    {
      return second->convert( first->get_converted(input), output );
    }

    virtual parameter_handle get_converted( const_parameter_handle const & input ) const
    {
      return second->get_converted( first->get_converted(input) );
    }

    virtual type_info_wrapper input_type() const
    { return first->input_type(); }

    virtual type_info_wrapper output_type() const
    { return second->output_type(); }

    shared_ptr<base_conversion_function> first;
    shared_ptr<base_conversion_function> second;
  };





  class converter
  {
  public:

    typedef std::map<type_info_wrapper, shared_ptr<base_conversion_function> > ConversionFunctionMapType;
    typedef std::map<type_info_wrapper, ConversionFunctionMapType> ConversionFunctionMapMapType;

    static type_info_wrapper get_type_id( const_parameter_handle const & tmp );

    shared_ptr<base_conversion_function> convert_function( const_parameter_handle const & input, const_parameter_handle const & output );
    template<typename ValueT>
    shared_ptr<base_conversion_function> convert_function( const_parameter_handle const & input );

    shared_ptr<base_conversion_function> best_convert_function( const_parameter_handle const & input, std::map<std::string, std::string> const & properties );
    shared_ptr<base_conversion_function> best_convert_function( const_parameter_handle const & input, std::string const & property_key, std::string const & property_value );


    bool is_convertable( const_parameter_handle const & input, const_parameter_handle const & output )
    { return (get_type_id(input) == get_type_id(output)) || convert_function(input, output); }
    template<typename ValueT>
    bool is_convertable( const_parameter_handle const & input )
    { return (get_type_id(input) == type_info_wrapper::make< parameter_wrapper<ValueT> >()) ||  convert_function<ValueT>(input); }

    bool convert( const_parameter_handle const & input, parameter_handle const & output )
    {
      shared_ptr<base_conversion_function> cf = convert_function(input, output);
      if (cf && cf->convert(input, output))
      {
        return true;
      }
      else
        return false;
    }
    template<typename ValueT>
    typename result_of::parameter_handle<ValueT>::type get_converted(const_parameter_handle const & input)
    {
      shared_ptr<base_conversion_function> cf = convert_function<ValueT>(input);
      if (cf)
        return static_pointer_cast< parameter_wrapper<ValueT> >(cf->get_converted(input));
      else
        return typename result_of::parameter_handle<ValueT>::type();
    }


    void print_conversions(const_parameter_handle const & input) const;


    template<typename InputValueT, typename OutputValueT>
    void register_conversion( bool (*fp)(InputValueT const &, OutputValueT &) )
    {
      typedef parameter_wrapper<InputValueT> InputParameterType;
      typedef parameter_wrapper<OutputValueT> OutputParameterType;

      type_info_wrapper input_type_id(typeid(InputParameterType));
      type_info_wrapper output_type_id(typeid(OutputParameterType));

      shared_ptr<base_conversion_function> current_conversion(new conversion_function<InputValueT, OutputValueT>(fp));
      simple_register_conversion( input_type_id, output_type_id, current_conversion );

      for (ConversionFunctionMapMapType::iterator imtit = conversions.begin();
           imtit != conversions.end();
           ++imtit)
      {
        ConversionFunctionMapType::iterator omtit = imtit->second.find(input_type_id);
        if (omtit != imtit->second.end())
        {
          simple_register_conversion(imtit->first, output_type_id, shared_ptr<base_conversion_function>(new dual_conversion_function(omtit->second, current_conversion)));
        }
      }

      ConversionFunctionMapMapType::iterator imtit = conversions.find(output_type_id);
      if (imtit != conversions.end())
      {
        for (ConversionFunctionMapType::iterator omtit = imtit->second.begin();
            omtit != imtit->second.end();
            ++omtit)
        {
          simple_register_conversion(input_type_id, omtit->first, shared_ptr<base_conversion_function>(new dual_conversion_function(current_conversion, omtit->second)));
        }
      }
    }

    static converter & get()
    {
      static converter converter_;
      return converter_;
    }

  private:

    void simple_register_conversion( type_info_wrapper const & input_type_id, type_info_wrapper const & output_type_id, shared_ptr<base_conversion_function> const & conversion )
    {
      shared_ptr<base_conversion_function> & entry = conversions[input_type_id][output_type_id];
      if (!entry)
        entry = conversion;
      else
      {
        if (entry->function_depth >= conversion->function_depth)
          entry = conversion;
      }
    }

    ConversionFunctionMapMapType conversions;
  };



  template<typename T>
  class parameter_wrapper;

  class base_parameter : public enable_shared_from_this<base_parameter>
  {
  public:

    virtual ~base_parameter() {}

    virtual parameter_handle unpack() = 0;
    virtual const_parameter_handle unpack() const = 0;
    virtual bool is_reference() const = 0;
    virtual std::string type_name() const = 0;


    std::pair<std::string, bool> get_property( std::string const & key ) const
    { return type_properties::get().get_property( typeid(*this), key ); }

    bool match_property( std::string const & key, std::string const & value ) const
    { return type_properties::get().match_property( typeid(*this), key, value ); }

    bool match_properties( std::map<std::string, std::string> const & properties ) const
    { return type_properties::get().match_properties( typeid(*this), properties ); }


    template<typename ValueT>
    bool is_type() const
    {
      return typeid(*unpack()) == typeid(parameter_wrapper<ValueT>);
    }

    template<typename ValueT>
    bool is_convertable_to() const
    {
      bool is_convertable = converter::get().is_convertable<ValueT>( shared_from_this() );

      info(10) << "is_convertable (" << (is_convertable ? "TRUE" : "FALSE")  << "): \"" << type_name() << "\" to \"" << type_information<ValueT>::name() << "\"" << std::endl;

      return is_convertable;
    }

    template<typename ValueT>
    typename result_of::parameter_handle<ValueT>::type get_converted() const
    {
      LoggingStack stack( std::string("get_converted") );
      info(5) << "Source type: " << type_name() << std::endl;
      info(5) << "Destination type: " << type_information<ValueT>::name() << std::endl;

      typename result_of::parameter_handle<ValueT>::type result_parameter = converter::get().get_converted<ValueT>( shared_from_this() );

      info(5) << "Success: " << std::boolalpha << static_cast<bool>(result_parameter) << std::endl;

      return result_parameter;

    }
  };

  inline bool is_convertable( const_parameter_handle const & source, parameter_handle & destination )
  { return converter::get().is_convertable( source, destination ); }

  inline bool convert( const_parameter_handle const & source, parameter_handle & destination )
  { return converter::get().convert( source, destination ); }



  template<typename InputValueT, typename OutputValueT>
  bool static_cast_convert( InputValueT const & input, OutputValueT & output )
  {
    output = static_cast<OutputValueT>(input);
    return true;
  }



  template<typename ValueT>
  class parameter_wrapper : public base_parameter
  {
  public:
    typedef ValueT value_type;

    parameter_wrapper() { static_init(); clear(); }
    parameter_wrapper(value_type const & val) { static_init(); set(val); }
    parameter_wrapper(value_type * val_ref) { static_init(); set(val_ref);}

    parameter_handle unpack() { return shared_from_this(); }
    const_parameter_handle unpack() const { return shared_from_this(); }

    bool is_reference() const { return value_ptr_ != 0; }

    std::string type_name() const { return type_information<ValueT>::name(); }

    static void static_init()
    {
      static bool to_init = true;
      if (to_init)
      {
        info(10) << "static_init< " << type_information<ValueT>::name() << " >::init" << std::endl;
        viennamesh::type_information<value_type>::init();
        to_init = false;
      }
    }

  private:

    void clear() { value_ptr_ = NULL; }
    void set( ValueT * val_ref ) { clear(); value_ptr_ = val_ref; }
    void set( ValueT const & val ) { clear(); value_ = val; }

    template<typename TypeT>
    friend class parameter_handle_t;

    template<typename InputValueT, typename OutputValueT>
    friend class conversion_function;

    ValueT & value() { if (value_ptr_) return *value_ptr_; else return value_; }
    ValueT const & value() const { if (value_ptr_) return *value_ptr_; else return value_; }

    ValueT value_;
    ValueT * value_ptr_;
  };



  class parameter_handle_reference : public base_parameter
  {
  public:

    parameter_handle_reference() {}
    parameter_handle_reference(parameter_handle const & handle_reference_) : handle_reference(handle_reference_) {}

    parameter_handle unpack() { return handle_reference; }
    const_parameter_handle unpack() const { return handle_reference; }

    bool is_reference() const { return true; }

    std::string type_name() const { return "parameter_handle_reference"; }

  private:
    parameter_handle handle_reference;
  };



  template<typename ValueT>
  shared_ptr<base_conversion_function> converter::convert_function( const_parameter_handle const & input )
  {
    type_information<ValueT>::init();

    ConversionFunctionMapMapType::iterator ipit = conversions.find(typeid(*input));
    if (ipit != conversions.end())
    {
      ConversionFunctionMapType::iterator opit = ipit->second.find(typeid(parameter_wrapper<ValueT>));
      if ( opit != ipit->second.end() )
      {
        return opit->second;
      }
    }

    return shared_ptr<base_conversion_function>();
  }


  class parameter_set
  {
  public:

    typedef std::map<std::string, parameter_handle> ParameterMapType;

    void set( std::string const & name, parameter_handle const & parameter )
    { internal_set(name, parameter); }

    void unset( std::string const & name )
    { parameters.erase(name); }

    const_parameter_handle get( std::string const & name ) const
    {
      ParameterMapType::const_iterator it = parameters.find(name);
      if (it == parameters.end())
        return parameter_handle();
      return it->second->unpack();
    }

    parameter_handle get( std::string const & name )
    {
      ParameterMapType::iterator it = parameters.find(name);
      if (it == parameters.end())
        return parameter_handle();
      return it->second->unpack();
    }

    template<typename ValueT>
    typename result_of::const_parameter_handle<ValueT>::type get( std::string const & name ) const
    {
      const_parameter_handle handle = get(name);
      if (!handle) return typename result_of::const_parameter_handle<ValueT>::type();
      typename result_of::const_parameter_handle<ValueT>::type result = dynamic_handle_cast<const ValueT>(handle);

      if (result)
        return result;

      return handle->template get_converted<ValueT>();
    }

    template<typename ValueT>
    typename result_of::parameter_handle<ValueT>::type get( std::string const & name )
    {
      parameter_handle handle = get(name);
      if (!handle) return typename result_of::parameter_handle<ValueT>::type();
      typename result_of::parameter_handle<ValueT>::type result = dynamic_handle_cast<ValueT>(handle);

      if (result)
        return result;

      return handle->template get_converted<ValueT>();
    }

    parameter_handle & get_create( std::string const & name )
    { return parameters[name]; }

    template<typename ValueT>
    bool copy_if_present( std::string const & name, ValueT & target ) const
    {
      typename result_of::const_parameter_handle<ValueT>::type ptr = get<ValueT>(name);
      if (ptr)
      {
        target = ptr->value;
        return true;
      }
      return false;
    }

    void clear()
    { parameters.clear(); }

    void clear_non_references()
    {
      ParameterMapType::iterator it = parameters.begin();
      while (it != parameters.end())
      {
        if (!it->second->is_reference())
          parameters.erase(it++);
        else
          ++it;
      }
    }

  private:
    void internal_set( std::string const & name, parameter_handle const & parameter )
    {
      if (parameter)
        parameters[name] = parameter;
      else
        unset(name);
    }

    ParameterMapType parameters;
  };


  class const_parameter_set
  {
  public:

    typedef std::map<std::string, const_parameter_handle> ParameterMapType;

    void set( std::string const & name, const_parameter_handle const & parameter )
    { internal_set(name, parameter); }

    void set( std::string const & name, parameter_handle const & parameter )
    { internal_set(name, parameter); }

    void unset( std::string const & name )
    { parameters.erase(name); }

    const_parameter_handle get( std::string const & name ) const
    {
      ParameterMapType::const_iterator it = parameters.find(name);
      if (it == parameters.end())
        return const_parameter_handle();
      return it->second->unpack();
    }

    template<typename ValueT>
    typename result_of::const_parameter_handle<ValueT>::type get( std::string const & name ) const
    {
      const_parameter_handle handle = get(name);
      if (!handle) return typename result_of::const_parameter_handle<ValueT>::type();
      typename result_of::const_parameter_handle<ValueT>::type result = dynamic_handle_cast<const ValueT>(handle);

      if (result)
        return result;

      return handle->template get_converted<ValueT>();
    }

    const_parameter_handle & get_create( std::string const & name )
    { return parameters[name]; }

    template<typename ValueT>
    bool copy_if_present( std::string const & name, ValueT & target ) const
    {
      typename result_of::const_parameter_handle<ValueT>::type handle = get<ValueT>(name);
      if (handle)
      {
        target = handle();
        return true;
      }
      return false;
    }

    void clear()
    { parameters.clear(); }

  private:

    void internal_set( std::string const & name, const_parameter_handle const & parameter )
    {
      if (parameter)
        parameters[name] = parameter;
      else
        unset(name);
    }

    ParameterMapType parameters;
  };

}


#endif
