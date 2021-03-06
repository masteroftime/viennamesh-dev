/* ************* :: Generic Scientific Simulation Environment :: *************
 **  http://www.gsse.at                                                     **

     Copyright (c) 2003-2008 Rene Heinzl                 rene@gsse.at
     Copyright (c) 2004-2008 Philipp Schwaha

     Use, modification and distribution is subject to the Boost Software
     License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
 **  http://www.boost.org/LICENSE_1_0.txt)                                  **
 *************************************************************************** */

#ifndef GSSE_LAMBDA_QUAN_ACCESS_LINEARIZED_HH
#define GSSE_LAMBDA_QUAN_ACCESS_LINEARIZED_HH

// *** system includes   
//
#include <iostream>

// *** BOOST includes   
//
// [RH] adaptation to boost >= 1.37
#include <boost/spirit/home/phoenix.hpp>   

// old boost < 1.37
// #include <boost/spirit/phoenix/actor.hpp>
// #include <boost/spirit/phoenix/actor.hpp>
// #include <boost/spirit/phoenix/core.hpp>
// #include <boost/spirit/phoenix/function.hpp>
// #include <boost/spirit/phoenix/operator.hpp>
// #include <boost/spirit/phoenix/scope.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>

// *** GSSE includes   
//
#include "gsse/math/linear_eqn.hpp"
#include "gsse/lambda/matrix_mapping.hpp"
#include "gsse/domain.hpp"


namespace gsse
{


// matrix insert functionality
//
template <typename LinearizedEquationT, typename InsertFunctor, typename RHSFunctor, typename ElementT, typename MappingT>
void matrix_line_insert(LinearizedEquationT const& eqn, ElementT const& elem, InsertFunctor& insert, RHSFunctor& rhs, MappingT& mapping)
{
   long index1 = mapping(elem);
   
   typename LinearizedEquationT::iterator iter;
   for(iter = eqn.begin(); iter != eqn.end(); ++iter)
   {
      long index2 = (*iter).first;
      double entry = (*iter).second;

      
      

#ifdef DEBUG
       std::cout << "  ## GSSE::mtx insert (i1/i2/entry): " << index1 << " / " << index2 << " / " << entry << std::endl;
#endif

      insert(index1,index2,entry);
   }

#ifdef DEBUG
   std::cout << "  rhs : " << index1 << " / " << eqn()  << std::endl;
#endif

   double rhs_entry = eqn();
   rhs(index1, rhs_entry);
}


// =============================================================================


struct scalar_quantity_entry {};              // used for scalar access

struct trivial_key_tag {};                    // takes the trivial (vertex only) key
struct object_and_quan_key_tag {};            // takes the vertex and the quantity key as key
struct domain_object_and_quan_key_tag {};    // takes the vertex and the quantity key and domain as key



////////////////////////////////////////////////////////////////////////
//
// base quan_access_linearized.. top of derivable types
//

template <typename DomainType, typename ObjectT, typename AccessTag, typename KeyTag>
class quan_access_linearized_impl {};


// ##############################################################################################
// ##############################################################################################
//
// [RH]  .. NEW gsse quan accessor/creator ..
//          modeled by a boost::phoenix::actor
//
// partial specialization of quan_access_linearized_impl with scalar_quan and trivial_key_tag
//
template <typename DomainType, typename ObjectT>
class quan_access_linearized_impl <DomainType, ObjectT, scalar_quantity_entry, trivial_key_tag>
{
   typedef quan_access_linearized_impl <DomainType, ObjectT, scalar_quantity_entry,  trivial_key_tag> self_type;

   typedef ObjectT                                                    object_t;
   typedef typename domain_traits<DomainType>::quan_key_t             quan_key_t;          
   typedef typename domain_traits<DomainType>::storage_type           storage_type; 
   typedef typename storage_type::numeric_t                           numeric_t; 
    
public:

   typedef gsse::matrix_mapping<object_t, self_type>                   mapping_type;
   typedef gsse::linear_equation<numeric_t>                            linear_equation_type;
//   typedef high_performance_matrix_mapping<object_t, self_type>        mapping_type;

public:
   // constructor
   //
   explicit quan_access_linearized_impl(DomainType& domain, const quan_key_t& quan_key, mapping_type& mapping)
      : domain(domain), quan_key(quan_key), mapping(mapping) {}

   // enable nullary/zero arguments for convenient lambda style 
   //
   typedef boost::mpl::false_ no_nullary;
  
   // Boost phoenix/fusion return type calculation 
   //
   // [RH].. new phoenixv2 / fusionv2 version
   //
   template <typename Env>
   struct result 
   { 
      typedef linear_equation_type type;
   };



   // boost::phoenix actor concept --> creates a simple/single linear equation
   //
   template <typename Env>
   typename result<Env>::type eval(Env const& arg) const
      {
         numeric_t res = domain(boost::fusion::at_c<0>(arg.args()), quan_key)(0,0);
         object_t obj = boost::fusion::at_c<0>(arg.args());
         
         return linear_equation_type(mapping(obj), res);
      }



   // matrix mapping functionality
   //
   template <typename ElementT>
   void register_entry(ElementT const& elem) 
   {	
      mapping.add(elem, this);
   }

   template <typename ElementT>
   long get_index(ElementT const& elem)
   {
      return mapping(elem);
   }

   // matrix/value set/get methods
   //
   template <typename ElementT>
   void set_value(ElementT const& elem, numeric_t const value)
   {
      domain(elem, quan_key)(0,0) = value;
   }
  template <typename ElementT>
  void update_value(ElementT const& elem, numeric_t const value)
  {
#ifdef DEBUG
     std::cout << ".. update elem/quan_key/value: " << elem.handle() << "/" << quan_key << "/" << value << std::endl;
#endif
    domain(elem, quan_key)(0,0) -= value;
  }

   // [RH][TODO] .. extracted into new method
   // .... remove this method
   //
   // OlD:: matrix insert functionality
   //
   template <typename InsertFunctor, typename RHSFunctor, typename ElementT>
   void insert(linear_equation_type const& eqn, ElementT const& elem, InsertFunctor& insert, RHSFunctor& rhs)
      {
         long index1 = get_index(elem);
         
         typename linear_equation_type::iterator iter;
         for(iter = eqn.begin(); iter != eqn.end(); ++iter)
         {
            long index2 = (*iter).first;
            double entry = (*iter).second;

	


#ifdef DEBUG
            std::cout << " ## GSSE::in_equ::  insert (i1/i2/entry): " << index1 << " / " << index2 << " / " << entry << std::endl;
#endif
            insert(index1,index2,entry);
         }
         double rhs_entry = eqn();
         rhs(index1, rhs_entry);
      }

      friend std::ostream& operator<<(std::ostream& ostr, self_type const& self)
      {
//         ostr << key.quan << ":" << key.v.handle();
         return ostr;
      }

   
private:
   DomainType&    domain;
   quan_key_t     quan_key;
   mapping_type&  mapping;      
};





// ##############################################################################################
// ##############################################################################################
//
// [RH]  .. NEW gsse quan accessor/creator ..
//          modeled by a boost::phoenix::actor
//
// partial specialization of quan_access_linearized_impl with scalar_quan and object-id/quan-key tag
//
template <typename DomainType, typename ObjectT>
class quan_access_linearized_impl <DomainType, ObjectT, scalar_quantity_entry, object_and_quan_key_tag>
{
   typedef quan_access_linearized_impl <DomainType, ObjectT, scalar_quantity_entry, object_and_quan_key_tag> self_type;

   // the object to be assembled into the matrix, e.g., vertex_type, edge_type, ...
   // 
   typedef ObjectT                                                    object_t;
   typedef typename domain_traits<DomainType>::quan_key_t             quan_key_t;          
   typedef typename domain_traits<DomainType>::storage_type           storage_type; 
   typedef typename storage_type::numeric_t                           numeric_t; 
    


  // -------------------------
  //
  // compound key_type for mapping: use object-id and quantity-key
  
   struct key_type
   {
    friend class quan_access_linearized_impl;

      key_type() {};
      key_type (quan_key_t const& quan, object_t const& object) : object(object), quan(quan) {}

      bool operator==(key_type const& other) const
      {
         return (quan == other.quan) && (object == other.object);
      }

      bool operator<(key_type const& other) const
      {
         if (quan != other.quan) 
            return (quan < other.quan);
         else 
            return (object < other.object);
	//if (v != other.v) return (v < other.v);
	//else return (quan < other.quan);      
      }

      friend std::ostream& operator<<(std::ostream& ostr, key_type const& key)
      {
         ostr << key.quan << ":" << key.object.handle();
         return ostr;
      }

   private:
      object_t   object;
      quan_key_t quan;

   };

public:

   // CLASS changes..
   //
   typedef gsse::matrix_mapping<key_type, self_type>                   mapping_type;
   typedef gsse::linear_equation<numeric_t>                            linear_equation_type;
//   typedef high_performance_matrix_mapping<key_type, self_type>        mapping_type;




public:
   // constructor
   //
   explicit quan_access_linearized_impl(DomainType& domain, const quan_key_t& quan_key, mapping_type& mapping)
      : domain(domain), quan_key(quan_key), mapping(mapping) {}

   // enable nullary/zero arguments for convenient lambda style 
   //
   typedef boost::mpl::false_ no_nullary;
  
   // Boost phoenix/fusion return type calculation 
   //
   template <typename Env>
   struct result
   { 
      typedef linear_equation_type type;
   };


   // boost::phoenix actor concept --> creates a simple/single linear equation
   //
   template <typename Env>
   typename result<Env>::type eval(Env const& arg) const
      {
         numeric_t res  = domain(boost::fusion::at_c<0>(arg.args()), quan_key)(0,0);
         object_t  elem = boost::fusion::at_c<0>(arg.args());
         key_type key(quan_key, elem);

#ifdef DEBUG
         std::cout << "  ##### quan_access__linearized.. build new one (vh, key, mapping): " << elem.handle() << " / " << quan_key << " / " << mapping(key) << std::endl;
#endif         
	
         return linear_equation_type(mapping(key), res);
      }



   // matrix mapping functionality
   //
   template <typename ElementT>
   void register_entry(ElementT const& elem) 
   {
      key_type key(quan_key, elem);
      mapping.add(key, this);
   }

   template <typename ElementT>
   long get_index(ElementT const& elem)
   {
      return mapping(key_type(quan_key, elem));
   }

   // matrix/value set/get methods
   //
   template <typename ElementT>
   void set_value(ElementT const& key, numeric_t const value)
   {
      domain(key.object, key.quan)(0,0) = value;
   }
  template <typename ElementT>
  void update_value(ElementT const& key, numeric_t const value)
  {
#ifdef DEBUG
     std::cout << ".. update elem/quan_key/value: " << key.object.handle() << "/" << key.quan << "/" << value << std::endl;
#endif
    domain(key.object, key.quan)(0,0) -= value;
  }

   
   quan_key_t get_key()
      {
         return quan_key;
      }
   
   template <typename LinearizedEquationT, typename InsertFunctor, typename RHSFunctor>
   void insert(LinearizedEquationT const& eqn, object_t const& elem, InsertFunctor& insert, RHSFunctor& rhs)
      {
         key_type key(quan_key, elem);
         long index1 = mapping(key);
         
         typename LinearizedEquationT::iterator iter;
         for(iter = eqn.begin(); iter != eqn.end(); ++iter)
         {
            long index2  = (*iter).first;
            double entry = (*iter).second;

	
            
#ifdef DEBUG
            std::cout << "  ## GSSE::mtx insert (i1/i2/entry): " << index1 << " / " << index2 << " / " << entry << std::endl;
#endif
            
            insert(index1,index2,entry);
         }
         
#ifdef DEBUG
         std::cout << "  rhs : " << index1 << " / " << eqn()  << std::endl;
#endif
         
         double rhs_entry = eqn();
         rhs(index1, rhs_entry);
}






      friend std::ostream& operator<<(std::ostream& ostr, self_type const& self)
      {
//         ostr << key.quan << ":" << key.v.handle();
         return ostr;
      }

   
private:
   DomainType&    domain;
   quan_key_t     quan_key;
   mapping_type&  mapping;      
};



////////////////////////////////////////////////////////////////////////
//
// partial specialization of quan_access_linearized_impl with scalar_quan and vertex_and_quantities 
//

/*
template <typename DomainType>
class quan_access_linearized_impl <DomainType, scalar_quantity_entry, object_and_quan_key_tag>
{
   typedef quan_access_linearized_impl <DomainType, scalar_quantity_entry, object_and_quan_key_tag> self_type;

   typedef typename domain_traits<DomainType>::quan_key_t      quan_key_t;          
   typedef typename domain_traits<DomainType>::storage_type  storage_type; 
   typedef typename domain_traits<DomainType>::vertex          vertex;
   typedef typename storage_type::numeric_t                     numeric_t; 
    
public:
   typedef typename gsse::linear_equation<numeric_t> linear_equation_type;

   template <typename TupleT>
   struct result { typedef linear_equation_type type;};

  // -------------------------
  //
  // internal data type key_type -> make most general!!! us tuples!!!
  
   struct key_type
   {
      friend class quan_access_linearized_impl;

      key_type() {};
      key_type (quan_key_t const& quan, vertex const& v) : v(v), quan(quan) {}

      bool operator==(key_type const& other) const
      {
         return (quan == other.quan) && (v == other.v);
      }

      bool operator<(key_type const& other) const
      {
	if (quan != other.quan) return (quan < other.quan);
	else return (v < other.v);
	//if (v != other.v) return (v < other.v);
	//else return (quan < other.quan);      
      }

      friend std::ostream& operator<<(std::ostream& ostr, key_type const& key)
      {
         ostr << key.quan << ":" << key.v.handle();
         return ostr;
      }

   private:
      vertex     v;
      quan_key_t quan;

   };
  //
  // -------------------------

   // end of key_type
    
   typedef matrix_mapping<key_type, self_type> mapping_type;

   quan_access_linearized_impl(DomainType& domain, const quan_key_t& quan_key, mapping_type & mapping)
      : domain(domain), quan_key(quan_key), mapping(mapping) {}

  
   template <typename TupleT>
   typename result<TupleT>::type
   eval(TupleT const& args) const
   {
      numeric_t res = domain(args[phoenix::tuple_index<TupleT::length - 1>()], quan_key)(0,0);
      key_type key(quan_key, args[phoenix::tuple_index<TupleT::length - 1>()]);
      return linear_equation_type(mapping(key), res);
   }

   template <typename ElementT>
   void register_entry(ElementT const& elem)
   {
      key_type key(quan_key, elem);
      mapping.add(key, this);
   }

   template <typename ElementT>
   long get_index(ElementT const& elem)
   {
      return mapping(key_type(quan_key, elem));
   }

   void set_value(key_type key, numeric_t  value) 
   {
      domain(key.v, key.quan)(0,0) = value;
   }    


   void update_value(key_type key, numeric_t  value) 
   {
      domain(key.v, key.quan)(0,0) -= value;
   }    


  template <typename InsertFunctor, typename RHSFunctor, typename ElementT>
  void insert(linear_equation_type & eqn, ElementT & elem, InsertFunctor & insert, RHSFunctor & rhs)
  {
    long index1 = get_index(elem);

    typename linear_equation_type::iterator iter;
    for(iter = eqn.begin(); iter != eqn.end(); ++iter)
      {
	long index2 = (*iter).first;
	double entry = (*iter).second;
	 //std::cout << index1 << "  " << index2 << "  " << entry << std::endl;
	insert(index1,index2,entry);
      }
    double rhs_entry = eqn();
    rhs(index1, rhs_entry);
  }


private:
   DomainType&    domain;
   quan_key_t     quan_key;
   mapping_type&  mapping;    
};
*/


// ===========================================================================================
//
// object generator 
//
template <typename KeyTag, typename ObjectT, typename DomainType>
inline boost::phoenix::actor<quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag> > const
scalar_quan_linearized(DomainType& domain, 
                       typename DomainType::quan_key_t const& key,
                       typename quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag>::mapping_type& mapping)  
{
   return boost::phoenix::actor<quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag> >
      (quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag>(domain, key, mapping));
}

template <typename DomainType, typename ObjectT, typename KeyTag >
struct scalar_quan_linearized_t
{
   typedef boost::phoenix::actor<quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag> >           type;

   typedef typename quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag>::mapping_type            mapping_type;
   typedef typename quan_access_linearized_impl<DomainType, ObjectT, scalar_quantity_entry, KeyTag>::linear_equation_type    linear_equation_type;
};





// template <unsigned int upper, typename BaseT>
// void scalar_quan_linearized(boost::phoenix::actor<BaseT> const& v);                         //  This is undefined and not allowed.


}  // end of namespace gsse
#endif
