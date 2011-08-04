/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_INSERT_10022005_1837)
#define FUSION_INSERT_10022005_1837

#include <boost/mpl/insert.hpp>
#include <boost/fusion/support/tag_of.hpp>
#include <boost/fusion/algorithm/transformation/insert.hpp>
#include <boost/fusion/sequence/intrinsic/mpl/detail/as.hpp>

namespace boost { namespace mpl
{
    template <typename Tag>
    struct insert_impl;

    template <>
    struct insert_impl<fusion::fusion_sequence_tag>
    {
        template <typename Sequence, typename Pos, typename T>
        struct apply
        {
            typedef typename 
                fusion::result_of::insert<Sequence, Pos, T>::type 
            result;
            
            typedef typename 
                fusion::detail::as<typename fusion::traits::tag_of<Sequence>::type, result>::type 
            type;
        };
    };
}}

#endif
