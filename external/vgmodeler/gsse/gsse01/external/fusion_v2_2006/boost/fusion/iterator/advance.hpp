/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_ADVANCE_09172005_1146)
#define FUSION_ADVANCE_09172005_1146

#include <boost/fusion/iterator/detail/advance.hpp>
#include <boost/fusion/support/category_of.hpp>

#include <boost/mpl/int.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/fusion/support/tag_of.hpp>

namespace boost { namespace fusion
{
    struct random_access_traversal_tag;

    namespace extension
    {
        template <typename Tag>
        struct advance_impl
        {
            // default implementation
            template <typename Iterator, typename N>
            struct apply :
                mpl::if_c<
                    (N::value > 0)
                  , advance_detail::forward<Iterator, N::value>
                  , advance_detail::backward<Iterator, N::value>
                >::type
            {
                typedef typename traits::category_of<Iterator>::type category;
                BOOST_MPL_ASSERT_NOT((is_same<category, random_access_traversal_tag>));
            };
        };
    }
    
    namespace result_of
    {
        template <typename Iterator, int N>
        struct advance_c
            : extension::advance_impl<typename traits::tag_of<Iterator>::type>::template apply<Iterator, mpl::int_<N> >
        {};

        template <typename Iterator, typename N>
        struct advance
            : extension::advance_impl<typename traits::tag_of<Iterator>::type>::template apply<Iterator, N>
        {};
    }

    template <int N, typename Iterator>
    inline typename result_of::advance_c<Iterator, N>::type
    advance_c(Iterator const& i)
    {
        return extension::advance_impl<typename Iterator::ftag>::
            template apply<Iterator, mpl::int_<N> >::call(i);
    }

    template<typename N, typename Iterator>
    inline typename result_of::advance<Iterator, N>::type
    advance(Iterator const& i)
    {
        return extension::advance_impl<typename Iterator::ftag>::
            template apply<Iterator, N>::call(i);
    }

}} // namespace boost::fusion

#endif
