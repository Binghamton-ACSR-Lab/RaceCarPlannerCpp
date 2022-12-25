//
// Created by acsr on 5/25/21.
//

#ifndef BOOST_RTREE_ARRAY_ADAPTOR_HPP
#define BOOST_RTREE_ARRAY_ADAPTOR_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <array>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


namespace boost::geometry::traits{

    template <typename CoordinateType, std::size_t DimensionCount>
    struct tag<std::array<CoordinateType, DimensionCount>>{
        using type = point_tag;
    };

    template <typename CoordinateType, std::size_t DimensionCount>
    struct coordinate_type<std::array<CoordinateType, DimensionCount>>{
        typedef CoordinateType type;
    };

    template <typename CoordinateType, std::size_t DimensionCount>
    struct dimension<std::array<CoordinateType, DimensionCount>>: boost::mpl::int_<DimensionCount> {

    };

    template <typename CoordinateType, size_t DimensionCount, size_t Dimension>
    struct access<std::array<CoordinateType, DimensionCount>, Dimension>{
        static inline CoordinateType get(std::array<CoordinateType, DimensionCount> const& a)
        {
            return a[Dimension];
        }

        static inline void set(std::array<CoordinateType, DimensionCount>& a,
                               CoordinateType const& value){
            a[Dimension] = value;
        }
    };
}

#define BOOST_GEOMETRY_REGISTER_ARRAY_CS(CoordinateSystem) \
namespace boost::geometry::traits{                     \
template <class T, size_t N>                              \
struct coordinate_system<std::array<T, N>>               \
{                                                              \
  typedef CoordinateSystem type;                               \
};                                                             \
}

#endif //BOOST_RTREE_EIGEN_ADAPTOR_HPP
