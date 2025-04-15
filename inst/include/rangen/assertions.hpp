

#ifndef ASSERTIONS_H
#define ASSERTIONS_H

namespace Assertion
{

#include <utility>

    using namespace std;

    template <typename, typename T = void>
    struct has_subscript_operator : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must support operator[].");
        }
    };

    template <typename T>
    struct has_subscript_operator<T, decltype(std::declval<T>()[std::declval<std::size_t>()], void())> : std::true_type
    {
        constexpr static void check_concept() {}
    };

    template <typename, typename T = void>
    struct has_size : false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!is_same<T, void>::value, "The template class must provide a function named 'size'.");
        }
    };

    template <typename T>
    struct has_size<T, void_t<decltype(declval<T>().size())>> : true_type
    {
        constexpr static void check_concept() {}
    };

    template <typename, typename T = void>
    struct has_value_type : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must have a 'value_type' member.");
        }
    };

    template <typename T>
    struct has_value_type<T, std::void_t<typename T::value_type>> : std::true_type
    {
        constexpr static void check_concept() {}
    };

    template <typename, typename T = std::void_t<>>
    struct has_ncol : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must define a function named 'ncol' for the number of the columns");
        }
    };

    template <typename T>
    struct has_ncol<T, std::void_t<decltype(std::declval<T>().ncol())>> : std::true_type
    {
    };

    template <typename, typename T = std::void_t<>>
    struct has_nrow : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must define a function named 'nrow' for the number of the rows");
        }
    };

    template <typename T>
    struct has_nrow<T, std::void_t<decltype(std::declval<T>().nrow())>> : std::true_type
    {
    };

    template <typename T>
    constexpr void is_value_int()
    {
        constexpr auto ok = std::is_same<typename T::value_type, int>::value;
        static_assert(!ok, "The template class value must be of type integer.");
    };
}

#endif