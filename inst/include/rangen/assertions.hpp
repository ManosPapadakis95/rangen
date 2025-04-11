

#ifndef ASSERTIONS_H
#define ASSERTIONS_H

namespace Assertion
{

#include <utility>
#include <string_view>

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

    struct has_col : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must define a function named 'col' for getting the i-th column.");
        }
    };

    template <typename T>
    struct has_row<T, decltype(std::declval<T>().row(0), void())> : std::true_type
    {
        constexpr static void check_concept() {}
    };

    struct has_row : std::false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!std::is_same<T, void>::value, "The template class must define a function named 'row' for getting the i-th row.");
        }
    };

    template <typename T>
    struct has_col<T, decltype(std::declval<T>().col(0), void())> : std::true_type
    {
        constexpr static void check_concept() {}
    };

    template <typename, typename = std::void_t<>>
    struct has_ncol : std::false_type
    {
    };

    template <typename T>
    struct has_ncol<T, std::void_t<decltype(std::declval<T>().ncol())>> : std::true_type
    {
    };

    template <typename, typename = std::void_t<>>
    struct has_nrow : std::false_type
    {
    };

    template <typename T>
    struct has_nrow<T, std::void_t<decltype(std::declval<T>().nrow())>> : std::true_type
    {
    };

    template <typename T>
    constexpr void has_value_int()
    {
        constexpr auto ok = std::is_same<typename T::value_type, int>::value;
        static_assert(!ok, "The template class value must be of type integer.");
    };

    template <typename T>
    constexpr void has_value_bool()
    {
        constexpr auto ok = std::is_same<typename T::value_type, bool>::value;
        static_assert(!ok, "The template class value must be of type bool.");
    };

}

#endif