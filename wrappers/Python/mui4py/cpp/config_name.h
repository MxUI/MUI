
#include <mui.h>
#include <string>

template <typename CONFIG>
std::string config_name()
{
    if (std::is_same<CONFIG, mui::mui_config_1d>())
        return "1d_f64_i32";
    if (std::is_same<CONFIG, mui::mui_config_2d>())
        return "2d_f64_i32";
    if (std::is_same<CONFIG, mui::mui_config_3d>())
        return "3d_f64_i32";
    if (std::is_same<CONFIG, mui::mui_config_1dx>())
        return "1d_f64_i64";
    if (std::is_same<CONFIG, mui::mui_config_2dx>())
        return "2d_f64_i64";
    if (std::is_same<CONFIG, mui::mui_config_3dx>())
        return "3d_f64_i64";
    if (std::is_same<CONFIG, mui::mui_config_1f>())
        return "1d_f32_i32";
    if (std::is_same<CONFIG, mui::mui_config_2f>())
        return "2d_f32_i32";
    if (std::is_same<CONFIG, mui::mui_config_3f>())
        return "3d_f32_i32";
    if (std::is_same<CONFIG, mui::mui_config_1fx>())
        return "1d_f32_i64";
    if (std::is_same<CONFIG, mui::mui_config_2fx>())
        return "2d_f32_i64";
    if (std::is_same<CONFIG, mui::mui_config_3fx>())
        return "3d_f32_i64";
    throw std::runtime_error("Invalid config type");
}

template <typename T>
std::string type_name()
{
    if (std::is_same<T, double>::value)
        return "double";
    if (std::is_same<T, float>::value)
        return "float";
    if (std::is_same<T, std::int32_t>::value)
        return "int32_t";
    if (std::is_same<T, std::int64_t>::value)
        return "int64_t";
    if (std::is_same<T, std::string>::value)
        return "string";
    throw std::runtime_error("Invalid type");
}
