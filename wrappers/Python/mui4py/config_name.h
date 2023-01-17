
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
}
