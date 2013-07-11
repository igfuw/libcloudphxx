#pragma once

#define libcloudphxx_decltype_return_macro(expr) -> decltype(expr) { return expr; }

#define libcloudphxx_declare_const_macro(name, value, unit) template <typename real_t> \
static constexpr auto name() libcloudphxx_decltype_return_macro(real_t(value) * unit)

#define libcloudphxx_derived_const_macro(name, value) template <typename real_t> \
static constexpr auto name() libcloudphxx_decltype_return_macro(value)

#define libcloudphxx_declare_funct_macro template <typename real_t>
