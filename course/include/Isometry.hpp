/**
 * @file Isometry.hpp
 * @author Gonzalo Cervetti (g.cervetti@creativa77.com)
 * @brief This Class provides implementation for vector operation.
 * @version 0.1
 * @date 2020-02-03
 *
 * @copyright Copyright (c) 2020
 *
 */

#pragma once

// Standard libraries
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace ekumen {
namespace math {

static constexpr int resolution{4};

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp) {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
           // unless the result is subnormal
           || std::fabs(x - y) < std::numeric_limits<T>::min();
}

class Vector3 {
   public:
    Vector3(const double& x, const double& y, const double& z);
    Vector3();
    /**
     * @brief Normalizes a function./
     *
     * @return double normalized value of the vector.
     */
    double norm() const { return std::sqrt(std::pow(x_, 2.0) + std::pow(y_, 2.0) + std::pow(z_, 2.0)); };
    double dot(const Vector3& v) const { return ((x_ * v.x_) + (y_ * v.y_) + (z_ * v.z_)); };
    Vector3 cross(const Vector3& v) const;

    // Getter for elements of the vector.
    const double& x() const { return x_; };
    const double& y() const { return y_; };
    const double& z() const { return z_; };

    double& x() { return x_; };
    double& y() { return y_; };
    double& z() { return z_; };

    // Operators overloading.
    Vector3 operator+(const Vector3& a) const;
    Vector3 operator-(const Vector3& a) const;
    Vector3 operator*(const double& a) const;
    Vector3 operator*(const Vector3& vect) const;
    Vector3 operator/(const Vector3& a) const;
    const double& operator[](int val) const;
    double& operator[](int val);
    bool operator==(const std::initializer_list<double>& rhs) const;
    bool operator==(const Vector3& a) const;
    bool operator!=(const std::initializer_list<double>& rhs) const;
    bool operator!=(const Vector3& a) const;

    // Static members
    static const Vector3 kUnitX;
    static const Vector3 kUnitY;
    static const Vector3 kUnitZ;
    static const Vector3 kZero;

   private:
    double x_, y_, z_;
};

// non-member binary operators

inline Vector3 operator*(const double scalar, const Vector3& vect) {
    Vector3 result(vect.x() * scalar, vect.y() * scalar, (vect.z() * scalar));
    return result;
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& vect) {
    return os << "(x: " << std::to_string(static_cast<int>(vect.x()))
              << ", y: " << std::to_string(static_cast<int>(vect.y()))
              << ", z: " << std::to_string(static_cast<int>(vect.z())) << ")";
}

}  // namespace math
}  // namespace ekumen
