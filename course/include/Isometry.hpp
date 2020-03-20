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
#include <vector>

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

class Matrix3 {
   public:
    // Constructors
    Matrix3(const double m00, const double m01, const double m02, const double m10, const double m11, const double m12,
            const double m20, const double m21, const double m22);
    Matrix3(const std::initializer_list<double>& row_1, const std::initializer_list<double>& row_2,
            const std::initializer_list<double>& row_3);
    Matrix3();
    Matrix3(const Matrix3& mat);
    Matrix3(Matrix3&& mat);
    // Operators
    Matrix3& operator=(Matrix3&& mat);
    Matrix3& operator=(const Matrix3& mat);
    bool operator==(const Matrix3& a) const;
    bool operator==(const std::initializer_list<double>& rhs) const;
    Matrix3 operator+(const Matrix3& a) const;
    Matrix3 operator-(const Matrix3& a) const;
    Matrix3 operator*(const double& a) const;
    Matrix3 operator*(const Matrix3& a) const;
    Matrix3 operator/(const Matrix3& a) const;
    const Vector3& operator[](int val) const;
    Vector3& operator[](int val);

    // Methods
    double det() const;
    Vector3 row(uint8_t row) const;
    Vector3 col(uint8_t col) const;

    // Static members
    static const Matrix3 kIdentity;
    static const Matrix3 kZero;
    static const Matrix3 kOnes;

   private:
    std::vector<Vector3> mat_;
};

// non-member binary operators

inline Matrix3 operator*(const double scalar, const Matrix3& rhs) {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
            result[idx] = rhs[idx] * scalar;
    }
    return result;
}

inline std::ostream& operator<<(std::ostream& os, const Matrix3& mat) {
    return os << "[[" << std::to_string(static_cast<int>(mat[0][0])) << ", "
              << std::to_string(static_cast<int>(mat[0][1])) << ", " << std::to_string(static_cast<int>(mat[0][2]))
              << "], "
              << "[" << std::to_string(static_cast<int>(mat[1][0])) << ", "
              << std::to_string(static_cast<int>(mat[1][1])) << ", " << std::to_string(static_cast<int>(mat[1][2]))
              << "], "
              << "[" << std::to_string(static_cast<int>(mat[2][0])) << ", "
              << std::to_string(static_cast<int>(mat[2][1])) << ", " << std::to_string(static_cast<int>(mat[2][2]))
              << "]]";
}

}  // namespace math
}  // namespace ekumen
