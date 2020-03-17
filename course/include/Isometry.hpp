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

// Gtest
#include "gtest/gtest.h"

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
    Vector3(std::initializer_list<double> elements);
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
    bool operator!=(const Matrix3& a) const;
    bool operator==(const std::initializer_list<double>& rhs) const;
    Matrix3 operator+(const Matrix3& a) const;
    Matrix3 operator-(const Matrix3& a) const;
    Matrix3 operator*(const double& a) const;
    Matrix3 operator*(const Matrix3& a) const;
    Vector3 operator*(const Vector3& a) const;
    Matrix3 operator/(const Matrix3& a) const;
    Matrix3 operator/(const double a) const;
    const Vector3& operator[](int val) const;
    Vector3& operator[](int val);

    // Methods
    double det() const;
    Vector3 row(uint8_t row) const;
    Vector3 col(uint8_t col) const;

    // Static members
    static const Matrix3 RotateAround(const Vector3& axis, double angle);
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

const inline Matrix3 matrix_multiplication(Matrix3 mat, Matrix3 rhs) {
    Matrix3 result;
    for (uint16_t row = 0; row < 3; ++row) {
        for (uint16_t column = 0; column < 3; ++column) {
            result[row][column] =
                mat[row][0] * rhs[0][column] + mat[row][1] * rhs[1][column] + mat[row][2] * rhs[2][column];
        }
    }
    return result;
}

class Isometry {
   public:
    Isometry(const Matrix3& rotation);
    Isometry(const Vector3& translation, const Matrix3& rotation);

    Isometry compose(const Isometry& isometry) const;
    Isometry inverse() const;
    Matrix3 rotation() const;
    Vector3 transform(Vector3 translation) const;
    Vector3 translation() const;

    static Isometry FromTranslation(const std::initializer_list<double>& values);
    static Isometry FromEulerAngles(double yaw, double pitch, double roll);
    static Matrix3 RotateAround(const Vector3& axis, double angle);

    // Operators
    bool operator==(const Isometry& a) const;
    Vector3 operator*(const Vector3& a) const;
    Isometry operator*(const Isometry& a) const;

   private:
    Matrix3 rotation_;
    Vector3 translation_;
};

inline ::testing::AssertionResult areAlmostEqual(const Isometry& lhs, const Isometry& rhs, const double tolerance) {
    auto lhs_rotation = lhs.rotation();
    auto lhs_translation = lhs.translation();

    auto rhs_rotation = rhs.rotation();
    auto rhs_translation = rhs.translation();

    for (uint8_t index_x = 0; index_x < 3; ++index_x) {
        for (uint8_t index_y = 0; index_y < 3; ++index_y) {
            if (!almost_equal(lhs_rotation[index_x][index_y], rhs_rotation[index_x][index_y], tolerance)) {
                return ::testing::AssertionFailure() << "Elements of rotational matrix are not equal";
            }
        }
    }
    for (uint8_t index = 0; index < 3; ++index) {
        if (!almost_equal(lhs_translation[index], rhs_translation[index], tolerance)) {
            return ::testing::AssertionFailure() << "Elements of translation vector are not equal";
            ;
        }
    }
    return ::testing::AssertionSuccess();
}

inline std::ostream& operator<<(std::ostream& os, const Isometry& rhs) {
    auto rhs_rotation = rhs.rotation();
    return os << "[T: (x: 0, y: 0, z: 0"
              << "), R:[[" << std::setprecision(9) << rhs_rotation[0][0] << ", " << rhs_rotation[0][1] << ", "
              << rhs_rotation[0][2] << "], "
              << "[" << rhs_rotation[1][0] << ", " << rhs_rotation[1][1] << ", " << rhs_rotation[1][2] << "], "
              << "[" << rhs_rotation[2][0] << ", " << rhs_rotation[2][1] << ", " << rhs_rotation[2][2] << "]]]";
}

}  // namespace math
}  // namespace ekumen
