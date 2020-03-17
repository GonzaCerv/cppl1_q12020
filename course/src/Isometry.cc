/**
 * @file Isometry.cc
 * @author Gonzalo Cervetti (g.cervetti@creativa77.com)
 * @brief
 * @version 0.1
 * @date 2020-02-03
 *
 * @copyright Copyright (c) 2020
 *
 */

// Challenge libraries.
#include <Isometry.hpp>
#include <array>

namespace ekumen {
namespace math {

// Vector3
Vector3::Vector3(const double& x, const double& y, const double& z) : x_(x), y_(y), z_(z) {}

Vector3::Vector3() : Vector3(0.0, 0.0, 0.0) {}

Vector3::Vector3(std::initializer_list<double> elements) {
    if (elements.size() < 3) {
        throw std::range_error("Not enough values in initializer list");
    }
    auto it = elements.begin();
    x_ = *it++;
    y_ = *it++;
    z_ = *it++;
}

Vector3 Vector3::cross(const Vector3& v) const {
    auto i = (y_ * v.z_) - (v.y_ * z_);
    auto j = (x_ * v.z_) - (v.x_ * z_);
    auto k = (x_ * v.y_) - (v.x_ * y_);

    return Vector3(i, j, k);
}

Vector3 Vector3::operator+(const Vector3& a) const { return Vector3(x_ + a.x_, y_ + a.y_, z_ + a.z_); }

Vector3 Vector3::operator-(const Vector3& a) const { return Vector3(x_ - a.x_, y_ - a.y_, z_ - a.z_); }

Vector3 Vector3::operator*(const double& a) const { return Vector3(x_ * a, y_ * a, z_ * a); }

Vector3 Vector3::operator*(const Vector3& vect) const { return Vector3(x_ * vect.x(), y_ * vect.y(), z_ * vect.z()); }

Vector3 Vector3::operator/(const Vector3& a) const { return Vector3(x_ / a.x_, y_ / a.y_, z_ / a.z_); }

const double& Vector3::operator[](int val) const {
    switch (val) {
        case 0:
            return x_;
        case 1:
            return y_;
        case 2:
            return z_;
        default:
            throw std::overflow_error("value too big");
    }
}

double& Vector3::operator[](int val) {
    switch (val) {
        case 0:
            return x_;
        case 1:
            return y_;
        case 2:
            return z_;
        default:
            throw std::overflow_error("value too big");
    }
}

bool Vector3::operator==(const std::initializer_list<double>& rhs) const {
    if (rhs.size() != 3) return false;
    auto it = rhs.begin();
    if (!almost_equal(x_, *it, resolution)) {
        return false;
    }
    ++it;
    if (!almost_equal(y_, *it, resolution)) {
        return false;
    }
    ++it;
    if (!almost_equal(z_, *it, resolution)) {
        return false;
    }
    return true;
}

bool Vector3::operator==(const Vector3& a) const {
    return (Vector3(x_, y_, z_) == std::initializer_list<double>({a.x_, a.y_, a.z_}));
}

bool Vector3::operator!=(const std::initializer_list<double>& rhs) const { return !(Vector3(x_, y_, z_) == rhs); }

bool Vector3::operator!=(const Vector3& a) const { return !(Vector3(x_, y_, z_) == a); }

const Vector3 Vector3::kUnitX(1.0, 0.0, 0.0);
const Vector3 Vector3::kUnitY(0.0, 1.0, 0.0);
const Vector3 Vector3::kUnitZ(0.0, 0.0, 1.0);
const Vector3 Vector3::kZero(0.0, 0.0, 0.0);

// Matrix3

Matrix3::Matrix3(const std::initializer_list<double>& row_1, const std::initializer_list<double>& row_2,
                 const std::initializer_list<double>& row_3) {
    mat_.resize(3);
    if ((row_1.size() < 3) || (row_2.size() < 3) || (row_3.size() < 3)) {
        throw "Size of the matrix isn't correct";
    }
    auto it = row_1.begin();
    mat_[0][0] = *it++;
    mat_[0][1] = *it++;
    mat_[0][2] = *it++;

    it = row_2.begin();
    mat_[1][0] = *it++;
    mat_[1][1] = *it++;
    mat_[1][2] = *it++;

    it = row_3.begin();
    mat_[2][0] = *it++;
    mat_[2][1] = *it++;
    mat_[2][2] = *it++;
}
Matrix3::Matrix3(const double m00, const double m01, const double m02, const double m10, const double m11,
                 const double m12, const double m20, const double m21, const double m22)
    : Matrix3{{m00, m01, m02}, {m10, m11, m12}, {m20, m21, m22}} {}

Matrix3::Matrix3() : Matrix3{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}} {}

Matrix3::Matrix3(const Matrix3& mat)
    : Matrix3{mat[0][0], mat[0][1], mat[0][2], mat[1][0], mat[1][1], mat[1][2], mat[2][0], mat[2][1], mat[2][2]} {}

Matrix3::Matrix3(Matrix3&& mat) : mat_(std::move(mat.mat_)) {}

Matrix3& Matrix3::operator=(Matrix3&& mat) {
    mat_ = std::move(mat.mat_);
    return *this;
}

Matrix3& Matrix3::operator=(const Matrix3& mat) {
    // check for self-assignment
    if (&mat == this) return *this;
    // reuse storage when possible
    mat_ = mat.mat_;
    return *this;
}

bool Matrix3::operator==(const Matrix3& a) const {
    for (auto idx = 0; idx < 3; ++idx) {
        if (mat_[idx] != a[idx]) {
            return false;
        }
    }
    return true;
}

bool Matrix3::operator!=(const Matrix3& a) const {
    for (auto idx = 0; idx < 3; ++idx) {
        if (mat_[idx] != a[idx]) {
            return true;
        }
    }
    return false;
}

bool Matrix3::operator==(const std::initializer_list<double>& rhs) const {
    if (rhs.size() < 9) {
        return false;
    }
    auto it = rhs.begin();

    for (auto idx = 0; idx < 3; ++idx) {
        auto val1 = *it++;
        auto val2 = *it++;
        auto val3 = *it++;
        if (mat_[idx] != std::initializer_list<double>({val1, val2, val3})) {
            return false;
        }
    }
    return true;
}

Matrix3 Matrix3::operator+(const Matrix3& a) const {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
        result[idx] = mat_[idx] + a.mat_[idx];
    }
    return result;
}

Matrix3 Matrix3::operator-(const Matrix3& a) const {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
        result[idx] = mat_[idx] - a.mat_[idx];
    }
    return result;
}

Matrix3 Matrix3::operator*(const Matrix3& a) const {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
        result[idx] = mat_[idx] * a.mat_[idx];
    }
    return result;
}

Matrix3 Matrix3::operator*(const double& a) const {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
        result[idx] = mat_[idx] * a;
    }
    return result;
}

Vector3 Matrix3::operator*(const Vector3& a) const {
    Vector3 result;
    result[0] = mat_[0][0] * a[0] + mat_[0][1] * a[1] + mat_[0][2] * a[2];
    result[1] = mat_[1][0] * a[0] + mat_[1][1] * a[1] + mat_[1][2] * a[2];
    result[2] = mat_[2][0] * a[0] + mat_[2][1] * a[1] + mat_[2][2] * a[2];
    return result;
}

Matrix3 Matrix3::operator/(const Matrix3& a) const {
    Matrix3 result;
    for (auto idx = 0; idx < 3; ++idx) {
        result[idx] = mat_[idx] / a.mat_[idx];
    }
    return result;
}

Matrix3 Matrix3::operator/(const double a) const {
    Matrix3 result;
    result.mat_ = mat_;
    for (auto row = 0; row < 2; ++row) {
        for (auto column = 0; column < 2; ++column) {
            result.mat_[row][column] /= a;
        }
    }
    return result;
}

const Vector3& Matrix3::operator[](int val) const {
    if (val > 2) {
        return mat_[2];
    }
    return mat_[val];
}

Vector3& Matrix3::operator[](int val) {
    if (val > 2) {
        return mat_[2];
    }
    return mat_[val];
}

double Matrix3::det() const {
    return ((mat_[0][0] * mat_[1][1] * mat_[2][2]) + (mat_[0][1] * mat_[1][2] * mat_[2][0]) +
            (mat_[0][2] * mat_[1][0] * mat_[2][1]) - (mat_[0][2] * mat_[1][1] * mat_[2][0]) -
            (mat_[0][0] * mat_[1][2] * mat_[2][1]) - (mat_[0][1] * mat_[1][0] * mat_[2][2]));
}

Vector3 Matrix3::row(uint8_t row) const {
    if (row < 2) {
        return mat_[row];
    }
    return mat_[2];
}

Vector3 Matrix3::col(uint8_t col) const {
    if (col < 2) {
        return Vector3(mat_[0][col], mat_[1][col], mat_[2][col]);
    }
    return Vector3(mat_[0][2], mat_[1][2], mat_[2][2]);
}

const Matrix3 Matrix3::kIdentity({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
const Matrix3 Matrix3::kZero({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
const Matrix3 Matrix3::kOnes({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0});

Isometry::Isometry(const Matrix3& rotation) {
    rotation_ = rotation;
    translation_ = Vector3(0.0, 0.0, 0.0);
}

Isometry::Isometry(const Vector3& translation, const Matrix3& rotation) {
    rotation_ = rotation;
    translation_ = translation;
}

Isometry Isometry::compose(const Isometry& isometry) const {
    Matrix3 rotation_result = matrix_multiplication(rotation_, isometry.rotation_);
    Vector3 translation_result = (rotation_ * isometry.translation_) + translation_;
    return Isometry(translation_result, rotation_result);
}

Isometry Isometry::inverse() const {
    auto det = rotation_.det();
    if (almost_equal(det, 0., 4)) {
        throw std::domain_error("Cannot get the inverse of the matrix");
    }
    Matrix3 rot_result;
    rot_result[0][0] = (rotation_[1][1] * rotation_[2][2]) - (rotation_[2][1] * rotation_[1][2]);
    rot_result[0][1] = (rotation_[0][2] * rotation_[2][1]) - (rotation_[2][2] * rotation_[0][1]);
    rot_result[0][2] = (rotation_[0][1] * rotation_[1][2]) - (rotation_[1][1] * rotation_[0][2]);

    rot_result[1][0] = (rotation_[1][2] * rotation_[2][0]) - (rotation_[2][2] * rotation_[1][0]);
    rot_result[1][1] = (rotation_[0][0] * rotation_[2][2]) - (rotation_[2][0] * rotation_[0][2]);
    rot_result[1][2] = (rotation_[0][2] * rotation_[1][0]) - (rotation_[1][2] * rotation_[0][0]);

    rot_result[2][0] = (rotation_[1][0] * rotation_[2][1]) - (rotation_[2][0] * rotation_[1][1]);
    rot_result[2][1] = (rotation_[0][1] * rotation_[2][0]) - (rotation_[2][1] * rotation_[0][0]);
    rot_result[2][2] = (rotation_[0][0] * rotation_[1][1]) - (rotation_[1][0] * rotation_[0][1]);

    rot_result = rot_result / det;
    auto vector_result = rot_result * translation_;
    vector_result = vector_result * (-1);

    return Isometry(vector_result, rot_result);
}

Matrix3 Isometry::rotation() const { return rotation_; }

Vector3 Isometry::transform(Vector3 translation) const { return (rotation_ * translation + translation_); }

Vector3 Isometry::translation() const { return translation_; }

Isometry Isometry::FromTranslation(const std::initializer_list<double>& values) {
    return Isometry(Vector3(values), Matrix3::kIdentity);
}

Isometry Isometry::FromEulerAngles(double yaw, double pitch, double roll) {
    Isometry result = Isometry(RotateAround(Vector3::kUnitX, yaw)) * Isometry(RotateAround(Vector3::kUnitY, pitch)) *
                      Isometry(RotateAround(Vector3::kUnitZ, roll));
    return result;
}

Matrix3 Isometry::RotateAround(const Vector3& axis, double angle) {
    Matrix3 result;
    auto cos = std::cos(angle);
    auto sin = std::sin(angle);
    result[0][0] = cos + (axis.x() * axis.x()) * (1 - cos);
    result[0][1] = axis.x() * axis.y() * (1 - cos) - axis.z() * sin;
    result[0][2] = axis.x() * axis.z() * (1 - cos) + axis.y() * sin;
    result[1][0] = axis.y() * axis.x() * (1 - cos) + axis.z() * sin;
    result[1][1] = cos + (axis.y() * axis.y()) * (1 - cos);
    result[1][2] = axis.y() * axis.z() * (1 - cos) - axis.x() * sin;
    result[2][0] = axis.z() * axis.x() * (1 - cos) + axis.y() * sin;
    result[2][1] = axis.z() * axis.y() * (1 - cos) + axis.x() * sin;
    result[2][2] = cos + (axis.z() * axis.z()) * (1 - cos);
    return result;
}

bool Isometry::operator==(const Isometry& rhs) const {
    if (rotation_ != rhs.rotation_) {
        return false;
    }
    if (translation_ != rhs.translation_) {
        return false;
    }
    return true;
}

Vector3 Isometry::operator*(const Vector3& rhs) const { return (rotation_ * rhs + translation_); }

Isometry Isometry::operator*(const Isometry& rhs) const {
    Matrix3 rotation_result = matrix_multiplication(rotation_, rhs.rotation_);
    Vector3 translation_result = (rotation_ * rhs.translation_) + translation_;
    return Isometry(translation_result, rotation_result);
}

}  // namespace math
}  // namespace ekumen
