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

namespace ekumen {
namespace math {
Vector3::Vector3(const double& x, const double& y, const double& z) : x_(x), y_(y), z_(z) {}

Vector3::Vector3() : Vector3(0.0, 0.0, 0.0) {}

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

}  // namespace math
}  // namespace ekumen
