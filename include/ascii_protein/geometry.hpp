#pragma once

#include <algorithm>
#include <cmath>

namespace ascii_protein {

constexpr double pi = 3.14159265358979323846;

struct Vec2 {
    double x{};
    double y{};
};

struct Vec3 {
    double x{};
    double y{};
    double z{};

    constexpr Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
    constexpr Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
    constexpr Vec3 operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    constexpr Vec3 operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }
    Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
};

inline double dot(const Vec3& left, const Vec3& right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

inline Vec3 cross(const Vec3& left, const Vec3& right) {
    return {left.y * right.z - left.z * right.y,
            left.z * right.x - left.x * right.z,
            left.x * right.y - left.y * right.x};
}

inline double length(const Vec3& value) { return std::sqrt(dot(value, value)); }

inline Vec3 normalized(const Vec3& value) {
    const double magnitude = length(value);
    return magnitude > 0.0 ? value / magnitude : Vec3{};
}

struct Mat3 {
    double m[3][3]{};

    static constexpr Mat3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }
    static Mat3 rotationX(double radians) {
        const double cosine = std::cos(radians);
        const double sine = std::sin(radians);
        return {{{1, 0, 0}, {0, cosine, -sine}, {0, sine, cosine}}};
    }
    static Mat3 rotationY(double radians) {
        const double cosine = std::cos(radians);
        const double sine = std::sin(radians);
        return {{{cosine, 0, sine}, {0, 1, 0}, {-sine, 0, cosine}}};
    }
    static Mat3 rotationZ(double radians) {
        const double cosine = std::cos(radians);
        const double sine = std::sin(radians);
        return {{{cosine, -sine, 0}, {sine, cosine, 0}, {0, 0, 1}}};
    }

    Vec3 operator*(const Vec3& value) const {
        return {m[0][0] * value.x + m[0][1] * value.y + m[0][2] * value.z,
                m[1][0] * value.x + m[1][1] * value.y + m[1][2] * value.z,
                m[2][0] * value.x + m[2][1] * value.y + m[2][2] * value.z};
    }

    Mat3 operator*(const Mat3& other) const {
        Mat3 result{};
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                for (int index = 0; index < 3; ++index) {
                    result.m[row][column] += m[row][index] * other.m[index][column];
                }
            }
        }
        return result;
    }
};

inline double radians(double degrees) { return degrees * pi / 180.0; }

} // namespace ascii_protein