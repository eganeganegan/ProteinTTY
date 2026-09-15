#pragma once

#include "ascii_protein/geometry.hpp"

namespace ascii_protein {

enum class Projection { perspective, orthographic };

struct Camera {
    double yaw{-25.0};
    double pitch{15.0};
    double roll{0.0};
    double zoom{1.0};
    Vec2 pan{};
    Projection projection{Projection::perspective};

    Mat3 rotation() const {
        return Mat3::rotationZ(radians(roll)) * Mat3::rotationX(radians(pitch)) *
               Mat3::rotationY(radians(yaw));
    }

    void reset() { *this = Camera{}; }
};

} // namespace ascii_protein