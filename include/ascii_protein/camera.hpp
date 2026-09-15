#pragma once

#include "ascii_protein/geometry.hpp"

#include <cmath>

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

    void rotate(double yawDelta, double pitchDelta, double rollDelta = 0.0) {
        yaw = wrapAngle(yaw + yawDelta);
        pitch = wrapAngle(pitch + pitchDelta);
        roll = wrapAngle(roll + rollDelta);
    }

    void reset() { *this = Camera{}; }

private:
    static double wrapAngle(double angle) {
        angle = std::fmod(angle + 180.0, 360.0);
        if (angle < 0.0) angle += 360.0;
        return angle - 180.0;
    }
};

} // namespace ascii_protein