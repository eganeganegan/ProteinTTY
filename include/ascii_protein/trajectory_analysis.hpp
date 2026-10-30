#pragma once
#include "ascii_protein/geometry.hpp"
#include <vector>
#include <cstddef>
namespace ascii_protein {
struct TrajectoryState;
struct RigidTransform {
    Mat3 rotation{Mat3::identity()};
    Vec3 translation;
    Vec3 apply(Vec3 p) const { return rotation * p + translation; }
};
RigidTransform fitCoordinates(const std::vector<Vec3>& moving, const std::vector<Vec3>& reference,
                              const std::vector<std::size_t>& selection);
std::vector<Vec3> trajectoryCoordinates(const TrajectoryState& state, std::size_t frame);
std::vector<double> trajectoryRmsd(const TrajectoryState& state, const std::vector<std::size_t>& selection,
                                   std::size_t reference);
std::vector<double> trajectoryRmsf(const TrajectoryState& state, const std::vector<std::size_t>& selection);
} // namespace ascii_protein
