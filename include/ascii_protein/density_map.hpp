#pragma once
#include "ascii_protein/geometry.hpp"
#include <array>
#include <filesystem>
#include <vector>
namespace ascii_protein {
struct DensityMap {
    std::filesystem::path path;
    std::array<int, 3> dimensions{};
    std::array<Vec3, 3> steps;
    Vec3 origin;
    std::vector<float> values;
    float minimum{}, maximum{};
    Vec3 position(int c, int r, int s) const { return origin + steps[0] * c + steps[1] * r + steps[2] * s; }
    float value(int c, int r, int s) const {
        return values[(static_cast<std::size_t>(s) * dimensions[1] + r) * dimensions[0] + c];
    }
};
DensityMap loadDensityMap(const std::filesystem::path& path);
} // namespace ascii_protein
