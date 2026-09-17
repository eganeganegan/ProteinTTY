#pragma once

#include "ascii_protein/geometry.hpp"

#include <cstddef>
#include <string_view>
#include <vector>

namespace ascii_protein {

class CoordinateFrameSource {
public:
    virtual ~CoordinateFrameSource() = default;
    virtual std::size_t atomCount() const = 0;
    virtual std::size_t frameCount() const = 0;
    virtual std::string_view atomLabel(std::size_t index) const = 0;
    virtual std::string_view frameLabel(std::size_t index) const = 0;
    virtual const std::vector<Vec3>& frame(std::size_t index) const = 0;
};

} // namespace ascii_protein