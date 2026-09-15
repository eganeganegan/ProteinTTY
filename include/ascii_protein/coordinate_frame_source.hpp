#pragma once

#include "ascii_protein/geometry.hpp"

#include <vector>

namespace ascii_protein {

class CoordinateFrameSource {
public:
    virtual ~CoordinateFrameSource() = default;
    virtual bool nextFrame(std::vector<Vec3>& positions) = 0;
};

} // namespace ascii_protein