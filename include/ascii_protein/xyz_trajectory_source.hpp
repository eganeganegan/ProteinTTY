#pragma once

#include "ascii_protein/coordinate_frame_source.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace ascii_protein {

class XyzTrajectorySource final : public CoordinateFrameSource {
public:
    explicit XyzTrajectorySource(const std::filesystem::path& path);

    std::size_t atomCount() const override { return atomCount_; }
    std::size_t frameCount() const override { return frames_.size(); }
    std::string_view atomLabel(std::size_t index) const override;
    std::string_view frameLabel(std::size_t index) const override;
    const std::vector<Vec3>& frame(std::size_t index) const override;

private:
    std::size_t atomCount_{};
    std::vector<std::string> atomLabels_;
    std::vector<std::string> labels_;
    std::vector<std::vector<Vec3>> frames_;
};

} // namespace ascii_protein