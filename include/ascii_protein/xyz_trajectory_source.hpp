#pragma once

#include "ascii_protein/coordinate_frame_source.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace ascii_protein {

class XyzTrajectorySource final : public CoordinateFrameSource {
public:
    explicit XyzTrajectorySource(const std::filesystem::path& path);

    std::size_t atomCount() const override { return atomCount_; }
    std::size_t frameCount() const override { return coordinateOffsets_.size(); }
    std::string_view atomLabel(std::size_t index) const override;
    std::string_view frameLabel(std::size_t index) const override;
    const std::vector<Vec3>& frame(std::size_t index) const override;

private:
    std::filesystem::path path_;
    std::size_t atomCount_{};
    std::vector<std::string> atomLabels_;
    std::vector<std::uint64_t> labelOffsets_;
    std::vector<std::uint64_t> coordinateOffsets_;
    mutable std::string labelBuffer_;
    mutable std::vector<Vec3> frameBuffer_;
    mutable std::optional<std::size_t> bufferedFrame_;
};

} // namespace ascii_protein
