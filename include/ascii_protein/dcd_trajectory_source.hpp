#pragma once

#include "ascii_protein/coordinate_frame_source.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

namespace ascii_protein {

class DcdTrajectorySource final : public CoordinateFrameSource {
public:
    explicit DcdTrajectorySource(const std::filesystem::path& path);

    std::size_t atomCount() const override { return atomCount_; }
    std::size_t frameCount() const override { return frames_.size(); }
    std::string_view atomLabel(std::size_t index) const override;
    std::string_view frameLabel(std::size_t index) const override;
    const std::vector<Vec3>& frame(std::size_t index) const override;

private:
    struct FrameIndex {
        std::array<std::uint64_t, 3> coordinateOffsets{};
        std::size_t coordinateCount{};
    };

    std::filesystem::path path_;
    std::size_t atomCount_{};
    bool swapBytes_{};
    std::vector<std::size_t> freeAtomIndices_;
    std::vector<FrameIndex> frames_;
    std::vector<Vec3> firstFrame_;
    mutable std::vector<Vec3> frameBuffer_;
    mutable std::optional<std::size_t> bufferedFrame_;
};

} // namespace ascii_protein
