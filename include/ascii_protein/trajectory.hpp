#pragma once

#include "ascii_protein/coordinate_frame_source.hpp"

#include <filesystem>
#include <memory>
#include <string>

namespace ascii_protein {

struct Structure;

struct TrajectoryState {
    std::unique_ptr<CoordinateFrameSource> source;
    std::filesystem::path path;
    std::size_t frameIndex{};
    std::size_t rangeStart{};
    std::size_t rangeEnd{};
    int stride{1};
    std::vector<std::size_t> alignmentAtoms;
    std::vector<Vec3> alignmentReference;
    int frameDelayMs{100};
    bool playing{};
    bool loop{};
};

void loadTrajectory(const std::filesystem::path& path, Structure& structure, TrajectoryState& state);
void seekTrajectory(Structure& structure, TrajectoryState& state, std::size_t frameIndex);
void setTrajectoryRange(Structure& structure, TrajectoryState& state, std::size_t first, std::size_t last);
bool stepTrajectory(Structure& structure, TrajectoryState& state, int offset = 1);
void stopTrajectory(Structure& structure, TrajectoryState& state);
std::string trajectoryStatus(const TrajectoryState& state);

} // namespace ascii_protein