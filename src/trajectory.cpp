#include "ascii_protein/trajectory.hpp"
#include "ascii_protein/trajectory_analysis.hpp"

#include "ascii_protein/dcd_trajectory_source.hpp"
#include "ascii_protein/structure.hpp"
#include "ascii_protein/xyz_trajectory_source.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace ascii_protein {
namespace {

std::string normalizedElement(std::string value) {
    for (auto& character : value) character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    return value;
}

void applyFrame(Structure& structure, const TrajectoryState& state, std::size_t frame) {
    const auto positions = trajectoryCoordinates(state, frame);
    if (state.alignmentAtoms.empty()) {
        structure.updatePositions(positions);
        return;
    }
    Vec3 center;
    for (auto p : state.alignmentReference)
        center += p;
    center = center / static_cast<double>(state.alignmentReference.size());
    structure.updatePositions(positions, center);
    double radius = 1;
    for (auto p : state.alignmentReference)
        radius = std::max(radius, length(p - center));
    structure.boundingRadius = radius;
}
} // namespace

void loadTrajectory(const std::filesystem::path& path, Structure& structure, TrajectoryState& state) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    std::unique_ptr<CoordinateFrameSource> source;
    if (extension == ".dcd") source = std::make_unique<DcdTrajectorySource>(path);
    else source = std::make_unique<XyzTrajectorySource>(path);
    if (source->atomCount() != structure.atoms.size()) {
        throw std::runtime_error("trajectory has " + std::to_string(source->atomCount()) +
                                 " atoms; structure has " + std::to_string(structure.atoms.size()));
    }
    for (std::size_t index = 0; index < structure.atoms.size(); ++index) {
        const auto sourceLabel = source->atomLabel(index);
        if (!sourceLabel.empty() && !structure.atoms[index].element.empty() &&
            sourceLabel != normalizedElement(structure.atoms[index].element)) {
            throw std::runtime_error("trajectory atom " + std::to_string(index + 1) + " is " +
                                     std::string(sourceLabel) + "; topology atom is " +
                                     structure.atoms[index].element);
        }
    }
    const auto firstFrame = source->frame(0);
    structure.updatePositions(firstFrame);
    state.source = std::move(source);
    state.path = path;
    state.frameIndex = 0;
    state.rangeStart = 0;
    state.rangeEnd = state.source->frameCount() - 1;
    state.stride = 1;
    state.alignmentAtoms.clear();
    state.alignmentReference.clear();
    state.playing = state.source->frameCount() > 1;
}

void seekTrajectory(Structure& structure, TrajectoryState& state, std::size_t frameIndex) {
    if (!state.source || frameIndex < state.rangeStart || frameIndex > state.rangeEnd)
        throw std::runtime_error("frame is outside the playback range");
    applyFrame(structure, state, frameIndex);
    state.frameIndex = frameIndex;
}

void setTrajectoryRange(Structure& structure, TrajectoryState& state, std::size_t first, std::size_t last) {
    if (!state.source || first > last || last >= state.source->frameCount())
        throw std::runtime_error("invalid trajectory range");
    // Read before changing state so an I/O failure preserves the current range.
    const auto next = std::clamp(state.frameIndex, first, last);
    applyFrame(structure, state, next);
    state.rangeStart = first;
    state.rangeEnd = last;
    state.frameIndex = next;
    state.playing = false;
}

bool stepTrajectory(Structure& structure, TrajectoryState& state, int offset) {
    if (!state.source || offset == 0) return false;
    const auto count = static_cast<long long>(state.rangeEnd - state.rangeStart + 1);
    long long next = static_cast<long long>(state.frameIndex - state.rangeStart) +
                     static_cast<long long>(offset) * state.stride;
    if (state.loop) next = (next % count + count) % count;
    else if (next < 0 || next >= count) return false;
    seekTrajectory(structure, state, state.rangeStart + static_cast<std::size_t>(next));
    return true;
}

void stopTrajectory(Structure& structure, TrajectoryState& state) {
    state.playing = false;
    if (!state.source) return;
    seekTrajectory(structure, state, state.rangeStart);
}

std::string trajectoryStatus(const TrajectoryState& state) {
    if (!state.source) return {};
    std::ostringstream output;
    output << "Frame " << state.frameIndex + 1 << '/' << state.source->frameCount() << "  "
           << (state.playing ? "playing" : "paused") << "  " << state.frameDelayMs << "ms";
    output << "  range " << state.rangeStart + 1 << '-' << state.rangeEnd + 1 << "  stride " << state.stride;
    if (state.loop) output << "  loop";
    return output.str();
}

} // namespace ascii_protein
