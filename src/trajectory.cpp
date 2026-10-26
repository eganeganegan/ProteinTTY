#include "ascii_protein/trajectory.hpp"

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
    state.source = std::move(source);
    state.path = path;
    state.frameIndex = 0;
    state.playing = state.source->frameCount() > 1;
    structure.updatePositions(state.source->frame(0));
}

bool stepTrajectory(Structure& structure, TrajectoryState& state, int offset) {
    if (!state.source || state.source->frameCount() == 0 || offset == 0) return false;
    const auto count = static_cast<long long>(state.source->frameCount());
    long long next = static_cast<long long>(state.frameIndex) + offset;
    if (state.loop) next = (next % count + count) % count;
    else if (next < 0 || next >= count) return false;
    state.frameIndex = static_cast<std::size_t>(next);
    structure.updatePositions(state.source->frame(state.frameIndex));
    return true;
}

void stopTrajectory(Structure& structure, TrajectoryState& state) {
    state.playing = false;
    if (!state.source) return;
    state.frameIndex = 0;
    structure.updatePositions(state.source->frame(0));
}

std::string trajectoryStatus(const TrajectoryState& state) {
    if (!state.source) return {};
    std::ostringstream output;
    output << "Frame " << state.frameIndex + 1 << '/' << state.source->frameCount() << "  "
           << (state.playing ? "playing" : "paused") << "  " << state.frameDelayMs << "ms";
    if (state.loop) output << "  loop";
    return output.str();
}

} // namespace ascii_protein
