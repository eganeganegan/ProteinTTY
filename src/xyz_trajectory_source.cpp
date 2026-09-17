#include "ascii_protein/xyz_trajectory_source.hpp"

#include <cmath>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ascii_protein {
namespace {

std::string normalizedLabel(std::string value) {
    for (auto& character : value) character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    return value;
}

} // namespace

XyzTrajectorySource::XyzTrajectorySource(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open XYZ trajectory: " + path.string());
    std::string line;
    std::size_t frameNumber = 0;
    while (std::getline(input, line)) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        ++frameNumber;
        std::size_t count = 0;
        try {
            std::size_t used = 0;
            count = std::stoul(line, &used);
            if (line.find_first_not_of(" \t\r", used) != std::string::npos || count == 0) throw std::invalid_argument("");
        } catch (...) {
            throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": invalid atom count");
        }
        if (atomCount_ == 0) atomCount_ = count;
        else if (count != atomCount_) {
            throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": atom count changed from " +
                                     std::to_string(atomCount_) + " to " + std::to_string(count));
        }
        if (!std::getline(input, line)) {
            throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": missing comment line");
        }
        labels_.push_back(line);
        std::vector<Vec3> positions;
        positions.reserve(count);
        for (std::size_t atom = 0; atom < count; ++atom) {
            if (!std::getline(input, line)) {
                throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": truncated atom list");
            }
            std::istringstream row(line);
            std::string element;
            Vec3 position;
            if (!(row >> element >> position.x >> position.y >> position.z) || !std::isfinite(position.x) ||
                !std::isfinite(position.y) || !std::isfinite(position.z)) {
                throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ", atom " +
                                         std::to_string(atom + 1) + ": invalid coordinate row");
            }
            element = normalizedLabel(element);
            if (frameNumber == 1) atomLabels_.push_back(element);
            else if (element != atomLabels_[atom]) {
                throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ", atom " +
                                         std::to_string(atom + 1) + ": atom label changed from " +
                                         atomLabels_[atom] + " to " + element);
            }
            positions.push_back(position);
        }
        frames_.push_back(std::move(positions));
    }
    if (frames_.empty()) throw std::runtime_error("XYZ trajectory contains no frames: " + path.string());
}

std::string_view XyzTrajectorySource::atomLabel(std::size_t index) const {
    if (index >= atomLabels_.size()) throw std::out_of_range("XYZ atom index out of range");
    return atomLabels_[index];
}

std::string_view XyzTrajectorySource::frameLabel(std::size_t index) const {
    if (index >= labels_.size()) throw std::out_of_range("XYZ frame index out of range");
    return labels_[index];
}

const std::vector<Vec3>& XyzTrajectorySource::frame(std::size_t index) const {
    if (index >= frames_.size()) throw std::out_of_range("XYZ frame index out of range");
    return frames_[index];
}

} // namespace ascii_protein