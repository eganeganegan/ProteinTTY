#include "ascii_protein/xyz_trajectory_source.hpp"

#include <charconv>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ascii_protein {
namespace {

std::string normalizedLabel(std::string value) {
    for (auto& character : value) character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    return value;
}

std::size_t parseAtomCount(const std::string& line, std::size_t frameNumber) {
    const auto first = line.find_first_not_of(" \t\r");
    const auto last = line.find_last_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '-') {
        throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": invalid atom count");
    }
    std::uint64_t value{};
    const char* begin = line.data() + first;
    const char* end = line.data() + last + 1;
    const auto parsed = std::from_chars(begin, end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end || value == 0 ||
        value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": invalid atom count");
    }
    return static_cast<std::size_t>(value);
}

Vec3 parseCoordinateRow(const std::string& line, std::size_t frameNumber, std::size_t atom,
                        std::string* label = nullptr) {
    std::istringstream row(line);
    std::string element;
    Vec3 position;
    if (!(row >> element >> position.x >> position.y >> position.z) || !std::isfinite(position.x) ||
        !std::isfinite(position.y) || !std::isfinite(position.z)) {
        throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ", atom " +
                                 std::to_string(atom + 1) + ": invalid coordinate row");
    }
    if (label) *label = normalizedLabel(std::move(element));
    return position;
}

std::uint64_t offsetOf(std::istream& input, const std::string& context) {
    const auto offset = input.tellg();
    if (offset < 0) throw std::runtime_error("cannot index XYZ trajectory " + context);
    return static_cast<std::uint64_t>(offset);
}

} // namespace

XyzTrajectorySource::XyzTrajectorySource(const std::filesystem::path& path) : path_(path) {
    std::ifstream input(path_, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open XYZ trajectory: " + path_.string());
    std::string line;
    std::size_t frameNumber = 0;
    while (std::getline(input, line)) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        ++frameNumber;
        const std::size_t count = parseAtomCount(line, frameNumber);
        if (atomCount_ == 0) atomCount_ = count;
        else if (count != atomCount_) {
            throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": atom count changed from " +
                                     std::to_string(atomCount_) + " to " + std::to_string(count));
        }
        labelOffsets_.push_back(offsetOf(input, "comment position"));
        if (!std::getline(input, line)) {
            throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": missing comment line");
        }
        coordinateOffsets_.push_back(offsetOf(input, "coordinate position"));
        for (std::size_t atom = 0; atom < count; ++atom) {
            if (!std::getline(input, line)) {
                throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ": truncated atom list");
            }
            std::string element;
            static_cast<void>(parseCoordinateRow(line, frameNumber, atom, &element));
            if (frameNumber == 1) atomLabels_.push_back(element);
            else if (element != atomLabels_[atom]) {
                throw std::runtime_error("XYZ frame " + std::to_string(frameNumber) + ", atom " +
                                         std::to_string(atom + 1) + ": atom label changed from " +
                                         atomLabels_[atom] + " to " + element);
            }
        }
    }
    if (coordinateOffsets_.empty()) throw std::runtime_error("XYZ trajectory contains no frames: " + path_.string());
}

std::string_view XyzTrajectorySource::atomLabel(std::size_t index) const {
    if (index >= atomLabels_.size()) throw std::out_of_range("XYZ atom index out of range");
    return atomLabels_[index];
}

std::string_view XyzTrajectorySource::frameLabel(std::size_t index) const {
    if (index >= labelOffsets_.size()) throw std::out_of_range("XYZ frame index out of range");
    std::ifstream input(path_, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(labelOffsets_[index]));
    if (!input || !std::getline(input, labelBuffer_)) {
        throw std::runtime_error("cannot read XYZ frame " + std::to_string(index + 1) + " comment");
    }
    if (!labelBuffer_.empty() && labelBuffer_.back() == '\r') labelBuffer_.pop_back();
    return labelBuffer_;
}

const std::vector<Vec3>& XyzTrajectorySource::frame(std::size_t index) const {
    if (index >= coordinateOffsets_.size()) throw std::out_of_range("XYZ frame index out of range");
    if (bufferedFrame_ == index) return frameBuffer_;
    std::ifstream input(path_, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(coordinateOffsets_[index]));
    if (!input) throw std::runtime_error("cannot seek to XYZ frame " + std::to_string(index + 1));
    frameBuffer_.clear();
    frameBuffer_.reserve(atomCount_);
    std::string line;
    for (std::size_t atom = 0; atom < atomCount_; ++atom) {
        if (!std::getline(input, line)) {
            throw std::runtime_error("XYZ frame " + std::to_string(index + 1) + ": truncated atom list");
        }
        frameBuffer_.push_back(parseCoordinateRow(line, index + 1, atom));
    }
    bufferedFrame_ = index;
    return frameBuffer_;
}

} // namespace ascii_protein
