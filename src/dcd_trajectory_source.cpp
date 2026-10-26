#include "ascii_protein/dcd_trajectory_source.hpp"

#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ascii_protein {
namespace {

std::uint32_t decodeU32(const char* bytes, bool swap) {
    const auto* value = reinterpret_cast<const unsigned char*>(bytes);
    if (swap) {
        return (static_cast<std::uint32_t>(value[0]) << 24) | (static_cast<std::uint32_t>(value[1]) << 16) |
               (static_cast<std::uint32_t>(value[2]) << 8) | static_cast<std::uint32_t>(value[3]);
    }
    return static_cast<std::uint32_t>(value[0]) | (static_cast<std::uint32_t>(value[1]) << 8) |
           (static_cast<std::uint32_t>(value[2]) << 16) | (static_cast<std::uint32_t>(value[3]) << 24);
}

std::uint32_t readU32(std::istream& input, bool swap, const std::string& context) {
    char bytes[4];
    if (!input.read(bytes, sizeof(bytes))) throw std::runtime_error("truncated DCD " + context);
    return decodeU32(bytes, swap);
}

std::uint64_t streamOffset(std::istream& input, const std::string& context) {
    const auto offset = input.tellg();
    if (offset < 0) throw std::runtime_error("cannot index DCD " + context);
    return static_cast<std::uint64_t>(offset);
}

std::vector<char> readRecord(std::istream& input, bool swap, std::size_t maximum, const std::string& context) {
    const std::uint32_t size = readU32(input, swap, context + " record marker");
    if (size > maximum) throw std::runtime_error("DCD " + context + " record is too large");
    std::vector<char> payload(size);
    if (size != 0 && !input.read(payload.data(), static_cast<std::streamsize>(size))) {
        throw std::runtime_error("truncated DCD " + context + " record");
    }
    if (readU32(input, swap, context + " closing marker") != size) {
        throw std::runtime_error("DCD " + context + " record markers do not match");
    }
    return payload;
}

std::uint64_t indexRecord(std::istream& input, bool swap, std::uint32_t expectedSize, const std::string& context) {
    const std::uint32_t size = readU32(input, swap, context + " record marker");
    if (size != expectedSize) {
        throw std::runtime_error("DCD " + context + " record has " + std::to_string(size) +
                                 " bytes; expected " + std::to_string(expectedSize));
    }
    const auto offset = streamOffset(input, context);
    input.seekg(static_cast<std::streamoff>(size), std::ios::cur);
    if (!input) throw std::runtime_error("truncated DCD " + context + " record");
    if (readU32(input, swap, context + " closing marker") != size) {
        throw std::runtime_error("DCD " + context + " record markers do not match");
    }
    return offset;
}

std::uint32_t peekRecordSize(std::istream& input, bool swap) {
    const auto position = input.tellg();
    const auto result = readU32(input, swap, "frame");
    input.seekg(position);
    if (!input) throw std::runtime_error("cannot seek while indexing DCD frame");
    return result;
}

std::vector<float> readCoordinates(std::istream& input, std::uint64_t offset, std::size_t count, bool swap,
                                   std::size_t frameNumber, char axis) {
    input.seekg(static_cast<std::streamoff>(offset));
    if (!input) throw std::runtime_error("cannot seek to DCD frame " + std::to_string(frameNumber));
    std::vector<float> values;
    values.reserve(count);
    char bytes[4];
    for (std::size_t index = 0; index < count; ++index) {
        if (!input.read(bytes, sizeof(bytes))) {
            throw std::runtime_error("truncated DCD frame " + std::to_string(frameNumber) + " " + axis +
                                     " coordinates");
        }
        const float value = std::bit_cast<float>(decodeU32(bytes, swap));
        if (!std::isfinite(value)) {
            throw std::runtime_error("DCD frame " + std::to_string(frameNumber) + " contains a non-finite coordinate");
        }
        values.push_back(value);
    }
    return values;
}

} // namespace

DcdTrajectorySource::DcdTrajectorySource(const std::filesystem::path& path) : path_(path) {
    std::ifstream input(path_, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open DCD trajectory: " + path_.string());

    char markerBytes[4];
    if (!input.read(markerBytes, sizeof(markerBytes))) throw std::runtime_error("truncated DCD header");
    const std::uint32_t littleMarker = decodeU32(markerBytes, false);
    if (littleMarker == 84) swapBytes_ = false;
    else if (decodeU32(markerBytes, true) == 84) swapBytes_ = true;
    else throw std::runtime_error("invalid DCD header record size");

    char header[84];
    if (!input.read(header, sizeof(header))) throw std::runtime_error("truncated DCD header");
    if (std::string_view(header, 4) != "CORD") throw std::runtime_error("DCD file is not a coordinate trajectory");
    if (readU32(input, swapBytes_, "header closing marker") != 84) {
        throw std::runtime_error("DCD header record markers do not match");
    }
    const auto control = [&](std::size_t index) { return decodeU32(header + 4 + index * 4, swapBytes_); };
    const std::uint32_t declaredFrames = control(0);
    const std::uint32_t fixedAtomCount = control(8);
    const bool hasUnitCell = control(10) != 0;
    const bool hasFourthDimension = control(11) != 0;

    const auto title = readRecord(input, swapBytes_, 16 * 1024 * 1024, "title");
    if (title.size() < 4) throw std::runtime_error("invalid DCD title record");
    const std::uint32_t titleCount = decodeU32(title.data(), swapBytes_);
    if (titleCount > (title.size() - 4) / 80) throw std::runtime_error("invalid DCD title count");

    const auto atomRecord = readRecord(input, swapBytes_, 4, "atom count");
    if (atomRecord.size() != 4) throw std::runtime_error("invalid DCD atom count record");
    atomCount_ = decodeU32(atomRecord.data(), swapBytes_);
    if (atomCount_ == 0 || atomCount_ > 100'000'000) throw std::runtime_error("invalid DCD atom count");
    if (fixedAtomCount >= atomCount_) throw std::runtime_error("invalid DCD fixed-atom count");

    if (fixedAtomCount != 0) {
        const std::size_t freeCount = atomCount_ - fixedAtomCount;
        if (freeCount > std::numeric_limits<std::uint32_t>::max() / 4) {
            throw std::runtime_error("DCD free-atom index is too large");
        }
        const auto indices = readRecord(input, swapBytes_, freeCount * 4, "free-atom index");
        if (indices.size() != freeCount * 4) throw std::runtime_error("invalid DCD free-atom index record");
        std::vector<bool> seen(atomCount_);
        freeAtomIndices_.reserve(freeCount);
        for (std::size_t index = 0; index < freeCount; ++index) {
            const std::uint32_t oneBased = decodeU32(indices.data() + index * 4, swapBytes_);
            if (oneBased == 0 || oneBased > atomCount_ || seen[oneBased - 1]) {
                throw std::runtime_error("invalid DCD free-atom index");
            }
            seen[oneBased - 1] = true;
            freeAtomIndices_.push_back(oneBased - 1);
        }
    }

    while (input.peek() != std::char_traits<char>::eof()) {
        const std::size_t frameNumber = frames_.size() + 1;
        const std::size_t coordinateCount = frames_.empty() ? atomCount_ :
                                            freeAtomIndices_.empty() ? atomCount_ : freeAtomIndices_.size();
        if (coordinateCount > std::numeric_limits<std::uint32_t>::max() / sizeof(float)) {
            throw std::runtime_error("DCD coordinate record is too large");
        }
        const std::uint32_t coordinateBytes = static_cast<std::uint32_t>(coordinateCount * sizeof(float));
        const std::uint32_t firstRecordSize = peekRecordSize(input, swapBytes_);
        if (hasUnitCell || firstRecordSize != coordinateBytes) {
            if (firstRecordSize != 24 && firstRecordSize != 48 && firstRecordSize != 56) {
                throw std::runtime_error("DCD frame " + std::to_string(frameNumber) +
                                         " has an invalid unit-cell or X-coordinate record");
            }
            static_cast<void>(readRecord(input, swapBytes_, 56,
                                         "frame " + std::to_string(frameNumber) + " unit cell"));
        }
        FrameIndex frame;
        frame.coordinateCount = coordinateCount;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            frame.coordinateOffsets[axis] = indexRecord(
                input, swapBytes_, coordinateBytes,
                "frame " + std::to_string(frameNumber) + " " + std::string(1, "XYZ"[axis]) + " coordinates");
        }
        if (hasFourthDimension) {
            static_cast<void>(indexRecord(input, swapBytes_, coordinateBytes,
                                          "frame " + std::to_string(frameNumber) + " fourth dimension"));
        }
        frames_.push_back(frame);
    }
    if (frames_.empty()) throw std::runtime_error("DCD trajectory contains no frames: " + path_.string());
    if (declaredFrames != 0 && frames_.size() != declaredFrames) {
        throw std::runtime_error("DCD header declares " + std::to_string(declaredFrames) + " frames; found " +
                                 std::to_string(frames_.size()));
    }
    if (!freeAtomIndices_.empty()) firstFrame_ = frame(0);
}

std::string_view DcdTrajectorySource::atomLabel(std::size_t index) const {
    if (index >= atomCount_) throw std::out_of_range("DCD atom index out of range");
    return {};
}

std::string_view DcdTrajectorySource::frameLabel(std::size_t index) const {
    if (index >= frames_.size()) throw std::out_of_range("DCD frame index out of range");
    return {};
}

const std::vector<Vec3>& DcdTrajectorySource::frame(std::size_t index) const {
    if (index >= frames_.size()) throw std::out_of_range("DCD frame index out of range");
    if (bufferedFrame_ == index) return frameBuffer_;
    const auto& frameIndex = frames_[index];
    std::ifstream input(path_, std::ios::binary);
    if (!input) throw std::runtime_error("cannot reopen DCD trajectory: " + path_.string());
    const auto x = readCoordinates(input, frameIndex.coordinateOffsets[0], frameIndex.coordinateCount, swapBytes_,
                                   index + 1, 'X');
    const auto y = readCoordinates(input, frameIndex.coordinateOffsets[1], frameIndex.coordinateCount, swapBytes_,
                                   index + 1, 'Y');
    const auto z = readCoordinates(input, frameIndex.coordinateOffsets[2], frameIndex.coordinateCount, swapBytes_,
                                   index + 1, 'Z');
    if (index != 0 && !freeAtomIndices_.empty()) frameBuffer_ = firstFrame_;
    else frameBuffer_.assign(atomCount_, {});
    for (std::size_t coordinate = 0; coordinate < frameIndex.coordinateCount; ++coordinate) {
        const std::size_t atom = index == 0 || freeAtomIndices_.empty() ? coordinate : freeAtomIndices_[coordinate];
        frameBuffer_[atom] = {x[coordinate], y[coordinate], z[coordinate]};
    }
    bufferedFrame_ = index;
    return frameBuffer_;
}

} // namespace ascii_protein
