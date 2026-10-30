#include "ascii_protein/density_map.hpp"
#include <bit>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <algorithm>
namespace ascii_protein {
DensityMap loadDensityMap(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open density map");
    std::array<unsigned char, 1024> header{};
    if (!file.read(reinterpret_cast<char*>(header.data()), 1024))
        throw std::runtime_error("truncated MRC header");
    bool big = header[212] == 0x11;
    auto word = [&](int offset) {
        std::uint32_t x = 0;
        for (int i = 0; i < 4; ++i)
            x |= static_cast<std::uint32_t>(header[offset + i]) << (big ? (3 - i) * 8 : i * 8);
        return x;
    };
    auto integer = [&](int offset) { return std::bit_cast<std::int32_t>(word(offset)); };
    auto real = [&](int offset) { return std::bit_cast<float>(word(offset)); };
    if (header[208] != 'M' || header[209] != 'A' || header[210] != 'P' || header[211] != ' ')
        throw std::runtime_error("invalid MRC MAP signature");
    DensityMap map;
    map.path = std::filesystem::absolute(path);
    std::size_t count = 1;
    for (int i = 0; i < 3; ++i) {
        const auto n = integer(i * 4);
        if (n <= 0 || n > 100000 || count > 134217728 / static_cast<std::size_t>(n))
            throw std::runtime_error("density map exceeds 512 MiB value limit or has invalid dimensions");
        map.dimensions[i] = n;
        count *= static_cast<std::size_t>(n);
    }
    const int mode = integer(12);
    int bytes = mode == 0 ? 1 : (mode == 1 || mode == 6 ? 2 : (mode == 2 ? 4 : 0));
    if (!bytes) throw std::runtime_error("unsupported MRC mode (use 0, 1, 2 or 6)");
    const int extended = integer(92);
    if (extended < 0) throw std::runtime_error("invalid MRC extended header size");
    file.seekg(0, std::ios::end);
    auto fileBytes = file.tellg();
    if (fileBytes < 0 ||
        static_cast<std::uint64_t>(fileBytes) < 1024ULL + static_cast<unsigned>(extended) + count * bytes)
        throw std::runtime_error("truncated MRC density data");
    const double a = real(40), b = real(44), c = real(48);
    const double alpha = radians(real(52)), beta = radians(real(56)), gamma = radians(real(60));
    if (!(a > 0 && b > 0 && c > 0) || !std::isfinite(a + b + c + alpha + beta + gamma) || alpha <= 0 ||
        alpha >= pi || beta <= 0 || beta >= pi || gamma <= 0 || gamma >= pi ||
        std::abs(std::sin(gamma)) < 1e-8)
        throw std::runtime_error("invalid MRC unit cell");
    const double cx = c * std::cos(beta),
                 cy = c * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    const double cz2 = c * c - cx * cx - cy * cy;
    if (cz2 <= 0) throw std::runtime_error("invalid MRC cell angles");
    std::array<Vec3, 3> basis{Vec3{a, 0, 0}, Vec3{b * std::cos(gamma), b * std::sin(gamma), 0},
                              Vec3{cx, cy, std::sqrt(cz2)}};
    std::array<bool, 3> used{};
    for (int i = 0; i < 3; ++i) {
        int axis = integer(64 + i * 4) - 1;
        if (axis < 0 || axis > 2 || used[axis]) throw std::runtime_error("invalid MRC axis mapping");
        used[axis] = true;
        int sampling = integer(28 + axis * 4);
        if (sampling <= 0) throw std::runtime_error("invalid MRC grid sampling");
        map.steps[i] = basis[axis] / sampling;
    }
    map.origin = {real(196), real(200), real(204)};
    if (!std::isfinite(map.origin.x + map.origin.y + map.origin.z))
        throw std::runtime_error("invalid MRC origin");
    // Explicit ORIGIN takes precedence; CCP4 start indices apply when ORIGIN is zero.
    if (length(map.origin) == 0)
        for (int i = 0; i < 3; ++i)
            map.origin += map.steps[i] * integer(16 + i * 4);
    map.values.resize(count);
    file.seekg(1024LL + extended);
    std::array<unsigned char, 65536> buffer{};
    for (std::size_t offset = 0; offset < count;) {
        const auto block = std::min(buffer.size() / static_cast<std::size_t>(bytes), count - offset);
        if (!file.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(block * bytes)))
            throw std::runtime_error("truncated MRC density data");
        for (std::size_t n = 0; n < block; ++n) {
            const auto* data = buffer.data() + n * bytes;
            std::uint32_t x = 0;
            for (int i = 0; i < bytes; ++i)
                x |= static_cast<std::uint32_t>(data[i]) << (big ? (bytes - i - 1) * 8 : i * 8);
            auto& value = map.values[offset + n];
            if (mode == 0)
                value = static_cast<float>(std::bit_cast<std::int8_t>(static_cast<std::uint8_t>(x)));
            else if (mode == 1)
                value = static_cast<float>(std::bit_cast<std::int16_t>(static_cast<std::uint16_t>(x)));
            else if (mode == 6)
                value = static_cast<float>(x);
            else
                value = std::bit_cast<float>(x);
            if (!std::isfinite(value)) throw std::runtime_error("non-finite MRC density value");
        }
        offset += block;
    }
    auto [lo, hi] = std::minmax_element(map.values.begin(), map.values.end());
    map.minimum = *lo;
    map.maximum = *hi;
    return map;
}
} // namespace ascii_protein
