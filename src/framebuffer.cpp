#include "ascii_protein/framebuffer.hpp"

#include <algorithm>
#include <stdexcept>

namespace ascii_protein {
namespace {

std::string utf8(char32_t codepoint) {
    std::string output;
    if (codepoint <= 0x7f) {
        output.push_back(static_cast<char>(codepoint));
    } else if (codepoint <= 0x7ff) {
        output.push_back(static_cast<char>(0xc0 | (codepoint >> 6)));
        output.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else if (codepoint <= 0xffff) {
        output.push_back(static_cast<char>(0xe0 | (codepoint >> 12)));
        output.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
        output.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else {
        output.push_back(static_cast<char>(0xf0 | (codepoint >> 18)));
        output.push_back(static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f)));
        output.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
        output.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    }
    return output;
}

} // namespace

Framebuffer::Framebuffer(int width, int height) { resize(width, height); }

void Framebuffer::resize(int width, int height) {
    width_ = std::max(1, width);
    height_ = std::max(1, height);
    cells_.resize(static_cast<std::size_t>(width_ * height_));
    clear();
}

void Framebuffer::clear() { std::fill(cells_.begin(), cells_.end(), Cell{}); }

void Framebuffer::put(int x, int y, double depth, char32_t glyph, int color) {
    if (x < 0 || y < 0 || x >= width_ || y >= height_) {
        return;
    }
    auto& cell = cells_[static_cast<std::size_t>(y * width_ + x)];
    if (depth <= cell.depth) {
        cell = {glyph, color, depth};
    }
}

void Framebuffer::text(int x, int y, const std::string& value, int color) {
    for (const char character : value) {
        put(x++, y, -std::numeric_limits<double>::infinity(), static_cast<unsigned char>(character), color);
    }
}

const Cell& Framebuffer::at(int x, int y) const {
    if (x < 0 || y < 0 || x >= width_ || y >= height_) {
        throw std::out_of_range("framebuffer coordinates out of range");
    }
    return cells_[static_cast<std::size_t>(y * width_ + x)];
}

std::string Framebuffer::toString(bool ansiColor, int backgroundColor) const {
    std::string output;
    output.reserve(cells_.size() * 2);
    int activeColor = -1;
    for (int y = 0; y < height_; ++y) {
        if (ansiColor && backgroundColor != 0) output += "\x1b[" + std::to_string(backgroundColor) + "m";
        for (int x = 0; x < width_; ++x) {
            const auto& cell = at(x, y);
            if (ansiColor && cell.color != activeColor) {
                output += cell.color == 0
                              ? "\x1b[0m" + (backgroundColor == 0 ? std::string{} :
                                               "\x1b[" + std::to_string(backgroundColor) + "m")
                              : "\x1b[" + std::to_string(cell.color) + "m";
                activeColor = cell.color;
            }
            output += utf8(cell.glyph);
        }
        if (ansiColor && activeColor != 0) {
            output += "\x1b[0m";
            activeColor = 0;
        }
        if (y + 1 < height_) {
            output.push_back('\n');
        }
    }
    return output;
}

std::string Framebuffer::toPlainText() const { return toString(false); }

} // namespace ascii_protein