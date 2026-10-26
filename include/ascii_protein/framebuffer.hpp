#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace ascii_protein {

constexpr int xterm256Color(int index) { return 256 + index; }

struct Cell {
    char32_t glyph{U' '};
    int color{0};
    double depth{std::numeric_limits<double>::infinity()};
};

class Framebuffer {
public:
    Framebuffer(int width, int height);

    void resize(int width, int height);
    void clear();
    void put(int x, int y, double depth, char32_t glyph, int color);
    void text(int x, int y, const std::string& value, int color = 0);

    int width() const { return width_; }
    int height() const { return height_; }
    const Cell& at(int x, int y) const;
    std::string toString(bool ansiColor, int backgroundColor = 0) const;
    std::string toPlainText() const;

private:
    int width_{};
    int height_{};
    std::vector<Cell> cells_;
};

} // namespace ascii_protein
