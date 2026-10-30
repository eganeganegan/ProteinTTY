#include "ascii_protein/renderer.hpp"
#include "ascii_protein/molecular.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <tuple>
#include <unordered_map>

namespace ascii_protein {
namespace {

double covalentRadius(const std::string& element) {
    static const std::unordered_map<std::string, double> radii{
        {"H", 0.31}, {"C", 0.76}, {"N", 0.71}, {"O", 0.66}, {"S", 1.05}, {"P", 1.07}};
    const auto found = radii.find(element);
    return found == radii.end() ? 0.77 : found->second;
}

double displayRadius(const std::string& element) {
    if (element == "H") return 0.65;
    if (element == "S" || element == "P") return 1.2;
    return 1.0;
}

char32_t shadeGlyph(Charset charset, double brightness) {
    brightness = std::clamp(brightness, 0.0, 0.999);
    if (charset == Charset::ascii) {
        constexpr std::string_view ramp = ".,:;irsXA253hMHGS#9B&@";
        return static_cast<unsigned char>(ramp[static_cast<std::size_t>(brightness * ramp.size())]);
    }
    constexpr std::array<char32_t, 6> ramp{U'\u00b7', U'\u2022', U'\u2591', U'\u2592', U'\u2593', U'\u2588'};
    return ramp[static_cast<std::size_t>(brightness * ramp.size())];
}

bool residuesAreAdjacent(const Atom& first, const Atom& second) {
    return first.chain == second.chain && std::abs(first.residueNumber - second.residueNumber) <= 1;
}

bool atomShown(const Atom& atom, RenderMode mode) {
    if (mode == RenderMode::atoms || mode == RenderMode::ballAndStick) return true;
    if (mode == RenderMode::backbone) return isBackboneAtom(atom);
    return atom.name == "CA";
}

std::string secondaryName(SecondaryStructure secondary) {
    if (secondary == SecondaryStructure::helix) return "helix";
    if (secondary == SecondaryStructure::sheet) return "sheet";
    return "coil";
}

std::string displayChain(const std::string& chain) { return chain.empty() ? "-" : chain; }

struct GridCell {
    int x{};
    int y{};
    int z{};

    bool operator==(const GridCell&) const = default;
};

struct GridCellHash {
    std::size_t operator()(const GridCell& cell) const {
        std::size_t hash = std::hash<int>{}(cell.x);
        hash ^= std::hash<int>{}(cell.y) + 0x9e3779b9U + (hash << 6U) + (hash >> 2U);
        hash ^= std::hash<int>{}(cell.z) + 0x9e3779b9U + (hash << 6U) + (hash >> 2U);
        return hash;
    }
};

GridCell gridCell(const Vec3& position, double cellSize) {
    return {static_cast<int>(std::floor(position.x / cellSize)),
            static_cast<int>(std::floor(position.y / cellSize)),
            static_cast<int>(std::floor(position.z / cellSize))};
}

Vec3 catmullRom(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, double amount) {
    const double squared = amount * amount;
    const double cubed = squared * amount;
    return (p1 * 2.0 + (p2 - p0) * amount +
            (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * squared +
            (p3 - p0 + (p1 - p2) * 3.0) * cubed) * 0.5;
}

struct DepthPalette {
    int near;
    int middle;
    int far;
};

DepthPalette depthPalette(int color) {
    switch (color) {
    case 30: return {240, 236, 232};
    case 31: case 91: return {203, 160, 88};
    case 32: case 92: return {83, 34, 22};
    case 33: case 93: return {227, 178, 58};
    case 34: case 94: return {75, 27, 17};
    case 35: case 95: return {207, 129, 53};
    case 36: case 96: return {87, 37, 23};
    case 37: case 97: return {255, 250, 244};
    case 90: return {250, 244, 238};
    default: return {255, 250, 244};
    }
}

} // namespace

std::vector<Bond> inferBonds(const Structure& structure) {
    constexpr double cellSize = 2.7;
    std::vector<Bond> bonds;
    std::vector<bool> connected(structure.atoms.size());
    for (auto [a, b] : depositedBonds(structure)) {
        bonds.push_back({a, b});
        connected[a] = connected[b] = true;
    }
    std::unordered_map<GridCell, std::vector<std::size_t>, GridCellHash> grid;
    grid.reserve(structure.atoms.size());
    for (std::size_t second = 0; second < structure.atoms.size(); ++second) {
        const auto& right = structure.atoms[second];
        const GridCell cell = gridCell(right.position, cellSize);
        for (int offsetZ = -1; offsetZ <= 1; ++offsetZ) {
            for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                    const auto found = grid.find({cell.x + offsetX, cell.y + offsetY, cell.z + offsetZ});
                    if (found == grid.end()) continue;
                    for (const std::size_t first : found->second) {
                        const auto& left = structure.atoms[first];
                        if (connected[first] || connected[second]) continue;
                        if (!residuesAreAdjacent(left, right)) continue;
                        const double threshold = covalentRadius(left.element) + covalentRadius(right.element) + 0.55;
                        const double distanceSquared = dot(left.position - right.position, left.position - right.position);
                        if (distanceSquared > 0.35 * 0.35 && distanceSquared <= threshold * threshold) {
                            bonds.push_back({first, second});
                        }
                    }
                }
            }
        }
        grid[cell].push_back(second);
    }
    std::sort(bonds.begin(), bonds.end(), [](const Bond& left, const Bond& right) {
        return std::tie(left.first, left.second) < std::tie(right.first, right.second);
    });
    return bonds;
}

Renderer::Renderer(const Structure& structure) : structure_(structure), bonds_(inferBonds(structure)) {
    if (!structure.atoms.empty()) {
        const auto [minimum, maximum] = std::minmax_element(
            structure.atoms.begin(), structure.atoms.end(),
            [](const Atom& left, const Atom& right) { return left.bFactor < right.bFactor; });
        minBFactor_ = minimum->bFactor;
        maxBFactor_ = maximum->bFactor;
    }
}

Renderer::RenderPoint Renderer::transform(const Vec3& position, const Camera& camera, int width, int height,
                                          RendererBackend backend) const {
    const Vec3 rotated = camera.rotation() * (position - camera.target);
    const double availableWidth = std::max(1, width - 2);
    const double availableHeight = std::max(1, height - 2);
    const double fitScale = std::min(availableWidth / (structure_.boundingRadius * 2.0),
                                     availableHeight / structure_.boundingRadius) * 0.82 * camera.zoom;
    const double centerX = width * 0.5 + camera.pan.x;
    const double centerY = height * 0.5 + camera.pan.y;
    const double cameraDistance = structure_.boundingRadius * 3.2;
    const double depth = cameraDistance - rotated.z;
    double scale = fitScale;
    bool visible = true;
    if (camera.projection == Projection::perspective) {
        visible = depth > structure_.boundingRadius * 0.1;
        scale *= cameraDistance / std::max(depth, 0.01);
    }
    const double verticalAspect = backend == RendererBackend::braille ? 1.0 : 0.5;
    return {{centerX + rotated.x * scale, centerY - rotated.y * scale * verticalAspect, depth, visible}, rotated};
}

ProjectedPoint Renderer::project(const Vec3& position, const Camera& camera, int width, int height) const {
    return transform(position, camera, width, height, RendererBackend::cells).projected;
}

int Renderer::atomColor(std::size_t atomIndex, const RenderSettings& settings) const {
    const auto& atom = structure_.atoms[atomIndex];
    if (settings.inspectedAtom && *settings.inspectedAtom < structure_.atoms.size() &&
        atomIndex == *settings.inspectedAtom) {
        return 96;
    }
    if (atomIndex < settings.atomColors.size() && settings.atomColors[atomIndex] != 0) {
        return settings.atomColors[atomIndex];
    }
    if (atomIndex < settings.selectedAtoms.size() && settings.selectedAtoms[atomIndex]) return 96;
    if (settings.selectedResidue && atom.residueIndex == *settings.selectedResidue) {
        return 97;
    }
    if (settings.colorMode == ColorMode::chain) {
        constexpr std::array<int, 8> colors{36, 33, 35, 32, 93, 34, 91, 96};
        return colors[std::hash<std::string>{}(atom.chain) % colors.size()];
    }
    if (settings.colorMode == ColorMode::bFactor) {
        const double range = std::max(0.001, maxBFactor_ - minBFactor_);
        const double value = (atom.bFactor - minBFactor_) / range;
        if (value < 0.25) return 34;
        if (value < 0.5) return 36;
        if (value < 0.75) return 93;
        return 91;
    }
    if (settings.colorMode == ColorMode::secondary) {
        const auto secondary = structure_.residues[atom.residueIndex].secondary;
        if (secondary == SecondaryStructure::helix) return 91;
        if (secondary == SecondaryStructure::sheet) return 93;
        return 36;
    }
    if (settings.colorMode == ColorMode::residue) {
        const std::string code = oneLetterCode(atom.residueName);
        if (std::string("AVILMFWY").find(code) != std::string::npos) return 33;
        if (std::string("STNQ").find(code) != std::string::npos) return 36;
        if (std::string("DE").find(code) != std::string::npos) return 91;
        if (std::string("RKH").find(code) != std::string::npos) return 94;
        return 95;
    }
    if (atom.element == "O") return 91;
    if (atom.element == "N") return 94;
    if (atom.element == "S") return 93;
    if (atom.element == "P") return 95;
    if (atom.element == "C") return 37;
    if (atom.element == "H") return 90;
    return 97;
}

int Renderer::depthColor(int color, double depth) const {
    if (color == 0 || color >= xterm256Color(0)) return color;
    const double radius = std::max(1.0, structure_.boundingRadius);
    const double distance = std::clamp((depth / radius - 2.15) / 2.1, 0.0, 1.0);
    const auto palette = depthPalette(color);
    const int index = distance < 0.33 ? palette.near : distance < 0.68 ? palette.middle : palette.far;
    return xterm256Color(index);
}

std::optional<std::size_t> Renderer::nearestAtom(const Camera& camera, const RenderSettings& settings, int width,
                                                 int height) const {
    const int panelWidth = settings.showInfo && width >= 72 ? 28 : 0;
    const int viewWidth = width - panelWidth;
    const int rasterWidth = settings.backend == RendererBackend::braille ? viewWidth * 2 : viewWidth;
    const int rasterHeight = settings.backend == RendererBackend::braille ? height * 4 : height;
    const double centerX = rasterWidth * 0.5;
    const double centerY = rasterHeight * 0.5;
    double nearestDistance = std::numeric_limits<double>::infinity();
    std::optional<std::size_t> nearest;
    for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
        if (index < settings.visibleAtoms.size() && !settings.visibleAtoms[index]) continue;
        if (!atomShown(structure_.atoms[index], settings.mode)) continue;
        const auto projected = transform(structure_.atoms[index].position, camera, rasterWidth, rasterHeight,
                                         settings.backend).projected;
        if (!projected.visible) continue;
        const double distance = std::pow(projected.x - centerX, 2.0) + std::pow(projected.y - centerY, 2.0);
        if (distance < nearestDistance) {
            nearestDistance = distance;
            nearest = index;
        }
    }
    return nearest;
}

void Renderer::drawAtom(std::size_t atomIndex, const Camera& camera, const RenderSettings& settings, int viewWidth,
                        int height) {
    if (atomIndex < settings.visibleAtoms.size() && !settings.visibleAtoms[atomIndex]) return;
    const auto& atom = structure_.atoms[atomIndex];
    const auto point = transform(atom.position, camera, viewWidth, height, settings.backend);
    if (!point.projected.visible) return;

    const double normalizedDepth = 1.0 - std::clamp(point.projected.depth / (structure_.boundingRadius * 6.4), 0.0, 1.0);
    const double resolutionScale = settings.backend == RendererBackend::braille ? 2.0 : 1.0;
    const int radius = std::clamp(
        static_cast<int>(std::round(displayRadius(atom.element) * camera.zoom * resolutionScale)), 0, 6);
    const int color = atomColor(atomIndex, settings);
    if (point.projected.x < -2 * radius - 1 || point.projected.x > viewWidth + 2 * radius ||
        point.projected.y < -radius - 1 || point.projected.y > height + radius)
        return;
    const int centerX = static_cast<int>(std::round(point.projected.x));
    const int centerY = static_cast<int>(std::round(point.projected.y));
    if (centerX + 2 * radius < 0 || centerX - 2 * radius >= viewWidth || centerY + radius < 0 ||
        centerY - radius >= height)
        return;
    struct SpritePixel {
        int x, y;
        double surface;
    };
    // Atom silhouette geometry depends only on the seven supported pixel radii.
    static const auto sprites = [] {
        std::array<std::vector<SpritePixel>, 7> result;
        for (int r = 0; r <= 6; ++r)
            for (int y = -r; y <= r; ++y)
                for (int x = -2 * r; x <= 2 * r; ++x) {
                    const double distance = std::sqrt(std::pow(x / 2.0, 2.0) + std::pow(y, 2.0));
                    if (distance > r + 0.25) continue;
                    const double surface =
                        r == 0 ? 1.0 : std::sqrt(std::max(0.0, 1.0 - std::pow(distance / (r + 0.25), 2.0)));
                    result[r].push_back({x, y, surface});
                }
        return result;
    }();
    for (const auto& pixel : sprites[radius]) {
        const double lighting =
            std::clamp(0.2 + normalizedDepth * 0.45 + pixel.surface * 0.35 - pixel.y * 0.06, 0.0, 1.0);
        const double depth = point.projected.depth - pixel.surface * 0.2;
        putPixel(centerX + pixel.x, centerY + pixel.y, depth, shadeGlyph(settings.charset, lighting),
                 depthColor(color, depth));
    }
}

void Renderer::drawLine(std::size_t first, std::size_t second, const Camera& camera, const RenderSettings& settings,
                        int viewWidth, int height, int thickness) {
    if ((first < settings.visibleAtoms.size() && !settings.visibleAtoms[first]) ||
        (second < settings.visibleAtoms.size() && !settings.visibleAtoms[second])) return;
    drawLine(structure_.atoms[first].position, structure_.atoms[second].position, first, camera, settings,
             viewWidth, height, thickness);
}

void Renderer::drawLine(const Vec3& startPosition, const Vec3& endPosition, std::size_t colorAtom,
                        const Camera& camera, const RenderSettings& settings, int viewWidth, int height,
                        int thickness) {
    const auto start = transform(startPosition, camera, viewWidth, height, settings.backend).projected;
    const auto end = transform(endPosition, camera, viewWidth, height, settings.backend).projected;
    if (!start.visible || !end.visible) return;
    const double deltaX = end.x - start.x;
    const double deltaY = end.y - start.y;
    const int steps = std::max(1, static_cast<int>(std::ceil(std::max(std::abs(deltaX), std::abs(deltaY)))));
    const int color = atomColor(colorAtom, settings);
    for (int step = 0; step <= steps; ++step) {
        const double amount = static_cast<double>(step) / steps;
        const int x = static_cast<int>(std::round(start.x + deltaX * amount));
        const int y = static_cast<int>(std::round(start.y + deltaY * amount));
        const double depth = start.depth + (end.depth - start.depth) * amount + 0.05;
        for (int offsetY = -thickness; offsetY <= thickness; ++offsetY) {
            for (int offsetX = -thickness; offsetX <= thickness; ++offsetX) {
                if (std::abs(offsetX) + std::abs(offsetY) > thickness) continue;
                putPixel(x + offsetX, y + offsetY, depth, settings.charset == Charset::ascii ? U'.' : U'\u00b7',
                         depthColor(color, depth));
            }
        }
    }
}

void Renderer::drawAlphaTrace(const Camera& camera, const RenderSettings& settings, int viewWidth, int height) {
    std::vector<std::size_t> run;
    const int thickness = renderingBraille_ ? 1 : 0;
    const auto drawRun = [&] {
        if (run.size() < 2) {
            run.clear();
            return;
        }
        for (std::size_t segment = 0; segment + 1 < run.size(); ++segment) {
            const Vec3& p0 = structure_.atoms[run[segment == 0 ? segment : segment - 1]].position;
            const Vec3& p1 = structure_.atoms[run[segment]].position;
            const Vec3& p2 = structure_.atoms[run[segment + 1]].position;
            const Vec3& p3 = structure_.atoms[run[segment + 2 < run.size() ? segment + 2 : segment + 1]].position;
            Vec3 previous = p1;
            constexpr int subdivisions = 6;
            for (int step = 1; step <= subdivisions; ++step) {
                const Vec3 current = catmullRom(p0, p1, p2, p3,
                                                static_cast<double>(step) / subdivisions);
                drawLine(previous, current, run[segment], camera, settings, viewWidth, height, thickness);
                previous = current;
            }
        }
        run.clear();
    };

    for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
        const auto& atom = structure_.atoms[index];
        if (atom.name != "CA") continue;
        const bool hidden = index < settings.visibleAtoms.size() && !settings.visibleAtoms[index];
        if (hidden || (!run.empty() && !residuesAreAdjacent(structure_.atoms[run.back()], atom))) drawRun();
        if (!hidden) run.push_back(index);
    }
    drawRun();
}

void Renderer::putPixel(int x, int y, double depth, char32_t glyph, int color) {
    if (!renderingBraille_) {
        framebuffer_.put(x, y, depth, glyph, color);
        return;
    }
    if (x < 0 || y < 0 || x >= brailleWidth_ * 2 || y >= brailleHeight_ * 4) return;
    static constexpr int dotBits[4][2]{{0, 3}, {1, 4}, {2, 5}, {6, 7}};
    const int bit = dotBits[y % 4][x % 2];
    auto& cell = brailleCells_[static_cast<std::size_t>((y / 4) * brailleWidth_ + x / 2)];
    if (depth <= cell.depths[static_cast<std::size_t>(bit)]) {
        cell.depths[static_cast<std::size_t>(bit)] = depth;
        cell.colors[static_cast<std::size_t>(bit)] = color;
    }
}

void Renderer::prepareBraille(int width, int height) {
    brailleWidth_ = width;
    brailleHeight_ = height;
    BrailleCell empty;
    empty.depths.fill(std::numeric_limits<double>::infinity());
    brailleCells_.assign(static_cast<std::size_t>(width * height), empty);
}

void Renderer::composeBraille(int viewWidth, int height) {
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < viewWidth; ++x) {
            const auto& cell = brailleCells_[static_cast<std::size_t>(y * viewWidth + x)];
            unsigned int mask = 0;
            double nearest = std::numeric_limits<double>::infinity();
            int color = 0;
            for (std::size_t bit = 0; bit < cell.depths.size(); ++bit) {
                if (!std::isfinite(cell.depths[bit])) continue;
                mask |= 1U << bit;
                if (cell.depths[bit] < nearest) {
                    nearest = cell.depths[bit];
                    color = cell.colors[bit];
                }
            }
            if (mask != 0) framebuffer_.put(x, y, nearest, static_cast<char32_t>(0x2800U + mask), color);
        }
    }
}

void Renderer::drawPanels(const Camera& camera, const RenderSettings& settings, int viewWidth) {
    if (settings.showInfo) {
        const int x = viewWidth + 2;
        int y = 1;
        framebuffer_.text(x, y++, "PROTEIN TTY", 96);
        framebuffer_.text(x, y++, std::string(18, '-'), 90);
        framebuffer_.text(x, y++, "File: " + structure_.sourceName, 37);
        framebuffer_.text(x, y++, "Atoms: " + std::to_string(structure_.atoms.size()), 37);
        framebuffer_.text(x, y++, "Residues: " + std::to_string(structure_.residues.size()), 37);
        framebuffer_.text(x, y++, "Chains: " + std::to_string(structure_.chains.size()), 37);
        ++y;
        framebuffer_.text(x, y++, "Mode: " + renderModeName(settings.mode), 37);
        framebuffer_.text(x, y++, "Color: " + colorModeName(settings.colorMode), 37);
        framebuffer_.text(x, y++, "View: " + projectionName(camera.projection), 37);
        std::ostringstream angles;
        angles << "Yaw " << std::fixed << std::setprecision(0) << camera.yaw << "  Pitch " << camera.pitch;
        framebuffer_.text(x, y++, angles.str(), 37);
        std::ostringstream zoom;
        zoom << "Zoom: " << std::fixed << std::setprecision(1) << camera.zoom << 'x';
        framebuffer_.text(x, y++, zoom.str(), 37);
        if (!settings.trajectoryInfo.empty()) framebuffer_.text(x, y++, settings.trajectoryInfo, 93);
        ++y;
        if (settings.inspectedAtom && *settings.inspectedAtom < structure_.atoms.size()) {
            const auto& atom = structure_.atoms[*settings.inspectedAtom];
            const Vec3 original = atom.position + structure_.originalCentroid;
            framebuffer_.text(x, y++, "INSPECTED ATOM", 96);
            framebuffer_.text(x, y++, atom.name + "  " + atom.element + "  " + atom.residueName +
                                             std::to_string(atom.residueNumber), 97);
            framebuffer_.text(x, y++, "Chain " + displayChain(atom.chain) +
                                             "  Serial " + std::to_string(atom.serial), 37);
            std::ostringstream position;
            position << std::fixed << std::setprecision(1) << "XYZ " << original.x << ' ' << original.y << ' '
                     << original.z;
            framebuffer_.text(x, y++, position.str(), 37);
            std::ostringstream factors;
            factors << std::fixed << std::setprecision(1) << "B " << atom.bFactor << "  Occ " << atom.occupancy;
            framebuffer_.text(x, y++, factors.str(), 37);
        } else if (settings.selectedResidue && *settings.selectedResidue < structure_.residues.size()) {
            const auto& residue = structure_.residues[*settings.selectedResidue];
            framebuffer_.text(x, y++, "SELECTED", 93);
            framebuffer_.text(x, y++, "Chain " + displayChain(residue.chain) + "  " +
                                             residue.name + std::to_string(residue.number), 97);
            framebuffer_.text(x, y++, "Atoms: " + std::to_string(residue.atomIndices.size()), 37);
            double average = 0.0;
            for (const auto atomIndex : residue.atomIndices) average += structure_.atoms[atomIndex].bFactor;
            average /= std::max<std::size_t>(1, residue.atomIndices.size());
            std::ostringstream factor;
            factor << "B-factor: " << std::fixed << std::setprecision(1) << average;
            framebuffer_.text(x, y++, factor.str(), 37);
            framebuffer_.text(x, y++, "Type: " + secondaryName(residue.secondary), 37);
        }
        framebuffer_.text(x, framebuffer_.height() - 2, "[h] help  [q] quit", 90);
    }

    if (settings.showSequence && settings.selectedResidue && *settings.selectedResidue < structure_.residues.size()) {
        const int y = framebuffer_.height() - 4;
        const std::size_t selected = *settings.selectedResidue;
        const std::size_t start = selected > 25 ? selected - 25 : 0;
        const std::size_t end = std::min(structure_.residues.size(), start + static_cast<std::size_t>(std::max(1, viewWidth - 4)));
        std::string sequence;
        for (std::size_t index = start; index < end; ++index) sequence += oneLetterCode(structure_.residues[index].name);
        framebuffer_.text(2, y, "Chain " + displayChain(structure_.residues[selected].chain) + "  " + sequence, 37);
        framebuffer_.text(9 + static_cast<int>(selected - start), y + 1, "^", 93);
    }
}

void Renderer::drawHelp() {
    const int boxWidth = std::min(54, framebuffer_.width() - 4);
    const int boxHeight = std::min(22, framebuffer_.height() - 2);
    const int left = std::max(0, (framebuffer_.width() - boxWidth) / 2);
    const int top = std::max(0, (framebuffer_.height() - boxHeight) / 2);
    for (int y = top; y < top + boxHeight; ++y) {
        framebuffer_.text(left, y, std::string(static_cast<std::size_t>(boxWidth), ' '), 0);
    }
    const std::array lines{
        "PROTEIN TTY CONTROLS", "", "Arrows / WASD   rotate", "Z / E           roll", "+ / -           zoom",
        "R               reset camera", "P               projection", "1..5            render mode",
        "C               color mode", "[ / ]           select residue", "O               inspect center atom",
        "I / S           info / sequence", "` / , / .       play / prev / next", ":               command console",
        "Space           autorotate", "X               save snapshot", "H / Escape      close help", "Q               quit"};
    int y = top + 2;
    for (const auto* line : lines) {
        if (y >= top + boxHeight - 1) break;
        const int color = y == top + 2 ? 96 : 37;
        framebuffer_.text(left + 3, y, line, color);
        ++y;
    }
}

void Renderer::drawCommandStatus(const RenderSettings& settings) {
    if (!settings.commandActive && settings.commandMessage.empty()) return;
    const int y = framebuffer_.height() - 1;
    framebuffer_.text(0, y, std::string(static_cast<std::size_t>(framebuffer_.width()), ' '));
    const std::string text = settings.commandActive ? ":" + settings.commandLine : settings.commandMessage;
    framebuffer_.text(0, y, text, settings.commandActive ? 96 : 93);
}

const Framebuffer& Renderer::render(const Camera& camera, const RenderSettings& settings, int width, int height) {
    framebuffer_.resize(width, height);
    const int panelWidth = settings.showInfo && width >= 72 ? 28 : 0;
    const int viewWidth = width - panelWidth;
    renderingBraille_ = settings.backend == RendererBackend::braille;
    if (renderingBraille_) prepareBraille(viewWidth, height);
    const int rasterWidth = renderingBraille_ ? viewWidth * 2 : viewWidth;
    const int rasterHeight = renderingBraille_ ? height * 4 : height;

    if (settings.densityMap) {
        const auto& map = *settings.densityMap;
        const int axis = settings.mapAxis, u = (axis + 1) % 3, v = (axis + 2) % 3;
        const int stride =
            std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(map.dimensions[u]) *
                                                             map.dimensions[v] / 200000.0))));
        for (int j = 0; j < map.dimensions[v]; j += stride)
            for (int i = 0; i < map.dimensions[u]; i += stride) {
                int index[3]{};
                index[axis] = settings.mapSlice;
                index[u] = i;
                index[v] = j;
                const auto value = map.value(index[0], index[1], index[2]);
                if (value < settings.mapLevel) continue;
                auto p = transform(map.position(index[0], index[1], index[2]) - structure_.originalCentroid,
                                   camera, rasterWidth, rasterHeight, settings.backend)
                             .projected;
                if (p.visible && p.x >= 0 && p.x < rasterWidth && p.y >= 0 && p.y < rasterHeight)
                    putPixel(static_cast<int>(p.x), static_cast<int>(p.y), p.depth, U'.', 36);
            }
    }

    if (settings.mode == RenderMode::ballAndStick) {
        for (const auto& bond : bonds_) drawLine(bond.first, bond.second, camera, settings, rasterWidth, rasterHeight);
    } else if (settings.mode == RenderMode::backbone) {
        std::optional<std::size_t> previous;
        for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
            if (!isBackboneAtom(structure_.atoms[index])) continue;
            if (previous && residuesAreAdjacent(structure_.atoms[*previous], structure_.atoms[index])) {
                drawLine(*previous, index, camera, settings, rasterWidth, rasterHeight);
            }
            previous = index;
        }
    } else if (settings.mode == RenderMode::alphaTrace) {
        drawAlphaTrace(camera, settings, rasterWidth, rasterHeight);
    }

    for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
        const auto& atom = structure_.atoms[index];
        if (atomShown(atom, settings.mode)) drawAtom(index, camera, settings, rasterWidth, rasterHeight);
    }
    if (renderingBraille_) composeBraille(viewWidth, height);
    renderingBraille_ = false;
    if (settings.residueLabels || settings.chainLabels) {
        std::vector<bool> occupied(static_cast<std::size_t>(viewWidth) * height);
        auto label = [&](std::size_t atomIndex, const std::string& text) {
            if (atomIndex < settings.visibleAtoms.size() && !settings.visibleAtoms[atomIndex]) return;
            auto p = transform(structure_.atoms[atomIndex].position, camera, rasterWidth, rasterHeight,
                               settings.backend)
                         .projected;
            if (settings.backend == RendererBackend::braille) {
                p.x /= 2;
                p.y /= 4;
            }
            if (!p.visible || p.x < 0 || p.x >= viewWidth || p.y < 0 || p.y >= height) return;
            const int px = static_cast<int>(p.x), py = static_cast<int>(p.y);
            const int size = static_cast<int>(text.size());
            const std::array<std::pair<int, int>, 8> candidates{{{px + 4, py},
                                                                 {px + 2, py - 2},
                                                                 {px + 2, py + 2},
                                                                 {px - size - 4, py},
                                                                 {px - size / 2, py - 3},
                                                                 {px - size / 2, py + 3},
                                                                 {px + 6, py - 1},
                                                                 {px - size - 6, py + 1}}};
            for (auto [x, y] : candidates) {
                if (x < 0 || y < 0 || y >= height - 2 || x + size > viewWidth) continue;
                bool clear = true;
                for (std::size_t i = 0; i < text.size(); ++i) {
                    auto index = static_cast<std::size_t>(y) * viewWidth + x + i;
                    if (occupied[index] || framebuffer_.at(x + static_cast<int>(i), y).glyph != U' ') {
                        clear = false;
                        break;
                    }
                }
                if (!clear) continue;
                framebuffer_.text(x, y, text, 37);
                for (std::size_t i = 0; i < text.size(); ++i)
                    occupied[static_cast<std::size_t>(y) * viewWidth + x + i] = true;
                break;
            }
        };
        if (settings.chainLabels)
            for (const auto& chain : structure_.chains)
                if (!chain.residueIndices.empty()) {
                    const auto& residue = structure_.residues[chain.residueIndices.front()];
                    if (!residue.atomIndices.empty()) label(residue.atomIndices.front(), "chain " + chain.id);
                }
        if (settings.residueLabels)
            for (const auto& residue : structure_.residues)
                if (!residue.atomIndices.empty())
                    label(residue.atomIndices.front(),
                          residue.chain + ":" + residue.name + std::to_string(residue.number) +
                              (residue.insertionCode == ' ' ? "" : std::string(1, residue.insertionCode)));
    }
    drawPanels(camera, settings, viewWidth);
    if (settings.showHelp) drawHelp();
    drawCommandStatus(settings);
    return framebuffer_;
}

std::string renderModeName(RenderMode mode) {
    constexpr std::array names{"atoms", "backbone", "C-alpha trace", "ball-and-stick", "residues"};
    return names[static_cast<std::size_t>(mode)];
}

std::string colorModeName(ColorMode mode) {
    constexpr std::array names{"element", "chain", "residue", "B-factor", "secondary"};
    return names[static_cast<std::size_t>(mode)];
}

std::string projectionName(Projection projection) {
    return projection == Projection::perspective ? "perspective" : "orthographic";
}

} // namespace ascii_protein
