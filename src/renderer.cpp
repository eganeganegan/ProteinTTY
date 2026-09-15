#include "ascii_protein/renderer.hpp"

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
    constexpr std::array<char32_t, 6> ramp{U'·', U'•', U'░', U'▒', U'▓', U'█'};
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

} // namespace

std::vector<Bond> inferBonds(const Structure& structure) {
    constexpr double cellSize = 2.7;
    std::vector<Bond> bonds;
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
    const Vec3 rotated = camera.rotation() * position;
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

int Renderer::atomColor(const Atom& atom, const RenderSettings& settings) const {
    if (settings.inspectedAtom && *settings.inspectedAtom < structure_.atoms.size() &&
        &atom == &structure_.atoms[*settings.inspectedAtom]) {
        return 96;
    }
    if (atom.residueIndex == settings.selectedResidue) {
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
    const auto& atom = structure_.atoms[atomIndex];
    const auto point = transform(atom.position, camera, viewWidth, height, settings.backend);
    if (!point.projected.visible) return;

    const double normalizedDepth = 1.0 - std::clamp(point.projected.depth / (structure_.boundingRadius * 6.4), 0.0, 1.0);
    const double resolutionScale = settings.backend == RendererBackend::braille ? 2.0 : 1.0;
    const int radius = std::clamp(
        static_cast<int>(std::round(displayRadius(atom.element) * camera.zoom * resolutionScale)), 0, 6);
    const int color = atomColor(atom, settings);
    for (int offsetY = -radius; offsetY <= radius; ++offsetY) {
        for (int offsetX = -radius * 2; offsetX <= radius * 2; ++offsetX) {
            const double distance = std::sqrt(std::pow(offsetX / 2.0, 2.0) + std::pow(offsetY, 2.0));
            if (distance > radius + 0.25) continue;
            const double surface = radius == 0 ? 1.0 : std::sqrt(std::max(0.0, 1.0 - std::pow(distance / (radius + 0.25), 2.0)));
            const double lighting = std::clamp(0.2 + normalizedDepth * 0.45 + surface * 0.35 - offsetY * 0.06, 0.0, 1.0);
            putPixel(static_cast<int>(std::round(point.projected.x)) + offsetX,
                     static_cast<int>(std::round(point.projected.y)) + offsetY,
                     point.projected.depth - surface * 0.2, shadeGlyph(settings.charset, lighting), color);
        }
    }
}

void Renderer::drawLine(std::size_t first, std::size_t second, const Camera& camera, const RenderSettings& settings,
                        int viewWidth, int height) {
    const auto start = transform(structure_.atoms[first].position, camera, viewWidth, height, settings.backend).projected;
    const auto end = transform(structure_.atoms[second].position, camera, viewWidth, height, settings.backend).projected;
    if (!start.visible || !end.visible) return;
    const double deltaX = end.x - start.x;
    const double deltaY = end.y - start.y;
    const int steps = std::max(1, static_cast<int>(std::ceil(std::max(std::abs(deltaX), std::abs(deltaY)))));
    const int color = atomColor(structure_.atoms[first], settings);
    for (int step = 0; step <= steps; ++step) {
        const double amount = static_cast<double>(step) / steps;
        const int x = static_cast<int>(std::round(start.x + deltaX * amount));
        const int y = static_cast<int>(std::round(start.y + deltaY * amount));
        const double depth = start.depth + (end.depth - start.depth) * amount + 0.05;
        putPixel(x, y, depth, settings.charset == Charset::ascii ? U'.' : U'·', color);
    }
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
        } else if (!structure_.residues.empty()) {
            const auto& residue = structure_.residues[std::min(settings.selectedResidue, structure_.residues.size() - 1)];
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

    if (settings.showSequence && !structure_.residues.empty()) {
        const int y = framebuffer_.height() - 4;
        const std::size_t selected = std::min(settings.selectedResidue, structure_.residues.size() - 1);
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
        "I               info panel", "S               sequence", ":               command console",
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
        std::optional<std::size_t> previous;
        for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
            if (structure_.atoms[index].name != "CA") continue;
            if (previous && residuesAreAdjacent(structure_.atoms[*previous], structure_.atoms[index])) {
                drawLine(*previous, index, camera, settings, rasterWidth, rasterHeight);
            }
            previous = index;
        }
    }

    for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
        const auto& atom = structure_.atoms[index];
        if (atomShown(atom, settings.mode)) drawAtom(index, camera, settings, rasterWidth, rasterHeight);
    }
    if (renderingBraille_) composeBraille(viewWidth, height);
    renderingBraille_ = false;
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