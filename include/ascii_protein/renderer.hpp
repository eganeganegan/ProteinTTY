#pragma once

#include "ascii_protein/camera.hpp"
#include "ascii_protein/framebuffer.hpp"
#include "ascii_protein/selection.hpp"
#include "ascii_protein/structure.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace ascii_protein {

enum class RenderMode { atoms, backbone, alphaTrace, ballAndStick, residues };
enum class ColorMode { element, chain, residue, bFactor, secondary };
enum class Charset { ascii, unicode };
enum class RendererBackend { cells, braille };

struct RenderSettings {
    RenderMode mode{RenderMode::ballAndStick};
    ColorMode colorMode{ColorMode::element};
    Charset charset{Charset::unicode};
    RendererBackend backend{RendererBackend::cells};
    bool color{true};
    bool showInfo{true};
    bool showSequence{false};
    bool showHelp{false};
    int backgroundColor{};
    std::optional<std::size_t> selectedResidue;
    std::optional<std::size_t> inspectedAtom;
    AtomMask visibleAtoms;
    AtomMask selectedAtoms;
    std::vector<int> atomColors;
    NamedSelections namedSelections;
    bool commandActive{};
    std::string commandLine;
    std::string commandMessage;
    std::string trajectoryInfo;
};

struct ProjectedPoint {
    double x{};
    double y{};
    double depth{};
    bool visible{};
};

struct Bond {
    std::size_t first{};
    std::size_t second{};
};

std::vector<Bond> inferBonds(const Structure& structure);

class Renderer {
public:
    explicit Renderer(const Structure& structure);

    const Framebuffer& render(const Camera& camera, const RenderSettings& settings, int width, int height);
    ProjectedPoint project(const Vec3& position, const Camera& camera, int width, int height) const;
    std::optional<std::size_t> nearestAtom(const Camera& camera, const RenderSettings& settings, int width,
                                           int height) const;

private:
    struct BrailleCell {
        std::array<double, 8> depths;
        std::array<int, 8> colors{};
    };

    struct RenderPoint {
        ProjectedPoint projected;
        Vec3 rotated;
    };

    RenderPoint transform(const Vec3& position, const Camera& camera, int width, int height,
                          RendererBackend backend) const;
    void drawAtom(std::size_t atomIndex, const Camera& camera, const RenderSettings& settings, int viewWidth,
                  int height);
    void drawLine(std::size_t first, std::size_t second, const Camera& camera, const RenderSettings& settings,
                  int viewWidth, int height);
    int atomColor(std::size_t atomIndex, const RenderSettings& settings) const;
    void putPixel(int x, int y, double depth, char32_t glyph, int color);
    void prepareBraille(int width, int height);
    void composeBraille(int viewWidth, int height);
    void drawPanels(const Camera& camera, const RenderSettings& settings, int viewWidth);
    void drawHelp();
    void drawCommandStatus(const RenderSettings& settings);

    const Structure& structure_;
    Framebuffer framebuffer_{1, 1};
    std::vector<Bond> bonds_;
    std::vector<BrailleCell> brailleCells_;
    int brailleWidth_{};
    int brailleHeight_{};
    bool renderingBraille_{};
    double minBFactor_{};
    double maxBFactor_{};
};

std::string renderModeName(RenderMode mode);
std::string colorModeName(ColorMode mode);
std::string projectionName(Projection projection);

} // namespace ascii_protein