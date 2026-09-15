#pragma once

#include "ascii_protein/renderer.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace ascii_protein {

struct CliOptions {
    std::optional<std::filesystem::path> input;
    std::optional<std::string> pdbId;
    RenderSettings render;
    Projection projection{Projection::perspective};
    bool spin{};
    bool snapshot{};
    bool help{};
    bool version{};
    int width{100};
    int height{32};
};

CliOptions parseCli(int argc, char** argv);
std::string helpText();

} // namespace ascii_protein