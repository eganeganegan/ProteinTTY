#pragma once

#include "ascii_protein/renderer.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace ascii_protein {

struct CliOptions {
    std::optional<std::filesystem::path> input;
    std::optional<std::string> pdbId;
    bool fetchCif{};
    RenderSettings render;
    Projection projection{Projection::perspective};
    bool spin{};
    bool snapshot{};
    bool help{};
    bool version{};
    int width{100};
    int height{32};
    std::vector<std::string> commands;
    std::vector<std::filesystem::path> scripts;
    std::optional<std::filesystem::path> trajectory;
    int frameDelayMs{100};
    bool trajectoryLoop{};
};

CliOptions parseCli(int argc, char** argv);
std::string helpText();

} // namespace ascii_protein