#pragma once

#include "ascii_protein/camera.hpp"
#include "ascii_protein/renderer.hpp"
#include "ascii_protein/structure.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace ascii_protein {

struct CommandState {
    Camera& camera;
    RenderSettings& render;
    const Structure& structure;
    bool& autorotate;
    bool& running;
};

struct CommandResult {
    bool success{true};
    std::vector<std::string> output;
};

CommandResult executeCommands(std::string_view source, CommandState state);
std::string commandHelp(std::string_view command = {});

} // namespace ascii_protein