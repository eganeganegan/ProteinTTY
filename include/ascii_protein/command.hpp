#pragma once

#include "ascii_protein/camera.hpp"
#include "ascii_protein/renderer.hpp"
#include "ascii_protein/structure.hpp"

#include <string>
#include <string_view>
#include <filesystem>
#include <functional>
#include <vector>
#include "ascii_protein/config.hpp"

namespace ascii_protein {

struct TrajectoryState;

struct CommandState {
    Camera& camera;
    RenderSettings& render;
    Structure& structure;
    bool& autorotate;
    bool& running;
    std::function<void()> structureChanged;
    std::function<void(const std::filesystem::path&)> saveFrame;
    TrajectoryState* trajectory{};
    CommandConfig* config{};
};

struct CommandResult {
    bool success{true};
    std::vector<std::string> output;
};

CommandResult executeCommands(std::string_view source, CommandState state);
std::string commandHelp(std::string_view command = {});
std::vector<std::string> commandCompletions(std::string_view input, CommandState state);

} // namespace ascii_protein