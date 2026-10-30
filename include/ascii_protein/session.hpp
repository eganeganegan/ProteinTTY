#pragma once
#include "ascii_protein/command.hpp"
namespace ascii_protein {
void saveSession(const std::filesystem::path& path, CommandState state);
void loadSession(const std::filesystem::path& path, CommandState state);
} // namespace ascii_protein
