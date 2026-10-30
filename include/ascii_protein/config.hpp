#pragma once
#include <filesystem>
#include <map>
#include <string>
namespace ascii_protein {
struct CommandConfig {
    std::map<std::string, std::string> macros;
    std::map<char32_t, std::string> bindings;
    int depth{};
};
CommandConfig loadCommandConfig(const std::filesystem::path& path);
} // namespace ascii_protein
