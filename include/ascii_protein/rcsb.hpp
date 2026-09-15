#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

namespace ascii_protein {

std::string normalizePdbId(std::string_view pdbId);
std::filesystem::path fetchPdb(const std::string& pdbId,
                               std::optional<std::filesystem::path> cacheDirectory = std::nullopt);
std::filesystem::path fetchMmcif(const std::string& pdbId,
                                 std::optional<std::filesystem::path> cacheDirectory = std::nullopt);

} // namespace ascii_protein