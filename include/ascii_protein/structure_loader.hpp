#pragma once

#include "ascii_protein/structure.hpp"

#include <filesystem>

namespace ascii_protein {

Structure loadStructure(const std::filesystem::path& path);

} // namespace ascii_protein