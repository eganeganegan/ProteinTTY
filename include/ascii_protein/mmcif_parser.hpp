#pragma once

#include "ascii_protein/structure.hpp"

#include <filesystem>
#include <istream>
#include <string>

namespace ascii_protein {

class MmcifParser {
public:
    static Structure parseFile(const std::filesystem::path& path);
    static Structure parse(std::istream& input, std::string sourceName = "<stream>");
};

} // namespace ascii_protein