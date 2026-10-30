#pragma once
#include "ascii_protein/structure.hpp"
#include <filesystem>
#include <map>
namespace ascii_protein {
using MetadataRow = std::map<std::string, std::string>;
void readMolecularRow(Structure& structure, const std::string& category, const MetadataRow& row);
std::vector<std::pair<std::size_t, std::size_t>> depositedBonds(const Structure& structure);
void selectModel(Structure& structure, int model);
void expandAssembly(Structure& structure, const std::string& id);
void exportStructure(const Structure& structure, const std::filesystem::path& path,
                     const std::vector<bool>& mask);
} // namespace ascii_protein
