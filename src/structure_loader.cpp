#include "ascii_protein/structure_loader.hpp"

#include "ascii_protein/mmcif_parser.hpp"
#include "ascii_protein/pdb_parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ascii_protein {

Structure loadStructure(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    if (extension == ".cif" || extension == ".mmcif") return MmcifParser::parseFile(path);
    if (extension == ".pdb" || extension == ".ent") return PdbParser::parseFile(path);

    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open structure file: " + path.string());
    std::string firstToken;
    input >> firstToken;
    input.clear();
    input.seekg(0);
    if (firstToken.starts_with("data_") || firstToken == "loop_" || firstToken.starts_with("_")) {
        return MmcifParser::parse(input, path.filename().string());
    }
    return PdbParser::parse(input, path.filename().string());
}

} // namespace ascii_protein