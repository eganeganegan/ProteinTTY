#include "ascii_protein/structure_loader.hpp"

#include "ascii_protein/binarycif_parser.hpp"
#include "ascii_protein/mmcif_parser.hpp"
#include "ascii_protein/pdb_parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ascii_protein {

static Structure readStructure(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    if (extension == ".bcif") return BinaryCifParser::parseFile(path);
    if (extension == ".cif" || extension == ".mmcif") return MmcifParser::parseFile(path);
    if (extension == ".pdb" || extension == ".ent") return PdbParser::parseFile(path);

    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open structure file: " + path.string());
    const int firstByte = input.peek();
    if (firstByte >= 0 && static_cast<unsigned char>(firstByte) >= 0x80) {
        return BinaryCifParser::parse(input, path.filename().string());
    }
    std::string firstToken;
    input >> firstToken;
    input.clear();
    input.seekg(0);
    if (firstToken.starts_with("data_") || firstToken == "loop_" || firstToken.starts_with("_")) {
        return MmcifParser::parse(input, path.filename().string());
    }
    return PdbParser::parse(input, path.filename().string());
}

Structure loadStructure(const std::filesystem::path& path) {
    auto structure = readStructure(path);
    structure.sourcePath = std::filesystem::absolute(path);
    return structure;
}

} // namespace ascii_protein