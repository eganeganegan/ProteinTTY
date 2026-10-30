#include "ascii_protein/pdb_parser.hpp"

#include <charconv>
#include <fstream>
#include <optional>
#include <sstream>
#include "ascii_protein/molecular.hpp"
#include <stdexcept>

namespace ascii_protein {
namespace {

std::string trim(std::string_view text) {
    const auto first = text.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return {};
    }
    return std::string(text.substr(first, text.find_last_not_of(' ') - first + 1));
}

std::string_view field(const std::string& line, std::size_t offset, std::size_t count) {
    if (offset >= line.size()) {
        return {};
    }
    return std::string_view(line).substr(offset, std::min(count, line.size() - offset));
}

char character(const std::string& line, std::size_t offset) {
    return offset < line.size() ? line[offset] : ' ';
}

std::string chainId(const std::string& line, std::size_t offset) {
    const char value = character(line, offset);
    return value == ' ' ? std::string{} : std::string(1, value);
}

std::optional<int> parseInteger(std::string_view text) {
    const std::string value = trim(text);
    if (value.empty()) {
        return std::nullopt;
    }
    int result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    return error == std::errc{} && end == value.data() + value.size() ? std::optional{result} : std::nullopt;
}

std::optional<double> parseDouble(std::string_view text) {
    const std::string value = trim(text);
    if (value.empty()) {
        return std::nullopt;
    }
    double result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    return error == std::errc{} && end == value.data() + value.size() ? std::optional{result} : std::nullopt;
}

std::string inferElement(const std::string& atomName) {
    for (const unsigned char value : atomName) {
        if (std::isalpha(value)) {
            return std::string(1, static_cast<char>(std::toupper(value)));
        }
    }
    return "?";
}

} // namespace

Structure PdbParser::parseFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open PDB file: " + path.string());
    }
    return parse(input, path.filename().string());
}

Structure PdbParser::parse(std::istream& input, std::string sourceName) {
    Structure structure;
    structure.sourceName = std::move(sourceName);
    std::string line;
    std::size_t lineNumber = 0;
    int model = 1;
    std::string assemblyId;
    std::vector<std::string> assemblyChains;
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string record = trim(field(line, 0, 6));
        if (record == "MODEL") {
            model = parseInteger(field(line, 10, 4)).value_or(1);
            continue;
        }
        if (record == "CONECT") {
            const auto first = parseInteger(field(line, 6, 5));
            if (first)
                for (std::size_t offset = 11; offset < 31 && offset < line.size(); offset += 5)
                    if (auto second = parseInteger(field(line, offset, 5)))
                        structure.explicitBonds.emplace_back(*first, *second);
            continue;
        }
        if (record == "REMARK" && field(line, 7, 3) == "350") {
            auto pos = line.find("BIOMOLECULE:");
            if (pos != std::string::npos) {
                assemblyId = trim(std::string_view(line).substr(pos + 12));
                assemblyChains.clear();
            }
            pos = line.find("CHAINS:");
            if (pos != std::string::npos) {
                if (line.find("APPLY") != std::string::npos) assemblyChains.clear();
                std::string names = line.substr(pos + 7);
                for (auto& c : names)
                    if (c == ',') c = ' ';
                std::istringstream ids(names);
                std::string id;
                while (ids >> id)
                    assemblyChains.push_back(id == "NULL" ? "" : id);
            }
            pos = line.find("BIOMT");
            if (pos != std::string::npos && pos + 5 < line.size()) {
                int row = line[pos + 5] - '1';
                std::istringstream values(line.substr(pos + 6));
                std::string id;
                double x, y, z, t;
                if (row >= 0 && row < 3 && values >> id >> x >> y >> z >> t) {
                    auto groupKey = assemblyId;
                    for (auto& c : groupKey)
                        if (c == ',' || c == ' ') c = '_';
                    const auto key = groupKey + ":" + id;
                    auto& op = structure.operators[key];
                    op.rotation.m[row][0] = x;
                    op.rotation.m[row][1] = y;
                    op.rotation.m[row][2] = z;
                    if (row == 0) op.translation.x = t;
                    if (row == 1) op.translation.y = t;
                    if (row == 2) op.translation.z = t;
                    if (row == 2) {
                        auto names = assemblyId;
                        for (auto& c : names)
                            if (c == ',') c = ' ';
                        std::istringstream ids(names);
                        std::string assembly;
                        while (ids >> assembly)
                            structure.assemblies.push_back({assembly, key, assemblyChains});
                    }
                }
            }
        }
        if (record == "HELIX") {
            const auto start = parseInteger(field(line, 21, 4));
            const auto end = parseInteger(field(line, 33, 4));
            if (start && end) {
                structure.secondaryRanges.push_back(
                    {SecondaryStructure::helix, chainId(line, 19), *start, chainId(line, 31), *end,
                     character(line, 25), character(line, 37)});
            }
            continue;
        }
        if (record == "SHEET") {
            const auto start = parseInteger(field(line, 22, 4));
            const auto end = parseInteger(field(line, 33, 4));
            if (start && end) {
                structure.secondaryRanges.push_back(
                    {SecondaryStructure::sheet, chainId(line, 21), *start, chainId(line, 32), *end,
                     character(line, 26), character(line, 37)});
            }
            continue;
        }
        if (record != "ATOM" && record != "HETATM") {
            continue;
        }

        const auto serial = parseInteger(field(line, 6, 5));
        const auto residueNumber = parseInteger(field(line, 22, 4));
        const auto x = parseDouble(field(line, 30, 8));
        const auto y = parseDouble(field(line, 38, 8));
        const auto z = parseDouble(field(line, 46, 8));
        if (!serial || !residueNumber || !x || !y || !z) {
            structure.warnings.push_back("line " + std::to_string(lineNumber) + ": malformed atom record skipped");
            continue;
        }

        Atom atom;
        atom.model = model;
        atom.serial = *serial;
        atom.name = trim(field(line, 12, 4));
        atom.alternateLocation = character(line, 16);
        atom.residueName = trim(field(line, 17, 3));
        atom.chain = chainId(line, 21);
        atom.residueNumber = *residueNumber;
        atom.insertionCode = character(line, 26);
        atom.position = {*x, *y, *z};
        atom.occupancy = parseDouble(field(line, 54, 6)).value_or(0.0);
        atom.bFactor = parseDouble(field(line, 60, 6)).value_or(0.0);
        atom.element = trim(field(line, 76, 2));
        if (atom.element.empty()) {
            atom.element = inferElement(atom.name);
        }
        atom.hetero = record == "HETATM";
        structure.atoms.push_back(std::move(atom));
    }

    if (structure.atoms.empty()) {
        throw std::runtime_error("no valid ATOM or HETATM records in " + structure.sourceName);
    }
    structure.finalize();
    return structure;
}

} // namespace ascii_protein