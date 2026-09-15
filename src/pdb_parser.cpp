#include "ascii_protein/pdb_parser.hpp"

#include <charconv>
#include <fstream>
#include <optional>
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
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string record = trim(field(line, 0, 6));
        if (record == "HELIX") {
            const auto start = parseInteger(field(line, 21, 4));
            const auto end = parseInteger(field(line, 33, 4));
            if (start && end) {
                structure.secondaryRanges.push_back(
                    {SecondaryStructure::helix, character(line, 19), *start, character(line, 31), *end});
            }
            continue;
        }
        if (record == "SHEET") {
            const auto start = parseInteger(field(line, 22, 4));
            const auto end = parseInteger(field(line, 33, 4));
            if (start && end) {
                structure.secondaryRanges.push_back(
                    {SecondaryStructure::sheet, character(line, 21), *start, character(line, 32), *end});
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
        atom.serial = *serial;
        atom.name = trim(field(line, 12, 4));
        atom.alternateLocation = character(line, 16);
        atom.residueName = trim(field(line, 17, 3));
        atom.chain = character(line, 21);
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