#include "ascii_protein/mmcif_parser.hpp"

#include <charconv>
#include <cctype>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ascii_protein {
namespace {

std::vector<std::string> tokenize(std::istream& input) {
    std::vector<std::string> tokens;
    std::string line;
    bool multiline = false;
    std::string block;
    while (std::getline(input, line)) {
        if (multiline) {
            if (!line.empty() && line.front() == ';') {
                tokens.push_back(std::move(block));
                block.clear();
                multiline = false;
            } else {
                if (!block.empty()) block.push_back('\n');
                block += line;
            }
            continue;
        }
        if (!line.empty() && line.front() == ';') {
            multiline = true;
            block = line.substr(1);
            continue;
        }

        std::size_t position = 0;
        while (position < line.size()) {
            while (position < line.size() && std::isspace(static_cast<unsigned char>(line[position]))) ++position;
            if (position >= line.size() || line[position] == '#') break;
            if (line[position] == '\'' || line[position] == '"') {
                const char quote = line[position++];
                const std::size_t start = position;
                while (position < line.size() &&
                       (line[position] != quote ||
                        (position + 1 < line.size() &&
                         !std::isspace(static_cast<unsigned char>(line[position + 1]))))) {
                    ++position;
                }
                if (position >= line.size()) throw std::runtime_error("unterminated quoted mmCIF value");
                tokens.push_back(line.substr(start, position - start));
                ++position;
            } else {
                const std::size_t start = position;
                while (position < line.size() && !std::isspace(static_cast<unsigned char>(line[position])) &&
                       line[position] != '#') {
                    ++position;
                }
                tokens.push_back(line.substr(start, position - start));
            }
        }
    }
    if (multiline) throw std::runtime_error("unterminated multiline mmCIF value");
    return tokens;
}

bool controlToken(const std::string& token) {
    return token == "loop_" || token == "stop_" || token.starts_with("data_") || token.starts_with("save_") ||
           token.starts_with("_");
}

std::optional<int> integer(std::string_view value) {
    if (value.empty() || value == "." || value == "?") return std::nullopt;
    int result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    return error == std::errc{} && end == value.data() + value.size() ? std::optional{result} : std::nullopt;
}

std::optional<double> number(std::string_view value) {
    if (value.empty() || value == "." || value == "?") return std::nullopt;
    double result{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    return error == std::errc{} && end == value.data() + value.size() ? std::optional{result} : std::nullopt;
}

std::string value(const std::vector<std::string>& row, const std::unordered_map<std::string, std::size_t>& columns,
                  std::initializer_list<std::string_view> names) {
    for (const auto name : names) {
        const auto found = columns.find(std::string(name));
        if (found != columns.end() && found->second < row.size()) {
            const auto& result = row[found->second];
            if (result != "." && result != "?") return result;
        }
    }
    return {};
}

char marker(const std::string& text) { return text.empty() ? ' ' : text.front(); }

void parseAtomRow(const std::vector<std::string>& row,
                  const std::unordered_map<std::string, std::size_t>& columns, Structure& structure,
                  std::size_t rowNumber) {
    const std::string record = value(row, columns, {"_atom_site.group_PDB"});
    if (record != "ATOM" && record != "HETATM") return;
    const auto serial = integer(value(row, columns, {"_atom_site.id"}));
    const auto residueNumber = integer(value(row, columns, {"_atom_site.auth_seq_id", "_atom_site.label_seq_id"}));
    const auto x = number(value(row, columns, {"_atom_site.Cartn_x"}));
    const auto y = number(value(row, columns, {"_atom_site.Cartn_y"}));
    const auto z = number(value(row, columns, {"_atom_site.Cartn_z"}));
    if (!serial || !residueNumber || !x || !y || !z) {
        structure.warnings.push_back("atom_site row " + std::to_string(rowNumber) + ": malformed atom skipped");
        return;
    }
    const auto model = integer(value(row, columns, {"_atom_site.pdbx_PDB_model_num"}));
    if (model && *model != 1) return;

    Atom atom;
    atom.serial = *serial;
    atom.name = value(row, columns, {"_atom_site.auth_atom_id", "_atom_site.label_atom_id"});
    atom.alternateLocation = marker(value(row, columns, {"_atom_site.label_alt_id"}));
    atom.element = value(row, columns, {"_atom_site.type_symbol"});
    atom.residueName = value(row, columns, {"_atom_site.auth_comp_id", "_atom_site.label_comp_id"});
    const std::string chain = value(row, columns, {"_atom_site.auth_asym_id", "_atom_site.label_asym_id"});
    atom.chain = chain;
    atom.residueNumber = *residueNumber;
    atom.insertionCode = marker(value(row, columns, {"_atom_site.pdbx_PDB_ins_code"}));
    atom.position = {*x, *y, *z};
    atom.occupancy = number(value(row, columns, {"_atom_site.occupancy"})).value_or(0.0);
    atom.bFactor = number(value(row, columns, {"_atom_site.B_iso_or_equiv"})).value_or(0.0);
    atom.hetero = record == "HETATM";
    structure.atoms.push_back(std::move(atom));
}

void parseSecondaryRow(const std::vector<std::string>& row,
                       const std::unordered_map<std::string, std::size_t>& columns, SecondaryStructure type,
                       Structure& structure) {
    const std::string prefix = type == SecondaryStructure::helix ? "_struct_conf." : "_struct_sheet_range.";
    if (type == SecondaryStructure::helix) {
        const std::string kind = value(row, columns, {"_struct_conf.conf_type_id"});
        if (!kind.starts_with("HELX")) return;
    }
    const std::string startChain = value(
        row, columns, {prefix + "beg_auth_asym_id", prefix + "beg_label_asym_id"});
    const std::string endChain = value(
        row, columns, {prefix + "end_auth_asym_id", prefix + "end_label_asym_id"});
    const auto startResidue = integer(value(
        row, columns, {prefix + "beg_auth_seq_id", prefix + "beg_label_seq_id"}));
    const auto endResidue = integer(value(
        row, columns, {prefix + "end_auth_seq_id", prefix + "end_label_seq_id"}));
    if (!startResidue || !endResidue || startChain.empty() || endChain.empty()) return;
    structure.secondaryRanges.push_back(
        {type, startChain, *startResidue, endChain, *endResidue,
         marker(value(row, columns, {prefix + "pdbx_beg_PDB_ins_code"})),
         marker(value(row, columns, {prefix + "pdbx_end_PDB_ins_code"}))});
}

} // namespace

Structure MmcifParser::parseFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open mmCIF file: " + path.string());
    return parse(input, path.filename().string());
}

Structure MmcifParser::parse(std::istream& input, std::string sourceName) {
    Structure structure;
    structure.sourceName = std::move(sourceName);
    const auto tokens = tokenize(input);
    std::size_t position = 0;
    std::size_t atomRow = 0;
    while (position < tokens.size()) {
        if (tokens[position++] != "loop_") continue;
        std::vector<std::string> tags;
        while (position < tokens.size() && tokens[position].starts_with("_")) tags.push_back(tokens[position++]);
        if (tags.empty()) continue;
        const bool atomLoop = tags.front().starts_with("_atom_site.");
        const bool helixLoop = tags.front().starts_with("_struct_conf.");
        const bool sheetLoop = tags.front().starts_with("_struct_sheet_range.");
        std::unordered_map<std::string, std::size_t> columns;
        for (std::size_t index = 0; index < tags.size(); ++index) columns.emplace(tags[index], index);
        while (position < tokens.size() && !controlToken(tokens[position])) {
            if (position + tags.size() > tokens.size()) {
                structure.warnings.push_back("incomplete mmCIF loop row skipped");
                position = tokens.size();
                break;
            }
            std::vector<std::string> row(tokens.begin() + static_cast<std::ptrdiff_t>(position),
                                         tokens.begin() + static_cast<std::ptrdiff_t>(position + tags.size()));
            position += tags.size();
            if (atomLoop) parseAtomRow(row, columns, structure, ++atomRow);
            else if (helixLoop) parseSecondaryRow(row, columns, SecondaryStructure::helix, structure);
            else if (sheetLoop) parseSecondaryRow(row, columns, SecondaryStructure::sheet, structure);
        }
    }
    if (structure.atoms.empty()) throw std::runtime_error("no valid atom_site records in " + structure.sourceName);
    structure.finalize();
    return structure;
}

} // namespace ascii_protein