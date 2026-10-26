#include "ascii_protein/command.hpp"
#include "ascii_protein/binarycif_parser.hpp"
#include "ascii_protein/cli.hpp"
#include "ascii_protein/dcd_trajectory_source.hpp"
#include "ascii_protein/geometry.hpp"
#include "ascii_protein/mmcif_parser.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/renderer.hpp"
#include "ascii_protein/selection.hpp"
#include "ascii_protein/structure.hpp"
#include "ascii_protein/structure_loader.hpp"
#include "ascii_protein/terminal_input.hpp"
#include "ascii_protein/trajectory.hpp"
#include "ascii_protein/xyz_trajectory_source.hpp"

#include <cmath>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool near(double left, double right) { return std::abs(left - right) < 1e-9; }

class MessagePackWriter {
public:
    void map(std::size_t size) { data_.push_back(static_cast<char>(0x80 | size)); }
    void array(std::size_t size) {
        if (size < 16) data_.push_back(static_cast<char>(0x90 | size));
        else {
            data_.push_back(static_cast<char>(0xdc));
            bigEndian(static_cast<std::uint16_t>(size));
        }
    }
    void string(std::string_view value) {
        require(value.size() < 32, "test MessagePack string length");
        data_.push_back(static_cast<char>(0xa0 | value.size()));
        data_.append(value);
    }
    void integer(int value) {
        if (value >= 0 && value < 128) data_.push_back(static_cast<char>(value));
        else {
            data_.push_back(static_cast<char>(0xd2));
            bigEndian(static_cast<std::uint32_t>(value));
        }
    }
    void boolean(bool value) { data_.push_back(static_cast<char>(value ? 0xc3 : 0xc2)); }
    void nil() { data_.push_back(static_cast<char>(0xc0)); }
    void binary(std::string_view value) {
        require(value.size() <= 0xffff, "test MessagePack binary length");
        data_.push_back(static_cast<char>(0xc5));
        bigEndian(static_cast<std::uint16_t>(value.size()));
        data_.append(value);
    }
    std::string finish() && { return std::move(data_); }

private:
    template <typename T>
    void bigEndian(T value) {
        for (std::size_t shift = sizeof(T); shift-- > 0;) {
            data_.push_back(static_cast<char>((value >> (shift * 8)) & 0xff));
        }
    }
    std::string data_;
};

std::string littleEndianIntegers(std::initializer_list<int> values) {
    std::string result;
    for (const int value : values) {
        const auto bits = static_cast<std::uint32_t>(value);
        for (int shift = 0; shift < 32; shift += 8) result.push_back(static_cast<char>((bits >> shift) & 0xff));
    }
    return result;
}

std::string littleEndianDoubles(std::initializer_list<double> values) {
    std::string result;
    for (const double value : values) {
        const auto bits = std::bit_cast<std::uint64_t>(value);
        for (int shift = 0; shift < 64; shift += 8) result.push_back(static_cast<char>((bits >> shift) & 0xff));
    }
    return result;
}

void appendDcdU32(std::string& output, std::uint32_t value, bool bigEndian) {
    if (bigEndian) {
        for (int shift = 24; shift >= 0; shift -= 8) output.push_back(static_cast<char>((value >> shift) & 0xff));
    } else {
        for (int shift = 0; shift <= 24; shift += 8) output.push_back(static_cast<char>((value >> shift) & 0xff));
    }
}

void setDcdU32(std::string& output, std::size_t offset, std::uint32_t value, bool bigEndian) {
    std::string encoded;
    appendDcdU32(encoded, value, bigEndian);
    output.replace(offset, encoded.size(), encoded);
}

void appendDcdRecord(std::string& output, const std::string& payload, bool bigEndian) {
    appendDcdU32(output, static_cast<std::uint32_t>(payload.size()), bigEndian);
    output += payload;
    appendDcdU32(output, static_cast<std::uint32_t>(payload.size()), bigEndian);
}

std::string dcdFixture(bool bigEndian = false, bool unitCell = false) {
    std::string output;
    std::string header(84, '\0');
    header.replace(0, 4, "CORD");
    setDcdU32(header, 4, 2, bigEndian);
    setDcdU32(header, 4 + 10 * 4, unitCell ? 1 : 0, bigEndian);
    setDcdU32(header, 4 + 19 * 4, 24, bigEndian);
    appendDcdRecord(output, header, bigEndian);
    std::string title;
    appendDcdU32(title, 1, bigEndian);
    title.append("ProteinTTY DCD test fixture", 26);
    title.resize(84, ' ');
    appendDcdRecord(output, title, bigEndian);
    std::string atomCount;
    appendDcdU32(atomCount, 2, bigEndian);
    appendDcdRecord(output, atomCount, bigEndian);

    const auto coordinates = [&](float first, float second) {
        std::string payload;
        const auto appendFloat = [&](float value) {
            appendDcdU32(payload, std::bit_cast<std::uint32_t>(value), bigEndian);
        };
        appendFloat(first);
        appendFloat(second);
        appendDcdRecord(output, payload, bigEndian);
    };
    for (int frame = 0; frame < 2; ++frame) {
        if (unitCell) appendDcdRecord(output, std::string(48, '\0'), bigEndian);
        coordinates(0.0F, frame == 0 ? 2.0F : 4.0F);
        coordinates(0.0F, 0.0F);
        coordinates(0.0F, 0.0F);
    }
    return output;
}

std::string fixedAtomDcdFixture() {
    std::string output;
    std::string header(84, '\0');
    header.replace(0, 4, "CORD");
    setDcdU32(header, 4, 2, false);
    setDcdU32(header, 4 + 8 * 4, 1, false);
    setDcdU32(header, 4 + 19 * 4, 24, false);
    appendDcdRecord(output, header, false);
    std::string title;
    appendDcdU32(title, 1, false);
    title.resize(84, ' ');
    appendDcdRecord(output, title, false);
    std::string atomCount;
    appendDcdU32(atomCount, 3, false);
    appendDcdRecord(output, atomCount, false);
    std::string freeAtoms;
    appendDcdU32(freeAtoms, 1, false);
    appendDcdU32(freeAtoms, 3, false);
    appendDcdRecord(output, freeAtoms, false);
    const auto coordinates = [&](std::initializer_list<float> values) {
        std::string payload;
        for (const float value : values) appendDcdU32(payload, std::bit_cast<std::uint32_t>(value), false);
        appendDcdRecord(output, payload, false);
    };
    coordinates({0.0F, 10.0F, 20.0F});
    coordinates({0.0F, 0.0F, 0.0F});
    coordinates({0.0F, 0.0F, 0.0F});
    coordinates({1.0F, 21.0F});
    coordinates({2.0F, 22.0F});
    coordinates({3.0F, 23.0F});
    return output;
}

void writeByteEncoding(MessagePackWriter& writer, int type) {
    writer.map(2); writer.string("kind"); writer.string("ByteArray"); writer.string("type"); writer.integer(type);
}

void writeEncodings(MessagePackWriter& writer, int type, std::size_t size, bool packed) {
    writer.array(packed ? 2 : 1);
    if (packed) {
        writer.map(4);
        writer.string("kind"); writer.string("IntegerPacking");
        writer.string("byteCount"); writer.integer(1);
        writer.string("srcSize"); writer.integer(static_cast<int>(size));
        writer.string("isUnsigned"); writer.boolean(true);
    }
    writeByteEncoding(writer, type);
}

void writeStringColumn(MessagePackWriter& writer, std::string_view name,
                       std::initializer_list<std::string_view> values,
                       std::initializer_list<int> mask = {}) {
    std::vector<std::string_view> unique;
    std::string indices;
    std::string text;
    std::string offsets(1, '\0');
    for (const auto value : values) {
        auto found = std::find(unique.begin(), unique.end(), value);
        if (found == unique.end()) {
            unique.push_back(value);
            text.append(value);
            offsets.push_back(static_cast<char>(text.size()));
            found = unique.end() - 1;
        }
        indices.push_back(static_cast<char>(found - unique.begin()));
    }
    writer.map(3);
    writer.string("name"); writer.string(name);
    writer.string("data"); writer.map(2);
    writer.string("data"); writer.binary(indices);
    writer.string("encoding"); writer.array(1); writer.map(5);
    writer.string("kind"); writer.string("StringArray");
    writer.string("dataEncoding"); writeEncodings(writer, 4, values.size(), true);
    writer.string("stringData"); writer.string(text);
    writer.string("offsetEncoding"); writeEncodings(writer, 4, offsets.size(), true);
    writer.string("offsets"); writer.binary(offsets);
    writer.string("mask");
    if (mask.size() == 0) writer.nil();
    else {
        std::string bytes;
        for (const int value : mask) bytes.push_back(static_cast<char>(value));
        writer.map(2); writer.string("data"); writer.binary(bytes);
        writer.string("encoding"); writeEncodings(writer, 4, mask.size(), false);
    }
}

void writeNumericColumn(MessagePackWriter& writer, std::string_view name, std::string bytes, int type) {
    writer.map(3);
    writer.string("name"); writer.string(name);
    writer.string("data"); writer.map(2); writer.string("data"); writer.binary(bytes);
    writer.string("encoding"); writeEncodings(writer, type, 0, false);
    writer.string("mask"); writer.nil();
}

std::string binaryCifFixture() {
    MessagePackWriter writer;
    writer.map(3);
    writer.string("version"); writer.string("0.3.0");
    writer.string("encoder"); writer.string("ProteinTTY tests");
    writer.string("dataBlocks"); writer.array(1); writer.map(2);
    writer.string("header"); writer.string("fixture");
    writer.string("categories"); writer.array(3);

    writer.map(3); writer.string("name"); writer.string("_struct_conf");
    writer.string("rowCount"); writer.integer(1); writer.string("columns"); writer.array(5);
    writeStringColumn(writer, "conf_type_id", {"HELX_P"});
    writeStringColumn(writer, "beg_auth_asym_id", {"ALPHA"});
    writeNumericColumn(writer, "beg_auth_seq_id", littleEndianIntegers({7}), 3);
    writeStringColumn(writer, "end_auth_asym_id", {"ALPHA"});
    writeNumericColumn(writer, "end_auth_seq_id", littleEndianIntegers({7}), 3);

    writer.map(3); writer.string("name"); writer.string("_struct_sheet_range");
    writer.string("rowCount"); writer.integer(1); writer.string("columns"); writer.array(6);
    writeStringColumn(writer, "beg_auth_asym_id", {"B"});
    writeNumericColumn(writer, "beg_auth_seq_id", littleEndianIntegers({8}), 3);
    writeStringColumn(writer, "pdbx_beg_PDB_ins_code", {"A"});
    writeStringColumn(writer, "end_auth_asym_id", {"B"});
    writeNumericColumn(writer, "end_auth_seq_id", littleEndianIntegers({8}), 3);
    writeStringColumn(writer, "pdbx_end_PDB_ins_code", {"A"});

    writer.map(3); writer.string("name"); writer.string("_atom_site");
    writer.string("rowCount"); writer.integer(3); writer.string("columns"); writer.array(14);
    writeStringColumn(writer, "group_PDB", {"ATOM", "HETATM", "ATOM"});
    writeNumericColumn(writer, "id", littleEndianIntegers({1, 2, 3}), 3);
    writeStringColumn(writer, "type_symbol", {"C", "O", "N"});
    writeStringColumn(writer, "auth_atom_id", {"CA", "O", "N"});
    writeStringColumn(writer, "auth_comp_id", {"ALA", "HOH", "GLY"});
    writeStringColumn(writer, "auth_asym_id", {"ALPHA", "B", "A"});
    writeNumericColumn(writer, "auth_seq_id", littleEndianIntegers({7, 8, 9}), 3);
    writeStringColumn(writer, "pdbx_PDB_ins_code", {"A", "A", "A"}, {1, 0, 1});
    writeNumericColumn(writer, "Cartn_x", littleEndianDoubles({1, 4, 7}), 33);
    writeNumericColumn(writer, "Cartn_y", littleEndianDoubles({2, 5, 8}), 33);
    writeNumericColumn(writer, "Cartn_z", littleEndianDoubles({3, 6, 9}), 33);
    writeNumericColumn(writer, "occupancy", littleEndianDoubles({1, 0.5, 1}), 33);
    writeNumericColumn(writer, "B_iso_or_equiv", littleEndianDoubles({12.5, 20, 10}), 33);
    writeNumericColumn(writer, "pdbx_PDB_model_num", littleEndianIntegers({1, 1, 2}), 3);
    return std::move(writer).finish();
}

std::string binaryCifEncodedColumnFixture(bool excessiveRunLength) {
    MessagePackWriter writer;
    writer.map(1); writer.string("dataBlocks"); writer.array(1); writer.map(1);
    writer.string("categories"); writer.array(1); writer.map(3);
    writer.string("name"); writer.string("_atom_site");
    writer.string("rowCount"); writer.integer(1);
    writer.string("columns"); writer.array(1); writer.map(3);
    writer.string("name"); writer.string("group_PDB");
    writer.string("data"); writer.map(2);
    writer.string("data");
    writer.binary(excessiveRunLength ? littleEndianIntegers({1, 10'000'001}) : std::string{});
    writer.string("encoding");
    if (excessiveRunLength) {
        writer.array(2);
        writer.map(3); writer.string("kind"); writer.string("RunLength");
        writer.string("srcType"); writer.integer(3); writer.string("srcSize"); writer.integer(1);
        writeByteEncoding(writer, 3);
    } else {
        writer.array(17);
        for (int index = 0; index < 17; ++index) writeByteEncoding(writer, 4);
    }
    writer.string("mask"); writer.nil();
    return std::move(writer).finish();
}

void requireBinaryCifRejected(std::string bytes, std::string_view expected) {
    std::istringstream input(std::move(bytes));
    try {
        static_cast<void>(ascii_protein::BinaryCifParser::parse(input, "malformed.bcif"));
        require(false, "malformed BinaryCIF rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find(expected) != std::string::npos,
                "BinaryCIF rejection reports " + std::string(expected));
    }
}

std::set<std::pair<std::size_t, std::size_t>> bruteForceBonds(const ascii_protein::Structure& structure) {
    std::set<std::pair<std::size_t, std::size_t>> result;
    const auto radius = [](const std::string& element) {
        if (element == "H") return 0.31;
        if (element == "C") return 0.76;
        if (element == "N") return 0.71;
        if (element == "O") return 0.66;
        if (element == "S") return 1.05;
        if (element == "P") return 1.07;
        return 0.77;
    };
    for (std::size_t first = 0; first < structure.atoms.size(); ++first) {
        for (std::size_t second = first + 1; second < structure.atoms.size(); ++second) {
            const auto& left = structure.atoms[first];
            const auto& right = structure.atoms[second];
            if (left.chain != right.chain || std::abs(left.residueNumber - right.residueNumber) > 1) continue;
            const double distance = ascii_protein::length(left.position - right.position);
            if (distance > 0.35 && distance <= radius(left.element) + radius(right.element) + 0.55) {
                result.emplace(first, second);
            }
        }
    }
    return result;
}

void testGeometry() {
    using namespace ascii_protein;
    const Vec3 x{1, 0, 0};
    const Vec3 y{0, 1, 0};
    require(dot(x, y) == 0.0, "orthogonal dot product");
    require(cross(x, y).z == 1.0, "cross product");
    require(near(length(Vec3{3, 4, 0}), 5.0), "vector length");
    const Vec3 rotated = Mat3::rotationZ(radians(90.0)) * x;
    require(near(rotated.x, 0.0) && near(rotated.y, 1.0), "z rotation");

    Camera camera;
    camera.rotate(725.0, -735.0, 540.0);
    require(near(camera.yaw, -20.0), "yaw wraps after full rotations");
    require(near(camera.pitch, 0.0), "pitch wraps after full rotations");
    require(near(camera.roll, -180.0), "roll wraps after full rotations");
}

void testTerminalInputDecoder() {
    using ascii_protein::Key;
    ascii_protein::TerminalInputDecoder decoder;
    decoder.appendBytes("paste");
    for (const char expected : std::string("paste")) {
        const auto event = decoder.next();
        require(event && event->key == Key::character && event->character == expected,
                "pasted terminal input preserves bytes");
    }
    require(!decoder.next(), "terminal input queue drains");

    decoder.appendBytes("\x1b[");
    require(!decoder.next(), "split arrow sequence waits for completion");
    decoder.appendBytes("A");
    require(decoder.next()->key == Key::up, "split arrow sequence decoding");
    decoder.appendBytes("\x1b[B\x1b[C\x1b[D");
    require(decoder.next()->key == Key::down, "down arrow decoding");
    require(decoder.next()->key == Key::right, "right arrow decoding");
    require(decoder.next()->key == Key::left, "left arrow decoding");

    decoder.appendBytes("\x1b");
    require(!decoder.next(), "standalone escape waits for flush");
    require(decoder.next(true)->key == Key::escape, "standalone escape decoding");

    decoder.appendEvent({Key::left}, 3);
    for (int repeat = 0; repeat < 3; ++repeat) require(decoder.next()->key == Key::left, "key repeat decoding");
    decoder.appendEvent({Key::resize});
    require(decoder.next()->key == Key::resize, "resize event decoding");

    decoder.appendCodepoint(U'\u00e9');
    decoder.appendCodepoint(U'\U0001f642');
    std::string unicode;
    while (const auto event = decoder.next()) {
        require(event->key == Key::character, "Unicode input produces character bytes");
        unicode.push_back(event->character);
    }
    require(unicode == "\xc3\xa9\xf0\x9f\x99\x82", "Unicode input is UTF-8 encoded");
}

void testParserAndModel() {
    std::istringstream pdb(
        "HELIX    1   1 ALA A    1  GLY A    2  1                                  2\n"
        "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 12.50           N  \n"
        "ATOM      2  CA AALA A   1      12.000  20.000  30.000  0.50 14.00           C  \n"
        "ATOM      3  CA  GLY A   2      14.000  20.000  30.000                      C  \n"
        "HETATM    4  O   HOH B   9      16.000  20.000  30.000  1.00 20.00           O  \n"
        "ATOM      x  BAD GLY A   2      nope      0.000   0.000\n");
    const auto structure = ascii_protein::PdbParser::parse(pdb, "fixture.pdb");
    require(structure.atoms.size() == 4, "parsed atom count");
    require(structure.residues.size() == 3, "residue grouping");
    require(structure.chains.size() == 2, "chain grouping");
    require(structure.atoms[1].alternateLocation == 'A', "alternate location");
    require(structure.atoms[3].hetero, "HETATM type");
    require(structure.warnings.size() == 1, "malformed line warning");
    require(near(structure.originalCentroid.x, 13.0), "original centroid");
    double centeredX = 0.0;
    for (const auto& atom : structure.atoms) {
        centeredX += atom.position.x;
    }
    require(near(centeredX, 0.0), "centered coordinates");
    require(structure.residues.front().secondary == ascii_protein::SecondaryStructure::helix,
            "HELIX association");
    require(ascii_protein::oneLetterCode("lys") == "K", "amino acid conversion");
    require(ascii_protein::oneLetterCode("UNK") == "X", "unknown amino acid conversion");

    const auto bonds = ascii_protein::inferBonds(structure);
    require(!bonds.empty(), "bond inference");

    ascii_protein::Renderer renderer(structure);
    ascii_protein::Camera camera;
    camera.projection = ascii_protein::Projection::orthographic;
    const auto center = renderer.project({}, camera, 80, 24);
    require(near(center.x, 40.0) && near(center.y, 12.0), "orthographic projection center");
    const auto& frame = renderer.render(camera, {}, 80, 24);
    require(frame.width() == 80 && frame.height() == 24, "framebuffer dimensions");
    require(frame.toPlainText().find("PROTEIN TTY") != std::string::npos, "information panel");

    ascii_protein::RenderSettings braille;
    braille.backend = ascii_protein::RendererBackend::braille;
    braille.mode = ascii_protein::RenderMode::alphaTrace;
    braille.showInfo = false;
    const auto& brailleFrame = renderer.render(camera, braille, 40, 16);
    const auto brailleText = brailleFrame.toPlainText();
    require(brailleText.find("\xe2\xa0") != std::string::npos || brailleText.find("\xe2\xa1") != std::string::npos ||
                brailleText.find("\xe2\xa2") != std::string::npos || brailleText.find("\xe2\xa3") != std::string::npos,
            "Braille subcell rendering");
    require(brailleFrame.toString(true).find("\x1b[38;5;") != std::string::npos,
            "Braille rendering uses depth-aware 256-color output");

    auto inspected = ascii_protein::RenderSettings{};
    inspected.inspectedAtom = renderer.nearestAtom(camera, inspected, 80, 24);
    require(inspected.inspectedAtom.has_value(), "nearest atom inspection");
    require(renderer.render(camera, inspected, 80, 24).toPlainText().find("INSPECTED ATOM") != std::string::npos,
            "inspected atom metadata");
}

void testAlternateConformers() {
    std::istringstream pdb(
        "ATOM      1  CA AALA A   1       1.000   0.000   0.000  0.30 10.00           C  \n"
        "ATOM      2  CA BALA A   1       2.000   0.000   0.000  0.70 11.00           C  \n"
        "ATOM      3  CB AALA A   1       3.000   0.000   0.000  0.50 12.00           C  \n"
        "ATOM      4  CB BALA A   1       4.000   0.000   0.000  0.50 13.00           C  \n"
        "ATOM      5  N   ALA A   1       5.000   0.000   0.000  1.00 14.00           N  \n");
    const auto structure = ascii_protein::PdbParser::parse(pdb, "alternates.pdb");
    require(structure.atoms.size() == 3, "alternate conformer filtering");
    require(structure.atoms[0].alternateLocation == 'B', "highest occupancy conformer");
    require(structure.atoms[1].alternateLocation == 'A', "A conformer tie preference");
    require(structure.atoms[2].alternateLocation == ' ', "ordinary atom preserved");
}

void testSecondaryInsertionCodes() {
    ascii_protein::Structure structure;
    for (const char insertionCode : {' ', 'A', 'B'}) {
        ascii_protein::Atom atom;
        atom.serial = static_cast<int>(structure.atoms.size() + 1);
        atom.name = "O";
        atom.element = "O";
        atom.residueName = "HOH";
        atom.chain = "A";
        atom.residueNumber = 100;
        atom.insertionCode = insertionCode;
        atom.position = {static_cast<double>(structure.atoms.size()), 0.0, 0.0};
        structure.atoms.push_back(atom);
    }
    structure.secondaryRanges.push_back(
        {ascii_protein::SecondaryStructure::sheet, "A", 100, "A", 100, ' ', 'A'});
    structure.finalize();
    require(structure.residues[0].secondary == ascii_protein::SecondaryStructure::sheet,
            "secondary range includes blank insertion code");
    require(structure.residues[1].secondary == ascii_protein::SecondaryStructure::sheet,
            "secondary range includes end insertion code");
    require(structure.residues[2].secondary == ascii_protein::SecondaryStructure::coil,
            "secondary range excludes later insertion code");

    std::string helix(80, ' ');
    helix.replace(0, 5, "HELIX");
    helix[19] = 'A';
    helix.replace(21, 4, " 100");
    helix[25] = 'A';
    helix[31] = 'A';
    helix.replace(33, 4, " 100");
    helix[37] = 'B';
    std::istringstream pdb(
        helix + "\nATOM      1  CA  ALA A 100A      1.000   0.000   0.000  1.00 10.00           C  \n");
    const auto parsed = ascii_protein::PdbParser::parse(pdb, "insertions.pdb");
    require(parsed.secondaryRanges.front().startInsertionCode == 'A' &&
                parsed.secondaryRanges.front().endInsertionCode == 'B',
            "PDB secondary insertion-code columns");
}

void testMmcifParser() {
    std::istringstream cif(
        "data_fixture\n"
        "_entry.details 'O'Brien model'\n"
        "loop_\n"
        "_struct_conf.conf_type_id\n"
        "_struct_conf.beg_auth_asym_id\n"
        "_struct_conf.beg_auth_seq_id\n"
        "_struct_conf.end_auth_asym_id\n"
        "_struct_conf.end_auth_seq_id\n"
        "HELX_P ALPHA 7 ALPHA 7\n"
        "loop_\n"
        "_struct_sheet_range.beg_auth_asym_id\n"
        "_struct_sheet_range.beg_auth_seq_id\n"
        "_struct_sheet_range.pdbx_beg_PDB_ins_code\n"
        "_struct_sheet_range.end_auth_asym_id\n"
        "_struct_sheet_range.end_auth_seq_id\n"
        "_struct_sheet_range.pdbx_end_PDB_ins_code\n"
        "B 8 A B 8 A\n"
        "loop_\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.group_PDB\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "3.0 ATOM 'CA' 1 C ALA ALPHA 7 ? 1.0 2.0 1.0 12.5 . 1\n"
        "6.0 HETATM O 2 O HOH B 8 A 4.0 5.0 ? ? . 1\n"
        "9.0 ATOM N 3 N GLY A bad ? 7.0 8.0 1.0 10.0 . 1\n"
        "12.0 ATOM CA 4 C GLY A 9 ? 10.0 11.0 1.0 10.0 . 2\n"
        "#\n");
    const auto structure = ascii_protein::MmcifParser::parse(cif, "fixture.cif");
    require(structure.atoms.size() == 2, "mmCIF atom and model filtering");
    require(structure.atoms[0].name == "CA" && structure.atoms[0].residueNumber == 7, "mmCIF shuffled columns");
    require(structure.atoms[0].chain == "ALPHA", "mmCIF multi-character chain ID");
    require(structure.atoms[1].hetero && structure.atoms[1].insertionCode == 'A', "mmCIF hetero and insertion");
    require(structure.atoms[1].occupancy == 0.0 && structure.atoms[1].bFactor == 0.0, "mmCIF missing values");
    require(structure.warnings.size() == 1, "mmCIF malformed row warning");
    require(structure.residues[0].secondary == ascii_protein::SecondaryStructure::helix, "mmCIF helix range");
    require(structure.residues[1].secondary == ascii_protein::SecondaryStructure::sheet, "mmCIF sheet range");
    require(structure.secondaryRanges[1].startInsertionCode == 'A' &&
                structure.secondaryRanges[1].endInsertionCode == 'A',
            "mmCIF secondary insertion codes");
}

void testBinaryCifParser() {
    std::istringstream input(binaryCifFixture());
    const auto structure = ascii_protein::BinaryCifParser::parse(input, "fixture.bcif");
    require(structure.atoms.size() == 2, "BinaryCIF atom and model filtering");
    require(structure.residues.size() == 2 && structure.chains.size() == 2, "BinaryCIF hierarchy");
    require(structure.atoms[0].name == "CA" && structure.atoms[0].chain == "ALPHA",
        "BinaryCIF string arrays");
    require(structure.atoms[1].hetero && structure.atoms[1].insertionCode == 'A',
        "BinaryCIF masks and HETATM");
    require(near(structure.atoms[1].occupancy, 0.5) && near(structure.atoms[0].bFactor, 12.5),
        "BinaryCIF floating-point columns");
    require(near(structure.originalCentroid.x, 2.5), "BinaryCIF coordinates");
    require(structure.residues[0].secondary == ascii_protein::SecondaryStructure::helix,
        "BinaryCIF helix range");
    require(structure.residues[1].secondary == ascii_protein::SecondaryStructure::sheet,
        "BinaryCIF sheet range");
    require(structure.secondaryRanges[1].startInsertionCode == 'A' &&
                structure.secondaryRanges[1].endInsertionCode == 'A',
            "BinaryCIF secondary insertion codes");

    std::istringstream malformed("not messagepack");
    try {
    static_cast<void>(ascii_protein::BinaryCifParser::parse(malformed, "bad.bcif"));
    require(false, "malformed BinaryCIF rejected");
    } catch (const std::runtime_error& error) {
    require(std::string(error.what()).find("BinaryCIF") != std::string::npos, "BinaryCIF parse error context");
    }
}

void testBinaryCifLimits() {
    requireBinaryCifRejected(std::string(66, static_cast<char>(0x91)) + static_cast<char>(0xc0),
                             "nesting limit exceeded");
    requireBinaryCifRejected(std::string{"\xdd\x00\x0f\x42\x41", 5}, "array size limit exceeded");
    requireBinaryCifRejected(std::string{"\xdf\x00\x0f\x42\x41", 5}, "map size limit exceeded");
    requireBinaryCifRejected(std::string{"\xdb\x00\x00\x00\x04x", 6}, "truncated string");
    requireBinaryCifRejected(binaryCifEncodedColumnFixture(false), "encoding count limit exceeded");
    requireBinaryCifRejected(binaryCifEncodedColumnFixture(true), "decoded value limit exceeded");
}

ascii_protein::Structure alphaCarbonTrace(const std::vector<ascii_protein::Vec3>& positions,
                                          const std::vector<int>& residueNumbers = {}) {
    ascii_protein::Structure structure;
    for (std::size_t index = 0; index < positions.size(); ++index) {
        ascii_protein::Atom atom;
        atom.serial = static_cast<int>(index + 1);
        atom.name = "CA";
        atom.element = "C";
        atom.residueName = "ALA";
        atom.residueNumber = residueNumbers.empty() ? static_cast<int>(index + 1) : residueNumbers[index];
        atom.chain = "A";
        atom.position = positions[index];
        atom.occupancy = 1.0;
        structure.atoms.push_back(std::move(atom));
    }
    return structure;
}

void testSecondaryStructureFallback() {
    std::vector<ascii_protein::Vec3> helixPositions;
    for (int index = 0; index < 6; ++index) {
        const double angle = ascii_protein::radians(index * 100.0);
        helixPositions.push_back({2.3 * std::cos(angle), 2.3 * std::sin(angle), index * 1.5});
    }
    auto helix = alphaCarbonTrace(helixPositions);
    helix.finalize();
    require(std::all_of(helix.residues.begin(), helix.residues.end(), [](const ascii_protein::Residue& residue) {
                return residue.secondary == ascii_protein::SecondaryStructure::helix;
            }),
            "geometry helix fallback");

    std::vector<ascii_protein::Vec3> sheetPositions;
    for (int index = 0; index < 6; ++index) {
        sheetPositions.push_back({index * 3.3, index % 2 == 0 ? -1.0 : 1.0, 0.0});
    }
    auto sheet = alphaCarbonTrace(sheetPositions);
    sheet.finalize();
    require(std::all_of(sheet.residues.begin(), sheet.residues.end(), [](const ascii_protein::Residue& residue) {
                return residue.secondary == ascii_protein::SecondaryStructure::sheet;
            }),
            "geometry sheet fallback");

    auto explicitSheet = alphaCarbonTrace(helixPositions);
    explicitSheet.secondaryRanges.push_back({ascii_protein::SecondaryStructure::sheet, "A", 2, "A", 4});
    explicitSheet.finalize();
    require(explicitSheet.residues[1].secondary == ascii_protein::SecondaryStructure::sheet &&
                explicitSheet.residues[3].secondary == ascii_protein::SecondaryStructure::sheet,
            "explicit secondary structure wins");

    auto gap = alphaCarbonTrace({helixPositions[0], helixPositions[1], helixPositions[2], helixPositions[3]},
                                {1, 2, 4, 5});
    gap.finalize();
    require(std::all_of(gap.residues.begin(), gap.residues.end(), [](const ascii_protein::Residue& residue) {
                return residue.secondary == ascii_protein::SecondaryStructure::coil;
            }),
            "numbering gap breaks geometry run");

    auto nucleic = alphaCarbonTrace(helixPositions);
    for (auto& atom : nucleic.atoms) atom.residueName = "DA";
    nucleic.finalize();
    require(std::all_of(nucleic.residues.begin(), nucleic.residues.end(), [](const ascii_protein::Residue& residue) {
                return residue.secondary == ascii_protein::SecondaryStructure::coil;
            }),
            "non-protein residues ignored by fallback");
}

void testRcsbCache() {
    require(ascii_protein::normalizePdbId("1CrN") == "1crn", "PDB ID normalization");
    bool rejected = false;
    try {
        (void)ascii_protein::normalizePdbId("../x");
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected, "invalid PDB ID rejection");

    const auto cache = std::filesystem::temp_directory_path() / "ascii-protein-test-cache";
    std::filesystem::remove_all(cache);
    std::filesystem::create_directories(cache);
    const auto expected = cache / "1crn.pdb";
    std::ofstream(expected) << "cached";
    require(ascii_protein::fetchPdb("1CRN", cache) == expected, "offline cache hit");
    const auto expectedCif = cache / "1crn.cif";
    std::ofstream(expectedCif) << "cached cif";
    require(ascii_protein::fetchMmcif("1CRN", cache) == expectedCif, "offline mmCIF cache hit");
    const auto expectedBinaryCif = cache / "1crn.bcif";
    std::ofstream(expectedBinaryCif) << "cached bcif";
    require(ascii_protein::fetchBinaryCif("1CRN", cache) == expectedBinaryCif, "offline BinaryCIF cache hit");
    std::filesystem::remove_all(cache);
}

void testSpatialBondInference() {
    ascii_protein::Structure structure;
    const std::vector<ascii_protein::Vec3> positions{
        {-5.45, 0.0, 0.0}, {-3.0, 0.0, 0.0}, {-2.65, 0.0, 0.0}, {-1.2, 0.0, 0.0},
        {0.0, 0.0, 0.0}, {1.45, 0.0, 0.0}, {2.75, 0.0, 0.0}, {8.0, 0.0, 0.0}};
    for (std::size_t index = 0; index < positions.size(); ++index) {
        ascii_protein::Atom atom;
        atom.serial = static_cast<int>(index + 1);
        atom.name = index % 2 == 0 ? "C" : "N";
        atom.element = atom.name;
        atom.residueName = "GLY";
        atom.residueNumber = static_cast<int>(index / 3) + 1;
        atom.chain = index == positions.size() - 1 ? "B" : "A";
        atom.position = positions[index];
        atom.occupancy = 1.0;
        structure.atoms.push_back(atom);
    }
    structure.finalize();
    const auto inferred = ascii_protein::inferBonds(structure);
    std::set<std::pair<std::size_t, std::size_t>> actual;
    for (const auto& bond : inferred) actual.emplace(bond.first, bond.second);
    require(actual == bruteForceBonds(structure), "spatial bond inference equivalence");
}

void testCommands() {
    ascii_protein::Structure structure;
    structure.sourceName = "command fixture";
    ascii_protein::Atom atom;
    atom.serial = 1;
    atom.name = "CA";
    atom.element = "C";
    atom.residueName = "GLY";
    atom.residueNumber = 1;
    atom.chain = "LONG CHAIN";
    atom.occupancy = 1.0;
    structure.atoms.push_back(atom);
    structure.finalize();
    ascii_protein::Camera camera;
    ascii_protein::RenderSettings render;
    bool autorotate = false;
    bool running = true;
    ascii_protein::CommandState state{camera, render, structure, autorotate, running, {}, {}};

    const auto result = ascii_protein::executeCommands(
        "show cartoon; color secondary; set renderer braille; set projection ortho; spin on; rotate yaw 10; "
        "zoom 2.5; sequence 'LONG CHAIN'; count atoms",
        state);
    require(result.success, "command chain execution");
    require(render.mode == ascii_protein::RenderMode::alphaTrace, "show command");
    require(render.colorMode == ascii_protein::ColorMode::secondary, "color command");
    require(render.backend == ascii_protein::RendererBackend::braille, "renderer command");
    require(camera.projection == ascii_protein::Projection::orthographic, "projection command");
    require(autorotate && near(camera.yaw, -15.0) && near(camera.zoom, 2.5), "numeric and toggle commands");
    require(result.output.size() == 2 && result.output[0] == "G" && result.output[1] == "1", "query commands");

    const auto wrapped = ascii_protein::executeCommands("rotate pitch 1080; rotate yaw -720", state);
    require(wrapped.success && near(camera.pitch, 15.0) && near(camera.yaw, -15.0), "command rotations wrap");

    const auto invalid = ascii_protein::executeCommands("show 'unterminated", state);
    require(!invalid.success && invalid.output.front().find("column 6") != std::string::npos,
            "positioned command error");

    const auto temporary = std::filesystem::temp_directory_path() / "protein-tty-command-load.pdb";
    {
        std::ofstream output(temporary);
        output << "ATOM      7  CA  ALA A   2       1.000   2.000   3.000  1.00 10.00           C  \n";
    }
    int reloads = 0;
    std::filesystem::path saved;
    state.structureChanged = [&] { ++reloads; };
    state.saveFrame = [&](const std::filesystem::path& path) { saved = path; };
    const auto io = ascii_protein::executeCommands(
        "load " + temporary.string() + "; set background black; save command-snapshot.txt", state);
    require(io.success && reloads == 1 && structure.atoms.size() == 1, "load command and renderer notification");
    require(saved == "command-snapshot.txt" && render.backgroundColor == 40, "save and background commands");
        const auto completions = ascii_protein::commandCompletions("set rend", state);
        require(completions.size() == 1 && completions.front() == "renderer", "setting completion");
        render.namedSelections["focus"] = ascii_protein::AtomMask(structure.atoms.size(), true);
        const auto selectionCompletions = ascii_protein::commandCompletions("hide fo", state);
        require(selectionCompletions.size() == 1 && selectionCompletions.front() == "focus", "named selection completion");
        const auto formatCompletions = ascii_protein::commandCompletions("fetch 1crn format b", state);
        require(formatCompletions.size() == 1 && formatCompletions.front() == "bcif", "BinaryCIF format completion");
    std::filesystem::remove(temporary);

        ascii_protein::Framebuffer backgroundFrame(2, 1);
        backgroundFrame.text(0, 0, "x", 91);
        const auto ansi = backgroundFrame.toString(true, 40);
        require(ansi.find("\x1b[40m") != std::string::npos && backgroundFrame.toPlainText() == "x ",
            "ANSI background and plain output");
}

void testSelections() {
    ascii_protein::Structure structure;
    const auto addAtom = [&](int serial, std::string name, std::string element, std::string residueName,
                             int residueNumber, std::string chain, bool hetero) {
        ascii_protein::Atom atom;
        atom.serial = serial;
        atom.name = std::move(name);
        atom.element = std::move(element);
        atom.residueName = std::move(residueName);
        atom.residueNumber = residueNumber;
        atom.chain = std::move(chain);
        atom.hetero = hetero;
        atom.occupancy = 1.0;
        if (serial == 2) atom.position = {1.0, 0.0, 0.0};
        if (serial == 3) atom.position = {1.0, 1.0, 0.0};
        if (serial == 4) atom.position = {1.0, 1.0, 1.0};
        structure.atoms.push_back(std::move(atom));
    };
    addAtom(1, "CA", "C", "ALA", 10, "A", false);
    addAtom(2, "OD1", "O", "ASP", 11, "A", false);
    addAtom(3, "CA", "C", "GLU", 12, "LONG CHAIN", false);
    addAtom(4, "O", "O", "HOH", 20, "B", true);
    structure.finalize();

    const auto acidic = ascii_protein::evaluateSelection("resn ASP+GLU", structure);
    require(ascii_protein::selectionSize(acidic) == 2, "multi-name residue selection");
    require(ascii_protein::selectionSize(ascii_protein::evaluateSelection("resi 10-11 and not name CA", structure)) == 1,
            "range and negation selection");
    require(ascii_protein::selectionSize(
                ascii_protein::evaluateSelection("(chain A or chain 'LONG CHAIN') and element C", structure)) == 2,
            "parenthesized quoted-chain selection");
    require(ascii_protein::selectionSize(ascii_protein::evaluateSelection("hetero and solvent", structure)) == 1,
            "solvent selection");
    ascii_protein::NamedSelections named{{"acidic", acidic}};
    require(ascii_protein::selectionSize(ascii_protein::evaluateSelection("acidic or serial 4", structure, named)) == 3,
            "named selection reference");

    ascii_protein::Camera camera;
    ascii_protein::RenderSettings render;
    bool autorotate = false;
    bool running = true;
    ascii_protein::CommandState state{camera, render, structure, autorotate, running, {}, {}};
    const auto commands = ascii_protein::executeCommands(
        "select acidic resn ASP+GLU; hide solvent; color acidic red; center acidic; zoom acidic; "
        "inspect serial 1; distance 1 2; angle 1 2 3; dihedral 1 2 3 4; count acidic",
        state);
    require(commands.success, "selection command chain");
    require(ascii_protein::selectionSize(render.selectedAtoms) == 2 && !render.visibleAtoms[3],
            "selection and visibility commands");
    require(render.atomColors[1] == 91 && render.atomColors[2] == 91, "selection color command");
    require(render.inspectedAtom == 0 && camera.zoom > 1.0, "selection navigation commands");
    require(commands.output.size() == 5 && commands.output[1] == "1.000 A" && commands.output[2] == "90.000 deg" &&
                commands.output[4] == "2",
            "measurement and selection count output");
    try {
        static_cast<void>(ascii_protein::evaluateSelection("chain A and", structure));
        require(false, "incomplete selection rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("column") != std::string::npos, "positioned selection error");
    }
    try {
        static_cast<void>(ascii_protein::evaluateSelection("resi 10oops", structure));
        require(false, "partial residue number rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("invalid residue range") != std::string::npos,
                "strict residue number error");
    }

    const auto undefinedAngle = ascii_protein::executeCommands("angle 1 1 2", state);
    require(!undefinedAngle.success && undefinedAngle.output.front().find("undefined") != std::string::npos,
            "degenerate angle rejected");

    ascii_protein::Structure repeatedNumbers;
    const auto addResidueAtom = [&](int serial, std::string name, std::string chain, int residueNumber) {
        ascii_protein::Atom atom;
        atom.serial = serial;
        atom.name = std::move(name);
        atom.element = "C";
        atom.residueName = "ALA";
        atom.residueNumber = residueNumber;
        atom.chain = std::move(chain);
        atom.occupancy = 1.0;
        repeatedNumbers.atoms.push_back(std::move(atom));
    };
    addResidueAtom(1, "N", "A", 10);
    addResidueAtom(2, "CA", "A", 10);
    addResidueAtom(3, "CA", "B", 10);
    addResidueAtom(4, "CA", "A", 11);
    repeatedNumbers.finalize();
    ascii_protein::Camera repeatedCamera;
    ascii_protein::RenderSettings repeatedRender;
    ascii_protein::CommandState repeatedState{
        repeatedCamera, repeatedRender, repeatedNumbers, autorotate, running, {}, {}};
    const auto ambiguous = ascii_protein::executeCommands("select ten resi 10", repeatedState);
    require(ambiguous.success && ambiguous.output.front() == "3 atoms selected (2 residues)",
            "residue number spans chains and reports both residues");
    const auto oneResidue = ascii_protein::executeCommands(
        "select one chain A and resi 10; show one", repeatedState);
    require(oneResidue.success && oneResidue.output.front() == "2 atoms selected (1 residue)",
            "chain-qualified residue feedback");
    require(ascii_protein::selectionSize(repeatedRender.visibleAtoms) == 2,
            "show selection hides nonmatching atoms");
    require(!repeatedRender.selectedResidue, "no phantom default residue highlight");
}

void testCommandCliOptions() {
    char executable[] = "ascii-protein";
    char input[] = "model.pdb";
    char commandOption[] = "--command";
    char firstCommand[] = "show cartoon";
    char scriptOption[] = "--script";
    char script[] = "demo.ptty";
    char secondCommand[] = "count atoms";
    char trajectoryOption[] = "--trajectory";
    char trajectory[] = "motion.xyz";
    char delayOption[] = "--frame-delay";
    char delay[] = "40";
    char loop[] = "--loop";
    char* arguments[] = {executable, input, commandOption, firstCommand, scriptOption, script, commandOption,
                         secondCommand, trajectoryOption, trajectory, delayOption, delay, loop};
    const auto options = ascii_protein::parseCli(13, arguments);
    require(options.commands.size() == 2 && options.commands[1] == "count atoms", "repeatable command options");
    require(options.scripts.size() == 1 && options.scripts[0] == "demo.ptty", "script option");
        require(options.trajectory == "motion.xyz" && options.frameDelayMs == 40 && options.trajectoryLoop,
            "trajectory CLI options");

    char pdbOption[] = "--pdb";
    char pdbId[] = "1CRN";
    char formatOption[] = "--format";
    char binaryCif[] = "bcif";
    char* binaryCifArguments[] = {executable, pdbOption, pdbId, formatOption, binaryCif};
    const auto binaryCifOptions = ascii_protein::parseCli(5, binaryCifArguments);
    require(binaryCifOptions.remoteFormat == ascii_protein::RemoteFormat::bcif, "BinaryCIF CLI format");

    char* missingValue[] = {executable, input, commandOption};
    try {
        static_cast<void>(ascii_protein::parseCli(3, missingValue));
        require(false, "missing command value rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("missing value") != std::string::npos, "missing command value error");
    }
    char invalidDelay[] = "9";
    char* invalidDelayArguments[] = {executable, input, delayOption, invalidDelay};
    try {
        static_cast<void>(ascii_protein::parseCli(4, invalidDelayArguments));
        require(false, "invalid frame delay rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("invalid frame delay") != std::string::npos,
                "invalid frame delay error");
    }
    char widthOption[] = "--width";
    char invalidWidth[] = "80columns";
    char* invalidWidthArguments[] = {executable, input, widthOption, invalidWidth};
    try {
        static_cast<void>(ascii_protein::parseCli(4, invalidWidthArguments));
        require(false, "dimension with trailing text rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("invalid value for --width") != std::string::npos,
                "invalid dimension error");
    }
}

void testXyzTrajectory() {
    const auto path = std::filesystem::temp_directory_path() / "protein-tty-trajectory.xyz";
    {
    std::ofstream output(path);
    output << "2\nframe zero\nC 0 0 0\nO 2 0 0\n\n"
          "2\nframe one\nC 0 0 0\nO 4 0 0\n";
    }
    ascii_protein::XyzTrajectorySource source(path);
    require(source.atomCount() == 2 && source.frameCount() == 2, "XYZ trajectory dimensions");
    require(source.frameLabel(1) == "frame one" && near(source.frame(1)[1].x, 4.0), "XYZ frame data");

    ascii_protein::Structure structure;
    for (int index = 0; index < 2; ++index) {
    ascii_protein::Atom atom;
    atom.serial = index + 1;
    atom.name = index == 0 ? "C" : "O";
    atom.element = atom.name;
    atom.residueName = "GLY";
    atom.residueNumber = 1;
    atom.chain = "A";
    atom.occupancy = 1.0;
    structure.atoms.push_back(std::move(atom));
    }
    structure.finalize();
    ascii_protein::TrajectoryState trajectory;
    ascii_protein::loadTrajectory(path, structure, trajectory);
    require(trajectory.playing && trajectory.frameIndex == 0 && near(structure.atoms[1].position.x, 1.0),
        "trajectory loads first centered frame");
    require(ascii_protein::stepTrajectory(structure, trajectory) && trajectory.frameIndex == 1 &&
        near(structure.atoms[1].position.x, 2.0),
        "trajectory advances frame");
    require(!ascii_protein::stepTrajectory(structure, trajectory), "trajectory stops at final frame");
    trajectory.loop = true;
    require(ascii_protein::stepTrajectory(structure, trajectory) && trajectory.frameIndex == 0,
        "trajectory loops to first frame");
    require(ascii_protein::stepTrajectory(structure, trajectory, -1) && trajectory.frameIndex == 1,
        "trajectory reverse step loops");
    ascii_protein::stopTrajectory(structure, trajectory);
    require(!trajectory.playing && trajectory.frameIndex == 0 && near(structure.atoms[1].position.x, 1.0),
        "trajectory stop resets first frame");
    require(ascii_protein::trajectoryStatus(trajectory).find("Frame 1/2") != std::string::npos,
        "trajectory status");

    ascii_protein::Camera camera;
    ascii_protein::RenderSettings render;
    bool autorotate = false;
    bool running = true;
    ascii_protein::CommandState commandState{
        camera, render, structure, autorotate, running, {}, {}, &trajectory};
    const auto commands = ascii_protein::executeCommands(
        "pause; step 1; set trajectory-loop on; step 1; set frame-delay 25; resume", commandState);
    require(commands.success && trajectory.playing && trajectory.loop && trajectory.frameDelayMs == 25 &&
                trajectory.frameIndex == 0,
            "trajectory command controls");
    const auto stopped = ascii_protein::executeCommands("stop", commandState);
    require(stopped.success && !trajectory.playing && trajectory.frameIndex == 0, "trajectory stop command");

    ascii_protein::Structure wrongSize;
    wrongSize.atoms.push_back(structure.atoms.front());
    wrongSize.finalize();
    try {
        ascii_protein::TrajectoryState wrongTrajectory;
        ascii_protein::loadTrajectory(path, wrongSize, wrongTrajectory);
        require(false, "trajectory atom mismatch rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("trajectory has 2 atoms; structure has 1") != std::string::npos,
                "trajectory atom mismatch error");
    }

    const auto malformed = std::filesystem::temp_directory_path() / "protein-tty-malformed.xyz";
    {
    std::ofstream output(malformed);
    output << "2\nbad\nC 0 0 0\nO nope 0 0\n";
    }
    try {
    static_cast<void>(ascii_protein::XyzTrajectorySource(malformed));
    require(false, "malformed XYZ rejected");
    } catch (const std::runtime_error& error) {
    require(std::string(error.what()).find("frame 1, atom 2") != std::string::npos,
        "positioned XYZ error");
    }
    const auto changedLabel = std::filesystem::temp_directory_path() / "protein-tty-label-change.xyz";
    {
        std::ofstream output(changedLabel);
        output << "1\nfirst\nC 0 0 0\n1\nsecond\nN 0 0 0\n";
    }
    try {
        static_cast<void>(ascii_protein::XyzTrajectorySource(changedLabel));
        require(false, "changed XYZ label rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("atom label changed") != std::string::npos,
                "changed XYZ label error");
    }
    const auto negativeCount = std::filesystem::temp_directory_path() / "protein-tty-negative-count.xyz";
    {
        std::ofstream output(negativeCount);
        output << "-1\ninvalid\n";
    }
    try {
        static_cast<void>(ascii_protein::XyzTrajectorySource(negativeCount));
        require(false, "negative XYZ atom count rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("invalid atom count") != std::string::npos,
                "negative XYZ atom count error");
    }
    auto mismatchedElements = structure;
    mismatchedElements.atoms[1].element = "N";
    try {
        ascii_protein::TrajectoryState mismatchedTrajectory;
        ascii_protein::loadTrajectory(path, mismatchedElements, mismatchedTrajectory);
        require(false, "trajectory topology label mismatch rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("topology atom") != std::string::npos,
                "trajectory topology label mismatch error");
    }
    std::filesystem::remove(path);
    std::filesystem::remove(malformed);
    std::filesystem::remove(changedLabel);
    std::filesystem::remove(negativeCount);
}

void testDcdTrajectory() {
    const auto path = std::filesystem::temp_directory_path() / "protein-tty-trajectory.dcd";
    {
        std::ofstream output(path, std::ios::binary);
        const auto fixture = dcdFixture(false, true);
        output.write(fixture.data(), static_cast<std::streamsize>(fixture.size()));
    }
    ascii_protein::DcdTrajectorySource source(path);
    require(source.atomCount() == 2 && source.frameCount() == 2, "DCD trajectory dimensions");
    require(source.atomLabel(0).empty() && source.frameLabel(1).empty(), "DCD label-free metadata");
    require(near(source.frame(1)[1].x, 4.0) && near(source.frame(0)[1].x, 2.0),
            "DCD random frame access");

    ascii_protein::Structure structure;
    for (int index = 0; index < 2; ++index) {
        ascii_protein::Atom atom;
        atom.serial = index + 1;
        atom.name = index == 0 ? "CA" : "N";
        atom.element = index == 0 ? "C" : "N";
        atom.residueName = "GLY";
        atom.residueNumber = 1;
        atom.chain = "A";
        atom.occupancy = 1.0;
        structure.atoms.push_back(std::move(atom));
    }
    structure.finalize();
    ascii_protein::TrajectoryState trajectory;
    ascii_protein::loadTrajectory(path, structure, trajectory);
    require(trajectory.frameIndex == 0 && near(structure.atoms[1].position.x, 1.0),
            "DCD trajectory dispatch and centering");
    require(ascii_protein::stepTrajectory(structure, trajectory) && near(structure.atoms[1].position.x, 2.0),
            "DCD trajectory playback");

    const auto bigEndianPath = std::filesystem::temp_directory_path() / "protein-tty-big-endian.dcd";
    {
        std::ofstream output(bigEndianPath, std::ios::binary);
        const auto fixture = dcdFixture(true);
        output.write(fixture.data(), static_cast<std::streamsize>(fixture.size()));
    }
    ascii_protein::DcdTrajectorySource bigEndian(bigEndianPath);
    require(near(bigEndian.frame(1)[1].x, 4.0), "big-endian DCD decoding");

    const auto fixedPath = std::filesystem::temp_directory_path() / "protein-tty-fixed-atoms.dcd";
    {
        std::ofstream output(fixedPath, std::ios::binary);
        const auto fixture = fixedAtomDcdFixture();
        output.write(fixture.data(), static_cast<std::streamsize>(fixture.size()));
    }
    ascii_protein::DcdTrajectorySource fixed(fixedPath);
    const auto& fixedFrame = fixed.frame(1);
    require(near(fixedFrame[0].x, 1.0) && near(fixedFrame[0].y, 2.0) && near(fixedFrame[0].z, 3.0) &&
                near(fixedFrame[1].x, 10.0) && near(fixedFrame[2].x, 21.0),
            "DCD fixed atoms retain first-frame coordinates");

    const auto truncatedPath = std::filesystem::temp_directory_path() / "protein-tty-truncated.dcd";
    {
        std::ofstream output(truncatedPath, std::ios::binary);
        auto fixture = dcdFixture();
        fixture.resize(fixture.size() - 2);
        output.write(fixture.data(), static_cast<std::streamsize>(fixture.size()));
    }
    try {
        static_cast<void>(ascii_protein::DcdTrajectorySource(truncatedPath));
        require(false, "truncated DCD rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("truncated DCD") != std::string::npos,
                "truncated DCD error");
    }
    std::filesystem::remove(path);
    std::filesystem::remove(bigEndianPath);
    std::filesystem::remove(fixedPath);
    std::filesystem::remove(truncatedPath);
}

} // namespace

int main() {
    try {
        testGeometry();
        testTerminalInputDecoder();
        testParserAndModel();
        testAlternateConformers();
        testSecondaryInsertionCodes();
        testMmcifParser();
        testBinaryCifParser();
        testBinaryCifLimits();
        testSecondaryStructureFallback();
        testRcsbCache();
        testSpatialBondInference();
        testCommands();
        testSelections();
        testCommandCliOptions();
        testXyzTrajectory();
        testDcdTrajectory();
        std::cout << "All tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}
