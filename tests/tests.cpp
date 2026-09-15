#include "ascii_protein/command.hpp"
#include "ascii_protein/cli.hpp"
#include "ascii_protein/geometry.hpp"
#include "ascii_protein/mmcif_parser.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/renderer.hpp"
#include "ascii_protein/selection.hpp"
#include "ascii_protein/structure.hpp"
#include "ascii_protein/structure_loader.hpp"

#include <cmath>
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
    braille.showInfo = false;
    const auto brailleText = renderer.render(camera, braille, 40, 16).toPlainText();
    require(brailleText.find("\xe2\xa0") != std::string::npos || brailleText.find("\xe2\xa1") != std::string::npos ||
                brailleText.find("\xe2\xa2") != std::string::npos || brailleText.find("\xe2\xa3") != std::string::npos,
            "Braille subcell rendering");

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
        "_struct_sheet_range.end_auth_asym_id\n"
        "_struct_sheet_range.end_auth_seq_id\n"
        "B 8 B 8\n"
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
    char* arguments[] = {
        executable, input, commandOption, firstCommand, scriptOption, script, commandOption, secondCommand};
    const auto options = ascii_protein::parseCli(8, arguments);
    require(options.commands.size() == 2 && options.commands[1] == "count atoms", "repeatable command options");
    require(options.scripts.size() == 1 && options.scripts[0] == "demo.ptty", "script option");

    char* missingValue[] = {executable, input, commandOption};
    try {
        static_cast<void>(ascii_protein::parseCli(3, missingValue));
        require(false, "missing command value rejected");
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("missing value") != std::string::npos, "missing command value error");
    }
}

} // namespace

int main() {
    try {
        testGeometry();
        testParserAndModel();
        testAlternateConformers();
        testMmcifParser();
        testRcsbCache();
        testSpatialBondInference();
        testCommands();
        testSelections();
        testCommandCliOptions();
        std::cout << "All tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}