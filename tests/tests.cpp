#include "ascii_protein/geometry.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/renderer.hpp"
#include "ascii_protein/structure.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
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

void testGeometry() {
    using namespace ascii_protein;
    const Vec3 x{1, 0, 0};
    const Vec3 y{0, 1, 0};
    require(dot(x, y) == 0.0, "orthogonal dot product");
    require(cross(x, y).z == 1.0, "cross product");
    require(near(length(Vec3{3, 4, 0}), 5.0), "vector length");
    const Vec3 rotated = Mat3::rotationZ(radians(90.0)) * x;
    require(near(rotated.x, 0.0) && near(rotated.y, 1.0), "z rotation");
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
    std::filesystem::remove_all(cache);
}

} // namespace

int main() {
    try {
        testGeometry();
        testParserAndModel();
        testAlternateConformers();
        testRcsbCache();
        std::cout << "All tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}