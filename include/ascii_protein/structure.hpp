#pragma once

#include "ascii_protein/geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <optional>
#include <utility>

namespace ascii_protein {

enum class SecondaryStructure { coil, helix, sheet };

struct Atom {
    int model{1};
    std::string labelChain, labelName;
    int labelResidue{};
    int serial{};
    std::string name;
    char alternateLocation{' '};
    std::string element;
    std::string residueName;
    int residueNumber{};
    char insertionCode{' '};
    std::string chain;
    Vec3 position;
    double occupancy{};
    double bFactor{};
    bool hetero{};
    std::size_t residueIndex{};
};

struct Residue {
    std::string name;
    int number{};
    char insertionCode{' '};
    std::string chain;
    std::vector<std::size_t> atomIndices;
    SecondaryStructure secondary{SecondaryStructure::coil};
};

struct Chain {
    std::string id;
    std::vector<std::size_t> residueIndices;
};

struct SecondaryRange {
    SecondaryStructure type{SecondaryStructure::coil};
    std::string startChain;
    int startResidue{};
    std::string endChain;
    int endResidue{};
    char startInsertionCode{' '};
    char endInsertionCode{' '};
};

struct AssemblyOperator {
    Mat3 rotation{Mat3::identity()};
    Vec3 translation;
};
struct AssemblyGeneration {
    std::string id, expression;
    std::vector<std::string> chains;
};
struct Structure {
    std::string sourceName;
    std::filesystem::path sourcePath;
    std::string assemblyApplied;
    std::vector<Atom> atoms;
    std::vector<Atom> depositedAtoms;
    int model{1};
    std::vector<std::pair<int, int>> explicitBonds;
    std::vector<std::map<std::string, std::string>> connections;
    std::map<std::string, AssemblyOperator> operators;
    std::vector<AssemblyGeneration> assemblies;
    std::vector<Residue> residues;
    std::vector<Chain> chains;
    std::vector<SecondaryRange> secondaryRanges;
    std::vector<std::string> warnings;
    Vec3 originalCentroid;
    Vec3 boundsMin;
    Vec3 boundsMax;
    double boundingRadius{1.0};

    void finalize();
    void updatePositions(const std::vector<Vec3>& positions, std::optional<Vec3> origin = std::nullopt);
};

std::string oneLetterCode(const std::string& residueName);
bool isBackboneAtom(const Atom& atom);

} // namespace ascii_protein