#pragma once

#include "ascii_protein/geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace ascii_protein {

enum class SecondaryStructure { coil, helix, sheet };

struct Atom {
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
};

struct Structure {
    std::string sourceName;
    std::vector<Atom> atoms;
    std::vector<Residue> residues;
    std::vector<Chain> chains;
    std::vector<SecondaryRange> secondaryRanges;
    std::vector<std::string> warnings;
    Vec3 originalCentroid;
    Vec3 boundsMin;
    Vec3 boundsMax;
    double boundingRadius{1.0};

    void finalize();
};

std::string oneLetterCode(const std::string& residueName);
bool isBackboneAtom(const Atom& atom);

} // namespace ascii_protein