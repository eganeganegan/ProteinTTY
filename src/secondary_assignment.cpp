#include "ascii_protein/secondary_assignment.hpp"

#include "ascii_protein/structure.hpp"

#include <array>
#include <optional>
#include <vector>

namespace ascii_protein {
namespace {

std::optional<Vec3> alphaCarbon(const Structure& structure, const Residue& residue) {
    if (oneLetterCode(residue.name) == "X") return std::nullopt;
    for (const auto atomIndex : residue.atomIndices) {
        if (structure.atoms[atomIndex].name == "CA") return structure.atoms[atomIndex].position;
    }
    return std::nullopt;
}

bool between(double value, double minimum, double maximum) {
    return value >= minimum && value <= maximum;
}

bool residuesAreSequential(const Residue& first, const Residue& second) {
    return second.number == first.number + 1 ||
           (second.number == first.number && first.insertionCode != second.insertionCode);
}

SecondaryStructure classifyWindow(const std::array<Vec3, 4>& positions) {
    for (std::size_t index = 1; index < positions.size(); ++index) {
        if (!between(length(positions[index] - positions[index - 1]), 3.4, 4.2)) {
            return SecondaryStructure::coil;
        }
    }
    const double secondNeighbor = length(positions[2] - positions[0]);
    const double thirdNeighbor = length(positions[3] - positions[0]);
    if (between(secondNeighbor, 4.8, 6.2) && between(thirdNeighbor, 4.5, 6.3)) {
        return SecondaryStructure::helix;
    }
    if (between(secondNeighbor, 6.0, 7.4) && between(thirdNeighbor, 8.5, 11.5)) {
        return SecondaryStructure::sheet;
    }
    return SecondaryStructure::coil;
}

} // namespace

void assignSecondaryStructureFromGeometry(Structure& structure) {
    std::vector<unsigned int> helixVotes(structure.residues.size());
    std::vector<unsigned int> sheetVotes(structure.residues.size());
    for (const auto& chain : structure.chains) {
        if (chain.residueIndices.size() < 4) continue;
        for (std::size_t start = 0; start + 3 < chain.residueIndices.size(); ++start) {
            std::array<Vec3, 4> positions;
            bool complete = true;
            for (std::size_t offset = 0; offset < positions.size(); ++offset) {
                const auto& residue = structure.residues[chain.residueIndices[start + offset]];
                if (offset > 0) {
                    const auto& previous = structure.residues[chain.residueIndices[start + offset - 1]];
                    if (!residuesAreSequential(previous, residue)) {
                        complete = false;
                        break;
                    }
                }
                const auto position = alphaCarbon(structure, residue);
                if (!position) {
                    complete = false;
                    break;
                }
                positions[offset] = *position;
            }
            if (!complete) continue;
            const auto type = classifyWindow(positions);
            if (type == SecondaryStructure::coil) continue;
            auto& votes = type == SecondaryStructure::helix ? helixVotes : sheetVotes;
            for (std::size_t offset = 0; offset < positions.size(); ++offset) {
                ++votes[chain.residueIndices[start + offset]];
            }
        }
    }

    for (std::size_t index = 0; index < structure.residues.size(); ++index) {
        auto& residue = structure.residues[index];
        if (residue.secondary != SecondaryStructure::coil) continue;
        if (helixVotes[index] > sheetVotes[index]) residue.secondary = SecondaryStructure::helix;
        else if (sheetVotes[index] > 0) residue.secondary = SecondaryStructure::sheet;
    }
}

} // namespace ascii_protein