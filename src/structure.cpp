#include "ascii_protein/structure.hpp"
#include "ascii_protein/secondary_assignment.hpp"

#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace ascii_protein {
namespace {

void selectAlternateConformers(std::vector<Atom>& source) {
    using AtomKey = std::tuple<bool, std::string, int, char, std::string, std::string>;
    std::map<AtomKey, std::size_t> selected;
    std::vector<Atom> atoms;
    atoms.reserve(source.size());
    for (auto& atom : source) {
        const AtomKey key{atom.hetero, atom.chain, atom.residueNumber, atom.insertionCode, atom.residueName, atom.name};
        const auto found = selected.find(key);
        if (found == selected.end()) {
            selected.emplace(key, atoms.size());
            atoms.push_back(std::move(atom));
            continue;
        }
        auto& current = atoms[found->second];
        const bool higherOccupancy = atom.occupancy > current.occupancy;
        const bool preferredTie = atom.occupancy == current.occupancy && current.alternateLocation != ' ' &&
                                  (atom.alternateLocation == ' ' || atom.alternateLocation == 'A');
        if (higherOccupancy || preferredTie) current = std::move(atom);
    }
    source = std::move(atoms);
}

void updateCoordinateBounds(Structure& structure) {
    structure.originalCentroid = {};
    for (const auto& atom : structure.atoms) structure.originalCentroid += atom.position;
    structure.originalCentroid = structure.originalCentroid / static_cast<double>(structure.atoms.size());
    structure.boundsMin = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::max()};
    structure.boundsMax = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                           std::numeric_limits<double>::lowest()};
    structure.boundingRadius = 0.0;
    for (auto& atom : structure.atoms) {
        atom.position = atom.position - structure.originalCentroid;
        structure.boundsMin.x = std::min(structure.boundsMin.x, atom.position.x);
        structure.boundsMin.y = std::min(structure.boundsMin.y, atom.position.y);
        structure.boundsMin.z = std::min(structure.boundsMin.z, atom.position.z);
        structure.boundsMax.x = std::max(structure.boundsMax.x, atom.position.x);
        structure.boundsMax.y = std::max(structure.boundsMax.y, atom.position.y);
        structure.boundsMax.z = std::max(structure.boundsMax.z, atom.position.z);
        structure.boundingRadius = std::max(structure.boundingRadius, length(atom.position));
    }
    structure.boundingRadius = std::max(1.0, structure.boundingRadius);
}

} // namespace

void Structure::finalize() {
    selectAlternateConformers(atoms);
    if (atoms.empty()) {
        boundingRadius = 1.0;
        return;
    }

    updateCoordinateBounds(*this);

    residues.clear();
    chains.clear();
    std::map<std::tuple<std::string, int, char>, std::size_t> residueLookup;
    std::map<std::string, std::size_t> chainLookup;
    for (std::size_t atomIndex = 0; atomIndex < atoms.size(); ++atomIndex) {
        auto& atom = atoms[atomIndex];
        const auto residueKey = std::make_tuple(atom.chain, atom.residueNumber, atom.insertionCode);
        auto [residueIt, inserted] = residueLookup.emplace(residueKey, residues.size());
        if (inserted) {
            residues.push_back(
                {atom.residueName, atom.residueNumber, atom.insertionCode, atom.chain, {}, SecondaryStructure::coil});
            auto [chainIt, chainInserted] = chainLookup.emplace(atom.chain, chains.size());
            if (chainInserted) {
                chains.push_back({atom.chain, {}});
            }
            chains[chainIt->second].residueIndices.push_back(residueIt->second);
        }
        atom.residueIndex = residueIt->second;
        residues[residueIt->second].atomIndices.push_back(atomIndex);
    }

    for (auto& residue : residues) {
        for (const auto& range : secondaryRanges) {
            if (residue.chain == range.startChain && residue.chain == range.endChain &&
                std::tie(residue.number, residue.insertionCode) >=
                    std::tie(range.startResidue, range.startInsertionCode) &&
                std::tie(residue.number, residue.insertionCode) <=
                    std::tie(range.endResidue, range.endInsertionCode)) {
                residue.secondary = range.type;
                break;
            }
        }
    }
    assignSecondaryStructureFromGeometry(*this);
}

void Structure::updatePositions(const std::vector<Vec3>& positions) {
    if (positions.size() != atoms.size()) throw std::runtime_error("coordinate count does not match structure");
    for (std::size_t index = 0; index < atoms.size(); ++index) atoms[index].position = positions[index];
    if (!atoms.empty()) updateCoordinateBounds(*this);
}

std::string oneLetterCode(const std::string& residueName) {
    static const std::unordered_map<std::string, char> codes{
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
        {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
        {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
        {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
        {"SEC", 'U'}, {"PYL", 'O'}, {"ASX", 'B'}, {"GLX", 'Z'}};
    std::string upper = residueName;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
    const auto found = codes.find(upper);
    return found == codes.end() ? "X" : std::string(1, found->second);
}

bool isBackboneAtom(const Atom& atom) {
    return atom.name == "N" || atom.name == "CA" || atom.name == "C" || atom.name == "O";
}

} // namespace ascii_protein
