#include "ascii_protein/structure.hpp"

#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>

namespace ascii_protein {

void Structure::finalize() {
    if (atoms.empty()) {
        boundingRadius = 1.0;
        return;
    }

    originalCentroid = {};
    for (const auto& atom : atoms) {
        originalCentroid += atom.position;
    }
    originalCentroid = originalCentroid / static_cast<double>(atoms.size());

    boundsMin = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max()};
    boundsMax = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                 std::numeric_limits<double>::lowest()};
    boundingRadius = 0.0;
    for (auto& atom : atoms) {
        atom.position = atom.position - originalCentroid;
        boundsMin.x = std::min(boundsMin.x, atom.position.x);
        boundsMin.y = std::min(boundsMin.y, atom.position.y);
        boundsMin.z = std::min(boundsMin.z, atom.position.z);
        boundsMax.x = std::max(boundsMax.x, atom.position.x);
        boundsMax.y = std::max(boundsMax.y, atom.position.y);
        boundsMax.z = std::max(boundsMax.z, atom.position.z);
        boundingRadius = std::max(boundingRadius, length(atom.position));
    }
    boundingRadius = std::max(1.0, boundingRadius);

    residues.clear();
    chains.clear();
    std::map<std::tuple<char, int, char>, std::size_t> residueLookup;
    std::map<char, std::size_t> chainLookup;
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
                residue.number >= range.startResidue && residue.number <= range.endResidue) {
                residue.secondary = range.type;
                break;
            }
        }
    }
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