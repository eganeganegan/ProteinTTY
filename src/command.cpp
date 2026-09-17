#include "ascii_protein/command.hpp"
#include "ascii_protein/selection.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/structure_loader.hpp"
#include "ascii_protein/trajectory.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace ascii_protein {
namespace {

struct Token {
    std::string text;
    std::size_t column{};
};

using Command = std::vector<Token>;

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return value;
}

std::vector<Command> tokenize(std::string_view source) {
    std::vector<Command> commands(1);
    std::size_t position = 0;
    while (position < source.size()) {
        while (position < source.size() && std::isspace(static_cast<unsigned char>(source[position]))) ++position;
        if (position >= source.size()) break;
        if (source[position] == ';') {
            if (!commands.back().empty()) commands.emplace_back();
            ++position;
            continue;
        }
        const std::size_t column = position + 1;
        std::string text;
        if (source[position] == '\'' || source[position] == '"') {
            const char quote = source[position++];
            while (position < source.size() && source[position] != quote) {
                if (source[position] == '\\' && position + 1 < source.size()) ++position;
                text.push_back(source[position++]);
            }
            if (position >= source.size()) {
                throw std::runtime_error("column " + std::to_string(column) + ": unterminated quoted string");
            }
            ++position;
        } else {
            while (position < source.size() && !std::isspace(static_cast<unsigned char>(source[position])) &&
                   source[position] != ';') {
                text.push_back(source[position++]);
            }
        }
        commands.back().push_back({std::move(text), column});
    }
    if (commands.back().empty()) commands.pop_back();
    return commands;
}

[[noreturn]] void fail(const Command& command, std::size_t token, const std::string& message) {
    const std::size_t column = token < command.size() ? command[token].column : command.back().column;
    throw std::runtime_error("column " + std::to_string(column) + ": " + message);
}

bool toggleValue(const Command& command, std::size_t index, bool current) {
    if (index >= command.size()) fail(command, command.size() - 1, "expected on, off, or toggle");
    const std::string value = lower(command[index].text);
    if (value == "on") return true;
    if (value == "off") return false;
    if (value == "toggle") return !current;
    fail(command, index, "expected on, off, or toggle");
}

double numericValue(const Command& command, std::size_t index) {
    if (index >= command.size()) fail(command, command.size() - 1, "expected a number");
    try {
        std::size_t used = 0;
        const double value = std::stod(command[index].text, &used);
        if (used != command[index].text.size()) fail(command, index, "invalid number: " + command[index].text);
        return value;
    } catch (const std::invalid_argument&) {
        fail(command, index, "invalid number: " + command[index].text);
    } catch (const std::out_of_range&) {
        fail(command, index, "number out of range: " + command[index].text);
    }
}

std::string sequence(const Structure& structure, std::string_view chain) {
    std::string result;
    for (const auto& residue : structure.residues) {
        if (chain.empty() || residue.chain == chain) result += oneLetterCode(residue.name);
    }
    return result;
}

void requireArguments(const Command& command, std::size_t count, std::string_view usage) {
    if (command.size() != count) fail(command, command.size() - 1, "usage: " + std::string(usage));
}

std::string selectionText(const Command& command, std::size_t begin, std::size_t end) {
    std::string result;
    for (std::size_t index = begin; index < end; ++index) {
        if (!result.empty()) result += ' ';
        if (command[index].text.find(' ') != std::string::npos) result += '"' + command[index].text + '"';
        else result += command[index].text;
    }
    return result;
}

AtomMask selectionFrom(const Command& command, std::size_t begin, std::size_t end, CommandState state) {
    if (begin >= end) fail(command, command.size() - 1, "expected a selection expression");
    try {
        return evaluateSelection(selectionText(command, begin, end), state.structure, state.render.namedSelections);
    } catch (const std::runtime_error& error) {
        fail(command, begin, error.what());
    }
}

int colorCode(const Command& command, std::size_t index) {
    const std::string name = lower(command[index].text);
    if (name == "default") return 0;
    if (name == "black") return 30;
    if (name == "red") return 91;
    if (name == "green") return 92;
    if (name == "yellow") return 93;
    if (name == "blue") return 94;
    if (name == "magenta") return 95;
    if (name == "cyan") return 96;
    if (name == "white") return 97;
    fail(command, index, "unknown color: " + command[index].text);
}

void ensureAtomState(RenderSettings& render, std::size_t atomCount) {
    if (render.visibleAtoms.size() != atomCount) render.visibleAtoms.assign(atomCount, true);
    if (render.selectedAtoms.size() != atomCount) render.selectedAtoms.assign(atomCount, false);
    if (render.atomColors.size() != atomCount) render.atomColors.assign(atomCount, 0);
}

void resetStructureState(CommandState state) {
    state.render.visibleAtoms.assign(state.structure.atoms.size(), true);
    state.render.selectedAtoms.assign(state.structure.atoms.size(), false);
    state.render.atomColors.assign(state.structure.atoms.size(), 0);
    state.render.namedSelections.clear();
    state.render.inspectedAtom.reset();
    state.render.selectedResidue.reset();
    state.camera.reset();
        if (state.trajectory) {
            state.trajectory->source.reset();
            state.trajectory->path.clear();
            state.trajectory->frameIndex = 0;
            state.trajectory->playing = false;
        }
    if (state.structureChanged) state.structureChanged();
}

std::vector<std::size_t> selectedIndices(const AtomMask& selection) {
    std::vector<std::size_t> result;
    for (std::size_t index = 0; index < selection.size(); ++index) {
        if (selection[index]) result.push_back(index);
    }
    return result;
}

std::size_t selectedResidueCount(const AtomMask& selection, const Structure& structure) {
    std::vector<bool> residues(structure.residues.size());
    for (std::size_t index = 0; index < selection.size(); ++index) {
        if (selection[index] && structure.atoms[index].residueIndex < residues.size()) {
            residues[structure.atoms[index].residueIndex] = true;
        }
    }
    return static_cast<std::size_t>(std::count(residues.begin(), residues.end(), true));
}

std::size_t atomOperand(const Command& command, std::size_t index, CommandState state) {
    if (index >= command.size()) fail(command, command.size() - 1, "expected an atom serial or named selection");
    try {
        std::size_t used = 0;
        const int serial = std::stoi(command[index].text, &used);
        if (used == command[index].text.size()) {
            const auto found = std::find_if(state.structure.atoms.begin(), state.structure.atoms.end(),
                                            [&](const Atom& atom) { return atom.serial == serial; });
            if (found == state.structure.atoms.end()) fail(command, index, "atom serial not found: " + command[index].text);
            return static_cast<std::size_t>(found - state.structure.atoms.begin());
        }
    } catch (const std::invalid_argument&) {
    } catch (const std::out_of_range&) {
        fail(command, index, "atom serial out of range: " + command[index].text);
    }
    const auto found = state.render.namedSelections.find(lower(command[index].text));
    if (found == state.render.namedSelections.end()) fail(command, index, "selection not found: " + command[index].text);
    const auto indices = selectedIndices(found->second);
    if (indices.size() != 1) fail(command, index, "atom selection must contain exactly one atom");
    return indices.front();
}

std::string formattedMeasurement(double value, std::string_view unit) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(3) << value << ' ' << unit;
    return output.str();
}

void execute(const Command& command, CommandState state, CommandResult& result) {
    if (command.empty()) return;
    const std::string verb = lower(command[0].text);
    if (verb == "show") {
        if (command.size() < 2) fail(command, 0, "usage: show <representation|selection>");
        const std::string mode = lower(command[1].text);
        if (command.size() == 2 && mode == "atoms") state.render.mode = RenderMode::atoms;
        else if (command.size() == 2 && mode == "backbone") state.render.mode = RenderMode::backbone;
        else if (command.size() == 2 && (mode == "cartoon" || mode == "ca" || mode == "trace")) {
            state.render.mode = RenderMode::alphaTrace;
        } else if (command.size() == 2 && (mode == "sticks" || mode == "ball-stick" || mode == "ball-and-stick")) {
            state.render.mode = RenderMode::ballAndStick;
        } else if (command.size() == 2 && (mode == "residue" || mode == "residues")) {
            state.render.mode = RenderMode::residues;
        }
        else {
            ensureAtomState(state.render, state.structure.atoms.size());
            state.render.visibleAtoms = selectionFrom(command, 1, command.size(), state);
        }
    } else if (verb == "hide") {
        ensureAtomState(state.render, state.structure.atoms.size());
        const auto selected = selectionFrom(command, 1, command.size(), state);
        for (std::size_t index = 0; index < selected.size(); ++index) {
            if (selected[index]) state.render.visibleAtoms[index] = false;
        }
    } else if (verb == "select") {
        if (command.size() < 3) fail(command, command.size() - 1, "usage: select <name> <expression>");
        const std::string name = lower(command[1].text);
        if (name == "all" || name == "none") fail(command, 1, "reserved selection name: " + name);
        auto selected = selectionFrom(command, 2, command.size(), state);
        state.render.namedSelections[name] = selected;
        state.render.selectedAtoms = std::move(selected);
        const std::size_t residueCount = selectedResidueCount(state.render.selectedAtoms, state.structure);
        result.output.push_back(std::to_string(selectionSize(state.render.selectedAtoms)) + " atoms selected (" +
                    std::to_string(residueCount) + (residueCount == 1 ? " residue)" : " residues)"));
    } else if (verb == "delete") {
        requireArguments(command, 2, "delete <selection|all>");
        const std::string name = lower(command[1].text);
        if (name == "all") state.render.namedSelections.clear();
        else if (state.render.namedSelections.erase(name) == 0) fail(command, 1, "selection not found: " + name);
        state.render.selectedAtoms.assign(state.structure.atoms.size(), false);
    } else if (verb == "center") {
        const auto indices = selectedIndices(selectionFrom(command, 1, command.size(), state));
        if (indices.empty()) fail(command, 1, "selection is empty");
        Vec3 center;
        for (const auto index : indices) center += state.structure.atoms[index].position;
        state.camera.target = center / static_cast<double>(indices.size());
    } else if (verb == "inspect") {
        const auto indices = selectedIndices(selectionFrom(command, 1, command.size(), state));
        if (indices.size() != 1) fail(command, 1, "inspection selection must contain exactly one atom");
        state.render.inspectedAtom = indices.front();
    } else if (verb == "color") {
        if (command.size() < 2) fail(command, 0, "usage: color <mode> | color <selection> <color>");
        const std::string mode = lower(command[1].text);
        if (command.size() == 2 && mode == "element") state.render.colorMode = ColorMode::element;
        else if (command.size() == 2 && mode == "chain") state.render.colorMode = ColorMode::chain;
        else if (command.size() == 2 && mode == "residue") state.render.colorMode = ColorMode::residue;
        else if (command.size() == 2 && (mode == "bfactor" || mode == "b-factor")) {
            state.render.colorMode = ColorMode::bFactor;
        } else if (command.size() == 2 && mode == "secondary") state.render.colorMode = ColorMode::secondary;
        else {
            if (command.size() < 3) fail(command, 1, "expected a selection and color");
            ensureAtomState(state.render, state.structure.atoms.size());
            const auto selected = selectionFrom(command, 1, command.size() - 1, state);
            const int color = colorCode(command, command.size() - 1);
            for (std::size_t index = 0; index < selected.size(); ++index) {
                if (selected[index]) state.render.atomColors[index] = color;
            }
        }
    } else if (verb == "set") {
        requireArguments(command, 3, "set projection|renderer|charset|info|sequence|background <value>");
        const std::string setting = lower(command[1].text);
        const std::string value = lower(command[2].text);
        if (setting == "projection") {
            if (value == "orthographic" || value == "ortho") state.camera.projection = Projection::orthographic;
            else if (value == "perspective") state.camera.projection = Projection::perspective;
            else fail(command, 2, "unknown projection: " + command[2].text);
        } else if (setting == "renderer") {
            if (value == "cells") state.render.backend = RendererBackend::cells;
            else if (value == "braille") state.render.backend = RendererBackend::braille;
            else fail(command, 2, "unknown renderer: " + command[2].text);
        } else if (setting == "charset") {
            if (value == "ascii") state.render.charset = Charset::ascii;
            else if (value == "unicode") state.render.charset = Charset::unicode;
            else fail(command, 2, "unknown charset: " + command[2].text);
        } else if (setting == "info") state.render.showInfo = toggleValue(command, 2, state.render.showInfo);
        else if (setting == "sequence") state.render.showSequence = toggleValue(command, 2, state.render.showSequence);
        else if (setting == "background") {
            const int foreground = colorCode(command, 2);
            state.render.backgroundColor = foreground == 0 ? 0 : foreground + 10;
        } else if (setting == "frame-delay") {
            if (!state.trajectory) fail(command, 1, "trajectory controls are unavailable");
            const double delay = numericValue(command, 2);
            if (delay < 10.0 || delay > 10000.0 || std::floor(delay) != delay) {
                fail(command, 2, "frame delay must be an integer from 10 to 10000 ms");
            }
            state.trajectory->frameDelayMs = static_cast<int>(delay);
        } else if (setting == "trajectory-loop") {
            if (!state.trajectory) fail(command, 1, "trajectory controls are unavailable");
            state.trajectory->loop = toggleValue(command, 2, state.trajectory->loop);
        }
        else fail(command, 1, "unknown setting: " + command[1].text);
    } else if (verb == "spin") {
        requireArguments(command, 2, "spin on|off|toggle");
        state.autorotate = toggleValue(command, 1, state.autorotate);
        } else if (verb == "play") {
            if (!state.trajectory) fail(command, 0, "trajectory controls are unavailable");
            if (command.size() > 2) fail(command, 2, "usage: play [trajectory.xyz]");
            if (command.size() == 2) loadTrajectory(command[1].text, state.structure, *state.trajectory);
            else if (!state.trajectory->source) fail(command, 0, "no trajectory loaded");
            state.trajectory->playing = state.trajectory->source->frameCount() > 1;
            result.output.push_back(trajectoryStatus(*state.trajectory));
        } else if (verb == "pause") {
            requireArguments(command, 1, "pause");
            if (!state.trajectory || !state.trajectory->source) fail(command, 0, "no trajectory loaded");
            state.trajectory->playing = false;
            result.output.push_back(trajectoryStatus(*state.trajectory));
        } else if (verb == "resume") {
            requireArguments(command, 1, "resume");
            if (!state.trajectory || !state.trajectory->source) fail(command, 0, "no trajectory loaded");
            state.trajectory->playing = state.trajectory->source->frameCount() > 1;
            result.output.push_back(trajectoryStatus(*state.trajectory));
        } else if (verb == "step") {
            if (command.size() > 2) fail(command, 2, "usage: step [offset]");
            if (!state.trajectory || !state.trajectory->source) fail(command, 0, "no trajectory loaded");
            int offset = 1;
            if (command.size() == 2) {
                const double value = numericValue(command, 1);
                if (std::floor(value) != value || value < -100000 || value > 100000) {
                    fail(command, 1, "step offset must be an integer");
                }
                offset = static_cast<int>(value);
            }
            state.trajectory->playing = false;
            if (!stepTrajectory(state.structure, *state.trajectory, offset)) fail(command, 0, "trajectory boundary reached");
            result.output.push_back(trajectoryStatus(*state.trajectory));
        } else if (verb == "stop") {
            requireArguments(command, 1, "stop");
            if (!state.trajectory || !state.trajectory->source) fail(command, 0, "no trajectory loaded");
            stopTrajectory(state.structure, *state.trajectory);
            result.output.push_back(trajectoryStatus(*state.trajectory));
    } else if (verb == "zoom") {
        if (command.size() == 2) {
            try {
                state.camera.zoom = std::clamp(numericValue(command, 1), 0.15, 8.0);
            } catch (const std::runtime_error&) {
                const auto indices = selectedIndices(selectionFrom(command, 1, command.size(), state));
                if (indices.empty()) fail(command, 1, "selection is empty");
                Vec3 center;
                for (const auto index : indices) center += state.structure.atoms[index].position;
                center = center / static_cast<double>(indices.size());
                double radius = 0.0;
                for (const auto index : indices) {
                    radius = std::max(radius, length(state.structure.atoms[index].position - center));
                }
                state.camera.target = center;
                state.camera.zoom = std::clamp(state.structure.boundingRadius / std::max(radius, 0.5), 0.15, 8.0);
            }
        } else {
            const auto indices = selectedIndices(selectionFrom(command, 1, command.size(), state));
            if (indices.empty()) fail(command, 1, "selection is empty");
            Vec3 center;
            for (const auto index : indices) center += state.structure.atoms[index].position;
            center = center / static_cast<double>(indices.size());
            double radius = 0.0;
            for (const auto index : indices) radius = std::max(radius, length(state.structure.atoms[index].position - center));
            state.camera.target = center;
            state.camera.zoom = std::clamp(state.structure.boundingRadius / std::max(radius, 0.5), 0.15, 8.0);
        }
    } else if (verb == "rotate") {
        requireArguments(command, 3, "rotate yaw|pitch|roll <degrees>");
        const double degrees = numericValue(command, 2);
        const std::string axis = lower(command[1].text);
        if (axis == "yaw") state.camera.rotate(degrees, 0.0);
        else if (axis == "pitch") state.camera.rotate(0.0, degrees);
        else if (axis == "roll") state.camera.rotate(0.0, 0.0, degrees);
        else fail(command, 1, "unknown rotation axis: " + command[1].text);
    } else if (verb == "reset") {
        requireArguments(command, 1, "reset");
        state.camera.reset();
    } else if (verb == "count") {
        if (command.size() < 2) fail(command, 0, "usage: count atoms|residues|chains|<selection>");
        const std::string kind = lower(command[1].text);
        if (command.size() == 2 && kind == "atoms") result.output.push_back(std::to_string(state.structure.atoms.size()));
        else if (kind == "residues") result.output.push_back(std::to_string(state.structure.residues.size()));
        else if (kind == "chains") result.output.push_back(std::to_string(state.structure.chains.size()));
        else result.output.push_back(std::to_string(selectionSize(selectionFrom(command, 1, command.size(), state))));
    } else if (verb == "list") {
        requireArguments(command, 2, "list chains|selections");
        const std::string target = lower(command[1].text);
        if (target == "chains") {
            std::string chains;
            for (const auto& chain : state.structure.chains) {
                if (!chains.empty()) chains += ' ';
                chains += chain.id.empty() ? "-" : chain.id;
            }
            result.output.push_back(chains);
        } else if (target == "selections") {
            if (state.render.namedSelections.empty()) result.output.push_back("No named selections");
            else {
                std::vector<std::string> names;
                for (const auto& [name, selection] : state.render.namedSelections) {
                    names.push_back(name + " (" + std::to_string(selectionSize(selection)) + ")");
                }
                std::sort(names.begin(), names.end());
                std::string selections;
                for (const auto& name : names) {
                    if (!selections.empty()) selections += ", ";
                    selections += name;
                }
                result.output.push_back(selections);
            }
        } else {
            fail(command, 1, "unknown list target: " + command[1].text);
        }
    } else if (verb == "sequence") {
        if (command.size() > 2) fail(command, 2, "usage: sequence [chain]");
        const std::string chain = command.size() == 2 ? command[1].text : std::string{};
        const std::string resultSequence = sequence(state.structure, chain);
        if (resultSequence.empty() && !chain.empty()) fail(command, 1, "chain not found: " + chain);
        result.output.push_back(resultSequence);
    } else if (verb == "load") {
        requireArguments(command, 2, "load <structure.pdb|structure.cif>");
        state.structure = loadStructure(command[1].text);
        resetStructureState(state);
        result.output.push_back("Loaded " + state.structure.sourceName);
    } else if (verb == "fetch") {
        if (command.size() != 2 && command.size() != 4) fail(command, command.size() - 1, "usage: fetch <PDB-ID> [format pdb|cif]");
        bool cif = false;
        if (command.size() == 4) {
            if (lower(command[2].text) != "format") fail(command, 2, "expected format pdb|cif");
            const std::string format = lower(command[3].text);
            if (format != "pdb" && format != "cif") fail(command, 3, "expected format pdb|cif");
            cif = format == "cif";
        }
        const auto path = cif ? fetchMmcif(command[1].text) : fetchPdb(command[1].text);
        state.structure = loadStructure(path);
        resetStructureState(state);
        result.output.push_back("Fetched " + normalizePdbId(command[1].text));
    } else if (verb == "save") {
        requireArguments(command, 2, "save <snapshot.txt>");
        if (!state.saveFrame) fail(command, 0, "saving is unavailable in this context");
        state.saveFrame(command[1].text);
        result.output.push_back("Saved " + command[1].text);
    } else if (verb == "distance") {
        requireArguments(command, 3, "distance <atom-serial|selection> <atom-serial|selection>");
        const Vec3 first = state.structure.atoms[atomOperand(command, 1, state)].position;
        const Vec3 second = state.structure.atoms[atomOperand(command, 2, state)].position;
        result.output.push_back(formattedMeasurement(length(first - second), "A"));
    } else if (verb == "angle") {
        requireArguments(command, 4, "angle <atom1> <atom2> <atom3>");
        const Vec3 first = state.structure.atoms[atomOperand(command, 1, state)].position;
        const Vec3 vertex = state.structure.atoms[atomOperand(command, 2, state)].position;
        const Vec3 third = state.structure.atoms[atomOperand(command, 3, state)].position;
        if (length(first - vertex) == 0.0 || length(third - vertex) == 0.0) {
            fail(command, 0, "angle is undefined for coincident atoms");
        }
        const double cosine = std::clamp(dot(normalized(first - vertex), normalized(third - vertex)), -1.0, 1.0);
        result.output.push_back(formattedMeasurement(std::acos(cosine) * 180.0 / std::numbers::pi, "deg"));
    } else if (verb == "dihedral") {
        requireArguments(command, 5, "dihedral <atom1> <atom2> <atom3> <atom4>");
        const Vec3 first = state.structure.atoms[atomOperand(command, 1, state)].position;
        const Vec3 second = state.structure.atoms[atomOperand(command, 2, state)].position;
        const Vec3 third = state.structure.atoms[atomOperand(command, 3, state)].position;
        const Vec3 fourth = state.structure.atoms[atomOperand(command, 4, state)].position;
        const Vec3 axis = normalized(third - second);
        const Vec3 firstCross = cross(second - first, third - second);
        const Vec3 secondCross = cross(third - second, fourth - third);
        if (length(third - second) == 0.0 || length(firstCross) == 0.0 || length(secondCross) == 0.0) {
            fail(command, 0, "dihedral is undefined for coincident or collinear atoms");
        }
        const Vec3 firstNormal = normalized(firstCross);
        const Vec3 secondNormal = normalized(secondCross);
        const double angle = std::atan2(dot(cross(firstNormal, secondNormal), axis), dot(firstNormal, secondNormal));
        result.output.push_back(formattedMeasurement(angle * 180.0 / std::numbers::pi, "deg"));
    } else if (verb == "info") {
        requireArguments(command, 1, "info");
        std::ostringstream info;
        info << state.structure.sourceName << ": " << state.structure.atoms.size() << " atoms, "
             << state.structure.residues.size() << " residues, " << state.structure.chains.size() << " chains";
        result.output.push_back(info.str());
    } else if (verb == "help") {
        if (command.size() > 2) fail(command, 2, "usage: help [command]");
        result.output.push_back(commandHelp(command.size() == 2 ? command[1].text : std::string_view{}));
    } else if (verb == "quit" || verb == "exit") {
        requireArguments(command, 1, "quit");
        state.running = false;
    } else {
        fail(command, 0, "unknown command: " + command[0].text);
    }
}

} // namespace

CommandResult executeCommands(std::string_view source, CommandState state) {
    CommandResult result;
    try {
        for (const auto& command : tokenize(source)) execute(command, state, result);
    } catch (const std::runtime_error& error) {
        result.success = false;
        result.output.push_back(error.what());
    }
    return result;
}

std::string commandHelp(std::string_view command) {
    const std::string topic = lower(std::string(command));
    if (topic == "show") return "show <representation> | show <selection> (only matching atoms)";
    if (topic == "hide") return "hide <selection>";
    if (topic == "select") return "select <name> <expression>";
    if (topic == "delete") return "delete <selection|all>";
    if (topic == "set") {
        return "set projection|renderer|charset|info|sequence|background|frame-delay|trajectory-loop <value>";
    }
    if (topic == "color") return "color <mode> | color <selection> <color>";
    if (topic == "rotate") return "rotate yaw|pitch|roll <degrees>";
    if (topic == "count") return "count atoms|residues|chains|<selection>";
    if (topic == "list") return "list chains|selections";
    if (topic == "center") return "center <selection>";
    if (topic == "zoom") return "zoom <factor|selection>";
    if (topic == "inspect") return "inspect <single-atom-selection>";
    if (topic == "load") return "load <structure.pdb|structure.cif>";
    if (topic == "fetch") return "fetch <PDB-ID> [format pdb|cif]";
    if (topic == "save") return "save <snapshot.txt>";
    if (topic == "distance") return "distance <atom-serial|selection> <atom-serial|selection>";
    if (topic == "angle") return "angle <atom1> <atom2> <atom3>";
    if (topic == "dihedral") return "dihedral <atom1> <atom2> <atom3> <atom4>";
    if (topic == "play") return "play [trajectory.xyz]";
    if (topic == "pause" || topic == "resume" || topic == "stop") return topic;
    if (topic == "step") return "step [signed-offset]";
    if (topic == "selection" || topic == "selections") {
        return "Selections: all, none, chain ID, resi N[-M], resn NAME[+NAME], name ATOM, element SYMBOL, "
               "serial N, hetero, solvent, named selections, and/or/not, and parentheses";
    }
    if (!topic.empty()) return "No help for command: " + std::string(command);
    return "Commands: show, hide, select, delete, color, set, spin, play, pause, resume, step, stop, zoom, center, "
           "inspect, rotate, reset, count, list, sequence, distance, angle, dihedral, info, load, fetch, save, help, quit";
}

std::vector<std::string> commandCompletions(std::string_view input, CommandState state) {
    const auto separator = input.find_last_of(';');
    std::string segment(input.substr(separator == std::string_view::npos ? 0 : separator + 1));
    const bool trailingSpace = !segment.empty() && std::isspace(static_cast<unsigned char>(segment.back()));
    std::istringstream stream(segment);
    std::vector<std::string> words;
    for (std::string word; stream >> word;) words.push_back(lower(word));
    if (trailingSpace || words.empty()) words.emplace_back();
    const std::string prefix = words.back();
    const std::size_t position = words.size() - 1;
    const std::string verb = position == 0 ? std::string{} : words.front();
    std::vector<std::string> candidates;
    const auto add = [&](std::initializer_list<std::string_view> values) {
        for (const auto value : values) candidates.emplace_back(value);
    };
    const auto addSelections = [&] {
        add({"all", "none", "chain", "resi", "resn", "name", "element", "serial", "hetero", "solvent",
             "and", "or", "not"});
        for (const auto& [name, selection] : state.render.namedSelections) {
            static_cast<void>(selection);
            candidates.push_back(name);
        }
    };
    if (position == 0) {
        add({"show", "hide", "select", "delete", "color", "set", "spin", "play", "pause", "resume", "step", "stop", "zoom", "center", "inspect",
             "rotate", "reset", "count", "list", "sequence", "distance", "angle", "dihedral", "info",
             "load", "fetch", "save", "help", "quit"});
    } else if (verb == "set" && position == 1) {
        add({"projection", "renderer", "charset", "info", "sequence", "background", "frame-delay", "trajectory-loop"});
    } else if (verb == "set" && words.size() >= 3) {
        if (words[1] == "projection") add({"perspective", "orthographic"});
        else if (words[1] == "renderer") add({"cells", "braille"});
        else if (words[1] == "charset") add({"ascii", "unicode"});
        else if (words[1] == "info" || words[1] == "sequence" || words[1] == "trajectory-loop") add({"on", "off", "toggle"});
        else if (words[1] == "background") add({"default", "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white"});
    } else if (verb == "rotate" && position == 1) {
        add({"yaw", "pitch", "roll"});
    } else if (verb == "spin") {
        add({"on", "off", "toggle"});
    } else if (verb == "list") {
        add({"chains", "selections"});
    } else if (verb == "show" && position == 1) {
        add({"atoms", "backbone", "cartoon", "sticks", "residues"});
        addSelections();
    } else if (verb == "color" && position == 1) {
        add({"element", "chain", "residue", "bfactor", "secondary"});
        addSelections();
    } else if (verb == "color" && position > 1) {
        add({"default", "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white"});
        addSelections();
    } else if (verb == "help") {
        add({"show", "hide", "select", "delete", "color", "set", "zoom", "center", "inspect", "rotate",
             "count", "load", "fetch", "save", "distance", "angle", "dihedral", "play", "pause", "resume",
             "step", "stop"});
    } else if (verb == "fetch" && position == 2) {
        add({"format"});
    } else if (verb == "fetch" && position == 3) {
        add({"pdb", "cif"});
    } else {
        addSelections();
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    std::erase_if(candidates, [&](const std::string& candidate) { return !candidate.starts_with(prefix); });
    return candidates;
}

} // namespace ascii_protein