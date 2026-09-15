#include "ascii_protein/command.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
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

void execute(const Command& command, CommandState state, CommandResult& result) {
    if (command.empty()) return;
    const std::string verb = lower(command[0].text);
    if (verb == "show") {
        requireArguments(command, 2, "show atoms|backbone|cartoon|sticks|residues");
        const std::string mode = lower(command[1].text);
        if (mode == "atoms") state.render.mode = RenderMode::atoms;
        else if (mode == "backbone") state.render.mode = RenderMode::backbone;
        else if (mode == "cartoon" || mode == "ca" || mode == "trace") state.render.mode = RenderMode::alphaTrace;
        else if (mode == "sticks" || mode == "ball-stick" || mode == "ball-and-stick") {
            state.render.mode = RenderMode::ballAndStick;
        } else if (mode == "residue" || mode == "residues") state.render.mode = RenderMode::residues;
        else fail(command, 1, "unknown representation: " + command[1].text);
    } else if (verb == "color") {
        requireArguments(command, 2, "color element|chain|residue|bfactor|secondary");
        const std::string mode = lower(command[1].text);
        if (mode == "element") state.render.colorMode = ColorMode::element;
        else if (mode == "chain") state.render.colorMode = ColorMode::chain;
        else if (mode == "residue") state.render.colorMode = ColorMode::residue;
        else if (mode == "bfactor" || mode == "b-factor") state.render.colorMode = ColorMode::bFactor;
        else if (mode == "secondary") state.render.colorMode = ColorMode::secondary;
        else fail(command, 1, "unknown color mode: " + command[1].text);
    } else if (verb == "set") {
        requireArguments(command, 3, "set projection|renderer|charset|info|sequence <value>");
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
        else fail(command, 1, "unknown setting: " + command[1].text);
    } else if (verb == "spin") {
        requireArguments(command, 2, "spin on|off|toggle");
        state.autorotate = toggleValue(command, 1, state.autorotate);
    } else if (verb == "zoom") {
        requireArguments(command, 2, "zoom <factor>");
        state.camera.zoom = std::clamp(numericValue(command, 1), 0.15, 8.0);
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
        requireArguments(command, 2, "count atoms|residues|chains");
        const std::string kind = lower(command[1].text);
        if (kind == "atoms") result.output.push_back(std::to_string(state.structure.atoms.size()));
        else if (kind == "residues") result.output.push_back(std::to_string(state.structure.residues.size()));
        else if (kind == "chains") result.output.push_back(std::to_string(state.structure.chains.size()));
        else fail(command, 1, "unknown count target: " + command[1].text);
    } else if (verb == "list") {
        requireArguments(command, 2, "list chains");
        if (lower(command[1].text) != "chains") fail(command, 1, "unknown list target: " + command[1].text);
        std::string chains;
        for (const auto& chain : state.structure.chains) {
            if (!chains.empty()) chains += ' ';
            chains += chain.id.empty() ? "-" : chain.id;
        }
        result.output.push_back(chains);
    } else if (verb == "sequence") {
        if (command.size() > 2) fail(command, 2, "usage: sequence [chain]");
        const std::string chain = command.size() == 2 ? command[1].text : std::string{};
        const std::string resultSequence = sequence(state.structure, chain);
        if (resultSequence.empty() && !chain.empty()) fail(command, 1, "chain not found: " + chain);
        result.output.push_back(resultSequence);
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
    if (topic == "show") return "show atoms|backbone|cartoon|sticks|residues";
    if (topic == "set") return "set projection|renderer|charset|info|sequence <value>";
    if (topic == "color") return "color element|chain|residue|bfactor|secondary";
    if (topic == "rotate") return "rotate yaw|pitch|roll <degrees>";
    if (topic == "count") return "count atoms|residues|chains";
    if (!topic.empty()) return "No help for command: " + std::string(command);
    return "Commands: show, color, set, spin, zoom, rotate, reset, count, list, sequence, info, help, quit";
}

} // namespace ascii_protein