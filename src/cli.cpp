#include "ascii_protein/cli.hpp"

#include <stdexcept>
#include <string_view>

namespace ascii_protein {
namespace {

std::string requireValue(int& index, int argc, char** argv, std::string_view option) {
    if (++index >= argc) throw std::runtime_error("missing value for " + std::string(option));
    return argv[index];
}

RenderMode parseMode(const std::string& value) {
    if (value == "atoms") return RenderMode::atoms;
    if (value == "backbone") return RenderMode::backbone;
    if (value == "ca") return RenderMode::alphaTrace;
    if (value == "ball-stick" || value == "ball-and-stick") return RenderMode::ballAndStick;
    if (value == "residue" || value == "residues") return RenderMode::residues;
    throw std::runtime_error("unknown render mode: " + value);
}

ColorMode parseColor(const std::string& value) {
    if (value == "element") return ColorMode::element;
    if (value == "chain") return ColorMode::chain;
    if (value == "residue") return ColorMode::residue;
    if (value == "bfactor" || value == "b-factor") return ColorMode::bFactor;
    if (value == "secondary") return ColorMode::secondary;
    throw std::runtime_error("unknown color mode: " + value);
}

int parseDimension(const std::string& value, std::string_view option) {
    try {
        const int result = std::stoi(value);
        if (result < 10 || result > 1000) throw std::runtime_error("");
        return result;
    } catch (...) {
        throw std::runtime_error("invalid value for " + std::string(option) + ": " + value);
    }
}

int parseFrameDelay(const std::string& value) {
    try {
        std::size_t used = 0;
        const int result = std::stoi(value, &used);
        if (used != value.size() || result < 10 || result > 10000) throw std::runtime_error("");
        return result;
    } catch (...) {
        throw std::runtime_error("invalid frame delay: " + value + " (expected 10..10000 ms)");
    }
}

} // namespace

CliOptions parseCli(int argc, char** argv) {
    CliOptions options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") options.help = true;
        else if (argument == "--version") options.version = true;
        else if (argument == "--command") options.commands.push_back(requireValue(index, argc, argv, argument));
        else if (argument == "--script") options.scripts.emplace_back(requireValue(index, argc, argv, argument));
            else if (argument == "--trajectory") options.trajectory = requireValue(index, argc, argv, argument);
            else if (argument == "--frame-delay") options.frameDelayMs = parseFrameDelay(requireValue(index, argc, argv, argument));
            else if (argument == "--loop") options.trajectoryLoop = true;
        else if (argument == "--pdb") options.pdbId = requireValue(index, argc, argv, argument);
        else if (argument == "--format") {
            const auto value = requireValue(index, argc, argv, argument);
            if (value != "pdb" && value != "cif") throw std::runtime_error("unknown structure format: " + value);
            options.fetchCif = value == "cif";
        }
        else if (argument == "--mode") options.render.mode = parseMode(requireValue(index, argc, argv, argument));
        else if (argument == "--color") options.render.colorMode = parseColor(requireValue(index, argc, argv, argument));
        else if (argument == "--charset") {
            const auto value = requireValue(index, argc, argv, argument);
            if (value != "ascii" && value != "unicode") throw std::runtime_error("unknown charset: " + value);
            options.render.charset = value == "ascii" ? Charset::ascii : Charset::unicode;
        } else if (argument == "--renderer") {
            const auto value = requireValue(index, argc, argv, argument);
            if (value != "cells" && value != "braille") throw std::runtime_error("unknown renderer: " + value);
            options.render.backend = value == "braille" ? RendererBackend::braille : RendererBackend::cells;
        } else if (argument == "--projection") {
            const auto value = requireValue(index, argc, argv, argument);
            if (value != "ortho" && value != "perspective") throw std::runtime_error("unknown projection: " + value);
            options.projection = value == "ortho" ? Projection::orthographic : Projection::perspective;
        } else if (argument == "--no-color") options.render.color = false;
        else if (argument == "--spin") {
            options.spin = true;
            options.render.showInfo = false;
        } else if (argument == "--snapshot") options.snapshot = true;
        else if (argument == "--width") options.width = parseDimension(requireValue(index, argc, argv, argument), argument);
        else if (argument == "--height") options.height = parseDimension(requireValue(index, argc, argv, argument), argument);
        else if (!argument.empty() && argument.front() == '-') throw std::runtime_error("unknown option: " + argument);
        else if (options.input) throw std::runtime_error("only one PDB file may be specified");
        else options.input = argument;
    }
    if (options.input && options.pdbId) throw std::runtime_error("use either a PDB file or --pdb, not both");
    if (options.fetchCif && !options.pdbId) throw std::runtime_error("--format is only used with --pdb");
    if (!options.help && !options.version && !options.input && !options.pdbId) {
        throw std::runtime_error("missing PDB file or --pdb ID (try --help)");
    }
    return options;
}

std::string helpText() {
    return
        "ProteinTTY - interactive terminal protein viewer\n\n"
        "Usage: ascii-protein <structure.pdb|structure.cif> [options]\n\n"
        "Options:\n"
        "  --command TEXT   Execute a viewer command (repeatable)\n"
        "  --script FILE    Execute commands from a .ptty script\n"
            "  --trajectory FILE  Play a multi-frame XYZ trajectory\n"
            "  --frame-delay MS   Trajectory frame delay, 10..10000 (default 100)\n"
            "  --loop             Loop trajectory playback\n"
        "  --pdb ID         Download and cache a structure from RCSB\n"
        "  --format pdb|cif Remote structure format (default pdb)\n"
        "  --mode atoms|backbone|ca|ball-stick|residue\n"
        "  --color element|chain|residue|bfactor|secondary\n"
        "  --charset ascii|unicode\n"
        "  --renderer cells|braille\n"
        "  --projection perspective|ortho\n"
        "  --no-color       Disable ANSI colors\n"
        "  --spin           Autorotate without the information panel\n"
        "  --snapshot       Render one plain-text frame and exit\n"
        "  --width N        Static snapshot width (default 100)\n"
        "  --height N       Static snapshot height (default 32)\n"
        "  --version        Print version\n"
        "  -h, --help       Print this help\n";
}

} // namespace ascii_protein