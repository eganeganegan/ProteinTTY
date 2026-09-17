#include "ascii_protein/cli.hpp"
#include "ascii_protein/command.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/structure_loader.hpp"
#include "ascii_protein/terminal.hpp"
#include "ascii_protein/trajectory.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <csignal>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

volatile std::sig_atomic_t interrupted = 0;

void handleInterrupt(int) { interrupted = 1; }

void saveSnapshot(const ascii_protein::Framebuffer& framebuffer, int& number) {
    std::string filename;
    do {
        const std::string padding = number < 10 ? "00" : number < 100 ? "0" : "";
        filename = "snapshot-" + padding + std::to_string(number) + ".txt";
        ++number;
    } while (std::ifstream(filename));
    std::ofstream output(filename);
    output << framebuffer.toPlainText() << '\n';
}

std::string joinOutput(const ascii_protein::CommandResult& result) {
    std::string output;
    for (const auto& line : result.output) {
        if (!output.empty()) output += " | ";
        output += line;
    }
    return output;
}

struct ConsoleState {
    bool active{};
    std::string input;
    std::string message;
    std::vector<std::string> history;
    std::size_t historyIndex{};
};

bool handleConsoleKey(const ascii_protein::KeyEvent& event, ConsoleState& console,
                      ascii_protein::CommandState commandState) {
    using ascii_protein::Key;
    if (event.key == Key::escape) {
        console.active = false;
        console.input.clear();
        return true;
    }
    if (event.key == Key::up && !console.history.empty()) {
        if (console.historyIndex > 0) --console.historyIndex;
        console.input = console.history[console.historyIndex];
        return true;
    }
    if (event.key == Key::down && !console.history.empty()) {
        if (console.historyIndex + 1 < console.history.size()) {
            console.input = console.history[++console.historyIndex];
        } else {
            console.historyIndex = console.history.size();
            console.input.clear();
        }
        return true;
    }
    if (event.key != Key::character) return true;
    if (event.character == '\r' || event.character == '\n') {
        if (!console.input.empty()) {
            const auto result = ascii_protein::executeCommands(console.input, commandState);
            console.message = result.success ? joinOutput(result) : "Error: " + joinOutput(result);
            if (console.message.empty()) console.message = "OK";
            if (console.history.empty() || console.history.back() != console.input) console.history.push_back(console.input);
        }
        console.historyIndex = console.history.size();
        console.input.clear();
        console.active = false;
        return true;
    }
    if (event.character == 127 || event.character == '\b') {
        if (!console.input.empty()) console.input.pop_back();
        return true;
    }
    if (event.character == '\t') {
        const auto completions = ascii_protein::commandCompletions(console.input, commandState);
        if (completions.empty()) return true;
        std::string common = completions.front();
        for (const auto& completion : completions) {
            while (!common.empty() && !completion.starts_with(common)) common.pop_back();
        }
        const auto tokenStart = console.input.find_last_of(" \t;");
        const std::size_t start = tokenStart == std::string::npos ? 0 : tokenStart + 1;
        if (common.size() > console.input.size() - start) console.input.replace(start, std::string::npos, common);
        if (completions.size() == 1) console.input += ' ';
        else {
            console.message.clear();
            for (const auto& completion : completions) {
                if (!console.message.empty()) console.message += ' ';
                console.message += completion;
            }
        }
        return true;
    }
    if (std::isprint(static_cast<unsigned char>(event.character))) console.input.push_back(event.character);
    return true;
}

ascii_protein::CommandResult runInitialCommands(const ascii_protein::CliOptions& options,
                                                ascii_protein::CommandState state) {
    ascii_protein::CommandResult combined;
    const auto run = [&](const std::string& command) {
        const auto result = ascii_protein::executeCommands(command, state);
        combined.output.insert(combined.output.end(), result.output.begin(), result.output.end());
        if (!result.success) combined.success = false;
        return result.success;
    };
    for (const auto& path : options.scripts) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot open command script: " + path.string());
        std::string line;
        std::size_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            const auto first = line.find_first_not_of(" \t");
            if (first == std::string::npos || line[first] == '#') continue;
            if (!run(line)) throw std::runtime_error(path.string() + ":" + std::to_string(lineNumber) + ": " +
                                                     combined.output.back());
        }
    }
    for (const auto& command : options.commands) {
        if (!run(command)) throw std::runtime_error("command: " + combined.output.back());
    }
    return combined;
}

void applyKey(const ascii_protein::KeyEvent& event, ascii_protein::Camera& camera,
              ascii_protein::RenderSettings& settings, std::size_t residueCount, bool& running, bool& autorotate,
              const ascii_protein::Framebuffer& framebuffer, int& snapshotNumber,
              std::optional<std::size_t> centerAtom, ascii_protein::Structure& structure,
              ascii_protein::TrajectoryState& trajectory) {
    using ascii_protein::ColorMode;
    using ascii_protein::Key;
    using ascii_protein::Projection;
    using ascii_protein::RenderMode;
    if (event.key == Key::left) camera.rotate(-5.0, 0.0);
    if (event.key == Key::right) camera.rotate(5.0, 0.0);
    if (event.key == Key::up) camera.rotate(0.0, -5.0);
    if (event.key == Key::down) camera.rotate(0.0, 5.0);
    if (event.key == Key::escape) settings.showHelp = false;
    if (event.key != Key::character) return;

    switch (event.character) {
    case 'q': case 'Q': running = false; break;
    case 'a': case 'A': camera.rotate(-5.0, 0.0); break;
    case 'd': case 'D': camera.rotate(5.0, 0.0); break;
    case 'w': case 'W': camera.rotate(0.0, -5.0); break;
    case 's': camera.rotate(0.0, 5.0); break;
    case 'z': case 'Z': camera.rotate(0.0, 0.0, -5.0); break;
    case 'e': case 'E': camera.rotate(0.0, 0.0, 5.0); break;
    case '+': case '=': camera.zoom = std::min(8.0, camera.zoom * 1.12); break;
    case '-': case '_': camera.zoom = std::max(0.15, camera.zoom / 1.12); break;
    case 'r': case 'R': camera.reset(); break;
    case 'p': case 'P':
        camera.projection = camera.projection == Projection::perspective ? Projection::orthographic
                                                                         : Projection::perspective;
        break;
    case '1': settings.mode = RenderMode::atoms; break;
    case '2': settings.mode = RenderMode::backbone; break;
    case '3': settings.mode = RenderMode::alphaTrace; break;
    case '4': settings.mode = RenderMode::ballAndStick; break;
    case '5': settings.mode = RenderMode::residues; break;
    case 'c': case 'C':
        settings.colorMode = static_cast<ColorMode>((static_cast<int>(settings.colorMode) + 1) % 5);
        break;
    case 'i': case 'I': settings.showInfo = !settings.showInfo; break;
    case 'S': settings.showSequence = !settings.showSequence; break;
    case 'h': case 'H': settings.showHelp = !settings.showHelp; break;
    case 'o': case 'O':
        if (settings.inspectedAtom) settings.inspectedAtom.reset();
        else settings.inspectedAtom = centerAtom;
        break;
    case '[':
        if (residueCount > 0) {
            settings.selectedResidue = !settings.selectedResidue || *settings.selectedResidue == 0
                                           ? residueCount - 1
                                           : *settings.selectedResidue - 1;
        }
        break;
    case ']':
        if (residueCount > 0) {
            settings.selectedResidue = settings.selectedResidue ? (*settings.selectedResidue + 1) % residueCount : 0;
        }
        break;
    case ' ': autorotate = !autorotate; break;
    case '`':
        if (trajectory.source) trajectory.playing = !trajectory.playing;
        break;
    case ',':
        if (trajectory.source) {
            trajectory.playing = false;
            ascii_protein::stepTrajectory(structure, trajectory, -1);
        }
        break;
    case '.':
        if (trajectory.source) {
            trajectory.playing = false;
            ascii_protein::stepTrajectory(structure, trajectory, 1);
        }
        break;
    case 'x': case 'X': saveSnapshot(framebuffer, snapshotNumber); break;
    default: break;
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        auto options = ascii_protein::parseCli(argc, argv);
        if (options.help) {
            std::cout << ascii_protein::helpText();
            return 0;
        }
        if (options.version) {
            std::cout << "ascii-protein " << ASCII_PROTEIN_VERSION << '\n';
            return 0;
        }

        const auto input = options.pdbId ? (options.fetchCif ? ascii_protein::fetchMmcif(*options.pdbId)
                                    : ascii_protein::fetchPdb(*options.pdbId))
                         : *options.input;
        auto structure = ascii_protein::loadStructure(input);
        for (const auto& warning : structure.warnings) std::cerr << "warning: " << warning << '\n';
        auto renderer = std::make_unique<ascii_protein::Renderer>(structure);
        ascii_protein::Camera camera;
        camera.projection = options.projection;
        bool running = true;
        bool autorotate = options.spin;
        ascii_protein::TrajectoryState trajectory;
        trajectory.frameDelayMs = options.frameDelayMs;
        trajectory.loop = options.trajectoryLoop;
        if (options.trajectory) ascii_protein::loadTrajectory(*options.trajectory, structure, trajectory);
        int currentWidth = options.width;
        int currentHeight = options.height;
        ascii_protein::CommandState commandState{
            camera, options.render, structure, autorotate, running,
            [&] { renderer = std::make_unique<ascii_protein::Renderer>(structure); },
            [&](const std::filesystem::path& path) {
                std::ofstream output(path);
                if (!output) throw std::runtime_error("cannot save snapshot: " + path.string());
                output << renderer->render(camera, options.render, currentWidth, currentHeight).toPlainText() << '\n';
            },
            &trajectory};
        const auto initialCommands = runInitialCommands(options, commandState);
        if (!running) {
            for (const auto& line : initialCommands.output) std::cout << line << '\n';
            return 0;
        }

        if (options.snapshot || !ascii_protein::Terminal::interactive()) {
            options.render.trajectoryInfo = ascii_protein::trajectoryStatus(trajectory);
            for (const auto& line : initialCommands.output) std::cerr << line << '\n';
            const auto& frame = renderer->render(camera, options.render, options.width, options.height);
            std::cout << frame.toString(options.render.color && !options.snapshot, options.render.backgroundColor) << '\n';
            return 0;
        }

        std::signal(SIGINT, handleInterrupt);
        std::signal(SIGTERM, handleInterrupt);
        ascii_protein::Terminal terminal;
        int snapshotNumber = 1;
        ConsoleState console;
        console.message = joinOutput(initialCommands);
        auto previous = std::chrono::steady_clock::now();
        double trajectoryElapsedMs = 0.0;
        while (running && !interrupted) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed = std::chrono::duration<double>(now - previous).count();
            previous = now;
            if (autorotate) camera.rotate(elapsed * 24.0, 0.0);
            if (trajectory.playing && trajectory.source) {
                trajectoryElapsedMs += elapsed * 1000.0;
                if (trajectoryElapsedMs >= trajectory.frameDelayMs) {
                    trajectoryElapsedMs = 0.0;
                    if (!ascii_protein::stepTrajectory(structure, trajectory)) {
                        trajectory.playing = false;
                        console.message = "Trajectory complete";
                    }
                }
            } else {
                trajectoryElapsedMs = 0.0;
            }
            const auto dimensions = terminal.size();
            currentWidth = dimensions.width;
            currentHeight = dimensions.height;
            options.render.commandActive = console.active;
            options.render.commandLine = console.input;
            options.render.commandMessage = console.message;
            options.render.trajectoryInfo = ascii_protein::trajectoryStatus(trajectory);
            const auto& frame = renderer->render(camera, options.render, dimensions.width, dimensions.height);
            if (dimensions.width < 20 || dimensions.height < 8) {
                terminal.present("Terminal too small (minimum 20x8). Resize or press q.");
            } else {
                terminal.present(frame.toString(options.render.color, options.render.backgroundColor));
            }
            const int timeout = trajectory.playing ? std::min(trajectory.frameDelayMs, 100) : autorotate ? 33 : 100;
            const auto event = terminal.readKey(timeout);
            if (options.spin && event.key != ascii_protein::Key::none) break;
            if (console.active) {
                handleConsoleKey(event, console, commandState);
                continue;
            }
            if (event.key == ascii_protein::Key::character && event.character == ':') {
                console.active = true;
                console.input.clear();
                console.message.clear();
                console.historyIndex = console.history.size();
                continue;
            }
            const auto centerAtom = renderer->nearestAtom(camera, options.render, dimensions.width, dimensions.height);
            applyKey(event, camera, options.render, structure.residues.size(), running, autorotate, frame,
                     snapshotNumber, centerAtom, structure, trajectory);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ascii-protein: " << error.what() << '\n';
        return 1;
    }
}