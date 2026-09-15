#include "ascii_protein/cli.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/rcsb.hpp"
#include "ascii_protein/terminal.hpp"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

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

void applyKey(const ascii_protein::KeyEvent& event, ascii_protein::Camera& camera,
              ascii_protein::RenderSettings& settings, std::size_t residueCount, bool& running, bool& autorotate,
              const ascii_protein::Framebuffer& framebuffer, int& snapshotNumber,
              std::optional<std::size_t> centerAtom) {
    using ascii_protein::ColorMode;
    using ascii_protein::Key;
    using ascii_protein::Projection;
    using ascii_protein::RenderMode;
    if (event.key == Key::left) camera.yaw -= 5.0;
    if (event.key == Key::right) camera.yaw += 5.0;
    if (event.key == Key::up) camera.pitch -= 5.0;
    if (event.key == Key::down) camera.pitch += 5.0;
    if (event.key == Key::escape) settings.showHelp = false;
    if (event.key != Key::character) return;

    switch (event.character) {
    case 'q': case 'Q': running = false; break;
    case 'a': case 'A': camera.yaw -= 5.0; break;
    case 'd': case 'D': camera.yaw += 5.0; break;
    case 'w': case 'W': camera.pitch -= 5.0; break;
    case 's': camera.pitch += 5.0; break;
    case 'z': case 'Z': camera.roll -= 5.0; break;
    case 'e': case 'E': camera.roll += 5.0; break;
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
        if (residueCount > 0) settings.selectedResidue = settings.selectedResidue == 0 ? residueCount - 1
                                                                                       : settings.selectedResidue - 1;
        break;
    case ']':
        if (residueCount > 0) settings.selectedResidue = (settings.selectedResidue + 1) % residueCount;
        break;
    case ' ': autorotate = !autorotate; break;
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

        const auto input = options.pdbId ? ascii_protein::fetchPdb(*options.pdbId) : *options.input;
        const auto structure = ascii_protein::PdbParser::parseFile(input);
        for (const auto& warning : structure.warnings) std::cerr << "warning: " << warning << '\n';
        ascii_protein::Renderer renderer(structure);
        ascii_protein::Camera camera;
        camera.projection = options.projection;

        if (options.snapshot || !ascii_protein::Terminal::interactive()) {
            const auto& frame = renderer.render(camera, options.render, options.width, options.height);
            std::cout << frame.toString(options.render.color && !options.snapshot) << '\n';
            return 0;
        }

        std::signal(SIGINT, handleInterrupt);
        std::signal(SIGTERM, handleInterrupt);
        ascii_protein::Terminal terminal;
        bool running = true;
        bool autorotate = options.spin;
        int snapshotNumber = 1;
        auto previous = std::chrono::steady_clock::now();
        while (running && !interrupted) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed = std::chrono::duration<double>(now - previous).count();
            previous = now;
            if (autorotate) camera.yaw += elapsed * 24.0;
            const auto dimensions = terminal.size();
            const auto& frame = renderer.render(camera, options.render, dimensions.width, dimensions.height);
            if (dimensions.width < 20 || dimensions.height < 8) {
                terminal.present("Terminal too small (minimum 20x8). Resize or press q.");
            } else {
                terminal.present(frame.toString(options.render.color));
            }
            const auto event = terminal.readKey(autorotate ? 33 : 100);
            if (options.spin && event.key != ascii_protein::Key::none) break;
            const auto centerAtom = renderer.nearestAtom(camera, options.render, dimensions.width, dimensions.height);
            applyKey(event, camera, options.render, structure.residues.size(), running, autorotate, frame,
                     snapshotNumber, centerAtom);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ascii-protein: " << error.what() << '\n';
        return 1;
    }
}