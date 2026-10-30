#include "ascii_protein/renderer.hpp"
#include <chrono>
#include <iostream>
#include <memory>
using namespace ascii_protein;
using Clock = std::chrono::steady_clock;
int main() {
    Structure structure;
    for (int i = 0; i < 100000; ++i) {
        Atom atom;
        atom.serial = i + 1;
        atom.name = "CA";
        atom.element = "C";
        atom.residueName = "GLY";
        atom.chain = std::to_string(i / 10000);
        atom.residueNumber = i % 10000;
        atom.position = {3. * (i % 100), 3. * ((i / 100) % 100), 3. * (i / 10000)};
        structure.atoms.push_back(std::move(atom));
    }
    auto start = Clock::now();
    structure.finalize();
    auto finalized = Clock::now();
    Renderer renderer(structure);
    auto indexed = Clock::now();
    Camera camera;
    RenderSettings settings;
    settings.showInfo = false;
    settings.mode = RenderMode::ballAndStick;
    constexpr int frames = 10;
    for (int i = 0; i < frames; ++i) {
        camera.rotate(3, 1);
        renderer.render(camera, settings, 120, 36);
    }
    auto rendered = Clock::now();
    auto map = std::make_shared<DensityMap>();
    map->dimensions = {256, 256, 256};
    map->steps = {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    map->values.resize(256 * 256 * 256, 1);
    settings.densityMap = map;
    settings.mapSlice = 128;
    auto mapped = Clock::now();
    for (int i = 0; i < frames; ++i) {
        camera.rotate(3, 1);
        renderer.render(camera, settings, 120, 36);
    }
    auto done = Clock::now();
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    std::cout << "100000 atoms, 120x36 Braille, " << frames << " frames\n"
              << "finalize_ms " << ms(start, finalized) << "\nbond_index_ms " << ms(finalized, indexed)
              << "\nstructure_frame_ms " << ms(indexed, rendered) / frames
              << "\nstructure_plus_256_cube_slice_frame_ms " << ms(mapped, done) / frames << '\n';
}
