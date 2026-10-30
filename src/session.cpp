#include "ascii_protein/session.hpp"
#include "ascii_protein/structure_loader.hpp"
#include "ascii_protein/molecular.hpp"
#include "ascii_protein/trajectory.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <map>
namespace ascii_protein {
namespace {
std::string bits(const AtomMask& mask) {
    std::string s;
    for (bool b : mask)
        s += b ? '1' : '0';
    return s;
}
AtomMask mask(const std::string& text, std::size_t count) {
    if (!text.empty() && text.size() != count) throw std::runtime_error("session atom count mismatch");
    AtomMask m;
    for (char c : text) {
        if (c != '0' && c != '1') throw std::runtime_error("invalid session mask");
        m.push_back(c == '1');
    }
    return m;
}
void finite(double x) {
    if (!std::isfinite(x)) throw std::runtime_error("invalid session number");
}
void checked(bool ok) {
    if (!ok) throw std::runtime_error("invalid or truncated session");
}
} // namespace
void saveSession(const std::filesystem::path& path, CommandState state) {
    if (state.structure.sourcePath.empty())
        throw std::runtime_error("session requires a structure loaded from a file");
    std::ostringstream out;
    out << std::setprecision(17) << "ProteinTTY-session 1\n";
    out << std::quoted(state.structure.sourcePath.string()) << ' ' << state.structure.model << ' '
        << std::quoted(state.structure.assemblyApplied) << '\n';
    const auto& c = state.camera;
    const auto& r = state.render;
    out << c.yaw << ' ' << c.pitch << ' ' << c.roll << ' ' << c.zoom << ' ' << c.pan.x << ' ' << c.pan.y
        << ' ' << c.target.x << ' ' << c.target.y << ' ' << c.target.z << ' '
        << static_cast<int>(c.projection) << '\n';
    out << static_cast<int>(r.mode) << ' ' << static_cast<int>(r.colorMode) << ' '
        << static_cast<int>(r.charset) << ' ' << static_cast<int>(r.backend) << ' ' << r.color << ' '
        << r.showInfo << ' ' << r.showSequence << ' ' << r.showHelp << ' ' << r.backgroundColor << ' '
        << r.residueLabels << ' ' << r.chainLabels << ' ' << state.autorotate << '\n';
    out << (r.selectedResidue ? static_cast<long long>(*r.selectedResidue) : -1) << ' '
        << (r.inspectedAtom ? static_cast<long long>(*r.inspectedAtom) : -1) << '\n';
    out << std::quoted(bits(r.visibleAtoms)) << ' ' << std::quoted(bits(r.selectedAtoms)) << '\n';
    out << r.atomColors.size();
    for (int color : r.atomColors)
        out << ' ' << color;
    out << '\n';
    std::map<std::string, AtomMask> sorted(r.namedSelections.begin(), r.namedSelections.end());
    out << sorted.size() << '\n';
    for (const auto& [name, m] : sorted)
        out << std::quoted(name) << ' ' << std::quoted(bits(m)) << '\n';
    auto* t = state.trajectory;
    out << static_cast<bool>(t && t->source) << '\n';
    if (t && t->source) {
        out << std::quoted(std::filesystem::absolute(t->path).string()) << ' ' << t->frameIndex << ' '
            << t->rangeStart << ' ' << t->rangeEnd << ' ' << t->stride << ' ' << t->frameDelayMs << ' '
            << t->playing << ' ' << t->loop << '\n';
        out << t->alignmentAtoms.size();
        for (auto i : t->alignmentAtoms)
            out << ' ' << i;
        out << '\n';
        out << t->alignmentReference.size();
        for (auto p : t->alignmentReference)
            out << ' ' << p.x << ' ' << p.y << ' ' << p.z;
        out << '\n';
    }
    out << std::quoted(r.densityMap ? r.densityMap->path.string() : "") << ' ' << r.mapAxis << ' '
        << r.mapSlice << ' ' << r.mapLevel << '\n';
    std::ofstream file(path);
    if (!file || !(file << out.str())) throw std::runtime_error("cannot save session");
}
void loadSession(const std::filesystem::path& path, CommandState state) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open session");
    if (std::filesystem::file_size(path) > 128ULL * 1024 * 1024)
        throw std::runtime_error("session exceeds size limit");
    std::string magic, source, assembly;
    int version, model;
    checked(static_cast<bool>(in >> magic >> version) && magic == "ProteinTTY-session" && version == 1);
    checked(static_cast<bool>(in >> std::quoted(source) >> model >> std::quoted(assembly)));
    auto structure = loadStructure(source);
    selectModel(structure, model);
    if (!assembly.empty()) expandAssembly(structure, assembly);
    Camera c;
    RenderSettings r;
    TrajectoryState trajectory;
    bool spin = false;
    int projection, mode, color, charset, backend;
    checked(static_cast<bool>(in >> c.yaw >> c.pitch >> c.roll >> c.zoom >> c.pan.x >> c.pan.y >>
                              c.target.x >> c.target.y >> c.target.z >> projection));
    for (auto x : {c.yaw, c.pitch, c.roll, c.zoom, c.pan.x, c.pan.y, c.target.x, c.target.y, c.target.z})
        finite(x);
    checked(c.zoom > 0 && projection >= 0 && projection <= 1);
    c.projection = static_cast<Projection>(projection);
    checked(static_cast<bool>(in >> mode >> color >> charset >> backend >> r.color >> r.showInfo >>
                              r.showSequence >> r.showHelp >> r.backgroundColor >> r.residueLabels >>
                              r.chainLabels >> spin));
    checked(mode >= 0 && mode <= 4 && color >= 0 && color <= 4 && charset >= 0 && charset <= 1 &&
            backend >= 0 && backend <= 1);
    r.mode = static_cast<RenderMode>(mode);
    r.colorMode = static_cast<ColorMode>(color);
    r.charset = static_cast<Charset>(charset);
    r.backend = static_cast<RendererBackend>(backend);
    long long residue, atom;
    checked(static_cast<bool>(in >> residue >> atom));
    checked(residue >= -1 && atom >= -1 &&
            (residue < 0 || static_cast<std::size_t>(residue) < structure.residues.size()) &&
            (atom < 0 || static_cast<std::size_t>(atom) < structure.atoms.size()));
    if (residue >= 0) r.selectedResidue = static_cast<std::size_t>(residue);
    if (atom >= 0) r.inspectedAtom = static_cast<std::size_t>(atom);
    std::string text;
    checked(static_cast<bool>(in >> std::quoted(text)));
    r.visibleAtoms = mask(text, structure.atoms.size());
    checked(static_cast<bool>(in >> std::quoted(text)));
    r.selectedAtoms = mask(text, structure.atoms.size());
    std::size_t count;
    checked(static_cast<bool>(in >> count) && (count == 0 || count == structure.atoms.size()));
    r.atomColors.resize(count);
    for (auto& x : r.atomColors)
        checked(static_cast<bool>(in >> x));
    checked(static_cast<bool>(in >> count) && count <= 10000);
    for (std::size_t i = 0; i < count; ++i) {
        std::string name;
        checked(static_cast<bool>(in >> std::quoted(name) >> std::quoted(text)));
        r.namedSelections[name] = mask(text, structure.atoms.size());
    }
    bool hasTrajectory;
    checked(static_cast<bool>(in >> hasTrajectory));
    if (hasTrajectory) {
        checked(state.trajectory != nullptr);
        std::string trajectoryPath;
        std::size_t frame, first, last;
        int stride, delay;
        bool playing, loop;
        checked(static_cast<bool>(in >> std::quoted(trajectoryPath) >> frame >> first >> last >> stride >>
                                  delay >> playing >> loop));
        checked(stride >= 1 && stride <= 100000 && delay >= 10 && delay <= 10000);
        loadTrajectory(trajectoryPath, structure, trajectory);
        setTrajectoryRange(structure, trajectory, first, last);
        trajectory.stride = stride;
        trajectory.frameDelayMs = delay;
        trajectory.loop = loop;
        checked(static_cast<bool>(in >> count) && count <= structure.atoms.size());
        trajectory.alignmentAtoms.resize(count);
        for (auto& i : trajectory.alignmentAtoms)
            checked(static_cast<bool>(in >> i) && i < structure.atoms.size());
        checked(static_cast<bool>(in >> count) && (count == 0 || count == structure.atoms.size()));
        trajectory.alignmentReference.resize(count);
        for (auto& p : trajectory.alignmentReference) {
            checked(static_cast<bool>(in >> p.x >> p.y >> p.z));
            finite(p.x);
            finite(p.y);
            finite(p.z);
        }
        checked(trajectory.alignmentAtoms.empty() || count == structure.atoms.size());
        seekTrajectory(structure, trajectory, frame);
        trajectory.playing = playing;
    }
    std::string mapPath;
    checked(static_cast<bool>(in >> std::quoted(mapPath) >> r.mapAxis >> r.mapSlice >> r.mapLevel));
    finite(r.mapLevel);
    checked(r.mapAxis >= 0 && r.mapAxis <= 2 && r.mapSlice >= 0);
    if (!mapPath.empty()) {
        r.densityMap = std::make_shared<DensityMap>(loadDensityMap(mapPath));
        checked(r.mapSlice < r.densityMap->dimensions[r.mapAxis]);
    }
    in >> std::ws;
    checked(in.eof());
    r.trajectoryInfo = trajectoryStatus(trajectory);
    state.structure = std::move(structure);
    state.camera = c;
    state.render = std::move(r);
    state.autorotate = spin;
    if (state.trajectory) *state.trajectory = std::move(trajectory);
    if (state.structureChanged) state.structureChanged();
}
} // namespace ascii_protein
