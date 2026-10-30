#include "ascii_protein/trajectory_analysis.hpp"
#include "ascii_protein/trajectory.hpp"
#include "ascii_protein/molecular.hpp"
#include "ascii_protein/structure_loader.hpp"
#include "ascii_protein/pdb_parser.hpp"
#include "ascii_protein/mmcif_parser.hpp"
#include "ascii_protein/session.hpp"
#include "ascii_protein/config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <bit>
#include <chrono>
#include <stdexcept>
using namespace ascii_protein;
namespace {
void check(bool value, const std::string& message) {
    if (!value)
        throw std::runtime_error(message);
}
void close(double a, double b, const std::string& message, double tolerance = 1e-8) {
    check(std::abs(a - b) < tolerance, message);
}
template <class F> void rejects(F f, const std::string& message) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return;
    }
    throw std::runtime_error(message);
}
struct Frames final : CoordinateFrameSource {
    std::vector<std::vector<Vec3>> data;
    std::size_t atomCount() const override { return data[0].size(); }
    std::size_t frameCount() const override { return data.size(); }
    std::string_view atomLabel(std::size_t) const override { return {}; }
    std::string_view frameLabel(std::size_t) const override { return {}; }
    const std::vector<Vec3>& frame(std::size_t i) const override { return data.at(i); }
};
void analysis() {
    std::vector<Vec3> ref{{0, 0, 0}, {1, 0, 0}, {0, 2, 0}, {0, 0, 3}}, moving;
    auto rotation = Mat3::rotationZ(radians(137)) * Mat3::rotationY(radians(71));
    for (auto p : ref)
        moving.push_back(rotation * p + Vec3{12, -8, 3});
    auto fit = fitCoordinates(moving, ref, {0, 1, 2, 3});
    for (std::size_t i = 0; i < ref.size(); ++i)
        close(length(fit.apply(moving[i]) - ref[i]), 0, "rigid fit");
    moving.clear();
    for (auto p : ref)
        moving.push_back(Mat3::rotationX(pi) * p + Vec3{2, 3, 4});
    fit = fitCoordinates(moving, ref, {0, 1, 2, 3});
    for (std::size_t i = 0; i < ref.size(); ++i)
        close(length(fit.apply(moving[i]) - ref[i]), 0, "180 degree fit");
    auto one = fitCoordinates({{3, 4, 5}}, {{1, 2, 3}}, {0});
    close(length(one.apply({3, 4, 5}) - Vec3{1, 2, 3}), 0, "one atom fit");
    rejects([&] { fitCoordinates(ref, ref, {}); }, "empty fit rejected");
    {
        Structure s;
        for (int i = 0; i < 4; ++i) {
            Atom a;
            a.serial = i + 1;
            a.name = std::to_string(i);
            s.atoms.push_back(a);
        }
        s.finalize();
        TrajectoryState playback;
        auto f = std::make_unique<Frames>();
        auto moved = ref;
        moved[3] = {100, 100, 100};
        f->data = {ref, moved};
        playback.source = std::move(f);
        playback.rangeEnd = 1;
        playback.alignmentReference = ref;
        playback.alignmentAtoms = {0, 1, 2};
        seekTrajectory(s, playback, 0);
        const auto anchor = s.atoms[0].position;
        const double radius = s.boundingRadius;
        seekTrajectory(s, playback, 1);
        close(length(s.atoms[0].position - anchor), 0,
              "selected alignment anchor stable despite unselected motion");
        close(s.boundingRadius, radius, "aligned playback keeps reference fit scale");
    }
    moving = ref;
    for (auto& p : moving)
        p.x = -p.x;
    fit = fitCoordinates(moving, ref, {0, 1, 2, 3});
    double residual = 0;
    for (std::size_t i = 0; i < ref.size(); ++i)
        residual += length(fit.apply(moving[i]) - ref[i]);
    check(residual > 0.1, "fit cannot reflect chiral coordinates");
    TrajectoryState t;
    auto frames = std::make_unique<Frames>();
    frames->data = {{{0, 0, 0}, {1, 0, 0}}, {{0, 0, 0}, {3, 0, 0}}, {{0, 0, 0}, {5, 0, 0}}};
    t.source = std::move(frames);
    t.rangeEnd = 2;
    auto rmsd = trajectoryRmsd(t, {0, 1}, 0);
    close(rmsd[1], std::sqrt(2.), "RMSD analytic fixture");
    auto rmsf = trajectoryRmsf(t, {0, 1});
    close(rmsf[0], 0, "fixed atom RMSF");
    close(rmsf[1], std::sqrt(8. / 3), "RMSF analytic fixture");
    t.stride = 2;
    close(trajectoryRmsf(t, {1})[0], 2, "RMSF stride fixture");
    auto rigid = std::make_unique<Frames>();
    rigid->data = {ref, moving};
    t.source = std::move(rigid);
    t.rangeEnd = 1;
    t.stride = 1;
    t.alignmentReference = ref;
    t.alignmentAtoms = {0, 1, 2, 3};
    check(trajectoryRmsd(t, {0, 1, 2, 3}, 0)[1] > 0.1, "aligned reflected RMSD remains nonzero");
    moving.clear();
    for (auto p : ref)
        moving.push_back(rotation * p + Vec3{12, -8, 3});
    static_cast<Frames*>(t.source.get())->data[1] = moving;
    close(trajectoryRmsd(t, {0, 1, 2, 3}, 0)[1], 0, "aligned RMSD");
    for (auto x : trajectoryRmsf(t, {0, 1, 2, 3}))
        close(x, 0, "aligned RMSF");
}
std::string atoms() {
    return "ATOM      1  CA AGLY A 100A     10.000  20.000  30.000  0.75 12.50           C  \n"
           "ATOM      2  N   GLY A 100A     20.000  20.000  30.000  1.00 13.50           N  \n"
           "ATOM      3  O   GLY A 100A     21.000  20.000  30.000  1.00 14.50           O  \n";
}
void molecular(const std::filesystem::path& dir) {
    std::istringstream input("MODEL        3\n" + atoms() + "ENDMDL\nMODEL        7\n" + atoms() +
                             "ENDMDL\nCONECT    1    2\n");
    auto s = PdbParser::parse(input, "models.pdb");
    check(s.model == 3 && s.atoms.size() == 3, "first deposited model selected");
    selectModel(s, 7);
    check(s.model == 7 && s.atoms.size() == 3, "model switching");
    rejects([&] { selectModel(s, 8); }, "missing model rejected");
    check(s.model == 7, "failed model preserves state");
    Renderer labelsRenderer(s);
    Camera labelsCamera;
    RenderSettings labels;
    labels.showInfo = false;
    labels.chainLabels = true;
    labels.residueLabels = true;
    const auto labeled = labelsRenderer.render(labelsCamera, labels, 100, 32).toPlainText();
    check(labeled.find("chain A") != std::string::npos && labeled.find("A:GLY100A") != std::string::npos,
          "chain and residue labels occupy separate empty areas");
    auto bonds = inferBonds(s);
    check(bonds.size() == 1 && bonds[0].first == 0 && bonds[0].second == 1,
          "explicit long bond overrides inference for connected atoms");
    for (auto ext : {".pdb", ".cif"}) {
        auto path = dir / (std::string("roundtrip") + ext);
        exportStructure(s, path, {true, true, true});
        auto copy = loadStructure(path);
        check(copy.atoms.size() == 3 && depositedBonds(copy).size() == 1, "export preserves topology");
        check(copy.atoms[0].name == "CA" && copy.atoms[0].alternateLocation == 'A' &&
                  copy.atoms[0].insertionCode == 'A',
              "export identifiers");
        close(copy.atoms[0].position.x + copy.originalCentroid.x, 10, "export original coordinates", 1e-3);
        exportStructure(s, path, {true, false, true});
        copy = loadStructure(path);
        check(copy.atoms.size() == 2 && depositedBonds(copy).empty(),
              "selection export excludes dangling bonds");
    }
    s.atoms[0].chain = "LONG";
    exportStructure(s, dir / "long.cif", {true, true, true});
    check(loadStructure(dir / "long.cif").atoms[0].chain == "LONG", "mmCIF long chain round trip");
    rejects([&] { exportStructure(s, dir / "long.pdb", {true, true, true}); }, "PDB width rejected");
    std::istringstream assembly("REMARK 350 BIOMOLECULE: 1\nREMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n"
                                "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
                                "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n"
                                "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n"
                                "REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000       10.00000\n"
                                "REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000\n"
                                "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000\n" +
                                atoms() + "CONECT    1    2\n");
    s = PdbParser::parse(assembly, "assembly.pdb");
    expandAssembly(s, "1");
    check(s.atoms.size() == 6 && s.chains.size() == 2 && depositedBonds(s).size() == 2,
          "assembly copies chains and bonds");
    close(s.atoms[3].position.x - s.atoms[0].position.x, 10, "assembly translation");
    exportStructure(s, dir / "assembly.cif", AtomMask(6, true));
    check(loadStructure(dir / "assembly.cif").atoms.size() == 6, "assembly export reload");
}
void cifMetadata(const std::filesystem::path& dir) {
    std::istringstream input(
        "data_metadata\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.auth_atom_id\n"
        "_atom_site.label_atom_id\n_atom_site.auth_asym_id\n_atom_site.label_asym_id\n"
        "_atom_site.auth_seq_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site."
        "Cartn_z\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 CA CA AUTH L 100 1 0 0 0 3\nATOM 2 N N AUTH L 100 1 10 0 0 3\nATOM 3 CA CA AUTH L 100 1 20 0 "
        "0 7\n#\n"
        "_struct_conn.conn_type_id covale\n_struct_conn.ptnr1_label_asym_id "
        "L\n_struct_conn.ptnr1_label_seq_id 1\n"
        "_struct_conn.ptnr1_label_atom_id CA\n_struct_conn.ptnr2_label_asym_id "
        "L\n_struct_conn.ptnr2_label_seq_id 1\n"
        "_struct_conn.ptnr2_label_atom_id N\n"
        "loop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list."
        "matrix[1][2]\n"
        "_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list."
        "matrix[2][2]\n"
        "_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list."
        "matrix[3][2]\n"
        "_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list."
        "vector[2]\n"
        "_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 1 0 0 0 1 0 0 0\n2 1 0 0 0 1 0 0 0 1 10 0 0\n"
        "_pdbx_struct_assembly_gen.assembly_id 1\n_pdbx_struct_assembly_gen.oper_expression '(1-2)(1)'\n"
        "_pdbx_struct_assembly_gen.asym_id_list L\n");
    auto s = MmcifParser::parse(input, "metadata.cif");
    check(s.atoms.size() == 2 && s.model == 3 && depositedBonds(s).size() == 1, "mmCIF singleton connection");
    exportStructure(s, dir / "metadata-roundtrip.cif", AtomMask(s.atoms.size(), true));
    const auto copy = loadStructure(dir / "metadata-roundtrip.cif");
    check(copy.atoms[0].chain == "AUTH" && copy.atoms[0].labelChain == "L" && copy.atoms[0].labelResidue == 1,
          "mmCIF preserves author and label identifiers");
    selectModel(s, 7);
    check(s.atoms.size() == 1, "mmCIF ensemble switching");
    selectModel(s, 3);
    expandAssembly(s, "1");
    check(s.atoms.size() == 4 && depositedBonds(s).size() == 2, "mmCIF assembly product expression");
    close(s.atoms[2].position.x - s.atoms[0].position.x, 10, "mmCIF assembly operators");
}
void put32(std::string& s, int offset, std::uint32_t value, bool big) {
    for (int i = 0; i < 4; ++i)
        s[offset + i] = static_cast<char>((value >> (big ? (3 - i) * 8 : i * 8)) & 255);
}
void makeMap(const std::filesystem::path& path, bool big, int mode = 2) {
    std::string data(1024 + 16 + (mode == 2 ? 4 : mode == 0 ? 1 : 2) * 8, '\0');
    auto integer = [&](int offset, int n) { put32(data, offset, static_cast<std::uint32_t>(n), big); };
    auto real = [&](int offset, float n) { put32(data, offset, std::bit_cast<std::uint32_t>(n), big); };
    integer(0, 2);
    integer(4, 2);
    integer(8, 2);
    integer(12, mode);
    integer(16, 3);
    integer(20, 4);
    integer(24, 5);
    integer(28, 2);
    integer(32, 2);
    integer(36, 2);
    real(40, 4);
    real(44, 6);
    real(48, 8);
    real(52, 90);
    real(56, 90);
    real(60, 90);
    integer(64, 2);
    integer(68, 1);
    integer(72, 3);
    integer(92, 16);
    data.replace(208, 4, "MAP ");
    data[212] = static_cast<char>(big ? 0x11 : 0x44);
    const int bytes = mode == 2 ? 4 : (mode == 0 ? 1 : 2);
    for (int i = 0; i < 8; ++i) {
        std::uint32_t x =
            mode == 2 ? std::bit_cast<std::uint32_t>(static_cast<float>(i)) : static_cast<std::uint32_t>(i);
        for (int j = 0; j < bytes; ++j)
            data[1040 + i * bytes + j] = static_cast<char>((x >> (big ? (bytes - 1 - j) * 8 : j * 8)) & 255);
    }
    std::ofstream(path, std::ios::binary).write(data.data(), static_cast<std::streamsize>(data.size()));
}
void maps(const std::filesystem::path& dir) {
    for (bool big : {false, true})
        for (int mode : {0, 1, 2, 6}) {
            auto path = dir / "map.mrc";
            makeMap(path, big, mode);
            auto map = loadDensityMap(path);
            close(map.value(1, 1, 1), 7, "MRC voxel ordering/endian/mode");
            auto p = map.position(1, 1, 1);
            close(p.x, 10, "MRC row permutation and starts");
            close(p.y, 12, "MRC column permutation and starts");
            close(p.z, 24, "MRC sections and starts");
            check(map.minimum == 0 && map.maximum == 7, "MRC computes stats");
        }
    auto path = dir / "map.mrc";
    makeMap(path, false);
    {
        std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
        std::string origin(12, '\0');
        put32(origin, 0, std::bit_cast<std::uint32_t>(100.f), false);
        f.seekp(196);
        f.write(origin.data(), 12);
    }
    close(loadDensityMap(path).position(0, 0, 0).x, 100, "MRC explicit origin");
    std::filesystem::resize_file(path, 1030);
    rejects([&] { loadDensityMap(path); }, "truncated map rejected");
    std::ofstream(path) << "bad";
    rejects([&] { loadDensityMap(path); }, "short map rejected");
    const auto corrupt = [&](int offset, std::uint32_t value) {
        makeMap(path, false);
        std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
        std::string bytes(4, '\0');
        put32(bytes, 0, value, false);
        f.seekp(offset);
        f.write(bytes.data(), 4);
        f.close();
        rejects([&] { loadDensityMap(path); }, "invalid MRC header/value rejected");
    };
    corrupt(0, 0);
    corrupt(0, 0x7fffffff);
    corrupt(12, 4);
    corrupt(64, 3);
    corrupt(92, 0xffffffff);
    corrupt(40, std::bit_cast<std::uint32_t>(-1.f));
    corrupt(1040, 0x7fc00000);
    makeMap(path, false);
}
void sessions(const std::filesystem::path& dir) {
    auto pdb = dir / "session.pdb";
    std::ofstream(pdb) << atoms();
    auto structure = loadStructure(pdb);
    Camera camera;
    RenderSettings render;
    TrajectoryState trajectory;
    CommandConfig config;
    bool spin = false, running = true;
    CommandState state{camera, render, structure, spin, running, {}, {}, &trajectory, &config};
    auto xyz = dir / "session.xyz";
    std::ofstream(xyz) << "3\nfirst\nC 10 20 30\nN 20 20 30\nO 21 20 30\n"
                          "3\nsecond\nC 11 20 30\nN 21 20 30\nO 22 20 30\n";
    loadTrajectory(xyz, structure, trajectory);
    check(executeCommands(
              "align all; frame 2; select test chain A; hide serial 3; color test red; set residue-labels on",
              state)
              .success,
          "workflow commands");
    render.densityMap = std::make_shared<DensityMap>(loadDensityMap(dir / "map.mrc"));
    render.mapSlice = 1;
    render.mapLevel = 3;
    camera.pan = {0.3, 0.7};
    camera.target = {1, 2, 3};
    camera.roll = 37;
    spin = true;
    render.trajectoryInfo = trajectoryStatus(trajectory);
    Renderer renderer(structure);
    const auto before = renderer.render(camera, render, 100, 32).toPlainText();
    saveSession(dir / "saved.ptty-session", state);
    camera.reset();
    render = {};
    trajectory = {};
    spin = false;
    loadSession(dir / "saved.ptty-session", state);
    Renderer restored(structure);
    check(restored.render(camera, render, 100, 32).toPlainText() == before, "session deterministic frame");
    check(trajectory.frameIndex == 1 && trajectory.alignmentAtoms.size() == 3 && spin &&
              render.namedSelections.contains("test"),
          "session analysis and selection state");
    std::ofstream(dir / "invalid.session") << "ProteinTTY-session 99\n";
    rejects([&] { loadSession(dir / "invalid.session", state); }, "invalid session rejected");
    close(camera.roll, 37, "failed session preserves state");
    auto configPath = dir / "config.ptty";
    std::ofstream(configPath) << "# configuration\nmacro home \"frame 1; pause\"\nbind g \"run home\"\nmacro "
                                 "recurse \"run recurse\"\n";
    config = loadCommandConfig(configPath);
    check(config.bindings.at('g') == "run home" && executeCommands("run home", state).success &&
              trajectory.frameIndex == 0,
          "macro dispatch");
    check(!executeCommands("run recurse", state).success && config.depth == 0, "macro recursion bounded");
    check(!executeCommands("map slice s 99", state).success, "invalid map slice rejected");
    check(executeCommands("rmsd all; rmsf color all; align off", state).success, "analysis command dispatch");
}
} // namespace
int main() {
    auto dir = std::filesystem::temp_directory_path() /
               ("protein-tty-workflows-" +
                std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    try {
        analysis();
        molecular(dir);
        cifMetadata(dir);
        maps(dir);
        sessions(dir);
        std::filesystem::remove_all(dir);
        std::cout << "Workflow tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << " (fixtures: " << dir << ")\n";
        return 1;
    }
}
