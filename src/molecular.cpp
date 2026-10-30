#include "ascii_protein/molecular.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <tuple>
#include <charconv>
#include <cctype>
namespace ascii_protein {
namespace {
std::string get(const MetadataRow& r, const std::string& k) {
    auto it = r.find(k);
    return it == r.end() || it->second == "." || it->second == "?" ? "" : it->second;
}
template <class T> T parseNumber(const std::string& text) {
    T number{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), number);
    if (error != std::errc{} || end != text.data() + text.size() ||
        !std::isfinite(static_cast<double>(number)))
        throw std::runtime_error("invalid assembly numeric value: " + text);
    return number;
}
std::vector<std::string> split(std::string s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string t;
    while (std::getline(in, t, ','))
        if (!t.empty())
            out.push_back(t);
    return out;
}
std::vector<std::string> operatorList(const std::string& expression) {
    std::vector<std::string> result;
    for (const auto& token : split(expression)) {
        auto dash = token.find('-');
        if (dash == std::string::npos)
            result.push_back(token);
        else {
            int a = parseNumber<int>(token.substr(0, dash)), b = parseNumber<int>(token.substr(dash + 1));
            if (a < 0 || b < a || b - a > 10000)
                throw std::runtime_error("invalid assembly operator range");
            for (long long i = a; i <= b; ++i)
                result.push_back(std::to_string(i));
        }
    }
    return result;
}
std::vector<std::vector<std::string>> expressionGroups(const std::string& expression) {
    if (expression.find('(') == std::string::npos)
        return {operatorList(expression)};
    std::vector<std::vector<std::string>> groups;
    std::size_t pos = 0;
    while (pos < expression.size()) {
        if (expression[pos] != '(')
            throw std::runtime_error("invalid assembly expression");
        auto end = expression.find(')', pos);
        if (end == std::string::npos)
            throw std::runtime_error("invalid assembly expression");
        groups.push_back(operatorList(expression.substr(pos + 1, end - pos - 1)));
        pos = end + 1;
    }
    return groups;
}
std::string cif(const std::string& value) {
    if (value.empty())
        return ".";
    if (value.find_first_of(" \t\r\n'\"#;") == std::string::npos && value.front() != '_' && value != "." &&
        value != "?")
        return value;
    if (value.find('\n') != std::string::npos)
        throw std::runtime_error("multiline identifiers cannot be exported");
    if (value.find('\'') == std::string::npos)
        return "'" + value + "'";
    if (value.find('"') == std::string::npos)
        return "\"" + value + "\"";
    throw std::runtime_error("identifier contains both quote characters");
}
bool endpoint(const Atom& a, const MetadataRow& row, int partner) {
    const std::string p = "ptnr" + std::to_string(partner) + "_";
    const auto chain = get(row, p + "auth_asym_id"), seq = get(row, p + "auth_seq_id"),
               name = get(row, p + "auth_atom_id");
    const bool auth = !seq.empty();
    const auto ins = get(row, "pdbx_" + p + "PDB_ins_code");
    const auto alt = get(row, "pdbx_" + p + "label_alt_id");
    return (auth ? a.chain == chain && std::to_string(a.residueNumber) == seq
                 : a.labelChain == get(row, p + "label_asym_id") &&
                       std::to_string(a.labelResidue) == get(row, p + "label_seq_id")) &&
           (name.empty() ? a.labelName == get(row, p + "label_atom_id") : a.name == name) &&
           (ins.empty() ? a.insertionCode == ' ' : a.insertionCode == ins.front()) &&
           (alt.empty() || a.alternateLocation == alt.front());
}
} // namespace
void readMolecularRow(Structure& s, const std::string& category, const MetadataRow& row) {
    if (category == "_struct_conn") {
        auto type = get(row, "conn_type_id");
        if (type == "covale" || type == "disulf" || type == "metalc")
            s.connections.push_back(row);
    } else if (category == "_pdbx_struct_oper_list") {
        AssemblyOperator op;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                auto x = get(row, "matrix[" + std::to_string(i + 1) + "][" + std::to_string(j + 1) + "]");
                if (x.empty())
                    throw std::runtime_error("incomplete assembly matrix");
                op.rotation.m[i][j] = parseNumber<double>(x);
            }
        op.translation = {parseNumber<double>(get(row, "vector[1]")),
                          parseNumber<double>(get(row, "vector[2]")),
                          parseNumber<double>(get(row, "vector[3]"))};
        s.operators[get(row, "id")] = op;
    } else if (category == "_pdbx_struct_assembly_gen") {
        s.assemblies.push_back(
            {get(row, "assembly_id"), get(row, "oper_expression"), split(get(row, "asym_id_list"))});
    }
}
// Resolve deposited connections only after model/conformer selection.
std::vector<std::pair<std::size_t, std::size_t>> depositedBonds(const Structure& s) {
    if (s.explicitBonds.empty() && s.connections.empty())
        return {};
    std::set<std::pair<std::size_t, std::size_t>> result;
    auto add = [&](std::size_t a, std::size_t b) {
        if (a != b)
            result.emplace(std::min(a, b), std::max(a, b));
    };
    std::map<int, std::size_t> serials;
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        serials[s.atoms[i].serial] = i;
    for (const auto& [a, b] : s.explicitBonds)
        if (serials.contains(a) && serials.contains(b))
            add(serials[a], serials[b]);
    if (s.connections.empty())
        return {result.begin(), result.end()};
    using Address = std::pair<std::string, std::string>;
    std::map<Address, std::vector<std::size_t>> authLookup, labelLookup;
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const auto& a = s.atoms[i];
        authLookup[{a.chain, std::to_string(a.residueNumber)}].push_back(i);
        labelLookup[{a.labelChain, std::to_string(a.labelResidue)}].push_back(i);
    }
    for (const auto& row : s.connections) {
        auto sym1 = get(row, "ptnr1_symmetry"), sym2 = get(row, "ptnr2_symmetry");
        if ((!sym1.empty() && sym1 != "1_555") || (!sym2.empty() && sym2 != "1_555"))
            continue;
        auto find = [&](int partner) {
            const auto prefix = "ptnr" + std::to_string(partner) + "_";
            const auto seq = get(row, prefix + "auth_seq_id");
            const bool auth = !seq.empty();
            const auto& lookup = auth ? authLookup : labelLookup;
            const Address key{get(row, prefix + (auth ? "auth_asym_id" : "label_asym_id")),
                              auth ? seq : get(row, prefix + "label_seq_id")};
            auto it = lookup.find(key);
            if (it != lookup.end())
                for (auto i : it->second)
                    if (endpoint(s.atoms[i], row, partner))
                        return i;
            return s.atoms.size();
        };
        const auto a = find(1), b = find(2);
        if (a < s.atoms.size() && b < s.atoms.size())
            add(a, b);
    }
    return {result.begin(), result.end()};
}
void selectModel(Structure& s, int model) {
    if (std::none_of(s.depositedAtoms.begin(), s.depositedAtoms.end(),
                     [&](const Atom& a) { return a.model == model; }))
        throw std::runtime_error("model does not exist");
    s.atoms = s.depositedAtoms;
    s.model = model;
    s.finalize();
}
void expandAssembly(Structure& s, const std::string& id) {
    std::vector<Atom> atoms;
    std::vector<std::pair<int, int>> bonds;
    std::vector<SecondaryRange> ranges;
    const auto originalBonds = depositedBonds(s);
    std::set<std::string> seen;
    for (const auto& gen : s.assemblies)
        if (gen.id == id) {
            auto groups = expressionGroups(gen.expression);
            std::vector<std::pair<std::string, AssemblyOperator>> transforms{{"", {}}};
            for (const auto& group : groups) {
                std::vector<std::pair<std::string, AssemblyOperator>> next;
                for (const auto& [key, t] : transforms)
                    for (const auto& opId : group) {
                        auto it = s.operators.find(opId);
                        if (it == s.operators.end())
                            throw std::runtime_error("assembly operator not found: " + opId);
                        const auto& op = it->second;
                        next.push_back(
                            {key + ":" + opId,
                             {t.rotation * op.rotation, t.rotation * op.translation + t.translation}});
                        if (next.size() > 10000)
                            throw std::runtime_error("assembly operator count limit exceeded");
                    }
                transforms = std::move(next);
            }
            for (const auto& [key, t] : transforms) {
                std::map<std::size_t, int> mapped;
                std::map<std::string, std::string> chainCopies;
                for (std::size_t i = 0; i < s.atoms.size(); ++i) {
                    const auto& a = s.atoms[i];
                    const auto& chain = a.labelChain.empty() ? a.chain : a.labelChain;
                    if (std::find(gen.chains.begin(), gen.chains.end(), chain) == gen.chains.end())
                        continue;
                    if (!seen.insert(key + "/" + std::to_string(i)).second)
                        continue;
                    if (atoms.size() >= 5000000)
                        throw std::runtime_error("assembly atom limit exceeded");
                    auto copy = a;
                    copy.position = t.rotation * (a.position + s.originalCentroid) + t.translation;
                    copy.chain = a.chain + "@" + key;
                    copy.labelChain = chain + "@" + key;
                    chainCopies[a.chain] = copy.chain;
                    copy.serial = static_cast<int>(atoms.size() + 1);
                    mapped[i] = copy.serial;
                    atoms.push_back(std::move(copy));
                }
                for (auto range : s.secondaryRanges)
                    if (chainCopies.contains(range.startChain) && chainCopies.contains(range.endChain)) {
                        range.startChain = chainCopies.at(range.startChain);
                        range.endChain = chainCopies.at(range.endChain);
                        ranges.push_back(std::move(range));
                    }
                for (auto [a, b] : originalBonds)
                    if (mapped.contains(a) && mapped.contains(b))
                        bonds.emplace_back(mapped[a], mapped[b]);
            }
        }
    if (atoms.empty())
        throw std::runtime_error("assembly not found or contains no atoms");
    s.assemblyApplied = id;
    s.atoms = std::move(atoms);
    s.depositedAtoms.clear();
    s.explicitBonds = std::move(bonds);
    s.connections.clear();
    s.assemblies.clear();
    s.operators.clear();
    s.secondaryRanges = std::move(ranges);
    s.finalize();
}
void exportStructure(const Structure& s, const std::filesystem::path& path, const std::vector<bool>& mask) {
    if (mask.size() != s.atoms.size())
        throw std::runtime_error("export selection size mismatch");
    auto extension = path.extension().string();
    for (auto& c : extension)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (extension != ".pdb" && extension != ".cif" && extension != ".mmcif")
        throw std::runtime_error("export requires .pdb or .cif");
    std::ostringstream out;
    out << std::fixed;
    const bool pdb = extension == ".pdb";
    if (pdb && (s.model < 1 || s.model > 9999))
        throw std::runtime_error("model number exceeds PDB field; export .cif instead");
    if (pdb && s.model != 1)
        out << "MODEL     " << std::setw(4) << s.model << '\n';
    if (!pdb)
        out << "data_export\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_"
               "site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_"
               "id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site."
               "Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site."
               "pdbx_PDB_model_num\n_atom_site.label_atom_id\n_atom_site.label_asym_id\n_atom_site.label_seq_"
               "id\n";
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (mask[i]) {
            const auto& a = s.atoms[i];
            auto p = a.position + s.originalCentroid;
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
                !std::isfinite(a.occupancy) || !std::isfinite(a.bFactor))
                throw std::runtime_error("cannot export non-finite atom values");
            if (pdb) {
                if (a.chain.size() > 1 || a.name.size() > 4 || a.residueName.size() > 3 ||
                    a.element.size() > 2 || a.serial < 0 || a.serial > 99999 || a.residueNumber < -999 ||
                    a.residueNumber > 9999 || p.x < -999.999 || p.y < -999.999 || p.z < -999.999 ||
                    p.x > 9999.999 || p.y > 9999.999 || p.z > 9999.999 || a.occupancy < -99.99 ||
                    a.occupancy > 999.99 || a.bFactor < -99.99 || a.bFactor > 999.99)
                    throw std::runtime_error(
                        "structure identifiers or coordinates exceed PDB fields; export .cif instead");
                out << std::left << std::setw(6) << (a.hetero ? "HETATM" : "ATOM") << std::right
                    << std::setw(5) << a.serial << ' ' << std::left << std::setw(4) << a.name
                    << a.alternateLocation << std::setw(3) << a.residueName << ' '
                    << (a.chain.empty() ? ' ' : a.chain.front()) << std::right << std::setw(4)
                    << a.residueNumber << a.insertionCode << "   " << std::setprecision(3) << std::setw(8)
                    << p.x << std::setw(8) << p.y << std::setw(8) << p.z << std::setprecision(2)
                    << std::setw(6) << a.occupancy << std::setw(6) << a.bFactor << "          "
                    << std::setw(2) << a.element << '\n';
            } else {
                out << (a.hetero ? "HETATM" : "ATOM") << ' ' << a.serial << ' ' << cif(a.element) << ' '
                    << cif(a.name) << ' '
                    << cif(a.alternateLocation == ' ' ? "" : std::string(1, a.alternateLocation)) << ' '
                    << cif(a.residueName) << ' ' << cif(a.chain) << ' ' << a.residueNumber << ' '
                    << cif(a.insertionCode == ' ' ? "" : std::string(1, a.insertionCode)) << ' '
                    << std::setprecision(6) << p.x << ' ' << p.y << ' ' << p.z << ' ' << a.occupancy << ' '
                    << a.bFactor << ' ' << a.model << ' ' << cif(a.labelName.empty() ? a.name : a.labelName)
                    << ' ' << cif(a.labelChain.empty() ? a.chain : a.labelChain) << ' '
                    << (a.labelResidue == 0 ? "." : std::to_string(a.labelResidue)) << '\n';
            }
        }
    if (pdb && s.model != 1)
        out << "ENDMDL\n";
    auto bonds = depositedBonds(s);
    bool header = false;
    for (auto [i, j] : bonds)
        if (mask[i] && mask[j]) {
            const auto& a = s.atoms[i];
            const auto& b = s.atoms[j];
            if (pdb)
                out << "CONECT" << std::setw(5) << a.serial << std::setw(5) << b.serial << '\n';
            else {
                if (!header) {
                    out << "#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n";
                    for (int partner = 1; partner <= 2; ++partner) {
                        const auto prefix = "_struct_conn.ptnr" + std::to_string(partner) + "_";
                        out << prefix << "auth_asym_id\n"
                            << prefix << "auth_seq_id\n"
                            << prefix << "auth_atom_id\n"
                            << "_struct_conn.pdbx_ptnr" << partner << "_PDB_ins_code\n";
                    }
                    header = true;
                }
                out << "bond" << i << '_' << j << " covale ";
                for (auto atom : {&a, &b})
                    out << cif(atom->chain) << ' ' << atom->residueNumber << ' ' << cif(atom->name) << ' '
                        << cif(atom->insertionCode == ' ' ? "" : std::string(1, atom->insertionCode)) << ' ';
                out << '\n';
            }
        }
    out << (pdb ? "END\n" : "#\n");
    std::ofstream file(path);
    if (!file || !(file << out.str()))
        throw std::runtime_error("cannot write structure export");
}
} // namespace ascii_protein
