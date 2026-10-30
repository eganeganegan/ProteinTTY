#include "ascii_protein/binarycif_parser.hpp"
#include "ascii_protein/molecular.hpp"
#include <sstream>
#include <iomanip>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

namespace ascii_protein {
namespace {

using Bytes = std::vector<std::uint8_t>;

constexpr std::size_t maxDocumentBytes = 512U * 1024U * 1024U;
constexpr std::size_t maxContainerEntries = 1'000'000U;
constexpr std::size_t maxDecodedValues = 10'000'000U;
constexpr std::size_t maxEncodingCount = 16U;

struct Value {
    using Array = std::vector<Value>;
    using Map = std::unordered_map<std::string, Value>;
    std::variant<std::monostate, bool, std::int64_t, double, std::string, Bytes, Array, Map> data;
};

class MessagePackReader {
public:
    explicit MessagePackReader(Bytes bytes) : bytes_(std::move(bytes)) {}

    Value read() {
        Value result = value(0);
        if (position_ != bytes_.size()) fail("trailing data");
        return result;
    }

private:
    [[noreturn]] void fail(std::string_view message) const {
        throw std::runtime_error("invalid BinaryCIF MessagePack at byte " + std::to_string(position_) + ": " +
                                 std::string(message));
    }

    std::uint8_t byte() {
        if (position_ >= bytes_.size()) fail("unexpected end of input");
        return bytes_[position_++];
    }

    template <typename T>
    T bigEndian() {
        static_assert(std::is_unsigned_v<T>);
        T result{};
        for (std::size_t index = 0; index < sizeof(T); ++index) result = static_cast<T>((result << 8) | byte());
        return result;
    }

    std::string string(std::size_t size) {
        if (size > maxDocumentBytes) fail("string size limit exceeded");
        if (size > bytes_.size() - position_) fail("truncated string");
        std::string result(reinterpret_cast<const char*>(bytes_.data() + position_), size);
        position_ += size;
        return result;
    }

    Bytes binary(std::size_t size) {
        if (size > maxDocumentBytes) fail("binary size limit exceeded");
        if (size > bytes_.size() - position_) fail("truncated binary value");
        Bytes result(bytes_.begin() + static_cast<std::ptrdiff_t>(position_),
                     bytes_.begin() + static_cast<std::ptrdiff_t>(position_ + size));
        position_ += size;
        return result;
    }

    Value::Array array(std::size_t size, int depth) {
        if (size > maxContainerEntries) fail("array size limit exceeded");
        Value::Array result;
        result.reserve(size);
        for (std::size_t index = 0; index < size; ++index) result.push_back(value(depth + 1));
        return result;
    }

    Value::Map map(std::size_t size, int depth) {
        if (size > maxContainerEntries) fail("map size limit exceeded");
        Value::Map result;
        for (std::size_t index = 0; index < size; ++index) {
            Value key = value(depth + 1);
            const auto* text = std::get_if<std::string>(&key.data);
            if (!text) fail("map key is not a string");
            result.emplace(*text, value(depth + 1));
        }
        return result;
    }

    Value value(int depth) {
        if (depth > 64) fail("nesting limit exceeded");
        const std::uint8_t marker = byte();
        if (marker <= 0x7f) return {static_cast<std::int64_t>(marker)};
        if (marker >= 0xe0) return {static_cast<std::int64_t>(static_cast<std::int8_t>(marker))};
        if ((marker & 0xe0) == 0xa0) return {string(marker & 0x1f)};
        if ((marker & 0xf0) == 0x90) return {array(marker & 0x0f, depth)};
        if ((marker & 0xf0) == 0x80) return {map(marker & 0x0f, depth)};
        switch (marker) {
        case 0xc0: return {};
        case 0xc2: return {false};
        case 0xc3: return {true};
        case 0xc4: return {binary(bigEndian<std::uint8_t>())};
        case 0xc5: return {binary(bigEndian<std::uint16_t>())};
        case 0xc6: return {binary(bigEndian<std::uint32_t>())};
        case 0xca: return {static_cast<double>(std::bit_cast<float>(bigEndian<std::uint32_t>()))};
        case 0xcb: return {std::bit_cast<double>(bigEndian<std::uint64_t>())};
        case 0xcc: return {static_cast<std::int64_t>(bigEndian<std::uint8_t>())};
        case 0xcd: return {static_cast<std::int64_t>(bigEndian<std::uint16_t>())};
        case 0xce: return {static_cast<std::int64_t>(bigEndian<std::uint32_t>())};
        case 0xcf: {
            const auto number = bigEndian<std::uint64_t>();
            if (number > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) fail("integer overflow");
            return {static_cast<std::int64_t>(number)};
        }
        case 0xd0: return {static_cast<std::int64_t>(static_cast<std::int8_t>(byte()))};
        case 0xd1: return {static_cast<std::int64_t>(static_cast<std::int16_t>(bigEndian<std::uint16_t>()))};
        case 0xd2: return {static_cast<std::int64_t>(static_cast<std::int32_t>(bigEndian<std::uint32_t>()))};
        case 0xd3: return {static_cast<std::int64_t>(bigEndian<std::uint64_t>())};
        case 0xd9: return {string(bigEndian<std::uint8_t>())};
        case 0xda: return {string(bigEndian<std::uint16_t>())};
        case 0xdb: return {string(bigEndian<std::uint32_t>())};
        case 0xdc: return {array(bigEndian<std::uint16_t>(), depth)};
        case 0xdd: return {array(bigEndian<std::uint32_t>(), depth)};
        case 0xde: return {map(bigEndian<std::uint16_t>(), depth)};
        case 0xdf: return {map(bigEndian<std::uint32_t>(), depth)};
        default: fail("unsupported value marker");
        }
    }

    Bytes bytes_;
    std::size_t position_{};
};

const Value::Map& asMap(const Value& value, std::string_view context) {
    const auto* result = std::get_if<Value::Map>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected map");
    return *result;
}

const Value::Array& asArray(const Value& value, std::string_view context) {
    const auto* result = std::get_if<Value::Array>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected array");
    return *result;
}

const Value& member(const Value::Map& map, std::string_view name, std::string_view context) {
    const auto found = map.find(std::string(name));
    if (found == map.end()) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": missing " + std::string(name));
    return found->second;
}

const std::string& asString(const Value& value, std::string_view context) {
    const auto* result = std::get_if<std::string>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected string");
    return *result;
}

const Bytes& asBytes(const Value& value, std::string_view context) {
    const auto* result = std::get_if<Bytes>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected bytes");
    return *result;
}

std::int64_t asInteger(const Value& value, std::string_view context) {
    const auto* result = std::get_if<std::int64_t>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected integer");
    return *result;
}

double asNumber(const Value& value, std::string_view context) {
    if (const auto* integer = std::get_if<std::int64_t>(&value.data)) return static_cast<double>(*integer);
    if (const auto* number = std::get_if<double>(&value.data)) return *number;
    throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected number");
}

bool asBool(const Value& value, std::string_view context) {
    const auto* result = std::get_if<bool>(&value.data);
    if (!result) throw std::runtime_error("invalid BinaryCIF " + std::string(context) + ": expected boolean");
    return *result;
}

template <typename T>
T littleEndian(const std::uint8_t* bytes) {
    using Bits = std::conditional_t<sizeof(T) == 8, std::uint64_t,
                 std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint16_t>>;
    Bits bits{};
    for (std::size_t index = 0; index < sizeof(T); ++index) bits |= static_cast<Bits>(bytes[index]) << (8 * index);
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) return std::bit_cast<T>(bits);
    return static_cast<T>(bits);
}

std::vector<double> decodeByteArray(const Bytes& bytes, int type) {
    std::vector<double> result;
    auto decode = [&](auto tag) {
        using T = decltype(tag);
        if (bytes.size() % sizeof(T) != 0) throw std::runtime_error("invalid BinaryCIF ByteArray length");
        const std::size_t count = bytes.size() / sizeof(T);
        if (count > maxDecodedValues) throw std::runtime_error("BinaryCIF decoded value limit exceeded");
        result.reserve(count);
        for (std::size_t offset = 0; offset < bytes.size(); offset += sizeof(T)) {
            if constexpr (sizeof(T) == 1) result.push_back(static_cast<double>(static_cast<T>(bytes[offset])));
            else result.push_back(static_cast<double>(littleEndian<T>(bytes.data() + offset)));
        }
    };
    switch (type) {
    case 1: decode(std::int8_t{}); break;
    case 2: decode(std::int16_t{}); break;
    case 3: decode(std::int32_t{}); break;
    case 4: decode(std::uint8_t{}); break;
    case 5: decode(std::uint16_t{}); break;
    case 6: decode(std::uint32_t{}); break;
    case 32: decode(float{}); break;
    case 33: decode(double{}); break;
    default: throw std::runtime_error("unsupported BinaryCIF ByteArray type " + std::to_string(type));
    }
    return result;
}

std::vector<double> decodeNumeric(const Bytes& bytes, const Value::Array& encodings);

std::vector<double> applyNumericEncoding(std::vector<double> values, const Value::Map& encoding) {
    const std::string& kind = asString(member(encoding, "kind", "encoding"), "encoding kind");
    if (kind == "IntegerPacking") {
        const int byteCount = static_cast<int>(asInteger(member(encoding, "byteCount", kind), kind));
        const bool isUnsigned = asBool(member(encoding, "isUnsigned", kind), kind);
        if (byteCount <= 0 || byteCount > 4) throw std::runtime_error("invalid BinaryCIF IntegerPacking byte count");
        const std::int64_t upper = isUnsigned ? (std::int64_t{1} << (byteCount * 8)) - 1
                                              : (std::int64_t{1} << (byteCount * 8 - 1)) - 1;
        const std::int64_t lower = isUnsigned ? 0 : -upper - 1;
        std::vector<double> unpacked;
        unpacked.reserve(std::min(values.size(), maxDecodedValues));
        std::int64_t accumulator = 0;
        for (const double raw : values) {
            const auto current = static_cast<std::int64_t>(raw);
            accumulator += current;
            if (current != upper && (isUnsigned || current != lower)) {
                if (unpacked.size() == maxDecodedValues) {
                    throw std::runtime_error("BinaryCIF decoded value limit exceeded");
                }
                unpacked.push_back(static_cast<double>(accumulator));
                accumulator = 0;
            }
        }
        if (accumulator != 0) throw std::runtime_error("invalid BinaryCIF IntegerPacking sequence");
        return unpacked;
    }
    if (kind == "RunLength") {
        if (values.size() % 2 != 0) throw std::runtime_error("invalid BinaryCIF RunLength sequence");
        std::vector<double> expanded;
        for (std::size_t index = 0; index < values.size(); index += 2) {
            const auto count = static_cast<std::int64_t>(values[index + 1]);
            if (count < 0 || static_cast<std::uint64_t>(count) > maxDecodedValues - expanded.size()) {
                throw std::runtime_error("BinaryCIF decoded value limit exceeded");
            }
            expanded.insert(expanded.end(), static_cast<std::size_t>(count), values[index]);
        }
        return expanded;
    }
    if (kind == "Delta") {
        double current = asNumber(member(encoding, "origin", kind), kind);
        for (double& value : values) {
            current += value;
            value = current;
        }
        return values;
    }
    if (kind == "FixedPoint") {
        const double factor = asNumber(member(encoding, "factor", kind), kind);
        if (factor == 0.0) throw std::runtime_error("invalid BinaryCIF FixedPoint factor");
        for (double& value : values) value /= factor;
        return values;
    }
    if (kind == "IntervalQuantization") {
        const double minimum = asNumber(member(encoding, "min", kind), kind);
        const double maximum = asNumber(member(encoding, "max", kind), kind);
        const auto steps = asInteger(member(encoding, "numSteps", kind), kind);
        if (steps <= 1) throw std::runtime_error("invalid BinaryCIF quantization steps");
        const double scale = (maximum - minimum) / static_cast<double>(steps - 1);
        for (double& value : values) value = minimum + value * scale;
        return values;
    }
    if (kind != "ByteArray") throw std::runtime_error("unsupported BinaryCIF encoding: " + kind);
    return values;
}

std::vector<double> decodeNumeric(const Bytes& bytes, const Value::Array& encodings) {
    if (encodings.empty()) throw std::runtime_error("BinaryCIF column has no encoding");
    if (encodings.size() > maxEncodingCount) throw std::runtime_error("BinaryCIF encoding count limit exceeded");
    const auto& byteEncoding = asMap(encodings.back(), "encoding");
    const std::string& kind = asString(member(byteEncoding, "kind", "encoding"), "encoding kind");
    if (kind != "ByteArray") throw std::runtime_error("BinaryCIF encoding chain does not end in ByteArray");
    auto values = decodeByteArray(bytes, static_cast<int>(asInteger(member(byteEncoding, "type", kind), kind)));
    for (std::size_t index = encodings.size() - 1; index-- > 0;) {
        values = applyNumericEncoding(std::move(values), asMap(encodings[index], "encoding"));
    }
    return values;
}

struct Column {
    std::vector<double> numbers;
    std::vector<std::string> strings;
    std::vector<double> mask;
    bool textual{};

    std::string text(std::size_t row) const {
        if (row >= size() || (row < mask.size() && mask[row] != 0.0)) return {};
        if (textual) return strings[row];
        const double number = numbers[row];
        if (number == static_cast<double>(static_cast<std::int64_t>(number))) {
            return std::to_string(static_cast<std::int64_t>(number));
        }
        std::ostringstream out;
        out << std::setprecision(17) << number;
        return out.str();
    }

    std::optional<double> number(std::size_t row) const {
        if (row >= size() || textual || (row < mask.size() && mask[row] != 0.0)) return std::nullopt;
        return numbers[row];
    }

    std::size_t size() const { return textual ? strings.size() : numbers.size(); }
};

Column decodeColumn(const Value::Map& columnMap) {
    Column column;
    const auto& dataMap = asMap(member(columnMap, "data", "column"), "column data");
    const auto& encodings = asArray(member(dataMap, "encoding", "column data"), "column encodings");
    const Bytes& bytes = asBytes(member(dataMap, "data", "column data"), "column bytes");
    if (!encodings.empty()) {
        const auto& first = asMap(encodings.front(), "encoding");
        const std::string& kind = asString(member(first, "kind", "encoding"), "encoding kind");
        if (kind == "StringArray") {
            const auto& dataEncodings = asArray(member(first, "dataEncoding", kind), "StringArray dataEncoding");
            const auto indices = decodeNumeric(bytes, dataEncodings);
            const auto& offsetEncodings = asArray(member(first, "offsetEncoding", kind), "StringArray offsetEncoding");
            const auto offsets = decodeNumeric(asBytes(member(first, "offsets", kind), "StringArray offsets"), offsetEncodings);
            const std::string& stringData = asString(member(first, "stringData", kind), "StringArray stringData");
            column.textual = true;
            column.strings.reserve(indices.size());
            for (const double rawIndex : indices) {
                const auto index = static_cast<std::int64_t>(rawIndex);
                if (index < 0) {
                    column.strings.emplace_back();
                    continue;
                }
                const auto offsetIndex = static_cast<std::size_t>(index);
                if (offsetIndex + 1 >= offsets.size()) throw std::runtime_error("invalid BinaryCIF string index");
                const auto start = static_cast<std::size_t>(offsets[offsetIndex]);
                const auto end = static_cast<std::size_t>(offsets[offsetIndex + 1]);
                if (start > end || end > stringData.size()) throw std::runtime_error("invalid BinaryCIF string offset");
                column.strings.push_back(stringData.substr(start, end - start));
            }
        } else {
            column.numbers = decodeNumeric(bytes, encodings);
        }
    }
    const auto maskFound = columnMap.find("mask");
    if (maskFound != columnMap.end() && !std::holds_alternative<std::monostate>(maskFound->second.data)) {
        const auto& maskMap = asMap(maskFound->second, "column mask");
        column.mask = decodeNumeric(asBytes(member(maskMap, "data", "column mask"), "column mask data"),
                                    asArray(member(maskMap, "encoding", "column mask"), "column mask encodings"));
    }
    return column;
}

using Columns = std::unordered_map<std::string, Column>;

std::string text(const Columns& columns, std::size_t row, std::initializer_list<std::string_view> names) {
    for (const auto name : names) {
        const auto found = columns.find(std::string(name));
        if (found != columns.end()) {
            const std::string value = found->second.text(row);
            if (!value.empty()) return value;
        }
    }
    return {};
}

std::optional<double> number(const Columns& columns, std::size_t row, std::initializer_list<std::string_view> names) {
    for (const auto name : names) {
        const auto found = columns.find(std::string(name));
        if (found != columns.end()) {
            if (const auto value = found->second.number(row)) return value;
        }
    }
    return std::nullopt;
}

char marker(const std::string& value) { return value.empty() ? ' ' : value.front(); }

void parseAtom(const Columns& columns, std::size_t row, Structure& structure) {
    const std::string record = text(columns, row, {"group_PDB"});
    if (record != "ATOM" && record != "HETATM") return;
    const auto serial = number(columns, row, {"id"});
    const auto residue = number(columns, row, {"auth_seq_id", "label_seq_id"});
    const auto x = number(columns, row, {"Cartn_x"});
    const auto y = number(columns, row, {"Cartn_y"});
    const auto z = number(columns, row, {"Cartn_z"});
    if (!serial || !residue || !x || !y || !z) {
        structure.warnings.push_back("atom_site row " + std::to_string(row + 1) + ": malformed atom skipped");
        return;
    }
    const auto model = number(columns, row, {"pdbx_PDB_model_num"});

    Atom atom;
    atom.model = model ? static_cast<int>(*model) : 1;
    atom.labelChain = text(columns, row, {"label_asym_id"});
    atom.labelName = text(columns, row, {"label_atom_id"});
    atom.labelResidue = static_cast<int>(number(columns, row, {"label_seq_id"}).value_or(0));
    atom.serial = static_cast<int>(*serial);
    atom.name = text(columns, row, {"auth_atom_id", "label_atom_id"});
    atom.alternateLocation = marker(text(columns, row, {"label_alt_id"}));
    atom.element = text(columns, row, {"type_symbol"});
    atom.residueName = text(columns, row, {"auth_comp_id", "label_comp_id"});
    atom.chain = text(columns, row, {"auth_asym_id", "label_asym_id"});
    atom.residueNumber = static_cast<int>(*residue);
    atom.insertionCode = marker(text(columns, row, {"pdbx_PDB_ins_code"}));
    atom.position = {*x, *y, *z};
    atom.occupancy = number(columns, row, {"occupancy"}).value_or(0.0);
    atom.bFactor = number(columns, row, {"B_iso_or_equiv"}).value_or(0.0);
    atom.hetero = record == "HETATM";
    structure.atoms.push_back(std::move(atom));
}

void parseSecondary(const Columns& columns, std::size_t row, SecondaryStructure type, Structure& structure) {
    if (type == SecondaryStructure::helix && !text(columns, row, {"conf_type_id"}).starts_with("HELX")) return;
    const std::string startChain = text(columns, row, {"beg_auth_asym_id", "beg_label_asym_id"});
    const std::string endChain = text(columns, row, {"end_auth_asym_id", "end_label_asym_id"});
    const auto start = number(columns, row, {"beg_auth_seq_id", "beg_label_seq_id"});
    const auto end = number(columns, row, {"end_auth_seq_id", "end_label_seq_id"});
    if (startChain.empty() || endChain.empty() || !start || !end) return;
    structure.secondaryRanges.push_back(
        {type, startChain, static_cast<int>(*start), endChain, static_cast<int>(*end),
         marker(text(columns, row, {"pdbx_beg_PDB_ins_code"})),
         marker(text(columns, row, {"pdbx_end_PDB_ins_code"}))});
}

} // namespace

Structure BinaryCifParser::parseFile(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open BinaryCIF file: " + path.string());
    return parse(input, path.filename().string());
}

Structure BinaryCifParser::parse(std::istream& input, std::string sourceName) {
    Bytes bytes;
    bytes.reserve(64U * 1024U);
    for (std::istreambuf_iterator<char> current(input), end; current != end; ++current) {
        if (bytes.size() == maxDocumentBytes) throw std::runtime_error("BinaryCIF document size limit exceeded");
        bytes.push_back(static_cast<std::uint8_t>(*current));
    }
    const Value root = MessagePackReader(std::move(bytes)).read();
    const auto& rootMap = asMap(root, "root");
    const auto& blocks = asArray(member(rootMap, "dataBlocks", "root"), "dataBlocks");
    if (blocks.empty()) throw std::runtime_error("BinaryCIF contains no data blocks");

    Structure structure;
    structure.sourceName = std::move(sourceName);
    const auto& block = asMap(blocks.front(), "data block");
    const auto& categories = asArray(member(block, "categories", "data block"), "categories");
    for (const Value& categoryValue : categories) {
        const auto& category = asMap(categoryValue, "category");
        const std::string& name = asString(member(category, "name", "category"), "category name");
        if (name != "_atom_site" && name != "_struct_conf" && name != "_struct_sheet_range" &&
            name != "_struct_conn" && name != "_pdbx_struct_oper_list" && name != "_pdbx_struct_assembly_gen")
            continue;
        const auto rowCount = asInteger(member(category, "rowCount", name), "category rowCount");
        if (rowCount < 0) throw std::runtime_error("invalid BinaryCIF category row count");
        Columns columns;
        for (const Value& columnValue : asArray(member(category, "columns", name), "category columns")) {
            const auto& columnMap = asMap(columnValue, "column");
            const std::string& columnName = asString(member(columnMap, "name", "column"), "column name");
            Column column = decodeColumn(columnMap);
            if (column.size() != static_cast<std::size_t>(rowCount)) {
                throw std::runtime_error("invalid BinaryCIF column length for " + name + "." + columnName +
                                         ": decoded " + std::to_string(column.size()) + ", expected " +
                                         std::to_string(rowCount));
            }
            columns.emplace(columnName, std::move(column));
        }
        for (std::size_t row = 0; row < static_cast<std::size_t>(rowCount); ++row) {
            if (name == "_atom_site") parseAtom(columns, row, structure);
            else if (name == "_struct_conf") parseSecondary(columns, row, SecondaryStructure::helix, structure);
            else if (name == "_struct_sheet_range")
                parseSecondary(columns, row, SecondaryStructure::sheet, structure);
            else {
                MetadataRow metadata;
                for (const auto& [key, column] : columns) {
                    static_cast<void>(column);
                    metadata[key] = text(columns, row, {key});
                }
                readMolecularRow(structure, name, metadata);
            }
        }
    }
    if (structure.atoms.empty()) throw std::runtime_error("no valid atom_site records in " + structure.sourceName);
    structure.finalize();
    return structure;
}

} // namespace ascii_protein
