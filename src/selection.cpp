#include "ascii_protein/selection.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace ascii_protein {
namespace {

struct Token {
    std::string text;
    std::size_t column{};
};

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return value;
}

std::vector<Token> tokenize(std::string_view source) {
    std::vector<Token> tokens;
    for (std::size_t position = 0; position < source.size();) {
        while (position < source.size() && std::isspace(static_cast<unsigned char>(source[position]))) ++position;
        if (position == source.size()) break;
        const std::size_t column = position + 1;
        if (source[position] == '(' || source[position] == ')') {
            tokens.push_back({std::string(1, source[position++]), column});
            continue;
        }
        std::string text;
        if (source[position] == '\'' || source[position] == '"') {
            const char quote = source[position++];
            while (position < source.size() && source[position] != quote) {
                if (source[position] == '\\' && position + 1 < source.size()) ++position;
                text.push_back(source[position++]);
            }
            if (position == source.size()) {
                throw std::runtime_error("column " + std::to_string(column) + ": unterminated quoted string");
            }
            ++position;
        } else {
            while (position < source.size() && !std::isspace(static_cast<unsigned char>(source[position])) &&
                   source[position] != '(' && source[position] != ')') {
                text.push_back(source[position++]);
            }
        }
        tokens.push_back({std::move(text), column});
    }
    return tokens;
}

class Parser {
public:
    Parser(std::string_view source, const Structure& structure, const NamedSelections& named)
        : tokens_(tokenize(source)), structure_(structure), named_(named) {}

    AtomMask parse() {
        if (tokens_.empty()) throw std::runtime_error("column 1: expected a selection expression");
        auto result = parseOr();
        if (position_ != tokens_.size()) fail("unexpected token: " + current().text);
        return result;
    }

private:
    const Token& current() const { return tokens_[position_]; }

    [[noreturn]] void fail(const std::string& message) const {
        const std::size_t column = position_ < tokens_.size() ? current().column : tokens_.back().column;
        throw std::runtime_error("column " + std::to_string(column) + ": " + message);
    }

    bool accept(std::string_view text) {
        if (position_ < tokens_.size() && lower(current().text) == text) {
            ++position_;
            return true;
        }
        return false;
    }

    AtomMask parseOr() {
        auto result = parseAnd();
        while (accept("or")) combine(result, parseAnd(), false);
        return result;
    }

    AtomMask parseAnd() {
        auto result = parseNot();
        while (accept("and")) combine(result, parseNot(), true);
        return result;
    }

    AtomMask parseNot() {
        if (!accept("not")) return parsePrimary();
        auto result = parseNot();
        for (std::size_t index = 0; index < result.size(); ++index) result[index] = !result[index];
        return result;
    }

    AtomMask parsePrimary() {
        if (accept("(")) {
            auto result = parseOr();
            if (!accept(")")) fail("expected ')'");
            return result;
        }
        if (position_ == tokens_.size()) fail("expected a selection expression");
        const std::string predicate = lower(current().text);
        ++position_;
        if (predicate == "all") return AtomMask(structure_.atoms.size(), true);
        if (predicate == "none") return AtomMask(structure_.atoms.size(), false);
        if (predicate == "hetero") return matching([](const Atom& atom) { return atom.hetero; });
        if (predicate == "solvent") {
            return matching([](const Atom& atom) {
                const std::string name = lower(atom.residueName);
                  return name == "hoh" || name == "wat" || name == "h2o" || name == "dod" || name == "sol" ||
                      name == "tip3" || name == "tip3p" || name == "spc" || name == "spce";
            });
        }
        if (predicate == "chain" || predicate == "resi" || predicate == "resn" || predicate == "name" ||
            predicate == "element" || predicate == "serial") {
            if (position_ == tokens_.size()) fail("expected a value after " + predicate);
            const std::string value = current().text;
            ++position_;
            return predicateMask(predicate, value);
        }
        const auto found = named_.find(predicate);
        if (found != named_.end()) {
            if (found->second.size() != structure_.atoms.size()) fail("selection is for a different structure");
            return found->second;
        }
        --position_;
        fail("unknown selection term: " + current().text);
    }

    AtomMask predicateMask(const std::string& predicate, const std::string& value) const {
        if (predicate == "resi") {
            const auto dash = value.find('-', 1);
            try {
                const int first = integer(value.substr(0, dash));
                const int last = dash == std::string::npos ? first : integer(value.substr(dash + 1));
                return matching([&](const Atom& atom) {
                    return atom.residueNumber >= std::min(first, last) && atom.residueNumber <= std::max(first, last);
                });
            } catch (...) {
                fail("invalid residue range: " + value);
            }
        }
        if (predicate == "serial") {
            try {
                const int serial = integer(value);
                return matching([&](const Atom& atom) { return atom.serial == serial; });
            } catch (...) {
                fail("invalid atom serial: " + value);
            }
        }
        const std::string expected = lower(value);
        return matching([&](const Atom& atom) {
            if (predicate == "chain") return lower(atom.chain) == expected;
            if (predicate == "name") return lower(atom.name) == expected;
            if (predicate == "element") return lower(atom.element) == expected;
            if (predicate == "resn") {
                const std::string residue = lower(atom.residueName);
                std::size_t start = 0;
                while (start <= expected.size()) {
                    const auto separator = expected.find('+', start);
                    if (residue == expected.substr(start, separator - start)) return true;
                    if (separator == std::string::npos) break;
                    start = separator + 1;
                }
            }
            return false;
        });
    }

    static int integer(const std::string& value) {
        std::size_t used = 0;
        const int result = std::stoi(value, &used);
        if (used != value.size()) throw std::invalid_argument("trailing characters");
        return result;
    }

    template <typename Predicate>
    AtomMask matching(Predicate predicate) const {
        AtomMask result(structure_.atoms.size());
        for (std::size_t index = 0; index < structure_.atoms.size(); ++index) {
            result[index] = predicate(structure_.atoms[index]);
        }
        return result;
    }

    static void combine(AtomMask& left, const AtomMask& right, bool conjunction) {
        for (std::size_t index = 0; index < left.size(); ++index) {
            left[index] = conjunction ? left[index] && right[index] : left[index] || right[index];
        }
    }

    std::vector<Token> tokens_;
    const Structure& structure_;
    const NamedSelections& named_;
    std::size_t position_{};
};

} // namespace

AtomMask evaluateSelection(std::string_view expression, const Structure& structure, const NamedSelections& named) {
    return Parser(expression, structure, named).parse();
}

std::size_t selectionSize(const AtomMask& selection) {
    return static_cast<std::size_t>(std::count(selection.begin(), selection.end(), true));
}

} // namespace ascii_protein