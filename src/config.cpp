#include "ascii_protein/config.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
namespace ascii_protein {
CommandConfig loadCommandConfig(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open configuration: " + path.string());
    CommandConfig config;
    std::string line;
    std::size_t number = 0;
    while (std::getline(input, line)) {
        ++number;
        auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream in(line);
        std::string kind, name, commands;
        if (!(in >> kind >> std::quoted(name) >> std::quoted(commands)) || commands.empty())
            throw std::runtime_error("invalid configuration line " + std::to_string(number));
        in >> std::ws;
        if (!in.eof()) throw std::runtime_error("unexpected configuration text");
        if (kind == "macro")
            config.macros[name] = commands;
        else if (kind == "bind" && name.size() == 1 && name.front() != ':' &&
                 static_cast<unsigned char>(name.front()) >= 32)
            config.bindings[static_cast<unsigned char>(name.front())] = commands;
        else
            throw std::runtime_error("expected macro or single ASCII key binding on line " +
                                     std::to_string(number));
    }
    return config;
}
} // namespace ascii_protein
