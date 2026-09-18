#pragma once

#include "ascii_protein/terminal.hpp"

#include <cstddef>
#include <deque>
#include <optional>
#include <string>
#include <string_view>

namespace ascii_protein {

class TerminalInputDecoder {
public:
    void appendBytes(std::string_view bytes);
    void appendEvent(KeyEvent event, std::size_t repeat = 1);
    void appendCodepoint(char32_t codepoint, std::size_t repeat = 1);
    std::optional<KeyEvent> next(bool flushIncompleteEscape = false);

private:
    std::string bytes_;
    std::deque<KeyEvent> events_;
};

} // namespace ascii_protein
