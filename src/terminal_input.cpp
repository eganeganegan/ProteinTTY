#include "ascii_protein/terminal_input.hpp"

#include <stdexcept>

namespace ascii_protein {
namespace {

constexpr std::size_t maxPendingEvents = 1'000'000;

std::string utf8(char32_t codepoint) {
    if (codepoint > 0x10ffff || (codepoint >= 0xd800 && codepoint <= 0xdfff)) {
        throw std::runtime_error("invalid Unicode terminal input");
    }
    std::string result;
    if (codepoint <= 0x7f) {
        result.push_back(static_cast<char>(codepoint));
    } else if (codepoint <= 0x7ff) {
        result.push_back(static_cast<char>(0xc0 | (codepoint >> 6)));
        result.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else if (codepoint <= 0xffff) {
        result.push_back(static_cast<char>(0xe0 | (codepoint >> 12)));
        result.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else {
        result.push_back(static_cast<char>(0xf0 | (codepoint >> 18)));
        result.push_back(static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    }
    return result;
}

} // namespace

void TerminalInputDecoder::appendBytes(std::string_view bytes) {
    if (bytes_.size() > maxPendingEvents - bytes.size()) throw std::runtime_error("terminal input buffer limit exceeded");
    bytes_.append(bytes);
}

void TerminalInputDecoder::appendEvent(KeyEvent event, std::size_t repeat) {
    if (repeat > maxPendingEvents - events_.size()) throw std::runtime_error("terminal input buffer limit exceeded");
    events_.insert(events_.end(), repeat, event);
}

void TerminalInputDecoder::appendCodepoint(char32_t codepoint, std::size_t repeat) {
    const std::string encoded = utf8(codepoint);
    if (repeat > maxPendingEvents / encoded.size() || encoded.size() * repeat > maxPendingEvents - bytes_.size()) {
        throw std::runtime_error("terminal input buffer limit exceeded");
    }
    for (std::size_t index = 0; index < repeat; ++index) bytes_.append(encoded);
}

std::optional<KeyEvent> TerminalInputDecoder::next(bool flushIncompleteEscape) {
    if (!events_.empty()) {
        const KeyEvent event = events_.front();
        events_.pop_front();
        return event;
    }
    if (bytes_.empty()) return std::nullopt;
    if (bytes_.front() == '\x1b') {
        if (bytes_.size() >= 3 && bytes_[1] == '[') {
            const char direction = bytes_[2];
            if (direction == 'A' || direction == 'B' || direction == 'C' || direction == 'D') {
                bytes_.erase(0, 3);
                if (direction == 'A') return KeyEvent{Key::up};
                if (direction == 'B') return KeyEvent{Key::down};
                if (direction == 'C') return KeyEvent{Key::right};
                return KeyEvent{Key::left};
            }
        }
        if (!flushIncompleteEscape && (bytes_.size() == 1 || (bytes_.size() == 2 && bytes_[1] == '['))) {
            return std::nullopt;
        }
        bytes_.erase(0, 1);
        return KeyEvent{Key::escape};
    }
    const char character = bytes_.front();
    bytes_.erase(0, 1);
    return KeyEvent{Key::character, character};
}

} // namespace ascii_protein
