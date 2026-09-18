#include "ascii_protein/terminal.hpp"
#include "ascii_protein/terminal_input.hpp"

#if defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace ascii_protein {
namespace {

void writeConsole(HANDLE output, std::string_view text) {
    const char* data = text.data();
    std::size_t remaining = text.size();
    while (remaining > 0) {
        const DWORD chunk = static_cast<DWORD>(std::min<std::size_t>(remaining, 0xffffffffUL));
        DWORD written = 0;
        if (!WriteFile(output, data, chunk, &written, nullptr) || written == 0) break;
        data += written;
        remaining -= written;
    }
}

} // namespace

struct Terminal::State {
    HANDLE input{INVALID_HANDLE_VALUE};
    HANDLE output{INVALID_HANDLE_VALUE};
    DWORD originalInputMode{};
    DWORD originalOutputMode{};
    UINT originalInputCodePage{};
    UINT originalOutputCodePage{};
    mutable TerminalInputDecoder decoder;
    mutable WCHAR highSurrogate{};
};

bool Terminal::interactive() {
    const HANDLE input = GetStdHandle(STD_INPUT_HANDLE);
    const HANDLE output = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    return input != INVALID_HANDLE_VALUE && output != INVALID_HANDLE_VALUE && GetConsoleMode(input, &mode) &&
           GetConsoleMode(output, &mode);
}

Terminal::Terminal() {
    if (!interactive()) return;
    state_ = new State;
    state_->input = GetStdHandle(STD_INPUT_HANDLE);
    state_->output = GetStdHandle(STD_OUTPUT_HANDLE);
    if (!GetConsoleMode(state_->input, &state_->originalInputMode) ||
        !GetConsoleMode(state_->output, &state_->originalOutputMode)) {
        delete state_;
        state_ = nullptr;
        throw std::runtime_error("failed to read Windows console settings");
    }
    state_->originalInputCodePage = GetConsoleCP();
    state_->originalOutputCodePage = GetConsoleOutputCP();
    const DWORD inputMode = (state_->originalInputMode | ENABLE_PROCESSED_INPUT | ENABLE_WINDOW_INPUT |
                             ENABLE_EXTENDED_FLAGS) &
                            ~(ENABLE_ECHO_INPUT | ENABLE_LINE_INPUT | ENABLE_QUICK_EDIT_MODE);
    const DWORD outputMode = state_->originalOutputMode | ENABLE_PROCESSED_OUTPUT | ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    if (!SetConsoleMode(state_->input, inputMode) || !SetConsoleMode(state_->output, outputMode)) {
        SetConsoleMode(state_->input, state_->originalInputMode);
        SetConsoleMode(state_->output, state_->originalOutputMode);
        delete state_;
        state_ = nullptr;
        throw std::runtime_error("Windows 10 virtual terminal support is required");
    }
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
    active_ = true;
    writeConsole(state_->output, "\x1b[?1049h\x1b[?25l\x1b[2J");
}

Terminal::~Terminal() {
    if (!active_) return;
    writeConsole(state_->output, "\x1b[0m\x1b[?25h\x1b[?1049l");
    SetConsoleMode(state_->input, state_->originalInputMode);
    SetConsoleMode(state_->output, state_->originalOutputMode);
    SetConsoleCP(state_->originalInputCodePage);
    SetConsoleOutputCP(state_->originalOutputCodePage);
    delete state_;
}

TerminalSize Terminal::size() const {
    CONSOLE_SCREEN_BUFFER_INFO information{};
    if (state_ && GetConsoleScreenBufferInfo(state_->output, &information)) {
        const auto& window = information.srWindow;
        return {static_cast<int>(window.Right - window.Left + 1), static_cast<int>(window.Bottom - window.Top + 1)};
    }
    return {};
}

KeyEvent Terminal::readKey(int timeoutMilliseconds) const {
    if (!state_) return {};
    if (const auto event = state_->decoder.next()) return *event;
    const DWORD timeout = timeoutMilliseconds < 0 ? INFINITE : static_cast<DWORD>(timeoutMilliseconds);
    while (WaitForSingleObject(state_->input, timeout) == WAIT_OBJECT_0) {
        INPUT_RECORD record{};
        DWORD count = 0;
        if (!ReadConsoleInputW(state_->input, &record, 1, &count) || count == 0) return {};
        if (record.EventType == WINDOW_BUFFER_SIZE_EVENT) return {Key::resize};
        if (record.EventType != KEY_EVENT || !record.Event.KeyEvent.bKeyDown) continue;
        const auto& key = record.Event.KeyEvent;
        KeyEvent event;
        if (key.wVirtualKeyCode == VK_LEFT) event = {Key::left};
        else if (key.wVirtualKeyCode == VK_RIGHT) event = {Key::right};
        else if (key.wVirtualKeyCode == VK_UP) event = {Key::up};
        else if (key.wVirtualKeyCode == VK_DOWN) event = {Key::down};
        else if (key.wVirtualKeyCode == VK_ESCAPE) event = {Key::escape};
        else if (key.uChar.UnicodeChar >= 0xd800 && key.uChar.UnicodeChar <= 0xdbff) {
            state_->highSurrogate = key.uChar.UnicodeChar;
            continue;
        } else if (key.uChar.UnicodeChar >= 0xdc00 && key.uChar.UnicodeChar <= 0xdfff && state_->highSurrogate) {
            const char32_t codepoint = 0x10000 +
                ((static_cast<char32_t>(state_->highSurrogate) - 0xd800) << 10) +
                (static_cast<char32_t>(key.uChar.UnicodeChar) - 0xdc00);
            state_->highSurrogate = 0;
            state_->decoder.appendCodepoint(codepoint, key.wRepeatCount);
            return state_->decoder.next().value_or(KeyEvent{});
        } else if (key.uChar.UnicodeChar >= 0xdc00 && key.uChar.UnicodeChar <= 0xdfff) {
            state_->highSurrogate = 0;
            continue;
        } else if (key.uChar.UnicodeChar > 0) {
            state_->highSurrogate = 0;
            state_->decoder.appendCodepoint(key.uChar.UnicodeChar, key.wRepeatCount);
            return state_->decoder.next().value_or(KeyEvent{});
        } else {
            continue;
        }
        state_->decoder.appendEvent(event, key.wRepeatCount);
        return state_->decoder.next().value_or(KeyEvent{});
    }
    return {};
}

void Terminal::present(std::string_view frame) const {
    const HANDLE output = state_ ? state_->output : GetStdHandle(STD_OUTPUT_HANDLE);
    if (output == INVALID_HANDLE_VALUE || output == nullptr) return;
    if (state_) writeConsole(output, "\x1b[H");
    writeConsole(output, frame);
}

} // namespace ascii_protein

#else

#include <cerrno>
#include <poll.h>
#include <stdexcept>
#include <string>
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

namespace ascii_protein {

struct Terminal::State {
    termios original{};
    mutable TerminalInputDecoder decoder;
};

bool Terminal::interactive() { return isatty(STDIN_FILENO) && isatty(STDOUT_FILENO); }

Terminal::Terminal() {
    if (!interactive()) return;
    state_ = new State;
    if (tcgetattr(STDIN_FILENO, &state_->original) != 0) {
        delete state_;
        state_ = nullptr;
        throw std::runtime_error("failed to read terminal settings");
    }
    termios raw = state_->original;
    raw.c_lflag &= static_cast<tcflag_t>(~(ICANON | ECHO));
    raw.c_iflag &= static_cast<tcflag_t>(~(IXON | ICRNL));
    raw.c_cc[VMIN] = 0;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw) != 0) {
        delete state_;
        state_ = nullptr;
        throw std::runtime_error("failed to enter raw terminal mode");
    }
    active_ = true;
    const std::string sequence = "\x1b[?1049h\x1b[?25l\x1b[2J";
    (void)::write(STDOUT_FILENO, sequence.data(), sequence.size());
}

Terminal::~Terminal() {
    if (!active_) return;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &state_->original);
    const std::string sequence = "\x1b[0m\x1b[?25h\x1b[?1049l";
    (void)::write(STDOUT_FILENO, sequence.data(), sequence.size());
    delete state_;
}

TerminalSize Terminal::size() const {
    winsize dimensions{};
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &dimensions) == 0 && dimensions.ws_col > 0 && dimensions.ws_row > 0) {
        return {dimensions.ws_col, dimensions.ws_row};
    }
    return {};
}

KeyEvent Terminal::readKey(int timeoutMilliseconds) const {
    if (!state_) return {};
    if (const auto event = state_->decoder.next()) return *event;
    pollfd descriptor{STDIN_FILENO, POLLIN, 0};
    const int result = poll(&descriptor, 1, timeoutMilliseconds);
    if (result <= 0 || !(descriptor.revents & POLLIN)) {
        return state_->decoder.next(true).value_or(KeyEvent{});
    }
    char bytes[64]{};
    const auto count = ::read(STDIN_FILENO, bytes, sizeof(bytes));
    if (count <= 0) return {};
    state_->decoder.appendBytes(std::string_view(bytes, static_cast<std::size_t>(count)));
    if (const auto event = state_->decoder.next()) return *event;

    descriptor.revents = 0;
    if (poll(&descriptor, 1, 5) > 0 && (descriptor.revents & POLLIN)) {
        const auto continuation = ::read(STDIN_FILENO, bytes, sizeof(bytes));
        if (continuation > 0) {
            state_->decoder.appendBytes(std::string_view(bytes, static_cast<std::size_t>(continuation)));
        }
    }
    return state_->decoder.next(true).value_or(KeyEvent{});
}

void Terminal::present(std::string_view frame) const {
    const std::string prefix = "\x1b[H";
    (void)::write(STDOUT_FILENO, prefix.data(), prefix.size());
    const char* data = frame.data();
    std::size_t remaining = frame.size();
    while (remaining > 0) {
        const auto written = ::write(STDOUT_FILENO, data, remaining);
        if (written > 0) {
            data += written;
            remaining -= static_cast<std::size_t>(written);
        } else if (errno != EINTR) break;
    }
}

} // namespace ascii_protein

#endif