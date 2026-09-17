#include "ascii_protein/terminal.hpp"

#if defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <algorithm>
#include <deque>
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

std::string utf8(WCHAR character) {
    char bytes[4]{};
    const int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, &character, 1, bytes, sizeof(bytes), nullptr, nullptr);
    return count > 0 ? std::string(bytes, static_cast<std::size_t>(count)) : std::string{};
}

} // namespace

struct Terminal::State {
    HANDLE input{INVALID_HANDLE_VALUE};
    HANDLE output{INVALID_HANDLE_VALUE};
    DWORD originalInputMode{};
    DWORD originalOutputMode{};
    UINT originalInputCodePage{};
    UINT originalOutputCodePage{};
    mutable std::deque<KeyEvent> pending;
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
    if (!state_->pending.empty()) {
        const auto event = state_->pending.front();
        state_->pending.pop_front();
        return event;
    }
    const DWORD timeout = timeoutMilliseconds < 0 ? INFINITE : static_cast<DWORD>(timeoutMilliseconds);
    while (WaitForSingleObject(state_->input, timeout) == WAIT_OBJECT_0) {
        INPUT_RECORD record{};
        DWORD count = 0;
        if (!ReadConsoleInputW(state_->input, &record, 1, &count) || count == 0) return {};
        if (record.EventType != KEY_EVENT || !record.Event.KeyEvent.bKeyDown) continue;
        const auto& key = record.Event.KeyEvent;
        KeyEvent event;
        std::string encoded;
        if (key.wVirtualKeyCode == VK_LEFT) event = {Key::left};
        else if (key.wVirtualKeyCode == VK_RIGHT) event = {Key::right};
        else if (key.wVirtualKeyCode == VK_UP) event = {Key::up};
        else if (key.wVirtualKeyCode == VK_DOWN) event = {Key::down};
        else if (key.wVirtualKeyCode == VK_ESCAPE) event = {Key::escape};
        else if (key.uChar.UnicodeChar > 0) {
            encoded = utf8(key.uChar.UnicodeChar);
            if (encoded.empty()) continue;
        } else {
            continue;
        }
        for (WORD repeat = 0; repeat < key.wRepeatCount; ++repeat) {
            if (encoded.empty()) {
                state_->pending.push_back(event);
            } else {
                for (const char byte : encoded) state_->pending.push_back({Key::character, byte});
            }
        }
        const auto decoded = state_->pending.front();
        state_->pending.pop_front();
        return decoded;
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
    std::string pending;
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
    if (state_->pending.empty()) {
        pollfd descriptor{STDIN_FILENO, POLLIN, 0};
        const int result = poll(&descriptor, 1, timeoutMilliseconds);
        if (result <= 0 || !(descriptor.revents & POLLIN)) return {};
        char bytes[64]{};
        const auto count = ::read(STDIN_FILENO, bytes, sizeof(bytes));
        if (count <= 0) return {};
        state_->pending.append(bytes, static_cast<std::size_t>(count));
    }
    if (state_->pending.front() == '\x1b') {
        if (state_->pending.size() >= 3 && state_->pending[1] == '[') {
            const char direction = state_->pending[2];
            state_->pending.erase(0, 3);
            if (direction == 'A') return {Key::up};
            if (direction == 'B') return {Key::down};
            if (direction == 'C') return {Key::right};
            if (direction == 'D') return {Key::left};
        }
        state_->pending.erase(0, 1);
        return {Key::escape};
    }
    const char character = state_->pending.front();
    state_->pending.erase(0, 1);
    return {Key::character, character};
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