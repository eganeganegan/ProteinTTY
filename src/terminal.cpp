#include "ascii_protein/terminal.hpp"

#include <cerrno>
#include <poll.h>
#include <stdexcept>
#include <string>
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

namespace ascii_protein {

struct Terminal::State { termios original{}; };

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
    pollfd descriptor{STDIN_FILENO, POLLIN, 0};
    const int result = poll(&descriptor, 1, timeoutMilliseconds);
    if (result <= 0 || !(descriptor.revents & POLLIN)) return {};
    char bytes[16]{};
    const auto count = ::read(STDIN_FILENO, bytes, sizeof(bytes));
    if (count <= 0) return {};
    if (bytes[0] == '\x1b') {
        if (count >= 3 && bytes[1] == '[') {
            if (bytes[2] == 'A') return {Key::up};
            if (bytes[2] == 'B') return {Key::down};
            if (bytes[2] == 'C') return {Key::right};
            if (bytes[2] == 'D') return {Key::left};
        }
        return {Key::escape};
    }
    return {Key::character, bytes[0]};
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