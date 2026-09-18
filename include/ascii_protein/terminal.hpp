#pragma once

#include <string_view>

namespace ascii_protein {

struct TerminalSize {
    int width{100};
    int height{32};
};

enum class Key { none, character, left, right, up, down, escape, resize };

struct KeyEvent {
    Key key{Key::none};
    char character{};
};

class Terminal {
public:
    Terminal();
    ~Terminal();
    Terminal(const Terminal&) = delete;
    Terminal& operator=(const Terminal&) = delete;

    static bool interactive();
    TerminalSize size() const;
    KeyEvent readKey(int timeoutMilliseconds) const;
    void present(std::string_view frame) const;

private:
    bool active_{};
    struct State;
    State* state_{};
};

} // namespace ascii_protein