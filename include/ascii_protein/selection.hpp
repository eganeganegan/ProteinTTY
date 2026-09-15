#pragma once

#include "ascii_protein/structure.hpp"

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ascii_protein {

using AtomMask = std::vector<bool>;
using NamedSelections = std::unordered_map<std::string, AtomMask>;

AtomMask evaluateSelection(std::string_view expression, const Structure& structure,
                           const NamedSelections& named = {});
std::size_t selectionSize(const AtomMask& selection);

} // namespace ascii_protein