#include "ascii_protein/rcsb.hpp"

#include <array>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <spawn.h>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace ascii_protein {
namespace {

std::filesystem::path defaultCacheDirectory() {
    if (const char* xdgCache = std::getenv("XDG_CACHE_HOME")) {
        if (*xdgCache != '\0') return std::filesystem::path(xdgCache) / "ascii-protein";
    }
    if (const char* home = std::getenv("HOME")) {
        if (*home != '\0') return std::filesystem::path(home) / ".cache" / "ascii-protein";
    }
    return std::filesystem::temp_directory_path() / "ascii-protein-cache";
}

void download(const std::string& url, const std::filesystem::path& destination) {
    std::string output = destination.string();
    std::string requestUrl = url;
    std::array<char*, 8> arguments{
        const_cast<char*>("curl"), const_cast<char*>("-fL"), const_cast<char*>("--silent"),
        const_cast<char*>("--show-error"), const_cast<char*>("--output"), output.data(), requestUrl.data(), nullptr};
    pid_t process{};
    const int spawnError = posix_spawnp(&process, "curl", nullptr, nullptr, arguments.data(), environ);
    if (spawnError != 0) throw std::runtime_error("could not start curl to download " + url);

    int status{};
    if (waitpid(process, &status, 0) < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        throw std::runtime_error("RCSB download failed for " + url);
    }
}

} // namespace

std::string normalizePdbId(std::string_view pdbId) {
    if (pdbId.size() != 4) throw std::runtime_error("PDB ID must contain exactly four letters or digits");
    std::string normalized;
    normalized.reserve(4);
    for (const unsigned char character : pdbId) {
        if (!std::isalnum(character)) throw std::runtime_error("PDB ID must contain only letters or digits");
        normalized.push_back(static_cast<char>(std::tolower(character)));
    }
    return normalized;
}

std::filesystem::path fetchFormat(const std::string& pdbId, std::string_view extension,
                                  std::optional<std::filesystem::path> cacheDirectory) {
    const std::string normalized = normalizePdbId(pdbId);
    const auto directory = cacheDirectory.value_or(defaultCacheDirectory());
    std::filesystem::create_directories(directory);
    const auto cached = directory / (normalized + std::string(extension));
    if (std::filesystem::is_regular_file(cached) && std::filesystem::file_size(cached) > 0) return cached;

    const auto temporary = directory / (normalized + std::string(extension) + ".part-" + std::to_string(getpid()));
    const std::string url = "https://files.rcsb.org/download/" + normalized + std::string(extension);
    try {
        download(url, temporary);
        if (!std::filesystem::is_regular_file(temporary) || std::filesystem::file_size(temporary) == 0) {
            throw std::runtime_error("RCSB returned an empty file for " + normalized);
        }
        std::filesystem::rename(temporary, cached);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
    return cached;
}

std::filesystem::path fetchPdb(const std::string& pdbId,
                               std::optional<std::filesystem::path> cacheDirectory) {
    return fetchFormat(pdbId, ".pdb", std::move(cacheDirectory));
}

std::filesystem::path fetchMmcif(const std::string& pdbId,
                                 std::optional<std::filesystem::path> cacheDirectory) {
    return fetchFormat(pdbId, ".cif", std::move(cacheDirectory));
}

} // namespace ascii_protein