#include "ascii_protein/trajectory_analysis.hpp"
#include "ascii_protein/trajectory.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace ascii_protein {
namespace {
void validate(const std::vector<std::size_t>& ids, std::size_t count) {
    if (ids.empty()) throw std::runtime_error("analysis selection is empty");
    for (auto id : ids)
        if (id >= count) throw std::runtime_error("analysis atom index out of range");
}
} // namespace
RigidTransform fitCoordinates(const std::vector<Vec3>& moving, const std::vector<Vec3>& reference,
                              const std::vector<std::size_t>& ids) {
    if (moving.size() != reference.size()) throw std::runtime_error("alignment coordinate count mismatch");
    validate(ids, moving.size());
    Vec3 a{}, b{};
    for (auto i : ids) {
        a += moving[i];
        b += reference[i];
    }
    a = a / static_cast<double>(ids.size());
    b = b / static_cast<double>(ids.size());
    double s[3][3]{};
    for (auto i : ids) {
        const auto p = moving[i] - a, q = reference[i] - b;
        const double x[]{p.x, p.y, p.z}, y[]{q.x, q.y, q.z};
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                s[r][c] += x[r] * y[c];
    }
    // Horn's symmetric quaternion matrix. Jacobi diagonalization also handles 180-degree rotations.
    double n[4][4]{{s[0][0] + s[1][1] + s[2][2], s[1][2] - s[2][1], s[2][0] - s[0][2], s[0][1] - s[1][0]},
                   {s[1][2] - s[2][1], s[0][0] - s[1][1] - s[2][2], s[0][1] + s[1][0], s[0][2] + s[2][0]},
                   {s[2][0] - s[0][2], s[0][1] + s[1][0], -s[0][0] + s[1][1] - s[2][2], s[1][2] + s[2][1]},
                   {s[0][1] - s[1][0], s[0][2] + s[2][0], s[1][2] + s[2][1], -s[0][0] - s[1][1] + s[2][2]}};
    double v[4][4]{};
    for (int i = 0; i < 4; ++i)
        v[i][i] = 1;
    for (int sweep = 0; sweep < 80; ++sweep) {
        int p = 0, q = 1;
        double largest = 0, scale = 0;
        for (int i = 0; i < 4; ++i) {
            scale = std::max(scale, std::abs(n[i][i]));
            for (int j = i + 1; j < 4; ++j)
                if (std::abs(n[i][j]) > largest) {
                    largest = std::abs(n[i][j]);
                    p = i;
                    q = j;
                }
        }
        if (largest <= 1e-14 * std::max(1.0, scale)) break;
        const double angle = 0.5 * std::atan2(2 * n[p][q], n[q][q] - n[p][p]);
        const double c = std::cos(angle), t = std::sin(angle);
        const double pp = n[p][p], qq = n[q][q], pq = n[p][q];
        n[p][p] = c * c * pp - 2 * c * t * pq + t * t * qq;
        n[q][q] = t * t * pp + 2 * c * t * pq + c * c * qq;
        n[p][q] = n[q][p] = 0;
        for (int k = 0; k < 4; ++k) {
            if (k != p && k != q) {
                const double kp = n[k][p], kq = n[k][q];
                n[k][p] = n[p][k] = c * kp - t * kq;
                n[k][q] = n[q][k] = t * kp + c * kq;
            }
            const double vp = v[k][p], vq = v[k][q];
            v[k][p] = c * vp - t * vq;
            v[k][q] = t * vp + c * vq;
        }
    }
    int best = 0;
    for (int i = 1; i < 4; ++i)
        if (n[i][i] > n[best][best]) best = i;
    const double w = v[0][best], x = v[1][best], y = v[2][best], z = v[3][best];
    Mat3 r{{{1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)},
            {2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)},
            {2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)}}};
    return {r, b - r * a};
}
std::vector<Vec3> trajectoryCoordinates(const TrajectoryState& state, std::size_t frame) {
    if (!state.source) throw std::runtime_error("no trajectory loaded");
    auto positions = state.source->frame(frame);
    if (!state.alignmentAtoms.empty()) {
        const auto transform = fitCoordinates(positions, state.alignmentReference, state.alignmentAtoms);
        for (auto& p : positions)
            p = transform.apply(p);
    }
    return positions;
}
std::vector<double> trajectoryRmsd(const TrajectoryState& state, const std::vector<std::size_t>& ids,
                                   std::size_t reference) {
    if (!state.source) throw std::runtime_error("no trajectory loaded");
    validate(ids, state.source->atomCount());
    const auto ref = trajectoryCoordinates(state, reference);
    std::vector<double> result;
    for (auto f = state.rangeStart; f <= state.rangeEnd;) {
        const auto positions = trajectoryCoordinates(state, f);
        double sum = 0;
        for (auto i : ids) {
            auto d = positions[i] - ref[i];
            sum += dot(d, d);
        }
        result.push_back(std::sqrt(sum / static_cast<double>(ids.size())));
        if (state.rangeEnd - f < static_cast<std::size_t>(state.stride)) break;
        f += static_cast<std::size_t>(state.stride);
    }
    return result;
}
std::vector<double> trajectoryRmsf(const TrajectoryState& state, const std::vector<std::size_t>& ids) {
    if (!state.source) throw std::runtime_error("no trajectory loaded");
    validate(ids, state.source->atomCount());
    std::vector<Vec3> mean(ids.size());
    std::vector<double> m2(ids.size());
    std::size_t count = 0;
    for (auto f = state.rangeStart; f <= state.rangeEnd;) {
        const auto positions = trajectoryCoordinates(state, f);
        ++count;
        for (std::size_t j = 0; j < ids.size(); ++j) {
            const auto delta = positions[ids[j]] - mean[j];
            mean[j] += delta / static_cast<double>(count);
            m2[j] += dot(delta, positions[ids[j]] - mean[j]);
        }
        if (state.rangeEnd - f < static_cast<std::size_t>(state.stride)) break;
        f += static_cast<std::size_t>(state.stride);
    }
    for (auto& x : m2)
        x = std::sqrt(std::max(0.0, x / static_cast<double>(count)));
    return m2;
}
} // namespace ascii_protein
