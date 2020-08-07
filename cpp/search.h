#ifndef _SEARCH_H
#define _SEARCH_H

#include "base.h"

namespace cube {
    template<u64 _size>
    struct array_u2 {
        std::array<u64, (_size + 31) / 32> a;

        constexpr u64 size() const {
            return _size;
        }

        u64 get(u64 i) const {
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            return (a[j] >> k) & u64(3);
        }

        void set(u64 i, u64 x) {
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            a[j] = (a[j] & ~(u64(3) << k)) | (x << k);
        }

        void fill(u64 x) {
            u64 y = x;
            for (u64 i = 1; i < 6; i++) {
                y = y | (y << (u64(1) << i));
            }
            a.fill(y);
        }
    };

    template<typename _solver>
    u64 set_alt(const _solver &s, array_u2<_solver::n_state> &distance_m3,
                const typename _solver::t_state &a, u64 i, u64 x) {
        auto alt = s.alt(a, i);
        for (u64 j: alt) {
            distance_m3.set(j, x);
        }
        return alt.size();
    }

    template<typename _solver>
    void bfs(const _solver &s, array_u2<_solver::n_state> &distance_m3) {
        constexpr std::array<u64, 3> modify{
                0xffffffffffffffff,
                0xaaaaaaaaaaaaaaaa,
                0x5555555555555555
        };
        constexpr u64 mask = 0x5555555555555555;
        distance_m3.fill(3);
        typename _solver::t_state a_start = s.cube_to_state(_solver::t_cube::i());
        u64 i_start = s.state_to_int(a_start);
        u64 count_start = set_alt<_solver>(s, distance_m3, a_start, i_start, 0);
        std::cout << "bfs: depth=" << 0 << ", count_distinct=" << 1 << ", count=" << count_start << std::endl;
        u64 depth = 1;
        u64 total_count_distinct = 1;
        u64 total_count = count_start;
        std::array<u64, 3> count_m3{count_start, 0, 0};
        auto t0 = std::chrono::steady_clock::now();
        auto t1 = t0;
        while (total_count != _solver::n_state) {
            u64 prev_depth_m3 = (depth - 1) % 3;
            u64 depth_m3 = depth % 3;
            u64 count_distinct = 0;
            u64 count = 0;
            if (count_m3[prev_depth_m3] <= _solver::n_state - total_count) {
                for (u64 i = 0; i < _solver::n_state; i += 32) {
                    u64 x = distance_m3.a[i / 32] ^modify[prev_depth_m3];
                    if (((x >> u64(1)) & x & mask) == 0) {
                        continue;
                    }
                    u64 end = std::min(i + 32, _solver::n_state);
                    for (u64 j = i; j < end; j++) {
                        if (distance_m3.get(j) == prev_depth_m3) {
                            typename _solver::t_state a = s.int_to_state(j);
                            for (const typename _solver::t_state &b: s.adj(a)) {
                                u64 k = s.state_to_int(b);
                                if (distance_m3.get(k) == 3) {
                                    count += set_alt<_solver>(s, distance_m3, b, k, depth_m3);
                                    count_distinct++;
                                }
                            }
                        }
                    }
                }
            } else {
                for (u64 i = 0; i < _solver::n_state; i += 32) {
                    u64 x = distance_m3.a[i / 32];
                    if (((x >> u64(1)) & x & mask) == 0) {
                        continue;
                    }
                    u64 end = std::min(i + 32, _solver::n_state);
                    for (u64 j = i; j < end; j++) {
                        if (distance_m3.get(j) == 3) {
                            typename _solver::t_state a = s.int_to_state(j);
                            for (const typename _solver::t_state &b: s.adj(a)) {
                                u64 k = s.state_to_int(b);
                                if (distance_m3.get(k) == prev_depth_m3) {
                                    count += set_alt<_solver>(s, distance_m3, a, j, depth_m3);
                                    count_distinct++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            auto t2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> d = t2 - t1;
            t1 = t2;
            std::cout << "bfs: depth=" << depth << ", count_distinct=" << count_distinct << ", count=" << count
                      << ", time=" << d.count() << "s" << std::endl;
            depth++;
            total_count_distinct += count_distinct;
            total_count += count;
            count_m3[depth_m3] += count;
        }
        std::chrono::duration<double> d = t1 - t0;
        std::cout << "bfs: total_count_distinct=" << total_count_distinct << ", total_count=" << total_count
                  << ", time=" << d.count() << "s" << std::endl;
    }

    template<typename _solver>
    std::tuple<u64, typename _solver::t_hint> get_distance(const _solver &s, const typename _solver::t_state &a) {
        typename _solver::t_state b = a;
        u64 i = s.state_to_int(b);
        u64 depth = 0;
        while (not s.is_start(b)) {
            u64 target = (s.distance_m3->get(i) + 2) % 3;
            bool found = false;
            for (const typename _solver::t_state &c: s.adj(b)) {
                u64 j = s.state_to_int(c);
                if (s.distance_m3->get(j) == target) {
                    b = c;
                    i = j;
                    depth++;
                    found = true;
                    break;
                }
            }
            assert (found);
        }
        return {depth, depth};
    }

    constexpr u64 computer_distance(u64 distance_m3, u64 distance_adj) {
        return distance_adj + (distance_m3 - distance_adj - 3) % 3 - 1;
    }

    template<typename _solver>
    std::tuple<u64, typename _solver::t_hint> get_distance_hint(
            const _solver &s, const typename _solver::t_state &a, const typename _solver::t_hint &hint) {
        u64 d = computer_distance(s.distance_m3->get(s.state_to_int(a)), hint);
        return {d, d};
    }
}

#endif
