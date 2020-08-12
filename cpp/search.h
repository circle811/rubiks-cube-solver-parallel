#ifndef _SEARCH_H
#define _SEARCH_H

#include "base.h"

namespace cube {
    template<u64 _size>
    struct array_u2 {
        std::array<std::atomic_uint64_t, (_size + 31) / 32> a;

        constexpr u64 size() const {
            return _size;
        }

        u64 get(u64 i) const {
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            return (a[j] >> k) & u64(3);
        }

        bool compare_and_set(u64 i, u64 old_x, u64 new_x) {
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            u64 old_y = a[j];
            while (((old_y >> k) & u64(3)) == old_x) {
                u64 new_y = (old_y & ~(u64(3) << k)) | (new_x << k);
                if (a[j].compare_exchange_strong(old_y, new_y)) {
                    return true;
                }
            }
            return false;
        }

        void fill(u64 x) {
            u64 y = x;
            for (u64 i = 1; i < 6; i++) {
                y = y | (y << (u64(1) << i));
            }
            for (u64 i = 0; i < a.size(); i++) {
                a[i] = y;
            }
        }
    };

    std::vector<u64> _split(u64 m, u64 n) {
        u64 mm = (m + 31) / 32;
        u64 div = mm / n;
        u64 mod = mm % n;
        std::vector<u64> sp{0};
        for (u64 i = 0; i < n; i++) {
            sp.push_back(std::min(sp[i] + (div + (i < mod ? 1 : 0)) * 32, m));
        }
        return sp;
    }

    template<u64 _size, typename C>
    u64 _set_multi(array_u2<_size> &distance_m3, const C &i_s, u64 old_x, u64 new_x) {
        auto it = i_s.begin();
        if (not distance_m3.compare_and_set(*it, old_x, new_x)) {
            return 0;
        }
        ++it;
        auto last = i_s.end();
        while (it != last) {
            assert(distance_m3.compare_and_set(*it, old_x, new_x));
            ++it;
        }
        return i_s.size();
    }

    template<typename _solver>
    std::tuple<u64, u64> _forward(const _solver &s, array_u2<_solver::n_state> &distance_m3,
                                  u64 start, u64 end, u64 prev_depth_m3, u64 depth_m3) {
        constexpr std::array<u64, 3> modify{
                0xffffffffffffffff,
                0xaaaaaaaaaaaaaaaa,
                0x5555555555555555
        };
        constexpr u64 mask = 0x5555555555555555;
        u64 count_distinct = 0;
        u64 count = 0;
        for (u64 i = start; i < end; i += 32) {
            u64 x = distance_m3.a[i / 32] ^modify[prev_depth_m3];
            if (((x >> u64(1)) & x & mask) == 0) {
                continue;
            }
            u64 j_end = std::min(i + 32, end);
            for (u64 j = i; j < j_end; j++) {
                if (distance_m3.get(j) == prev_depth_m3) {
                    typename _solver::t_state a = s.int_to_state(j);
                    for (typename _solver::t_state b: s.adj(a)) {
                        u64 k = s.state_to_int(b);
                        if (distance_m3.get(k) == 3) {
                            u64 c = _set_multi<_solver::n_state>(distance_m3, s.alt(b, k), 3, depth_m3);
                            if (c > 0) {
                                count_distinct++;
                                count += c;
                            }
                        }
                    }
                }
            }
        }
        return {count_distinct, count};
    }

    template<typename _solver>
    std::tuple<u64, u64> _backward(const _solver &s, array_u2<_solver::n_state> &distance_m3,
                                   u64 start, u64 end, u64 prev_depth_m3, u64 depth_m3) {
        constexpr u64 mask = 0x5555555555555555;
        u64 count_distinct = 0;
        u64 count = 0;
        for (u64 i = start; i < end; i += 32) {
            u64 x = distance_m3.a[i / 32];
            if (((x >> u64(1)) & x & mask) == 0) {
                continue;
            }
            u64 j_end = std::min(i + 32, end);
            for (u64 j = i; j < j_end; j++) {
                if (distance_m3.get(j) == 3) {
                    typename _solver::t_state a = s.int_to_state(j);
                    for (typename _solver::t_state b: s.adj(a)) {
                        u64 k = s.state_to_int(b);
                        if (distance_m3.get(k) == prev_depth_m3) {
                            u64 c = _set_multi<_solver::n_state>(distance_m3, s.alt(a, j), 3, depth_m3);
                            if (c > 0) {
                                count_distinct++;
                                count += c;
                            }
                        }
                    }
                }
            }
        }
        return {count_distinct, count};
    }

    template<typename _solver>
    void bfs(const _solver &s, array_u2<_solver::n_state> &distance_m3, u64 n_thread = 1) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<u64> sp = _split(_solver::n_state, n_thread);
        distance_m3.fill(3);
        u64 count_start;
        {
            typename _solver::t_state a_start = s.cube_to_state(_solver::t_cube::i());
            u64 i_start = s.state_to_int(a_start);
            count_start = _set_multi<_solver::n_state>(distance_m3, s.alt(a_start, i_start), 3, 0);
            assert(count_start > 0);
        }
        std::cout << "bfs: n_state=" << _solver::n_state << ", n_thread=" << n_thread << std::endl;
        std::cout << "bfs: depth=" << 0 << ", count_distinct=" << 1 << ", count=" << count_start << std::endl;
        u64 total_count_distinct = 1;
        u64 total_count = count_start;
        std::array<u64, 3> count_m3{count_start, 0, 0};
        for (u64 depth = 1; total_count != _solver::n_state; depth++) {
            auto t1 = std::chrono::steady_clock::now();
            u64 count_distinct = 0;
            u64 count = 0;
            std::vector<std::future<std::tuple<u64, u64>>> result{};
            if (count_m3[(depth - 1) % 3] <= _solver::n_state - total_count) {
                for (u64 i = 0; i < n_thread; i++) {
                    result.push_back(std::async(
                            std::launch::async,
                            &_forward<_solver>,
                            std::cref(s), std::ref(distance_m3), sp[i], sp[i + 1], (depth - 1) % 3, depth % 3));
                }
            } else {
                for (u64 i = 0; i < n_thread; i++) {
                    result.push_back(std::async(
                            std::launch::async,
                            &_backward<_solver>,
                            std::cref(s), std::ref(distance_m3), sp[i], sp[i + 1], (depth - 1) % 3, depth % 3));
                }
            }
            for (u64 i = 0; i < n_thread; i++) {
                auto[cd, c] = result[i].get();
                count_distinct += cd;
                count += c;
            }
            auto t2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> d = t2 - t1;
            std::cout << "bfs: depth=" << depth << ", count_distinct=" << count_distinct << ", count=" << count
                      << ", time=" << d.count() << "s" << std::endl;
            total_count_distinct += count_distinct;
            total_count += count;
            count_m3[depth % 3] += count;
        }
        auto t3 = std::chrono::steady_clock::now();
        std::chrono::duration<double> d = t3 - t0;
        std::cout << "bfs: total_count_distinct=" << total_count_distinct << ", total_count=" << total_count
                  << ", total_time=" << d.count() << "s" << std::endl;
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
