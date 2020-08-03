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
            assert (x < 4);
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            a[j] = (a[j] & ~(u64(3) << k)) | (x << k);
        }

        void fill(u64 x) {
            assert (x < 4);
            u64 y = x;
            for (u64 i = 1; i < 6; i++) {
                y = y | (y << (u64(1) << i));
            }
            a.fill(y);
        }
    };

    template<typename _solver>
    void bfs_m3(const _solver &s, array_u2<_solver::size> &table_dist_m3) {
        constexpr std::array<u64, 3> alt{
                0xffffffffffffffff,
                0xaaaaaaaaaaaaaaaa,
                0x5555555555555555,
        };
        constexpr u64 mask = 0x5555555555555555;
        table_dist_m3.fill(3);
        table_dist_m3.set(s.v_to_int(_solver::start), 0);
        std::cout << "bfs_m3: depth=" << 0 << ", count=" << 1 << std::endl;
        u64 depth = 1;
        u64 total = 1;
        std::array<u64, 3> count_m3{1, 0, 0};
        while (total != _solver::size) {
            u64 prev_depth_m3 = (depth - 1) % 3;
            u64 depth_m3 = depth % 3;
            u64 count = 0;
            if (count_m3[prev_depth_m3] <= _solver::size - total) {
                for (u64 i = 0; i < _solver::size; i += 32) {
                    u64 x = table_dist_m3.a[i / 32] ^alt[prev_depth_m3];
                    if (((x >> u64(1)) & x & mask) == 0) {
                        continue;
                    }
                    u64 end = i + 32 < _solver::size ? i + 32 : _solver::size;
                    for (u64 i_a = i; i_a < end; i_a++) {
                        if (table_dist_m3.get(i_a) == prev_depth_m3) {
                            const typename _solver::v_type &a = s.int_to_v(i_a);
                            for (const typename _solver::v_type &b: s.adj(a)) {
                                u64 i_b = s.v_to_int(b);
                                if (table_dist_m3.get(i_b) == 3) {
                                    table_dist_m3.set(i_b, depth_m3);
                                    count++;
                                }
                            }
                        }
                    }
                }
            } else {
                for (u64 i = 0; i < _solver::size; i += 32) {
                    u64 x = table_dist_m3.a[i / 32];
                    if (((x >> u64(1)) & x & mask) == 0) {
                        continue;
                    }
                    u64 end = i + 32 < _solver::size ? i + 32 : _solver::size;
                    for (u64 i_b = i; i_b < end; i_b++) {
                        if (table_dist_m3.get(i_b) == 3) {
                            const typename _solver::v_type &b = s.int_to_v(i_b);
                            for (const typename _solver::v_type &a: s.adj(b)) {
                                u64 i_a = s.v_to_int(a);
                                if (table_dist_m3.get(i_a) == prev_depth_m3) {
                                    table_dist_m3.set(i_b, depth_m3);
                                    count++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            std::cout << "bfs_m3: depth=" << depth << ", count=" << count << std::endl;
            depth++;
            total += count;
            count_m3[depth_m3] += count;
        }
    }

    template<typename _solver>
    u64 get_dist_m3(const _solver &s, const typename _solver::v_type &a) {
        return s.table_dist_m3->get(s.v_to_int(a));
    }

    template<typename _solver>
    u64 _dfs_dist(const _solver &s, const typename _solver::v_type &a, u64 depth) {
        if (a == _solver::start) {
            return depth;
        } else {
            u64 t = (get_dist_m3<_solver>(s, a) + 2) % 3;
            for (const typename _solver::v_type &b: s.adj(a)) {
                if (get_dist_m3<_solver>(s, b) == t) {
                    return _dfs_dist<_solver>(s, b, depth + 1);
                }
            }
            assert (false);
        }
    }

    template<typename _solver>
    u64 get_dist(const _solver &s, const typename _solver::v_type &a) {
        return _dfs_dist<_solver>(s, a, 0);
    }

    constexpr u64 compute_dist(u64 dist_adj, u64 dist_m3) {
        return dist_adj + (dist_m3 - dist_adj - 3) % 3 - 1;
    }

    typedef std::function<void(const std::vector<u64> &)> action_type;

    template<typename _solver>
    void _dfs_m3(const _solver &s, const typename _solver::v_type &a, u64 dist_a, u64 max_depth,
                 std::vector<u64> &moves, std::array<u64, 3> &count, const action_type &action) {
        count[0]++;
        if (moves.size() == max_depth) {
            count[1]++;
            if (a == _solver::start) {
                count[2]++;
                action(moves);
            }
        } else {
            const std::array<typename _solver::v_type, _solver::base_size> &bs = s.adj(a);
            for (u64 i = 0; i < _solver::base_size; i++) {
                const typename _solver::v_type &b = bs[i];
                u64 dist_b = compute_dist(dist_a, get_dist_m3<_solver>(s, b));
                if (moves.size() + 1 + dist_b <= max_depth) {
                    moves.push_back(i);
                    _dfs_m3<_solver>(s, b, dist_b, max_depth, moves, count, action);
                    moves.pop_back();
                }
            }
        }
    }

    template<typename _solver>
    void ida_star_m3(const _solver &s, const typename _solver::v_type &a, u64 addition, const action_type &action) {
        u64 dist_a = get_dist<_solver>(s, a);
        for (u64 max_depth = dist_a; max_depth <= dist_a + addition; max_depth++) {
            std::cout << "ida_star_m3: max_depth=" << max_depth << std::endl;
            std::vector<u64> moves{};
            std::array<u64, 3> count{0, 0, 0};
            _dfs_m3<_solver>(s, a, dist_a, max_depth, moves, count, action);
            std::cout << "ida_star_m3: max_depth=" << max_depth
                      << ", count=(" << count[0] << " " << count[1] << " " << count[2] << ")"
                      << std::endl;
        }
    }

    template<typename _solver>
    std::vector<std::vector<u64>> solve(const _solver &s, const typename _solver::v_type &a) {
        std::vector<std::vector<u64>> vector_moves{};
        ida_star_m3<_solver>(
                s, a, 0,
                [&vector_moves](const std::vector<u64> &moves) -> void {
                    vector_moves.push_back(moves);
                }
        );
        return vector_moves;
    }

    template<typename cube, u64 base_size>
    cube moves_to_cube(const std::array<cube, base_size> &base, const std::vector<u64> &moves) {
        cube c = cube::i();
        for (u64 m: moves) {
            c = c * base[m];
        }
        return c;
    }

    template<u64 base_size>
    std::string moves_to_string(const std::array<std::string, base_size> &base_name, const std::vector<u64> &moves) {
        std::string s{};
        for (u64 m: moves) {
            s += base_name[m];
            s += " ";
        }
        return s;
    }

    struct random_moves {
        std::default_random_engine e;
        std::uniform_int_distribution<u64> d;

        random_moves(u64 base_size, u64 seed) : e(seed), d(0, base_size - 1) {
        }

        std::vector<u64> operator()(u64 size) {
            std::vector<u64> moves{};
            for (u64 i = 0; i < size; i++) {
                moves.push_back(d(e));
            }
            return moves;
        }
    };
}

#endif
