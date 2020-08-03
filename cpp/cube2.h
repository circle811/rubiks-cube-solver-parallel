#ifndef _CUBE2_H
#define _CUBE2_H

#include "base.h"
#include "group.h"
#include "search.h"

//  number of corner blocks
//
//     3    1
//  2    0
//
//     *    5
//  6    4

namespace cube::_2 {
    constexpr u64 n_cp = number_p<7>();
    constexpr u64 n_co = number_o<7, 3>();

    struct cube2a {
        array_u8<7> cp;
        array_u8<7> co;

        static constexpr cube2a i() {
            return cube2a{i_p<7>(), i_o<7, 3>()};
        }

        constexpr cube2a inv() const {
            const po<7, 3> &r = inv_po<7, 3>(cp, co);
            return cube2a{std::get<0>(r), std::get<1>(r)};
        }

        constexpr cube2a operator*(const cube2a &other) const {
            const po<7, 3> &r = mul_po<7, 3>(cp, co, other.cp, other.co);
            return cube2a{std::get<0>(r), std::get<1>(r)};
        }

        constexpr bool operator==(const cube2a &other) const {
            return (cp == other.cp) and (co == other.co);
        }
    };

    struct cube2b {
        u16 cp;
        u16 co;

        constexpr bool operator==(const cube2b &other) const {
            return (cp == other.cp) and (co == other.co);
        }
    };

    constexpr cube2b cube2a_to_b(const cube2a &a) {
        return cube2b{u16(p_to_int<7>(a.cp)), u16(o_to_int<7, 3>(a.co))};
    }

    constexpr cube2a cube2b_to_a(const cube2b &b) {
        return cube2a{int_to_p<7>(b.cp), int_to_o<7, 3>(b.co)};
    }

    constexpr u64 cube2b_to_int(const cube2b &b) {
        return b.cp * n_co + b.co;
    }

    constexpr cube2b int_to_cube2b(u64 x) {
        return cube2b{u16(x / n_co), u16(x % n_co)};
    }

    constexpr cube2a U1{{1, 3, 0, 2, 4, 5, 6},
                        {0, 0, 0, 0, 0, 0, 0}};
    constexpr cube2a R1{{4, 0, 2, 3, 5, 1, 6},
                        {2, 1, 0, 0, 1, 2, 0}};
    constexpr cube2a F1{{2, 1, 6, 3, 0, 5, 4},
                        {1, 0, 2, 0, 2, 0, 1}};
    constexpr cube2a U2 = U1 * U1;
    constexpr cube2a R2 = R1 * R1;
    constexpr cube2a F2 = F1 * F1;
    constexpr cube2a U3 = U1.inv();
    constexpr cube2a R3 = R1.inv();
    constexpr cube2a F3 = F1.inv();

    constexpr std::array<cube2a, 9> cube2_base{
            U1, U2, U3, R1, R2, R3, F1, F2, F3
    };

    const std::array<std::string, cube2_base.size()> cube2_base_name{
            "U", "U2", "U'", "R", "R2", "R'", "F", "F2", "F'"
    };

    struct cube2_solver {
        typedef cube2b v_type;

        static constexpr u64 size = n_co * n_cp;

        static constexpr u64 base_size = cube2_base.size();

        static constexpr v_type start = cube2a_to_b(cube2a::i());

        static constexpr u64 v_to_int(const v_type &a) {
            return cube2b_to_int(a);
        }

        static constexpr v_type int_to_v(u64 x) {
            return int_to_cube2b(x);
        }

        std::unique_ptr<array_2d < u16, n_cp, base_size>> table_mul_cp;

        std::unique_ptr<array_2d < u16, n_co, base_size>> table_mul_co;

        std::unique_ptr<array_u2 < size>> table_dist_m3;

        cube2_solver() {
            table_mul_cp = cache_data<array_2d<u16, n_cp, base_size>>(
                    "cube2.table_mul_cp",
                    [](array_2d<u16, n_cp, base_size> &t) -> void {
                        for (u64 i = 0; i < n_cp; i++) {
                            const cube2a &a = cube2b_to_a(cube2b{u16(i), 0});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = cube2a_to_b(a * cube2_base[j]).cp;
                            }
                        }
                    }
            );
            table_mul_co = cache_data<array_2d<u16, n_co, base_size>>(
                    "cube2.table_mul_co",
                    [](array_2d<u16, n_co, base_size> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            const cube2a &a = cube2b_to_a(cube2b{0, u16(i)});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = cube2a_to_b(a * cube2_base[j]).co;
                            }
                        }
                    }
            );
            table_dist_m3 = cache_data<array_u2<size>>(
                    "cube2.table_dist_m3",
                    [this](array_u2<size> &t) {
                        bfs_m3<cube2_solver>(*this, t);
                    }
            );
        }

        std::array<v_type, base_size> adj(const v_type &a) const {
            std::array<v_type, base_size> bs{};
            for (u64 i = 0; i < base_size; i++) {
                bs[i] = v_type{(*table_mul_cp)[a.cp][i], (*table_mul_co)[a.co][i]};
            }
            return bs;
        }
    };
}

#endif
