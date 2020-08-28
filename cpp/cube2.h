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
    struct cube2 {
        array_u8<7> cp;
        array_u8<7> co;

        static constexpr cube2 i() {
            return cube2{
                    array_u8<7>::i(),
                    i_o<7, 3>()
            };
        }

        constexpr cube2 inv() const {
            array_u8<7> cp_inv = cp.inv();
            return cube2{
                    cp_inv,
                    inv_o<7, 3>(co) * cp_inv
            };
        }

        constexpr cube2 operator*(const cube2 &other) const {
            return cube2{
                    cp * other.cp,
                    mul_o<7, 3>(co * other.cp, other.co)
            };
        }

        constexpr bool operator==(const cube2 &other) const {
            return cp == other.cp and co == other.co;
        }

        std::string to_string() const {
            return "cube2(" + cp.to_string() + " " + co.to_string() + ")";
        }
    };

    constexpr cube2 U1{{1, 3, 0, 2, 4, 5, 6},
                       {0, 0, 0, 0, 0, 0, 0}};
    constexpr cube2 R1{{4, 0, 2, 3, 5, 1, 6},
                       {2, 1, 0, 0, 1, 2, 0}};
    constexpr cube2 F1{{2, 1, 6, 3, 0, 5, 4},
                       {1, 0, 2, 0, 2, 0, 1}};
    constexpr cube2 U2 = U1 * U1;
    constexpr cube2 R2 = R1 * R1;
    constexpr cube2 F2 = F1 * F1;
    constexpr cube2 U3 = U1.inv();
    constexpr cube2 R3 = R1.inv();
    constexpr cube2 F3 = F1.inv();

    constexpr u64 n_cube2_base = 9;

    constexpr std::array<cube2, n_cube2_base> cube2_base{
            U1, U2, U3,
            R1, R2, R3,
            F1, F2, F3
    };

    constexpr std::array<const char *, n_cube2_base> cube2_base_name{
            "U", "U2", "U'",
            "R", "R2", "R'",
            "F", "F2", "F'"
    };

    struct cube2_solver {
        typedef cube2 t_cube;
        typedef u64 t_state;
        typedef u64 t_hint;

        static constexpr u64 n_cp = number_p<7>();
        static constexpr u64 n_co = number_o<7, 3>();
        static constexpr u64 n_state = n_cp * n_co;
        static constexpr u64 n_base = n_cube2_base;

        static constexpr std::array<t_cube, n_base> base = cube2_base;
        static constexpr std::array<const char *, n_base> base_name = cube2_base_name;
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr t_state cube_to_state(const t_cube &a) {
            return p_to_int<7>(a.cp) * n_co + o_to_int<7, 3>(a.co);
        }

        static constexpr u64 state_to_int(const t_state &a) {
            return a;
        }

        static constexpr t_state int_to_state(u64 x) {
            return x;
        }

        static constexpr bool is_start(const t_state &a) {
            constexpr t_state _start = cube_to_state(t_cube::i());
            return a == _start;
        }

        static constexpr std::array<u64, 1> alt(const t_state &a, u64 i) {
            return std::array<u64, 1>{i};
        }

        u64 n_thread;
        std::unique_ptr<array_2d < u16, n_cp, n_base>> mul_cp;
        std::unique_ptr<array_2d < u16, n_co, n_base>> mul_co;
        std::unique_ptr<array_u2 < n_state>> distance_m3;

        explicit cube2_solver(u64 _n_thread) : n_thread(_n_thread) {
            mul_cp = cache_data<array_2d<u16, n_cp, n_base>>(
                    "cube2.mul_cp",
                    [](array_2d<u16, n_cp, n_base> &t) -> void {
                        for (u64 i = 0; i < n_cp; i++) {
                            array_u8<7> cp = int_to_p<7>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(p_to_int<7>(cp * base[j].cp));
                            }
                        }
                    }
            );
            mul_co = cache_data<array_2d<u16, n_co, n_base>>(
                    "cube2.mul_co",
                    [](array_2d<u16, n_co, n_base> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            array_u8<7> co = int_to_o<7, 3>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                mul_o<7, 3>(co * base[j].cp, base[j].co);
                                t[i][j] = u16(o_to_int<7, 3>(mul_o<7, 3>(co * base[j].cp, base[j].co)));
                            }
                        }
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube2.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<cube2_solver>(*this, t, n_thread);
                    }
            );
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_cp = (*mul_cp)[a / n_co];
            const std::array<u16, n_base> &a_co = (*mul_co)[a % n_co];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                a_s[i] = a_cp[i] * n_co + a_co[i];
            }
            return a_s;
        }

        template<u64 capacity>
        ida_star <cube2_solver, capacity> solve(const t_cube &a, u64 max_n_moves = capacity) const {
            return ida_star<cube2_solver, capacity>(*this, a, max_n_moves);
        }
    };
}

#endif
