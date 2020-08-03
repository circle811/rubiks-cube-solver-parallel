#ifndef _CUBE3_2P_H
#define _CUBE3_2P_H

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube3.h"

//  G = <U, D, R, L, F, B>
//  H = <U, D, R2, L2, F2, B2>
//  I = <>
//
//  phase 0: G -> H
//  phase 1: H -> I

namespace cube::_3::_2p {
    typedef orbits <orbit<0, 1, 2, 3>, orbit<4, 5, 6, 7, 8, 9, 10, 11>> os_e;
    typedef orbits <orbit<0>, orbit<1>, orbit<2>, orbit<3>, orbit<4>, orbit<5>, orbit<6, 7>> os_c;

    // phase 0
    constexpr u64 n_co = number_o<8, 3>();
    constexpr u64 n_egp = number_gp<os_e>();
    constexpr u64 n_eo = number_o<12, 2>();

    struct p0a {
        array_u8<8> co;
        array_u8<12> egp;
        array_u8<12> eo;

        constexpr p0a operator*(const cube3a &other) const {
            const po<8, 3> &rc = mul_po<8, 3>(i_p<8>(), co, other.cp, other.co);
            const po<12, 2> &re = mul_po<12, 2>(egp, eo, other.ep, other.eo);
            return p0a{std::get<1>(rc),
                       std::get<0>(re),
                       std::get<1>(re)};
        }

        constexpr bool operator==(const p0a &other) const {
            return (co == other.co) and (egp == other.egp) and (eo == other.eo);
        }
    };

    struct p0b {
        u32 co;
        u32 egp_eo;

        constexpr bool operator==(const p0b &other) const {
            return (co == other.co) and (egp_eo == other.egp_eo);
        }
    };

    constexpr p0a cube3a_to_p0a(const cube3a &a) {
        return p0a{a.co,
                   p_to_gp<os_e>(a.ep),
                   a.eo};
    }

    constexpr p0b p0a_to_b(const p0a &a) {
        return p0b{u32(o_to_int<8, 3>(a.co)),
                   u32(gp_to_int<os_e>(a.egp) * n_eo + o_to_int<12, 2>(a.eo))};
    }

    constexpr p0a p0b_to_a(const p0b &b) {
        return p0a{int_to_o<8, 3>(b.co),
                   int_to_gp<os_e>(b.egp_eo / n_eo),
                   int_to_o<12, 2>(b.egp_eo % n_eo)};
    }

    constexpr u64 p0b_to_int(const p0b &b) {
        return b.co * (n_egp * n_eo) + b.egp_eo;
    }

    constexpr p0b int_to_p0b(u64 x) {
        return p0b{u32(x / (n_egp * n_eo)),
                   u32(x % (n_egp * n_eo))};
    }

    constexpr std::array<cube3a, 18> p0_base = cube3_base;

    const std::array<std::string, cube3_base.size()> p0_base_name = cube3_base_name;

    struct p0_solver {
        typedef p0b v_type;

        static constexpr u64 size = n_co * n_egp * n_eo;

        static constexpr u64 base_size = p0_base.size();

        static constexpr v_type start = p0a_to_b(cube3a_to_p0a(cube3a::i()));

        static constexpr u64 v_to_int(const v_type &a) {
            return p0b_to_int(a);
        }

        static constexpr v_type int_to_v(u64 x) {
            return int_to_p0b(x);
        }

        std::unique_ptr<array_2d < u32, n_co, base_size>> table_mul_co;

        std::unique_ptr<array_2d < u32, n_egp * n_eo, base_size>> table_mul_egp_eo;

        std::unique_ptr<array_u2 < size>> table_dist_m3;

        p0_solver() {
            table_mul_co = cache_data<array_2d<u32, n_co, base_size>>(
                    "p0.table_mul_co",
                    [](array_2d<u32, n_co, base_size> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            const p0a &a = p0b_to_a(p0b{u32(i), 0});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = p0a_to_b(a * p0_base[j]).co;
                            }
                        }
                    }
            );
            table_mul_egp_eo = cache_data<array_2d<u32, n_egp * n_eo, base_size>>(
                    "p0.table_mul_egp_eo",
                    [](array_2d<u32, n_egp * n_eo, base_size> &t) -> void {
                        for (u64 i = 0; i < (n_egp * n_eo); i++) {
                            const p0a &a = p0b_to_a(p0b{0, u32(i)});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = p0a_to_b(a * p0_base[j]).egp_eo;
                            }
                        }
                    }
            );
            table_dist_m3 = cache_data<array_u2<size>>(
                    "p0.table_dist_m3",
                    [this](array_u2<size> &t) {
                        bfs_m3<p0_solver>(*this, t);
                    }
            );
        }

        std::array<v_type, base_size> adj(const v_type &a) const {
            std::array<v_type, base_size> bs{};
            for (u64 i = 0; i < base_size; i++) {
                bs[i] = v_type{(*table_mul_co)[a.co][i],
                               (*table_mul_egp_eo)[a.egp_eo][i]};
            }
            return bs;
        }
    };

    // phase 1
    constexpr u64 n_cgp = number_gp<os_c>();
    constexpr u64 n_ep0 = number_p<4>();
    constexpr u64 n_ep1 = number_p<8>();

    struct p1a {
        array_u8<8> cgp;
        array_u8<4> ep0;
        array_u8<8> ep1;

        constexpr p1a operator*(const cube3a &other) const {
            const os_e::pp_type &pp = p_to_pp<os_e>(other.ep);
            return p1a{mul_p<8>(cgp, other.cp),
                       mul_p<4>(ep0, std::get<0>(pp)),
                       mul_p<8>(ep1, std::get<1>(pp))};
        }

        constexpr bool operator==(const p1a &other) const {
            return cgp == other.cgp and ep0 == other.ep0 and ep1 == other.ep1;
        }
    };

    struct p1b {
        u32 cgp;
        u32 ep0_ep1;

        constexpr bool operator==(const p1b &other) const {
            return (cgp == other.cgp) and (ep0_ep1 == other.ep0_ep1);
        }
    };

    constexpr p1a cube3a_to_p1a(const cube3a &a) {
        const os_e::pp_type &pp = p_to_pp<os_e>(a.ep);
        return p1a{p_to_gp<os_c>(a.cp),
                   std::get<0>(pp),
                   std::get<1>(pp)};
    }

    constexpr p1b p1a_to_b(const p1a &a) {
        return p1b{u32(gp_to_int<os_c>(a.cgp)),
                   u32(p_to_int<4>(a.ep0) * n_ep1 + p_to_int<8>(a.ep1))};
    }

    constexpr p1a p1b_to_a(const p1b &b) {
        return p1a{int_to_gp<os_c>(b.cgp),
                   int_to_p<4>(b.ep0_ep1 / n_ep1),
                   int_to_p<8>(b.ep0_ep1 % n_ep1)};
    }

    constexpr u64 p1b_to_int(const p1b &b) {
        return b.cgp * (n_ep0 * n_ep1) + b.ep0_ep1;
    }

    constexpr p1b int_to_p1b(u64 x) {
        return p1b{u32(x / (n_ep0 * n_ep1)),
                   u32(x % (n_ep0 * n_ep1))};
    }

    constexpr std::array<cube3a, 10> p1_base{
            U1, U2, U3, D1, D2, D3,
            R2, L2, F2, B2
    };

    const std::array<std::string, p1_base.size()> p1_base_name{
            "U", "U2", "U'", "D", "D2", "D'",
            "R2", "L2", "F2", "B2"
    };

    struct p1_solver {
        typedef p1b v_type;

        static constexpr u64 size = n_cgp * n_ep0 * n_ep1;

        static constexpr u64 base_size = p1_base.size();

        static constexpr v_type start = p1a_to_b(cube3a_to_p1a(cube3a::i()));

        static constexpr u64 v_to_int(const v_type &a) {
            return p1b_to_int(a);
        }

        static constexpr v_type int_to_v(u64 x) {
            return int_to_p1b(x);
        }

        std::unique_ptr<array_2d < u32, n_cgp, base_size>> table_mul_cgp;

        std::unique_ptr<array_2d < u32, n_ep0 * n_ep1, base_size>> table_mul_ep0_ep1;

        std::unique_ptr<array_u2 < size>> table_dist_m3;

        p1_solver() {
            table_mul_cgp = cache_data<array_2d<u32, n_cgp, base_size>>(
                    "p1.table_mul_cgp",
                    [](array_2d<u32, n_cgp, base_size> &t) -> void {
                        for (u64 i = 0; i < n_cgp; i++) {
                            const p1a &a = p1b_to_a(p1b{u32(i), 0});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = p1a_to_b(a * p1_base[j]).cgp;
                            }
                        }
                    }
            );
            table_mul_ep0_ep1 = cache_data<array_2d<u32, n_ep0 * n_ep1, base_size>>(
                    "p1.table_mul_ep0_ep1",
                    [](array_2d<u32, n_ep0 * n_ep1, base_size> &t) -> void {
                        for (u64 i = 0; i < (n_ep0 * n_ep1); i++) {
                            const p1a &a = p1b_to_a(p1b{0, u32(i)});
                            for (u64 j = 0; j < base_size; j++) {
                                t[i][j] = p1a_to_b(a * p1_base[j]).ep0_ep1;
                            }
                        }
                    }
            );
            table_dist_m3 = cache_data<array_u2<size>>(
                    "p1.table_dist_m3",
                    [this](array_u2<size> &t) {
                        bfs_m3<p1_solver>(*this, t);
                    }
            );
        }

        std::array<v_type, base_size> adj(const v_type &a) const {
            std::array<v_type, base_size> bs{};
            for (u64 i = 0; i < base_size; i++) {
                bs[i] = v_type{(*table_mul_cgp)[a.cgp][i],
                               (*table_mul_ep0_ep1)[a.ep0_ep1][i]};
            }
            return bs;
        }
    };
}

#endif
