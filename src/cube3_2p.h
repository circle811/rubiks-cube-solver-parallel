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
    typedef orbits <orbit<0, 1, 2, 3>, orbit<4, 5, 6, 7, 8, 9, 10, 11>> os_4_8;
    typedef orbits <orbit<0>, orbit<1>, orbit<2>, orbit<3>, orbit<4>, orbit<5>, orbit<6, 7>> os_1x6_2;

    // phase 0
    struct p0_solver {
        typedef cube3 t_cube;
        typedef u64 t_state;
        typedef u64 t_hint;

        static constexpr u64 n_co = number_o<8, 3>();
        static constexpr u64 n_egp = number_gp<os_4_8>();
        static constexpr u64 n_eo = number_o<12, 2>();
        static constexpr u64 n_state = n_co * n_egp * n_eo;
        static constexpr u64 n_super_base = n_cube3_base;
        static constexpr u64 n_base = n_cube3_base;

        static constexpr std::array<t_cube, n_super_base> super_base = cube3_base;
        static constexpr std::array<const char *, n_super_base> super_base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_index{
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
        };
        static constexpr std::array<t_cube, n_base> base = array_sub<t_cube, n_super_base, n_base>(
                super_base, base_index);
        static constexpr std::array<const char *, n_base> base_name = array_sub<const char *, n_super_base, n_base>(
                super_base_name, base_index);
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr t_state cube_to_state(const t_cube &a) {
            return o_to_int<8, 3>(a.co) * (n_egp * n_eo)
                   + gp_to_int<os_4_8>(p_to_gp<os_4_8>(a.ep)) * n_eo
                   + o_to_int<12, 2>(a.eo);
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
        std::unique_ptr<array_2d < u16, n_co, n_base>> mul_co;
        std::unique_ptr<array_2d < u32, n_egp * n_eo, n_base>> mul_egp_eo;
        std::unique_ptr<array_u2 < n_state>> distance_m3;

        explicit p0_solver(u64 _n_thread) : n_thread(_n_thread) {
            mul_co = cache_data<array_2d<u16, n_co, n_base>>(
                    "cube3.p0.mul_co",
                    [](array_2d<u16, n_co, n_base> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            array_u8<8> co = int_to_o<8, 3>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(o_to_int<8, 3>(mul_o<8, 3>(co * base[j].cp, base[j].co)));
                            }
                        }
                    }
            );
            mul_egp_eo = cache_data<array_2d<u32, n_egp * n_eo, n_base>>(
                    "cube3.p0.mul_egp_eo",
                    [](array_2d<u32, n_egp * n_eo, n_base> &t) -> void {
                        for (u64 i = 0; i < n_egp; i++) {
                            array_u8<12> egp = int_to_gp<os_4_8>(i);
                            for (u64 j = 0; j < n_eo; j++) {
                                array_u8<12> eo = int_to_o<12, 2>(j);
                                for (u64 k = 0; k < n_base; k++) {
                                    t[i * n_eo + j][k] =
                                            u32(gp_to_int<os_4_8>(egp * base[k].ep) * n_eo
                                                + o_to_int<12, 2>(mul_o<12, 2>(eo * base[k].ep, base[k].eo)));
                                }
                            }
                        }
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube3.p0.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<p0_solver>(*this, t, n_thread);
                    }
            );
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_co = (*mul_co)[a / (n_egp * n_eo)];
            const std::array<u32, n_base> &a_egp_eo = (*mul_egp_eo)[a % (n_egp * n_eo)];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                a_s[i] = a_co[i] * (n_egp * n_eo) + a_egp_eo[i];
            }
            return a_s;
        }

        template<u64 capacity>
        ida_star <p0_solver, capacity> solve(const t_cube &a, u64 max_n_moves = capacity) const {
            return ida_star<p0_solver, capacity>(*this, a, max_n_moves);
        }
    };

    // phase 0 symmetry
    template<typename _os_e, typename U_SC, u64 _n_sc_egp_eo, const char *name>
    struct g_p0s_solver {
        typedef cube3 t_cube;

        struct t_state {
            u8 sym;
            u16 co;
            U_SC sc_egp_eo;
        };

        typedef u64 t_hint;

        static constexpr u64 n_co = number_o<8, 3>();
        static constexpr u64 n_egp = number_gp<_os_e>();
        static constexpr u64 n_eo = number_o<12, 2>();
        static constexpr u64 n_sc_egp_eo = _n_sc_egp_eo;
        static constexpr u64 n_state = n_co * n_sc_egp_eo;
        static constexpr u64 n_super_base = n_cube3_base;
        static constexpr u64 n_base = n_cube3_base;

        static constexpr std::array<t_cube, n_super_base> super_base = cube3_base;
        static constexpr std::array<const char *, n_super_base> super_base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_index{
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
        };
        static constexpr std::array<t_cube, n_base> base = array_sub<t_cube, n_super_base, n_base>(
                super_base, base_index);
        static constexpr std::array<const char *, n_base> base_name = array_sub<const char *, n_super_base, n_base>(
                super_base_name, base_index);
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr array_2d <u8, n_base, n_s16> conj_base =
                generate_table_conj_base<cube3, n_base, n_s16>(base, elements_s16);

        u64 n_thread;
        std::unique_ptr<array_2d < u16, n_co, n_s16>> conj_co;
        std::unique_ptr<array_2d < u16, n_co, n_base>> mul_co;
        std::unique_ptr<table_conj_mul < u32, U_SC, u16, n_egp * n_eo, n_sc_egp_eo, n_s16, n_base>> conj_mul_egp_eo;
        std::unique_ptr<array_u2 < n_state>> distance_m3;
        t_state _start;

        explicit g_p0s_solver(u64 _n_thread) : n_thread(_n_thread) {
            conj_co = cache_data<array_2d<u16, n_co, n_s16>>(
                    std::string("cube3.") + name + ".conj_co",
                    [](array_2d<u16, n_co, n_s16> &t) -> void {
                        cube3 a = cube3::i();
                        for (u64 i = 0; i < n_co; i++) {
                            a.co = int_to_o<8, 3>(i);
                            for (u64 s = 0; s < n_s16; s++) {
                                cube3 b = elements_s16[inv_s16[s]] * a * elements_s16[s];
                                t[i][s] = u16(o_to_int<8, 3>(b.co));
                            }
                        }
                    }
            );
            mul_co = cache_data<array_2d<u16, n_co, n_base>>(
                    std::string("cube3.") + name + ".mul_co",
                    [](array_2d<u16, n_co, n_base> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            array_u8<8> co = int_to_o<8, 3>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(o_to_int<8, 3>(mul_o<8, 3>(co * base[j].cp, base[j].co)));
                            }
                        }
                    }
            );
            conj_mul_egp_eo = cache_data<table_conj_mul<u32, U_SC, u16, n_egp * n_eo, n_sc_egp_eo, n_s16, n_base>>(
                    std::string("cube3.") + name + ".conj_mul_egp_eo",
                    [](table_conj_mul<u32, U_SC, u16, n_egp * n_eo, n_sc_egp_eo, n_s16, n_base> &cm) -> void {
                        cm.init(
                                [](u32 i) -> std::array<u32, n_s16> {
                                    cube3 a = cube3::i();
                                    a.ep = gp_to_p<_os_e>(int_to_gp<_os_e>(i / n_eo));
                                    a.eo = int_to_o<12, 2>(i % n_eo);
                                    std::array<u32, n_s16> conj_i{};
                                    for (u64 s = 0; s < n_s16; s++) {
                                        cube3 b = elements_s16[inv_s16[s]] * a * elements_s16[s];
                                        conj_i[s] = u32(gp_to_int<_os_e>(p_to_gp<_os_e>(b.ep)) * n_eo
                                                        + o_to_int<12, 2>(b.eo));
                                    }
                                    return conj_i;
                                },
                                [](u32 i) -> std::array<u32, n_base> {
                                    array_u8<12> egp = int_to_gp<_os_e>(i / n_eo);
                                    array_u8<12> eo = int_to_o<12, 2>(i % n_eo);
                                    std::array<u32, n_base> mul_i{};
                                    for (u64 j = 0; j < n_base; j++) {
                                        mul_i[j] = u32(gp_to_int<_os_e>(egp * base[j].ep) * n_eo
                                                       + o_to_int<12, 2>(mul_o<12, 2>(eo * base[j].ep, base[j].eo)));
                                    }
                                    return mul_i;
                                }
                        );
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    std::string("cube3.") + name + ".distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<g_p0s_solver<_os_e, U_SC, _n_sc_egp_eo, name>>(*this, t, n_thread);
                    }
            );
            _start = cube_to_state(t_cube::i());
        }

        t_state cube_to_state(const t_cube &a) const {
            u64 egp_eo = gp_to_int<_os_e>(p_to_gp<_os_e>(a.ep)) * n_eo + o_to_int<12, 2>(a.eo);
            auto[sym, sc] = conj_mul_egp_eo->g_to_sym_sc(egp_eo);
            return t_state{
                    sym,
                    (*conj_co)[o_to_int<8, 3>(a.co)][inv_s16[sym]],
                    sc
            };
        }

        u64 state_to_int(const t_state &a) const {
            return a.co * n_sc_egp_eo + a.sc_egp_eo;
        }

        t_state int_to_state(u64 x) const {
            return t_state{
                    0,
                    u16(x / n_sc_egp_eo),
                    U_SC(x % n_sc_egp_eo)
            };
        }

        bool is_start(const t_state &a) const {
            return a.co == _start.co and a.sc_egp_eo == _start.sc_egp_eo;
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_co = (*mul_co)[a.co];
            const std::array<u8, n_base> &a_sym = conj_mul_egp_eo->mul_sc_sym[a.sc_egp_eo];
            const std::array<U_SC, n_base> &a_sc = conj_mul_egp_eo->mul_sc_sc[a.sc_egp_eo];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                u64 conj_i = conj_base[i][inv_s16[a.sym]];
                u64 sym = a_sym[conj_i];
                a_s[i] = t_state{
                        mul_s16[sym][a.sym],
                        (*conj_co)[a_co[conj_i]][inv_s16[sym]],
                        a_sc[conj_i]
                };
            }
            return a_s;
        }

        std::set<u64> alt(const t_state &a, u64 i) const {
            u64 ss = conj_mul_egp_eo->sc_to_ss[a.sc_egp_eo];
            if (ss == 1) {
                return std::set<u64>{i};
            }
            std::set<u64> a_i{};
            for (u64 s = 0; s < n_s16; s++) {
                if ((ss >> s) & u64(1)) {
                    a_i.insert(state_to_int(t_state{
                            0,
                            (*conj_co)[a.co][s],
                            a.sc_egp_eo
                    }));
                }
            }
            return a_i;
        }

        template<u64 capacity>
        ida_star <g_p0s_solver<_os_e, U_SC, _n_sc_egp_eo, name>, capacity> solve(
                const t_cube &a, u64 max_n_moves = capacity) const {
            return ida_star<g_p0s_solver<_os_e, U_SC, _n_sc_egp_eo, name>, capacity>(*this, a, max_n_moves);
        }
    };

    constexpr char _p0s[] = "p0s";
    typedef g_p0s_solver<os_4_8, u16, 64430, _p0s> p0s_solver;

    // phase 1
    struct p1_solver {
        typedef cube3 t_cube;
        typedef u64 t_state;
        typedef u64 t_hint;

        static constexpr u64 n_cgp = number_gp<os_1x6_2>();
        static constexpr u64 n_ep4 = number_p<4>();
        static constexpr u64 n_ep8 = number_p<8>();
        static constexpr u64 n_state = n_cgp * n_ep4 * n_ep8;
        static constexpr u64 n_super_base = n_cube3_base;
        static constexpr u64 n_base = 10;

        static constexpr std::array<t_cube, n_super_base> super_base = cube3_base;
        static constexpr std::array<const char *, n_super_base> super_base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_index{
                0, 1, 2, 3, 4, 5, 7, 10, 13, 16
        };
        static constexpr std::array<t_cube, n_base> base = array_sub<t_cube, n_super_base, n_base>(
                super_base, base_index);
        static constexpr std::array<const char *, n_base> base_name = array_sub<const char *, n_super_base, n_base>(
                super_base_name, base_index);
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr t_state cube_to_state(const t_cube &a) {
            array_u8<8> cgp = p_to_gp<os_1x6_2>(a.cp);
            auto[ep4, ep8] = p_to_pp<os_4_8>(a.ep);
            return gp_to_int<os_1x6_2>(cgp) * (n_ep4 * n_ep8)
                   + p_to_int<4>(ep4) * n_ep8
                   + p_to_int<8>(ep8);
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
        std::unique_ptr<array_2d < u16, n_cgp, n_base>> mul_cgp;
        std::unique_ptr<array_2d < u32, n_ep4 * n_ep8, n_base>> mul_ep4_ep8;
        std::unique_ptr<array_u2 < n_state>> distance_m3;

        explicit p1_solver(u64 _n_thread) : n_thread(_n_thread) {
            mul_cgp = cache_data<array_2d<u16, n_cgp, n_base>>(
                    "cube3.p1.mul_cgp",
                    [](array_2d<u16, n_cgp, n_base> &t) -> void {
                        for (u64 i = 0; i < n_cgp; i++) {
                            array_u8<8> cgp = int_to_gp<os_1x6_2>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(gp_to_int<os_1x6_2>(cgp * base[j].cp));
                            }
                        }
                    }
            );
            mul_ep4_ep8 = cache_data<array_2d<u32, n_ep4 * n_ep8, n_base>>(
                    "cube3.p1.mul_ep4_ep8",
                    [](array_2d<u32, n_ep4 * n_ep8, n_base> &t) -> void {
                        for (u64 i = 0; i < n_ep4; i++) {
                            array_u8<4> ep4 = int_to_p<4>(i);
                            for (u64 j = 0; j < n_ep8; j++) {
                                array_u8<8> ep8 = int_to_p<8>(j);
                                for (u64 k = 0; k < n_base; k++) {
                                    auto[b_ep4, b_ep8] = p_to_pp<os_4_8>(base[k].ep);
                                    t[i * n_ep8 + j][k] =
                                            u32(p_to_int<4>(ep4 * b_ep4) * n_ep8 + p_to_int<8>(ep8 * b_ep8));
                                }
                            }
                        }
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube3.p1.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<p1_solver>(*this, t, n_thread);
                    }
            );
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_cgp = (*mul_cgp)[a / (n_ep4 * n_ep8)];
            const std::array<u32, n_base> &a_ep4_ep8 = (*mul_ep4_ep8)[a % (n_ep4 * n_ep8)];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                a_s[i] = a_cgp[i] * (n_ep4 * n_ep8) + a_ep4_ep8[i];
            }
            return a_s;
        }

        template<u64 capacity>
        ida_star <p1_solver, capacity> solve(const t_cube &a, u64 max_n_moves = capacity) const {
            return ida_star<p1_solver, capacity>(*this, a, max_n_moves);
        }
    };

    // phase 1 symmetry
    struct p1s_solver {
        typedef cube3 t_cube;

        struct t_state {
            u8 sym;
            u16 cp;
            u16 sc_ep4_ep8;
        };

        typedef u64 t_hint;

        static constexpr u64 n_cp = number_p<8>();
        static constexpr u64 n_ep4 = number_p<4>();
        static constexpr u64 n_ep8 = number_p<8>();
        static constexpr u64 n_sc_ep4_ep8 = 62432;
        static constexpr u64 n_state = n_cp * n_sc_ep4_ep8 / 2;
        static constexpr u64 n_super_base = n_cube3_base;
        static constexpr u64 n_base = 10;

        static constexpr std::array<t_cube, n_super_base> super_base = cube3_base;
        static constexpr std::array<const char *, n_super_base> super_base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_index{
                0, 1, 2, 3, 4, 5, 7, 10, 13, 16
        };
        static constexpr std::array<t_cube, n_base> base = array_sub<t_cube, n_super_base, n_base>(
                super_base, base_index);
        static constexpr std::array<const char *, n_base> base_name = array_sub<const char *, n_super_base, n_base>(
                super_base_name, base_index);
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr array_2d <u8, n_base, n_s16> conj_base =
                generate_table_conj_base<cube3, n_base, n_s16>(base, elements_s16);

        u64 n_thread;
        std::unique_ptr<array_2d < u16, n_cp, n_s16>> conj_cp;
        std::unique_ptr<array_2d < u16, n_cp, n_base>> mul_cp;
        std::unique_ptr<table_conj_mul < u32, u16, u16, n_ep4 * n_ep8, n_sc_ep4_ep8, n_s16, n_base>> conj_mul_ep4_ep8;
        std::unique_ptr<std::array<u8, n_cp>> parity_p8;
        std::unique_ptr<std::array<u8, n_sc_ep4_ep8>> parity_sc_ep4_ep8;
        std::unique_ptr<array_u2 < n_state>> distance_m3;
        t_state _start;

        explicit p1s_solver(u64 _n_thread) : n_thread(_n_thread) {
            conj_cp = cache_data<array_2d<u16, n_cp, n_s16>>(
                    "cube3.p1s.conj_cp",
                    [](array_2d<u16, n_cp, n_s16> &t) -> void {
                        cube3 a = cube3::i();
                        for (u64 i = 0; i < n_cp; i++) {
                            a.cp = int_to_p<8>(i);
                            for (u64 s = 0; s < n_s16; s++) {
                                cube3 b = elements_s16[inv_s16[s]] * a * elements_s16[s];
                                t[i][s] = u16(p_to_int<8>(b.cp));
                            }
                        }
                    }
            );
            mul_cp = cache_data<array_2d<u16, n_cp, n_base>>(
                    "cube3.p1s.mul_cp",
                    [](array_2d<u16, n_cp, n_base> &t) -> void {
                        for (u64 i = 0; i < n_cp; i++) {
                            array_u8<8> cp = int_to_p<8>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(p_to_int<8>(cp * base[j].cp));
                            }
                        }
                    }
            );
            conj_mul_ep4_ep8 = cache_data<table_conj_mul<u32, u16, u16, n_ep4 * n_ep8, n_sc_ep4_ep8, n_s16, n_base>>(
                    "cube3.p1s.conj_mul_ep4_ep8",
                    [](table_conj_mul<u32, u16, u16, n_ep4 * n_ep8, n_sc_ep4_ep8, n_s16, n_base> &cm) -> void {
                        cm.init(
                                [](u32 i) -> std::array<u32, n_s16> {
                                    cube3 a = cube3::i();
                                    a.ep = pp_to_p<os_4_8>({int_to_p<4>(i / n_ep8), int_to_p<8>(i % n_ep8)});
                                    std::array<u32, n_s16> conj_i{};
                                    for (u64 s = 0; s < n_s16; s++) {
                                        cube3 b = elements_s16[inv_s16[s]] * a * elements_s16[s];
                                        auto[ep4, ep8] = p_to_pp<os_4_8>(b.ep);
                                        conj_i[s] = u32(p_to_int<4>(ep4) * n_ep8 + p_to_int<8>(ep8));
                                    }
                                    return conj_i;
                                },
                                [](u32 i) -> std::array<u32, n_base> {
                                    array_u8<4> ep4 = int_to_p<4>(i / n_ep8);
                                    array_u8<8> ep8 = int_to_p<8>(i % n_ep8);
                                    std::array<u32, n_base> mul_i{};
                                    for (u64 j = 0; j < n_base; j++) {
                                        auto[b_ep4, b_ep8] = p_to_pp<os_4_8>(base[j].ep);
                                        mul_i[j] = u32(p_to_int<4>(ep4 * b_ep4) * n_ep8 + p_to_int<8>(ep8 * b_ep8));
                                    }
                                    return mul_i;
                                }
                        );
                    }
            );
            parity_p8 = cache_data<std::array<u8, n_cp>>(
                    "cube3.p1s.parity_p8",
                    [](std::array<u8, n_cp> &t) -> void {
                        for (u64 i = 0; i < n_cp; i++) {
                            t[i] = u8(parity_p<8>(int_to_p<8>(i)));
                        }
                    }
            );
            parity_sc_ep4_ep8 = cache_data<std::array<u8, n_sc_ep4_ep8>>(
                    "cube3.p1s.parity_sc_ep4_ep8",
                    [this](std::array<u8, n_sc_ep4_ep8> &t) -> void {
                        for (u64 i = 0; i < n_sc_ep4_ep8; i++) {
                            u64 ep4_ep8 = conj_mul_ep4_ep8->sc_to_g[i];
                            t[i] = u8((*parity_p8)[ep4_ep8 / n_ep8] ^ (*parity_p8)[ep4_ep8 % n_ep8]);
                        }
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube3.p1s.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<p1s_solver>(*this, t, n_thread);
                    }
            );
            _start = cube_to_state(t_cube::i());
        }

        t_state cube_to_state(const t_cube &a) const {
            auto[ep4, ep8] = p_to_pp<os_4_8>(a.ep);
            u64 ep4_ep8 = p_to_int<4>(ep4) * n_ep8 + p_to_int<8>(ep8);
            auto[sym, sc] = conj_mul_ep4_ep8->g_to_sym_sc(ep4_ep8);
            return t_state{
                    sym,
                    (*conj_cp)[p_to_int<8>(a.cp)][inv_s16[sym]],
                    sc
            };
        }

        u64 state_to_int(const t_state &a) const {
            return (a.cp / 2) * n_sc_ep4_ep8 + a.sc_ep4_ep8;
        }

        t_state int_to_state(u64 x) const {
            u64 cp_a = (x / n_sc_ep4_ep8) * 2;
            u64 sc = x % n_sc_ep4_ep8;
            u64 cp = cp_a ^(*parity_p8)[cp_a] ^(*parity_sc_ep4_ep8)[sc];
            return t_state{
                    0,
                    u16(cp),
                    u16(sc)
            };
        }

        bool is_start(const t_state &a) const {
            return a.cp == _start.cp and a.sc_ep4_ep8 == _start.sc_ep4_ep8;
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_cp = (*mul_cp)[a.cp];
            const std::array<u8, n_base> &a_sym = conj_mul_ep4_ep8->mul_sc_sym[a.sc_ep4_ep8];
            const std::array<u16, n_base> &a_sc = conj_mul_ep4_ep8->mul_sc_sc[a.sc_ep4_ep8];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                u64 conj_i = conj_base[i][inv_s16[a.sym]];
                u64 sym = a_sym[conj_i];
                a_s[i] = t_state{
                        mul_s16[sym][a.sym],
                        (*conj_cp)[a_cp[conj_i]][inv_s16[sym]],
                        a_sc[conj_i]
                };
            }
            return a_s;
        }

        std::set<u64> alt(const t_state &a, u64 i) const {
            u64 ss = conj_mul_ep4_ep8->sc_to_ss[a.sc_ep4_ep8];
            if (ss == 1) {
                return std::set<u64>{i};
            }
            std::set<u64> a_i{};
            for (u64 s = 0; s < n_s16; s++) {
                if ((ss >> s) & u64(1)) {
                    a_i.insert(state_to_int(t_state{
                            0,
                            (*conj_cp)[a.cp][s],
                            a.sc_ep4_ep8
                    }));
                }
            }
            return a_i;
        }

        template<u64 capacity>
        ida_star <p1s_solver, capacity> solve(const t_cube &a, u64 max_n_moves = capacity) const {
            return ida_star<p1s_solver, capacity>(*this, a, max_n_moves);
        }
    };

    typedef combine_solver <p0_solver, p1_solver> _2p_solver;
    typedef combine_solver <p0s_solver, p1s_solver> _2ps_solver;
}

#endif
