#ifndef _CUBE3_OPT_H
#define _CUBE3_OPT_H

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube3.h"
#include "cube3_2p.h"

namespace cube::_3::opt {
    constexpr u64 n_s3 = 3;

    constexpr std::array<cube3, n_s3> elements_s3{
            cube3::i(), AU1 * AF1, (AU1 * AF1).inv()
    };

    typedef orbits <orbit<0>, orbit<1>, orbit<2>, orbit<3>, orbit<4, 5, 6, 7, 8, 9, 10, 11>> os_1x4_8;
    constexpr char _p0es[] = "p0es";
    typedef _2p::g_p0s_solver<os_1x4_8, u32, 1523864, _p0es> p0es_solver;

    struct c8_solver {
        typedef cube3 t_cube;
        typedef u64 t_state;
        typedef u64 t_hint;

        static constexpr u64 n_cp = number_p<8>();
        static constexpr u64 n_co = number_o<8, 3>();
        static constexpr u64 n_state = n_cp * n_co;
        static constexpr u64 n_base = n_cube3_base;

        static constexpr std::array<t_cube, n_base> base = cube3_base;
        static constexpr std::array<const char *, n_base> base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr array_2d <u8, n_base, n_s48> conj_base =
                generate_table_conj_base<cube3, n_base, n_s48>(base, elements_s48);

        static constexpr t_state cube_to_state(const t_cube &a) {
            return p_to_int<8>(a.cp) * n_co + o_to_int<8, 3>(a.co);
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

        static void init_self_sym_subgroup(
                std::array<u8, n_state> &self_sym_subgroup, const std::vector<u64> &subgroups_s48) {
            constexpr u64 n_sc_state = 1841970;
            typedef table_conj_mul<u32, u32, u64, n_state, n_sc_state, n_s48, n_base> t_conj_mul;
            std::unique_ptr<t_conj_mul> conj_mul = std::make_unique<t_conj_mul>();
            conj_mul->init(
                    [](u32 i) -> std::array<u32, n_s48> {
                        cube3 a = cube3::i();
                        a.cp = int_to_p<8>(i / n_co);
                        a.co = int_to_o<8, 3>(i % n_co);
                        std::array<u32, n_s48> conj_i{};
                        for (u64 s = 0; s < n_s48; s++) {
                            cube3 b = elements_s48[inv_s48[s]] * a * elements_s48[s];
                            conj_i[s] = u32(p_to_int<8>(b.cp) * n_co + o_to_int<8, 3>(b.co));
                        }
                        return conj_i;
                    },
                    [](u32 i) -> std::array<u32, n_base> {
                        return {};
                    }
            );

            std::map<u64, std::array<u64, n_s48>> conj_subgroup_s48 =
                    generate_table_conj_subgroup<u64, n_s48>(subgroups_s48, inv_s48, mul_s48);
            std::map<u64, u8> subgroup_to_int{};
            for (u64 i = 0; i < subgroups_s48.size(); i++) {
                subgroup_to_int[subgroups_s48[i]] = i;
            }
            for (u64 i = 0; i < n_state; i++) {
                auto[sym, sc] = conj_mul->g_to_sym_sc(i);
                u64 subgroup0 = conj_mul->sc_to_ss[sc];
                u64 subgroup1 = conj_subgroup_s48.at(subgroup0)[sym];
                self_sym_subgroup[i] = subgroup_to_int[subgroup1];
            }
        }

        u64 n_thread;
        std::vector<u64> subgroups_s48;
        std::map<u64, u64> sym_mask;
        std::unique_ptr<std::array<u8, n_state>> self_sym_subgroup;
        std::unique_ptr<array_2d < u16, n_cp, n_base>> mul_cp;
        std::unique_ptr<array_2d < u16, n_co, n_base>> mul_co;
        std::unique_ptr<array_u2 < n_state>> distance_m3;

        explicit c8_solver(u64 _n_thread) : n_thread(_n_thread) {
            subgroups_s48 = generate_table_subgroups<u64, n_s48>(mul_s48);
            sym_mask = generate_table_sym_mask<u64, n_s48, n_base>(subgroups_s48, conj_base);
            self_sym_subgroup = cache_data<std::array<u8, n_state>>(
                    "cube3.c8.self_sym_subgroup",
                    std::bind(&init_self_sym_subgroup, std::placeholders::_1, std::cref(subgroups_s48))
            );
            mul_cp = cache_data<array_2d<u16, n_cp, n_base>>(
                    "cube3.c8.mul_cp",
                    [](array_2d<u16, n_cp, n_base> &t) -> void {
                        for (u64 i = 0; i < n_cp; i++) {
                            array_u8<8> cp = int_to_p<8>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(p_to_int<8>(cp * base[j].cp));
                            }
                        }
                    }
            );
            mul_co = cache_data<array_2d<u16, n_co, n_base>>(
                    "cube3.c8.mul_co",
                    [](array_2d<u16, n_co, n_base> &t) -> void {
                        for (u64 i = 0; i < n_co; i++) {
                            array_u8<8> co = int_to_o<8, 3>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                mul_o<8, 3>(co * base[j].cp, base[j].co);
                                t[i][j] = u16(o_to_int<8, 3>(mul_o<8, 3>(co * base[j].cp, base[j].co)));
                            }
                        }
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube3.c8.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<c8_solver>(*this, t, n_thread);
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

        u64 get_self_sym_subgroup(const t_state &a) const {
            return subgroups_s48[(*self_sym_subgroup)[a]];
        }

        template<u64 capacity>
        ida_star <c8_solver, capacity> solve(
                const t_cube &a, u64 max_n_moves = capacity, u64 sym_mask_n_moves = 0) const {
            return ida_star<c8_solver, capacity>(*this, a, max_n_moves, sym_mask_n_moves);
        }
    };

    struct opt_solver {
        typedef cube3 t_cube;

        struct t_state {
            std::array<p0es_solver::t_state, n_s3> p0es;
            c8_solver::t_state c8;
        };

        typedef std::array<u8, n_s3 + 1> t_hint;

        static constexpr u64 n_base = n_cube3_base;

        static constexpr std::array<t_cube, n_base> base = cube3_base;
        static constexpr std::array<const char *, n_base> base_name = cube3_base_name;
        static constexpr std::array<u64, n_base> base_mask =
                generate_table_base_mask<t_cube, n_base>(base, t_cube::i());

        static constexpr array_2d <u8, n_base, n_s3> conj_base =
                generate_table_conj_base<cube3, n_base, n_s3>(base, elements_s3);

        u64 n_thread;
        p0es_solver p0es_s;
        c8_solver c8_s;

        explicit opt_solver(u64 _n_thread) : n_thread(_n_thread), p0es_s(_n_thread), c8_s(_n_thread) {
        }

        t_state cube_to_state(const t_cube &a) const {
            t_state b{};
            for (u64 i = 0; i < n_s3; i++) {
                b.p0es[i] = p0es_s.cube_to_state(elements_s3[i].inv() * a * elements_s3[i]);
            }
            b.c8 = c8_s.cube_to_state(a);
            return b;
        }

        bool is_start(const t_state &a) const {
            for (u64 i = 0; i < n_s3; i++) {
                if (not p0es_s.is_start(a.p0es[i])) {
                    return false;
                }
            }
            return c8_s.is_start(a.c8);
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            std::array<t_state, n_base> a_s{};
            for (u64 j = 0; j < n_s3; j++) {
                std::array<p0es_solver::t_state, n_base> adj_p = p0es_s.adj(a.p0es[j]);
                for (u64 i = 0; i < n_base; i++) {
                    a_s[i].p0es[j] = adj_p[conj_base[i][j]];
                }
            }
            std::array<c8_solver::t_state, n_base> adj_c = c8_s.adj(a.c8);
            for (u64 i = 0; i < n_base; i++) {
                a_s[i].c8 = adj_c[i];
            }
            return a_s;
        }

        template<u64 capacity>
        ida_star <opt_solver, capacity> solve(
                const t_cube &a, u64 max_n_moves = capacity, u64 sym_mask_n_moves = 0) const {
            return ida_star<opt_solver, capacity>(*this, a, max_n_moves, sym_mask_n_moves);
        }
    };
}

namespace cube {
    using cube::_3::n_s48;
    using cube::_3::elements_s48;
    using cube::_3::inv_s48;
    using cube::_3::opt::n_s3;
    using cube::_3::opt::p0es_solver;
    using cube::_3::opt::c8_solver;
    using cube::_3::opt::opt_solver;

    template<>
    std::tuple<u64, typename opt_solver::t_hint> get_distance<opt_solver>(
            const opt_solver &s, const typename opt_solver::t_state &a) {
        u64 min_d = u64(-1);
        u64 max_d = 0;
        opt_solver::t_hint r_hint{};
        for (u64 i = 0; i < n_s3; i++) {
            auto[d, h] = get_distance<p0es_solver>(s.p0es_s, a.p0es[i]);
            min_d = std::min(min_d, d);
            max_d = std::max(max_d, d);
            r_hint[i] = u8(h);
        }
        if (max_d == min_d and max_d > 0) {
            max_d++;
        }
        {
            auto[d, h] = get_distance<c8_solver>(s.c8_s, a.c8);
            max_d = std::max(max_d, d);
            r_hint[n_s3] = u8(h);
        }
        return {max_d, r_hint};
    }

    template<>
    std::tuple<u64, typename opt_solver::t_hint> get_distance_hint<opt_solver>(
            const opt_solver &s, const typename opt_solver::t_state &a, const typename opt_solver::t_hint &hint) {
        u64 min_d = u64(-1);
        u64 max_d = 0;
        opt_solver::t_hint r_hint{};
        for (u64 i = 0; i < n_s3; i++) {
            auto[d, h] = get_distance_hint<p0es_solver>(s.p0es_s, a.p0es[i], hint[i]);
            min_d = std::min(min_d, d);
            max_d = std::max(max_d, d);
            r_hint[i] = u8(h);
        }
        if (max_d == min_d and max_d > 0) {
            max_d++;
        }
        {
            auto[d, h] = get_distance_hint<c8_solver>(s.c8_s, a.c8, hint[n_s3]);
            max_d = std::max(max_d, d);
            r_hint[n_s3] = u8(h);
        }
        return {max_d, r_hint};
    }

    template<u64 capacity>
    struct get_sym_mask<c8_solver, capacity> {
        static u64 call(
                const c8_solver &s, const typename c8_solver::t_cube &a, const ida_star_node <c8_solver, capacity> &b) {
            u64 subgroup = s.get_self_sym_subgroup(b.state);
            return s.sym_mask.at(subgroup);
        }
    };

    template<u64 capacity>
    struct get_sym_mask<opt_solver, capacity> {
        static u64 call(
                const opt_solver &s, const typename opt_solver::t_cube &a,
                const ida_star_node <opt_solver, capacity> &b) {
            u64 subgroup0 = s.c8_s.get_self_sym_subgroup(b.state.c8);
            if (subgroup0 == 1) {
                return s.c8_s.sym_mask.at(subgroup0);
            }
            typename opt_solver::t_cube bb = a * moves_to_cube<opt_solver, capacity>(b.moves);
            u64 subgroup1 = 0;
            for (u64 s = 0; s < n_s48; s++) {
                if ((subgroup0 >> s) & u64(1)) {
                    if (elements_s48[inv_s48[s]] * bb * elements_s48[s] == bb) {
                        subgroup1 = subgroup1 | (u64(1) << s);
                    }
                }
            }
            return s.c8_s.sym_mask.at(subgroup1);
        }
    };
}

#endif
