#ifndef _CUBE3_E12_H
#define _CUBE3_E12_H

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube3.h"

namespace cube::_3::e12 {
    template<typename U_G, typename U_SC, typename U_SS, typename U_EXT, u64 n_g, u64 n_sc, u64 n_sym, u64 n_base>
    struct table_conj_mul_ext {
        std::array<u8, n_g> g_to_sym;
        std::array<U_SC, n_g> g_to_sc;
        std::array<U_EXT, n_g> g_to_ext;
        std::array<U_G, n_sc> sc_to_g;
        std::array<U_SS, n_sc> sc_to_ss;
        array_2d <u8, n_sc, n_base> mul_sc_sym;
        array_2d <U_SC, n_sc, n_base> mul_sc_sc;
        array_2d <U_EXT, n_sc, n_base> mul_sc_ext;

        void init(const std::function<std::array<std::tuple<U_G, U_EXT>, n_sym>(U_G)> &conj,
                  const std::function<std::array<U_G, n_base>(U_G)> &mul) {
            std::vector<bool> mark(n_g, false);
            u64 k = 0;
            for (u64 i = 0; i < n_g; i++) {
                if (not mark[i]) {
                    u64 ss = 0;
                    std::array<std::tuple<U_G, U_EXT>, n_sym> conj_i = conj(i);
                    for (u64 s = 0; s < n_sym; s++) {
                        auto[j, ext] = conj_i[s];
                        if (not mark[j]) {
                            mark[j] = true;
                            g_to_sym[j] = u8(s);
                            g_to_sc[j] = U_SC(k);
                            g_to_ext[j] = ext;
                        }
                        if (j == i) {
                            ss = ss | (u64(1) << s);
                        }
                    }
                    sc_to_g[k] = U_G(i);
                    sc_to_ss[k] = U_SS(ss);
                    k++;
                }
            }
            assert (k == n_sc);
            for (u64 i = 0; i < n_sc; i++) {
                std::array<U_G, n_base> mul_i = mul(sc_to_g[i]);
                for (u64 j = 0; j < n_base; j++) {
                    mul_sc_sym[i][j] = g_to_sym[mul_i[j]];
                    mul_sc_sc[i][j] = g_to_sc[mul_i[j]];
                    mul_sc_ext[i][j] = g_to_ext[mul_i[j]];
                }
            }
        }
    };

    constexpr u64 inv_eo_fast(u64 eo) {
        return eo;
    }

    constexpr u64 mul_eo_fast(u64 eo0, u64 eo1) {
        return eo0 ^ eo1;
    }

    struct e12s_solver {
        typedef cube3 t_cube;

        struct t_state {
            u8 sym;
            u16 eo;
            u32 sc_ep;
        };

        typedef u64 t_hint;

        static constexpr u64 n_eo = number_o<12, 2>();
        static constexpr u64 n_ep = number_p<12>();
        static constexpr u64 n_sc_ep = 9985968;
        static constexpr u64 n_state = n_eo * n_sc_ep;
        static constexpr u64 n_base = n_cube3_base;

        static constexpr std::array<t_cube, n_base> base = cube3_base;
        static constexpr std::array<const char *, n_base> base_name = cube3_base_name;

        static constexpr array_2d <u8, n_base, n_s48> conj_base = generate_table_conj<cube3, n_base, n_s48>(
                base, elements_s48);

        u64 n_thread;
        std::unique_ptr<array_2d < u16, n_eo, n_s48>> conj_eo;
        std::unique_ptr<array_2d < u16, n_eo, n_base>> mul_eo;
        std::unique_ptr<table_conj_mul_ext<u32, u32, u64, u16, n_ep, n_sc_ep, n_s48, n_base>> conj_mul_ep;
        std::unique_ptr<array_u2 < n_state>> distance_m3;
        t_state _start;

        explicit e12s_solver(u64 _n_thread) : n_thread(_n_thread) {
            conj_eo = cache_data<array_2d<u16, n_eo, n_s48>>(
                    "cube3.e12s.conj_eo",
                    [this](array_2d<u16, n_eo, n_s48> &t) -> void {
                        cube3 a = cube3::i();
                        for (u64 i = 0; i < n_eo; i++) {
                            a.eo = int_to_o<12, 2>(i);
                            for (u64 s = 0; s < n_s48; s++) {
                                cube3 b = elements_s48[inv_s48[s]] * a * elements_s48[s];
                                t[i][s] = u16(o_to_int<12, 2>(b.eo));
                            }
                        }
                    }
            );
            mul_eo = cache_data<array_2d<u16, n_eo, n_base>>(
                    "cube3.e12s.mul_eo",
                    [](array_2d<u16, n_eo, n_base> &t) -> void {
                        for (u64 i = 0; i < n_eo; i++) {
                            array_u8<12> eo = int_to_o<12, 2>(i);
                            for (u64 j = 0; j < n_base; j++) {
                                t[i][j] = u16(o_to_int<12, 2>(mul_o<12, 2>(eo * base[j].ep, base[j].eo)));
                            }
                        }
                    }
            );
            conj_mul_ep = cache_data<table_conj_mul_ext<u32, u32, u64, u16, n_ep, n_sc_ep, n_s48, n_base>>(
                    "cube3.e12s.conj_mul_ep",
                    [this](table_conj_mul_ext<u32, u32, u64, u16, n_ep, n_sc_ep, n_s48, n_base> &cm) -> void {
                        cm.init(
                                [this](u32 i) -> std::array<std::tuple<u32, u16>, n_s48> {
                                    cube3 a = cube3::i();
                                    a.ep = int_to_p<12>(i);
                                    std::array<std::tuple<u32, u16>, n_s48> conj_i{};
                                    for (u64 s = 0; s < n_s48; s++) {
                                        cube3 b = elements_s48[inv_s48[s]] * a * elements_s48[s];
                                        conj_i[s] = {u32(p_to_int<12>(b.ep)),
                                                     u16(o_to_int<12, 2>(b.eo))};
                                    }
                                    return conj_i;
                                },
                                [](u32 i) -> std::array<u32, n_base> {
                                    array_u8<12> ep = int_to_p<12>(i);
                                    std::array<u32, n_base> mul_i{};
                                    for (u64 j = 0; j < n_base; j++) {
                                        mul_i[j] = p_to_int<12>(ep * base[j].ep);
                                    }
                                    return mul_i;
                                }
                        );
                    }
            );
            distance_m3 = cache_data<array_u2<n_state>>(
                    "cube3.e12s.distance_m3",
                    [this](array_u2<n_state> &t) -> void {
                        bfs<e12s_solver>(*this, t, n_thread);
                    }
            );
            _start = cube_to_state(t_cube::i());
        }

        t_state cube_to_state(const t_cube &a) const {
            u64 ep = p_to_int<12>(a.ep);
            u64 eo = o_to_int<12, 2>(a.eo);
            u64 sym = conj_mul_ep->g_to_sym[ep];
            u64 eo_ext = conj_mul_ep->g_to_ext[ep];
            return t_state{
                    u8(sym),
                    (*conj_eo)[mul_eo_fast(inv_eo_fast(eo_ext), eo)][inv_s48[sym]],
                    conj_mul_ep->g_to_sc[ep]
            };
        }

        u64 state_to_int(const t_state &a) const {
            return a.eo * n_sc_ep + a.sc_ep;
        }

        t_state int_to_state(u64 x) const {
            return t_state{
                    0,
                    u16(x / n_sc_ep),
                    u32(x % n_sc_ep)
            };
        }

        bool is_start(const t_state &a) const {
            return a.eo == _start.eo and a.sc_ep == _start.sc_ep;
        }

        std::array<t_state, n_base> adj(const t_state &a) const {
            const std::array<u16, n_base> &a_eo = (*mul_eo)[a.eo];
            const std::array<u8, n_base> &a_sym = conj_mul_ep->mul_sc_sym[a.sc_ep];
            const std::array<u32, n_base> &a_sc = conj_mul_ep->mul_sc_sc[a.sc_ep];
            const std::array<u16, n_base> &a_ext = conj_mul_ep->mul_sc_ext[a.sc_ep];
            std::array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                u64 conj_i = conj_base[i][inv_s48[a.sym]];
                u64 sym = a_sym[conj_i];
                a_s[i] = t_state{
                        mul_s48[sym][a.sym],
                        (*conj_eo)[mul_eo_fast(inv_eo_fast(a_ext[conj_i]), a_eo[conj_i])][inv_s48[sym]],
                        a_sc[conj_i]
                };
            }
            return a_s;
        }

        std::set<u64> alt(const t_state &a, u64 i) const {
            u64 ss = conj_mul_ep->sc_to_ss[a.sc_ep];
            if (ss == 1) {
                return std::set<u64>{i};
            }
            std::set<u64> a_i{};
            for (u64 s = 0; s < n_s48; s++) {
                if ((ss >> s) & u64(1)) {
                    a_i.insert(state_to_int(t_state{
                            0,
                            (*conj_eo)[a.eo][s],
                            a.sc_ep
                    }));
                }
            }
            return a_i;
        }
    };
}

#endif
