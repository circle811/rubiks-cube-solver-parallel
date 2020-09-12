#ifndef _GROUP_H
#define _GROUP_H

#include "base.h"

namespace cube {
    template<typename G, u64 n>
    constexpr std::array<u8, n> generate_table_inv(const std::array<G, n> &elements) {
        std::array<u8, n> t{};
        for (u64 i = 0; i < n; i++) {
            t[i] = array_find<G, n>(elements, elements[i].inv());
        }
        return t;
    }

    template<typename G, u64 n>
    constexpr array_2d <u8, n, n> generate_table_mul(const std::array<G, n> &elements) {
        array_2d <u8, n, n> t{};
        for (u64 i = 0; i < n; i++) {
            for (u64 j = 0; j < n; j++) {
                t[i][j] = array_find<G, n>(elements, elements[i] * elements[j]);
            }
        }
        return t;
    }

    template<typename G, u64 n_base>
    constexpr std::array<u64, n_base> generate_table_base_mask(const std::array<G, n_base> &base, const G &identity) {
        std::array<u64, n_base> t{};
        for (u64 i = 0; i < n_base; i++) {
            for (u64 j = i + 1; j < n_base; j++) {
                G a = base[i] * base[j];
                if (array_find<G, n_base>(base, a) == u64(-1) and not(a == identity)) {
                    t[i] = t[i] | (u64(1) << j);
                    if (not(a == base[j] * base[i])) {
                        t[j] = t[j] | (u64(1) << i);
                    }
                }
            }
        }
        return t;
    }

    template<typename G, u64 n_base, u64 n_sym>
    constexpr array_2d <u8, n_base, n_sym> generate_table_conj_base(
            const std::array<G, n_base> &base, const std::array<G, n_sym> &elements_sym) {
        array_2d <u8, n_base, n_sym> t{};
        for (u64 i = 0; i < n_base; i++) {
            for (u64 s = 0; s < n_sym; s++) {
                t[i][s] = array_find<G, n_base>(base, elements_sym[s].inv() * base[i] * elements_sym[s]);
            }
        }
        return t;
    }

    template<u64 m>
    struct array_u8 {
        std::array<u8, m> _a;

        constexpr const u8 &operator[](u64 i) const {
            return _a[i];
        }

        constexpr u8 &operator[](u64 i) {
            return _a[i];
        }

        static constexpr array_u8<m> i() {
            array_u8<m> p{};
            for (u64 i = 0; i < m; i++) {
                p[i] = i;
            }
            return p;
        }

        constexpr array_u8<m> inv() const {
            array_u8<m> r{};
            for (u64 i = 0; i < m; i++) {
                r[_a[i]] = i;
            }
            return r;
        }

        constexpr array_u8<m> operator*(const array_u8<m> &other) const {
            array_u8<m> r{};
            for (u64 i = 0; i < m; i++) {
                r[i] = _a[other._a[i]];
            }
            return r;
        }

        constexpr bool operator==(const array_u8<m> &other) const {
            return array_eq<u8, m>(_a, other._a);
        }

        std::string to_string() const {
            std::string s = "(";
            for (u64 i = 0; i < m; i++) {
                s += std::to_string(_a[i]);
                if (i < m - 1) {
                    s += " ";
                }
            }
            s += ")";
            return s;
        }
    };

    template<u64 ... _e>
    struct orbit {
        static constexpr u64 n = sizeof...(_e);

        template<u64 m>
        static constexpr std::array<u64, m> e{_e ...};
    };

    template<typename ... _orbits>
    struct orbits {
        typedef std::tuple<array_u8<_orbits::n>...> pp_type;

        static constexpr u64 l = sizeof...(_orbits);
        static constexpr u64 m = array_sum<u64, l>({_orbits::n...});
        static constexpr std::array<u64, l> ns{_orbits::n...};
        static constexpr array_2d <u64, l, m> es{_orbits::template e<m>...};

        static constexpr array_u8<m> x = []() -> array_u8<m> {
            array_u8<m> _x{};
            for (u64 i = 0; i < l; i++) {
                for (u64 j = 0; j < ns[i]; j++) {
                    _x[es[i][j]] = i;
                }
            }
            return _x;
        }();

        static constexpr array_u8<m> y = []() -> array_u8<m> {
            array_u8<m> _y{};
            for (u64 i = 0; i < l; i++) {
                for (u64 j = 0; j < ns[i]; j++) {
                    _y[es[i][j]] = j;
                }
            }
            return _y;
        }();
    };

    template<u64 m>
    constexpr u64 parity_p(const array_u8<m> &p) {
        u64 x = 0;
        for (u64 i = 0; i < m; i++) {
            for (u64 j = i + 1; j < m; j++) {
                if (p[i] > p[j]) {
                    x++;
                }
            }
        }
        return x % 2;
    }

    template<u64 m>
    constexpr u64 number_p() {
        u64 number = 1;
        for (u64 i = 0; i < m; i++) {
            number = number * (m - i);
        }
        return number;
    }

    template<u64 m>
    constexpr array_u8<m> int_to_p(u64 x) {
        u64 x_c = x;
        array_u8<m> p{};
        for (u64 i = m - 1; i < m; i--) {
            p[i] = x_c % (m - i);
            x_c = x_c / (m - i);
        }
        for (u64 i = m - 1; i < m; i--) {
            for (u64 j = m - 1; j > i; j--) {
                if (p[j] >= p[i]) {
                    p[j]++;
                }
            }
        }
        return p;
    }

    template<u64 m>
    constexpr u64 p_to_int(const array_u8<m> &p) {
        array_u8<m> p_c = p;
        for (u64 i = 0; i < m; i++) {
            for (u64 j = i + 1; j < m; j++) {
                if (p_c[j] > p_c[i]) {
                    p_c[j]--;
                }
            }
        }
        u64 x = 0;
        for (u64 i = 0; i < m; i++) {
            x = x * (m - i) + p_c[i];
        }
        return x;
    }

    template<typename os>
    constexpr u64 number_gp() {
        u64 number = 1;
        u64 i = 0;
        for (u64 n: os::ns) {
            for (u64 j = 0; j < n; j++) {
                number = number * (os::m - i) / (j + 1);
                i++;
            }
        }
        return number;
    }

    template<typename os>
    constexpr array_u8<os::m> int_to_gp(u64 x) {
        u64 number = number_gp<os>();
        std::array<u64, os::l> ns_c = os::ns;
        u64 x_c = x;
        array_u8<os::m> gp{};
        for (u64 i = 0; i < os::m; i++) {
            for (u64 j = 0; j < os::l; j++) {
                if (ns_c[j] > 0) {
                    u64 number_s = number * ns_c[j] / (os::m - i);
                    if (x_c < number_s) {
                        number = number_s;
                        ns_c[j]--;
                        gp[i] = j;
                        break;
                    } else {
                        x_c -= number_s;
                    }
                }
            }
        }
        return gp;
    }

    template<typename os>
    constexpr u64 gp_to_int(const array_u8<os::m> &gp) {
        u64 number = number_gp<os>();
        std::array<u64, os::l> ns_c = os::ns;
        u64 x = 0;
        for (u64 i = 0; i < os::m; i++) {
            u64 j = gp[i];
            for (u64 k = 0; k < j; k++) {
                x += number * ns_c[k] / (os::m - i);
            }
            number = number * ns_c[j] / (os::m - i);
            ns_c[j]--;
        }
        return x;
    }

    template<typename os>
    constexpr array_u8<os::m> p_to_gp(const array_u8<os::m> &p) {
        return os::x * p;
    }

    template<typename os>
    constexpr array_u8<os::m> gp_to_p(const array_u8<os::m> &gp) {
        array_u8<os::m> p{};
        array_u8<os::l> c{};
        for (u64 i = 0; i < os::m; i++) {
            u64 j = gp[i];
            p[i] = os::es[j][c[j]];
            c[j]++;
        }
        return p;
    }

    template<typename os, u64 l, u64 i>
    struct _loop_p {
        static constexpr void run(const array_u8<os::m> &p, typename os::pp_type &pp) {
            for (u64 j = 0; j < os::ns[i]; j++) {
                std::get<i>(pp)[j] = os::y[p[os::es[i][j]]];
            }
            _loop_p<os, l, i + 1>::run(p, pp);
        }
    };

    template<typename os, u64 l>
    struct _loop_p<os, l, l> {
        static constexpr void run(const array_u8<os::m> &p, typename os::pp_type &pp) {
        }
    };

    template<typename os>
    constexpr typename os::pp_type p_to_pp(const array_u8<os::m> &p) {
        typename os::pp_type pp{};
        _loop_p<os, os::l, 0>::run(p, pp);
        return pp;
    }

    template<typename os, u64 l, u64 i>
    struct _loop_pp {
        static constexpr void run(const typename os::pp_type &pp, array_u8<os::m> &p) {
            for (u64 j = 0; j < os::ns[i]; j++) {
                p[os::es[i][j]] = os::es[i][std::get<i>(pp)[j]];
            }
            _loop_pp<os, l, i + 1>::run(pp, p);
        }
    };

    template<typename os, u64 l>
    struct _loop_pp<os, l, l> {
        static constexpr void run(const typename os::pp_type &pp, array_u8<os::m> &p) {
        }
    };

    template<typename os>
    constexpr array_u8<os::m> pp_to_p(const typename os::pp_type &pp) {
        array_u8<os::m> p{};
        _loop_pp<os, os::l, 0>::run(pp, p);
        return p;
    }

    template<u64 n>
    constexpr u64 _n_d = number_p<n>();

    template<u64 n>
    using _t_ed = std::array<array_u8<n>, _n_d<n>>;

    template<u64 n>
    constexpr _t_ed<n> _elements_d{};

    template<>
    constexpr _t_ed<2> _elements_d<2>{
            array_u8<2>{0, 1},
            array_u8<2>{1, 0}
    };

    template<>
    constexpr _t_ed<3> _elements_d<3>{
            array_u8<3>{0, 1, 2},
            array_u8<3>{2, 0, 1},
            array_u8<3>{1, 2, 0},
            array_u8<3>{0, 2, 1},
            array_u8<3>{2, 1, 0},
            array_u8<3>{1, 0, 2}
    };

    template<u64 n>
    constexpr std::array<u8, _n_d<n>> _inv_d = generate_table_inv<array_u8<n>, _n_d<n>>(_elements_d<n>);

    template<u64 n>
    constexpr array_2d <u8, _n_d<n>, _n_d<n>> _mul_d = generate_table_mul<array_u8<n>, _n_d<n>>(_elements_d<n>);

    template<u64 m, u64 n>
    constexpr array_u8<m> i_o() {
        return array_u8<m>{};
    }

    template<u64 m, u64 n>
    constexpr array_u8<m> inv_o(const array_u8<m> &o) {
        array_u8<m> r{};
        for (u64 i = 0; i < m; i++) {
            r[i] = _inv_d<n>[o[i]];
        }
        return r;
    }

    template<u64 m, u64 n>
    constexpr array_u8<m> mul_o(const array_u8<m> &o0, const array_u8<m> &o1) {
        array_u8<m> r{};
        for (u64 i = 0; i < m; i++) {
            r[i] = _mul_d<n>[o0[i]][o1[i]];
        }
        return r;
    }

    template<u64 m, u64 n>
    constexpr u64 number_o() {
        u64 number = 1;
        for (u64 i = 0; i < m - 1; i++) {
            number = number * n;
        }
        return number;
    }

    template<u64 m, u64 n>
    constexpr array_u8<m> int_to_o(u64 x) {
        u64 x_c = x;
        array_u8<m> o{};
        u64 y = 0;
        for (u64 i = 1; i < m; i++) {
            u64 d = x_c % n;
            x_c = x_c / n;
            o[m - i - 1] = d;
            y += n - d;
        }
        o[m - 1] = y % n;
        return o;
    }

    template<u64 m, u64 n>
    constexpr u64 o_to_int(const array_u8<m> &o) {
        u64 x = 0;
        for (u64 i = 0; i < m - 1; i++) {
            x = x * n + o[i];
        }
        return x;
    }

    template<typename U_SS, u64 n>
    U_SS _generate_subgroup(const array_2d <u8, n, n> &mul, U_SS subset) {
        while (true) {
            U_SS next_subset = 0;
            for (u64 i = 0; i < n; i++) {
                if ((subset >> i) & U_SS(1)) {
                    for (u64 j = 0; j < n; j++) {
                        if ((subset >> j) & U_SS(1)) {
                            next_subset = next_subset | (U_SS(1) << mul[i][j]);
                        }
                    }
                }
            }
            if (subset != next_subset) {
                subset = next_subset;
            } else {
                break;
            }
        }
        return subset;
    }

    template<typename U_SS, u64 n>
    std::vector<U_SS> generate_table_subgroups(const array_2d <u8, n, n> &mul) {
        std::set<U_SS> t{1};
        for (u64 i = 0; i < n; i++) {
            std::set<U_SS> tt(t);
            for (u64 subgroup: t) {
                t.insert(_generate_subgroup<U_SS, n>(mul, subgroup | (U_SS(1) << i)));
            }
        }
        return std::vector<U_SS>(t.begin(), t.end());
    }

    template<typename U_SS, u64 n>
    std::map<U_SS, std::array<U_SS, n>> generate_table_conj_subgroup(
            const std::vector<U_SS> &subgroups, const std::array<u8, n> &inv, const array_2d <u8, n, n> &mul) {
        std::map<U_SS, std::array<U_SS, n>> t{};
        for (U_SS subgroup: subgroups) {
            for (u64 s = 0; s < n; s++) {
                U_SS subgroup1 = 0;
                for (u64 i = 0; i < n; i++) {
                    if ((subgroup >> i) & U_SS(1)) {
                        subgroup1 = subgroup1 | (U_SS(1) << mul[inv[s]][mul[i][s]]);
                    }
                }
                t[subgroup][s] = subgroup1;
            }
        }
        return t;
    }

    template<typename U_SS, u64 n, u64 n_base>
    std::map<U_SS, u64> generate_table_sym_mask(
            const std::vector<U_SS> &subgroups, const array_2d <u8, n_base, n> &conj_base) {
        std::map<U_SS, u64> t{};
        for (U_SS subgroup: subgroups) {
            u64 mask = 0;
            u64 mark = 0;
            for (u64 i = 0; i < n_base; i++) {
                if (not((mark >> i) & u64(1))) {
                    mask = mask | (u64(1) << i);
                    for (u64 s = 0; s < n; s++) {
                        if ((subgroup >> s) & U_SS(1)) {
                            mark = mark | (u64(1) << conj_base[i][s]);
                        }
                    }
                }

            }
            t[subgroup] = mask;
        }
        return t;
    }

    template<typename U_G, typename U_SC, typename U_SS, u64 n_g, u64 n_sc, u64 n_sym, u64 n_base>
    struct table_conj_mul {
        std::array<u8, n_g> g_to_sym;
        std::array<U_SC, n_g> g_to_sc;
        std::array<U_G, n_sc> sc_to_g;
        std::array<U_SS, n_sc> sc_to_ss;
        array_2d <u8, n_sc, n_base> mul_sc_sym;
        array_2d <U_SC, n_sc, n_base> mul_sc_sc;

        void init(const std::function<std::array<U_G, n_sym>(U_G)> &conj,
                  const std::function<std::array<U_G, n_base>(U_G)> &mul) {
            std::vector<bool> mark(n_g, false);
            u64 k = 0;
            for (u64 i = 0; i < n_g; i++) {
                if (not mark[i]) {
                    u64 ss = 0;
                    std::array<U_G, n_sym> conj_i = conj(i);
                    for (u64 s = 0; s < n_sym; s++) {
                        u64 j = conj_i[s];
                        if (not mark[j]) {
                            mark[j] = true;
                            g_to_sym[j] = u8(s);
                            g_to_sc[j] = U_SC(k);
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
                }
            }
        }

        std::tuple<u8, U_SC> g_to_sym_sc(U_G g) {
            return {g_to_sym[g], g_to_sc[g]};
        }
    };
}

#endif
