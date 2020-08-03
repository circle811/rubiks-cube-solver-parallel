#ifndef _GROUP_H
#define _GROUP_H

#include "base.h"

namespace cube {
    template<u64 _size>
    constexpr u64 sum_a(const std::array<u64, _size> &a) {
        u64 s = 0;
        for (u64 x: a) {
            s += x;
        }
        return s;
    }

    template<u64 ... _e>
    struct orbit {
        static constexpr u64 n = sizeof...(_e);

        template<u64 m>
        static constexpr std::array<u64, m> e{_e ...};
    };

    template<typename ... _orbits>
    struct orbits {
        typedef std::tuple<array_u8 < _orbits::n>...> pp_type;

        static constexpr u64 l = sizeof...(_orbits);
        static constexpr u64 m = sum_a<l>({_orbits::n...});
        static constexpr std::array<u64, l> ns{_orbits::n...};
        static constexpr array_2d <u64, l, m> es{_orbits::template e<m>...};

        static constexpr array_u8 <m> x = []() -> array_u8 <m> {
            array_u8 <m> _x{};
            for (u64 i = 0; i < l; i++) {
                for (u64 j = 0; j < ns[i]; j++) {
                    _x[es[i][j]] = i;
                }
            }
            return _x;
        }();

        static constexpr array_u8 <m> y = []() -> array_u8 <m> {
            array_u8 <m> _y{};
            for (u64 i = 0; i < l; i++) {
                for (u64 j = 0; j < ns[i]; j++) {
                    _y[es[i][j]] = j;
                }
            }
            return _y;
        }();
    };

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
    constexpr array_u8 <os::m> int_to_gp(u64 x) {
        u64 number = number_gp<os>();
        std::array<u64, os::l> ns_c = os::ns;
        u64 x_c = x;
        array_u8 <os::m> gp{};
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
    constexpr u64 gp_to_int(const array_u8 <os::m> &gp) {
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
    constexpr array_u8 <os::m> p_to_gp(const array_u8 <os::m> &p) {
        array_u8 <os::m> gp{};
        for (u64 i = 0; i < os::m; i++) {
            gp[i] = os::x[p[i]];
        }
        return gp;
    }

    template<typename os>
    constexpr array_u8 <os::m> gp_tor_p(const array_u8 <os::m> &gp) {
        array_u8 <os::m> p{};
        array_u8 <os::l> c{};
        for (u64 i = 0; i < os::m; i++) {
            u64 j = gp[i];
            p[i] = os::es[j][c[j]];
            c[j]++;
        }
        return p;
    }

    template<typename os, u64 l, u64 i>
    struct _loop_pp {
        static constexpr void run(const array_u8 <os::m> &p, typename os::pp_type &pp) {
            for (u64 j = 0; j < os::ns[i]; j++) {
                std::get<i>(pp)[j] = os::y[p[os::es[i][j]]];
            }
            _loop_pp<os, l, i + 1>::run(p, pp);
        }
    };

    template<typename os, u64 l>
    struct _loop_pp<os, l, l> {
        static constexpr void run(const array_u8 <os::m> &p, typename os::pp_type &pp) {
        }
    };

    template<typename os>
    constexpr typename os::pp_type p_to_pp(const array_u8 <os::m> &p) {
        typename os::pp_type pp{};
        _loop_pp<os, os::l, 0>::run(p, pp);
        return pp;
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
    constexpr array_u8 <m> int_to_p(u64 x) {
        u64 x_c = x;
        array_u8 <m> p{};
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
    constexpr u64 p_to_int(const array_u8 <m> &p) {
        array_u8 <m> p_c = p;
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

    template<u64 m>
    constexpr u64 parity_p(const array_u8 <m> &p) {
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

    template<u64 m, u64 n>
    constexpr u64 number_o() {
        u64 number = 1;
        for (u64 i = 0; i < m - 1; i++) {
            number = number * n;
        }
        return number;
    }

    template<u64 m, u64 n>
    constexpr array_u8 <m> int_to_o(u64 x) {
        u64 x_c = x;
        array_u8 <m> o{};
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
    constexpr u64 o_to_int(const array_u8 <m> &o) {
        u64 x = 0;
        for (u64 i = 0; i < m - 1; i++) {
            x = x * n + o[i];
        }
        return x;
    }

    template<u64 m>
    constexpr array_u8 <m> i_p() {
        array_u8 <m> p{};
        for (u64 i = 0; i < m; i++) {
            p[i] = i;
        }
        return p;
    }

    template<u64 m>
    constexpr array_u8 <m> inv_p(const array_u8 <m> &p) {
        array_u8 <m> q{};
        for (u64 i = 0; i < m; i++) {
            q[p[i]] = i;
        }
        return q;
    }

    template<u64 m>
    constexpr array_u8 <m> mul_p(const array_u8 <m> &p, const array_u8 <m> &q) {
        array_u8 <m> r{};
        for (u64 i = 0; i < m; i++) {
            r[i] = p[q[i]];
        }
        return r;
    }

    template<u64 n>
    constexpr u64 _np = number_p<n>();

    template<u64 n>
    using _ed = std::array<array_u8 < n>, _np<n>>;

    template<u64 n>
    constexpr _ed<n> _elements_d{};

    template<>
    constexpr _ed<2> _elements_d<2>{
            array_u8 < 2 > {0, 1},
            array_u8 < 2 > {1, 0}
    };

    template<>
    constexpr _ed<3> _elements_d<3>{
            array_u8 < 3 > {0, 1, 2},
            array_u8 < 3 > {2, 0, 1},
            array_u8 < 3 > {1, 2, 0},
            array_u8 < 3 > {0, 2, 1},
            array_u8 < 3 > {2, 1, 0},
            array_u8 < 3 > {1, 0, 2}
    };

    template<u64 n>
    constexpr std::array<u8, _np<n>> _inv_d = generate_table<u8, _np<n>>(
            [](u64 i) -> u8 {
                return u8(find<array_u8<n>, _np<n>>(
                        _elements_d<n>,
                        inv_p<n>(_elements_d<n>[i])));
            }
    );

    template<u64 n>
    constexpr array_2d <u8, _np<n>, _np<n>> _mul_d = generate_table<u8, _np<n>, _np<n>>(
            [](u64 i, u64 j) -> u8 {
                return u8(find<array_u8<n>, _np<n>>(
                        _elements_d<n>,
                        mul_p<n>(_elements_d<n>[i], _elements_d<n>[j])));
            }
    );

    template<u64 m, u64 n>
    constexpr array_u8 <m> i_o() {
        return array_u8 < m > {};
    }

    template<u64 m, u64 n>
    using po = std::tuple<array_u8 < m>, array_u8 <m >>;

    template<u64 m, u64 n>
    constexpr po<m, n> inv_po(const array_u8 <m> &p, const array_u8 <m> &o) {
        po<m, n> r{};
        for (u64 i = 0; i < m; i++) {
            std::get<0>(r)[p[i]] = i;
            std::get<1>(r)[p[i]] = _inv_d<n>[o[i]];
        }
        return r;
    }

    template<u64 m, u64 n>
    constexpr po<m, n> mul_po(const array_u8 <m> &p0, const array_u8 <m> &o0,
                              const array_u8 <m> &p1, const array_u8 <m> &o1) {
        po<m, n> r{};
        for (u64 i = 0; i < m; i++) {
            std::get<0>(r)[i] = p0[p1[i]];
            std::get<1>(r)[i] = _mul_d<n>[o0[p1[i]]][o1[i]];
        }
        return r;
    }
}

#endif
