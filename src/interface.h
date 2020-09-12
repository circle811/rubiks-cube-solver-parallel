#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "base.h"
#include "group.h"
#include "cube3.h"

namespace cube::_3 {
    struct cube3_interface {
        static constexpr std::array<char, 6> default_scheme = {'U', 'L', 'F', 'R', 'B', 'D'};

//        static constexpr array_u8<54> facelet = {
//                             0,  1,  2,
//                             3,  4,  5,
//                             6,  7,  8,
//                 9, 10, 11, 18, 19, 20, 27, 28, 29, 36, 37, 38,
//                12, 13, 14, 21, 22, 23, 30, 31, 32, 39, 40, 41,
//                15, 16, 17, 24, 25, 26, 33, 34, 35, 42, 43, 44,
//                            45, 46, 47,
//                            48, 49, 50,
//                            51, 52, 53
//        };

        static constexpr array_u8<54> facelet = {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11, 18, 19, 20, 27, 28, 29, 36, 37, 38,
                12, 13, 14, 21, 22, 23, 30, 31, 32, 39, 40, 41,
                15, 16, 17, 24, 25, 26, 33, 34, 35, 42, 43, 44,
                45, 46, 47,
                48, 49, 50,
                51, 52, 53
        };

        static constexpr array_u8<54> init_colors = []() -> array_u8<54> {
            array_u8<54> _init_colors{};
            for (u64 i = 0; i < 54; i++) {
                _init_colors[i] = i / 9;
            }
            return _init_colors;
        }();

        static constexpr array_2d<u8, 8, 3> corner_facelet{
                std::array<u8, 3>{8, 27, 20},
                std::array<u8, 3>{2, 36, 29},
                std::array<u8, 3>{6, 18, 11},
                std::array<u8, 3>{0, 9, 38},
                std::array<u8, 3>{47, 26, 33},
                std::array<u8, 3>{53, 35, 42},
                std::array<u8, 3>{45, 17, 24},
                std::array<u8, 3>{51, 44, 15},
        };

        static constexpr array_2d<u8, 12, 2> edge_facelet{
                std::array<u8, 2>{23, 30},
                std::array<u8, 2>{21, 14},
                std::array<u8, 2>{41, 12},
                std::array<u8, 2>{39, 32},
                std::array<u8, 2>{7, 19},
                std::array<u8, 2>{1, 37},
                std::array<u8, 2>{52, 43},
                std::array<u8, 2>{46, 25},
                std::array<u8, 2>{5, 28},
                std::array<u8, 2>{50, 34},
                std::array<u8, 2>{48, 16},
                std::array<u8, 2>{3, 10}
        };

        static constexpr std::array<u8, 6> center_facelet{
                4, 13, 22, 31, 40, 49
        };

        static constexpr std::array<std::tuple<u8, u8, u8>, 54> facelet_block = []()
                -> std::array<std::tuple<u8, u8, u8>, 54> {
            std::array<std::tuple<u8, u8, u8>, 54> _facelet_block{};
            for (u64 i = 0; i < 8; i++) {
                for (u64 j = 0; j < 3; j++) {
                    std::get<0>(_facelet_block[corner_facelet[i][j]]) = u8(i);
                    std::get<1>(_facelet_block[corner_facelet[i][j]]) = u8(j);
                    std::get<2>(_facelet_block[corner_facelet[i][j]]) = u8(0);
                }
            }
            for (u64 i = 0; i < 12; i++) {
                for (u64 j = 0; j < 2; j++) {
                    std::get<0>(_facelet_block[edge_facelet[i][j]]) = u8(i);
                    std::get<1>(_facelet_block[edge_facelet[i][j]]) = u8(j);
                    std::get<2>(_facelet_block[edge_facelet[i][j]]) = u8(1);
                }
            }
            for (u64 i = 0; i < 6; i++) {
                std::get<0>(_facelet_block[center_facelet[i]]) = u8(i);
                std::get<2>(_facelet_block[center_facelet[i]]) = u8(2);
            }
            return _facelet_block;
        }();

        static constexpr array_u8<54> cube3_to_p(const cube3 &a) {
            array_u8<54> p{};
            for (u64 i = 0; i < 8; i++) {
                for (u64 j = 0; j < 3; j++) {
                    p[corner_facelet[i][j]] = corner_facelet[a.cp[i]][(j + 3 - a.co[i]) % 3];
                }
            }
            for (u64 i = 0; i < 12; i++) {
                for (u64 j = 0; j < 2; j++) {
                    p[edge_facelet[i][j]] = edge_facelet[a.ep[i]][(j + 2 - a.eo[i]) % 2];
                }
            }
            for (u64 i = 0; i < 6; i++) {
                p[center_facelet[i]] = center_facelet[i];
            }
            return p;
        }

        static constexpr cube3 p_to_cube3(const array_u8<54> &p) {
            cube3 a{};
            for (u64 i = 0; i < 8; i++) {
                auto[ii, jj, b] = facelet_block[p[corner_facelet[i][0]]];
                a.cp[i] = ii;
                a.co[i] = (3 - jj) % 3;
            }
            for (u64 i = 0; i < 12; i++) {
                auto[ii, jj, b] = facelet_block[p[edge_facelet[i][0]]];
                a.ep[i] = ii;
                a.eo[i] = (2 - jj) % 2;
            }
            return a;
        }

        static std::string cube3_to_string(const cube3 &a, const std::array<char, 6> &scheme = default_scheme) {
            array_u8<54> colors = init_colors * cube3_to_p(a);
            std::string s = "";
            u64 k = 0;
            for (u64 i = 0; i < 3; i++) {
                s += "    ";
                for (u64 j = 0; j < 3; j++) {
                    s += scheme[colors[facelet[k]]];
                    k++;
                }
                s += "\n";
            }
            s += "\n";
            for (u64 i = 0; i < 3; i++) {
                for (u64 j = 0; j < 12; j++) {
                    if (j > 0 and j % 3 == 0){
                        s += " ";
                    }
                    s += scheme[colors[facelet[k]]];
                    k++;
                }
                s += "\n";
            }
            s += "\n";
            for (u64 i = 0; i < 3; i++) {
                s += "    ";
                for (u64 j = 0; j < 3; j++) {
                    s += scheme[colors[facelet[k]]];
                    k++;
                }
                s += "\n";
            }
            return s;
        }

        static std::string corner_name(u64 i, u64 j) {
            std::string s = "";
            s += default_scheme[init_colors[corner_facelet[i][j]]];
            s += default_scheme[init_colors[corner_facelet[i][(j + 1) % 3]]];
            s += default_scheme[init_colors[corner_facelet[i][(j + 2) % 3]]];
            return s;
        }

        static std::string edge_name(u64 i, u64 j) {
            std::string s = "";
            s += default_scheme[init_colors[edge_facelet[i][j]]];
            s += default_scheme[init_colors[edge_facelet[i][(j + 1) % 2]]];
            return s;
        }

        static std::string center_name(u64 i) {
            std::string s = "";
            s += default_scheme[init_colors[center_facelet[i]]];
            return s;
        }

        std::map<std::tuple<u8, u8, u8>, std::tuple<u8, u8>> color_to_corner;
        std::map<std::tuple<u8, u8>, std::tuple<u8, u8>> color_to_edge;

        cube3_interface() {
            for (u64 i = 0; i < 8; i++) {
                for (u64 j = 0; j < 3; j++) {
                    color_to_corner[{
                            init_colors[corner_facelet[i][j]],
                            init_colors[corner_facelet[i][(j + 1) % 3]],
                            init_colors[corner_facelet[i][(j + 2) % 3]]
                    }] = {u8(i), u8(j)};
                }
            }
            for (u64 i = 0; i < 12; i++) {
                for (u64 j = 0; j < 2; j++) {
                    color_to_edge[{
                            init_colors[edge_facelet[i][j]],
                            init_colors[edge_facelet[i][(j + 1) % 2]]
                    }] = {u8(i), u8(j)};
                }
            }
        }

        std::tuple<bool, std::string, std::array<char, 6>, cube3> parse(const std::string &s) const {
            std::array<char, 54> names{};
            {
                u64 k = 0;
                for (char c: s) {
                    if (std::isalnum(c)) {
                        if (k < 54) {
                            names[facelet[k]] = c;
                            k++;
                        } else {
                            return {false, "char number more than 54", std::array<char, 6>{}, cube3::i()};
                        }
                    } else if (std::isspace(c)) {
                        // ignore
                    } else {
                        return {false, "unknown char", std::array<char, 6>{}, cube3::i()};
                    };
                }
                if (k < 54) {
                    return {false, "char number less than 54", std::array<char, 6>{}, cube3::i()};
                }
            }

            std::array<char, 6> scheme{};
            std::map<char, u8> inv_scheme{};
            for (u64 i = 0; i < 6; i++) {
                char n = names[center_facelet[i]];
                auto it = inv_scheme.find(n);
                if (it != inv_scheme.end()) {
                    std::string e = "center at " + center_name(it->second) + " and " + center_name(i) + " conflict";
                    return {false, e, std::array<char, 6>{}, cube3::i()};
                } else {
                    scheme[i] = n;
                    inv_scheme[n] = i;
                }
            }

            array_u8<54> colors{};
            for (u64 i = 0; i < 54; i++) {
                auto it = inv_scheme.find(names[i]);
                if (it != inv_scheme.end()) {
                    colors[i] = it->second;
                } else {
                    auto[ii, jj, b] = facelet_block[i];
                    if (b == 0) {
                        std::string e = "corner facelet at " + corner_name(ii, jj) + " error";
                        return {false, e, std::array<char, 6>{}, cube3::i()};
                    } else if (b == 1) {
                        std::string e = "edge facelet at " + edge_name(ii, jj) + " error";
                        return {false, e, std::array<char, 6>{}, cube3::i()};
                    } else {
                        assert(0);
                    }
                }
            }

            cube3 a{};

            for (u64 i = 0; i < 8; i++) {
                auto it = color_to_corner.find(
                        {
                                colors[corner_facelet[i][0]],
                                colors[corner_facelet[i][1]],
                                colors[corner_facelet[i][2]]
                        });
                if (it != color_to_corner.end()) {
                    auto[ii, jj] = it->second;
                    a.cp[i] = ii;
                    a.co[i] = (3 - jj) % 3;
                } else {
                    std::string e = "corner at " + corner_name(i, 0) + " error";
                    return {false, e, std::array<char, 6>{}, cube3::i()};
                }
            }

            for (u64 i = 0; i < 12; i++) {
                auto it = color_to_edge.find(
                        {
                                colors[edge_facelet[i][0]],
                                colors[edge_facelet[i][1]]
                        });
                if (it != color_to_edge.end()) {
                    auto[ii, jj] = it->second;
                    a.ep[i] = ii;
                    a.eo[i] = (2 - jj) % 2;
                } else {
                    std::string e = "edge at " + edge_name(i, 0) + " error";
                    return {false, e, std::array<char, 6>{}, cube3::i()};
                }
            }

            {
                std::array<u64, 8> inv{};
                inv.fill(u64(-1));
                for (u64 i = 0; i < 8; i++) {
                    u64 j = a.cp[i];
                    if (inv[j] != u64(-1)) {
                        std::string e = "corner at " + corner_name(inv[j], 0)
                                        + " and " + corner_name(i, 0) + " conflict";
                        return {false, e, std::array<char, 6>{}, cube3::i()};
                    } else {
                        inv[j] = i;
                    }
                }
            }

            {
                std::array<u64, 12> inv{};
                inv.fill(u64(-1));
                for (u64 i = 0; i < 12; i++) {
                    u64 j = a.ep[i];
                    if (inv[j] != u64(-1)) {
                        std::string e = "edge at " + edge_name(inv[j], 0)
                                        + " and " + edge_name(i, 0) + " conflict";
                        return {false, e, std::array<char, 6>{}, cube3::i()};
                    } else {
                        inv[j] = i;
                    }
                }
            }

            if (parity_p<8>(a.cp) ^ parity_p<12>(a.ep)) {
                return {false, "parity error", std::array<char, 6>{}, cube3::i()};
            }

            if (not(a.co == int_to_o<8, 3>(o_to_int<8, 3>(a.co)))) {
                return {false, "corner orientation error", std::array<char, 6>{}, cube3::i()};
            }

            if (not(a.eo == int_to_o<12, 2>(o_to_int<12, 2>(a.eo)))) {
                return {false, "edge orientation error", std::array<char, 6>{}, cube3::i()};
            }

            return {true, "", scheme, a};
        }
    };
}

#endif
