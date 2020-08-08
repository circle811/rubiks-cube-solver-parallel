#ifndef _CUBE3_H
#define _CUBE3_H

#include "base.h"
#include "group.h"
#include "search.h"

//  number of corner and edge blocks
//
//        *3      5     *1
//     11             8
//  *2      4     *0
//
//         2             3
//
//   1             0
//
//        *7      6     *5
//     10             9
//  *6      7     *4

namespace cube::_3 {
    struct cube3 {
        array_u8<8> cp;
        array_u8<8> co;
        array_u8<12> ep;
        array_u8<12> eo;

        static constexpr cube3 i() {
            return cube3{
                    array_u8<8>::i(),
                    i_o<8, 3>(),
                    array_u8<12>::i(),
                    i_o<12, 2>()
            };
        }

        constexpr cube3 inv() const {
            array_u8<8> cp_inv = cp.inv();
            array_u8<12> ep_inv = ep.inv();
            return cube3{
                    cp_inv,
                    inv_o<8, 3>(co) * cp_inv,
                    ep_inv,
                    inv_o<12, 2>(eo) * ep_inv
            };
        }

        constexpr cube3 operator*(const cube3 &other) const {
            return cube3{
                    cp * other.cp,
                    mul_o<8, 3>(co * other.cp, other.co),
                    ep * other.ep,
                    mul_o<12, 2>(eo * other.ep, other.eo)
            };
        }

        constexpr bool operator==(const cube3 &other) const {
            return cp == other.cp and co == other.co and ep == other.ep and eo == other.eo;
        }

        std::string to_string() const {
            return "cube3(" + cp.to_string() + " " + co.to_string() + " "
                   + ep.to_string() + " " + eo.to_string() + ")";
        }
    };

    constexpr cube3 U1{{1, 3, 0, 2, 4, 5,  6, 7},
                       {0, 0, 0, 0, 0, 0,  0, 0},
                       {0, 1, 2, 3, 8, 11, 6, 7, 5, 9, 10, 4},
                       {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0}};
    constexpr cube3 D1{{0, 1, 2, 3, 6, 4, 7, 5},
                       {0, 0, 0, 0, 0, 0, 0, 0},
                       {0, 1, 2, 3, 4, 5, 9, 10, 8, 7, 6, 11},
                       {0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0}};
    constexpr cube3 R1{{4, 0, 2, 3, 5, 1, 6, 7},
                       {2, 1, 0, 0, 1, 2, 0, 0},
                       {9, 1, 2, 8, 4, 5, 6, 7, 0, 3, 10, 11},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0}};
    constexpr cube3 L1{{0, 1,  3,  7, 4, 5, 2, 6},
                       {0, 0,  1,  2, 0, 0, 2, 1},
                       {0, 11, 10, 3, 4, 5, 6, 7, 8, 9, 1, 2},
                       {0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0}};
    constexpr cube3 F1{{2, 1, 6, 3, 0, 5, 4, 7},
                       {1, 0, 2, 0, 2, 0, 1, 0},
                       {4, 7, 2, 3, 1, 5, 6, 0, 8, 9, 10, 11},
                       {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0}};
    constexpr cube3 B1{{0, 5, 2, 1, 4, 7, 6, 3},
                       {0, 2, 0, 1, 0, 1, 0, 2},
                       {0, 1, 5, 6, 4, 3, 2, 7, 8, 9, 10, 11},
                       {0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,  0}};
    constexpr cube3 AU1{{1, 3, 0, 2, 5, 7,  4,  6},
                        {0, 0, 0, 0, 0, 0,  0,  0},
                        {3, 0, 1, 2, 8, 11, 10, 9, 5, 6, 7, 4},
                        {1, 1, 1, 1, 0, 0,  0,  0, 0, 0, 0, 0}};
    constexpr cube3 AR1{{4, 0,  6,  2, 5, 1, 7, 3},
                        {2, 1,  1,  2, 1, 2, 2, 1},
                        {9, 10, 11, 8, 7, 4, 5, 6, 0, 3, 2, 1},
                        {0, 0,  0,  0, 1, 1, 1, 1, 0, 0, 0, 0}};
    constexpr cube3 AF1{{2, 3, 6, 7, 0, 1, 4, 5},
                        {1, 2, 2, 1, 2, 1, 1, 2},
                        {4, 7, 6, 5, 1, 2, 3, 0, 11, 8, 9, 10},
                        {1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1}};
    constexpr cube3 REF_UD{{4, 5, 6, 7, 0, 1, 2, 3},
                           {3, 3, 3, 3, 3, 3, 3, 3},
                           {0, 1, 2, 3, 7, 6, 5, 4, 9, 8, 11, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0}};

    constexpr cube3 U2 = U1 * U1;
    constexpr cube3 D2 = D1 * D1;
    constexpr cube3 R2 = R1 * R1;
    constexpr cube3 L2 = L1 * L1;
    constexpr cube3 F2 = F1 * F1;
    constexpr cube3 B2 = B1 * B1;

    constexpr cube3 U3 = U1.inv();
    constexpr cube3 D3 = D1.inv();
    constexpr cube3 R3 = R1.inv();
    constexpr cube3 L3 = L1.inv();
    constexpr cube3 F3 = F1.inv();
    constexpr cube3 B3 = B1.inv();

    constexpr u64 n_cube3_base = 18;

    constexpr std::array<cube3, n_cube3_base> cube3_base{
            U1, U2, U3, D1, D2, D3,
            R1, R2, R3, L1, L2, L3,
            F1, F2, F3, B1, B2, B3
    };

    constexpr std::array<const char *, n_cube3_base> cube3_base_name{
            "U", "U2", "U'", "D", "D2", "D'",
            "R", "R2", "R'", "L", "L2", "L'",
            "F", "F2", "F'", "B", "B2", "B'"
    };

    constexpr u64 n_s48 = 48;

    constexpr std::array<cube3, n_s48> elements_s48 = []() -> std::array<cube3, n_s48> {
        std::array<cube3, 2> ref_ud{cube3::i(), REF_UD};
        std::array<cube3, 2> ar2{cube3::i(), AR1 * AR1};
        std::array<cube3, 4> au1{cube3::i(), AU1, AU1 * AU1, AU1.inv()};
        std::array<cube3, 3> au1_af1{cube3::i(), AU1 * AF1, (AU1 * AF1).inv()};
        std::array<cube3, n_s48> _elements_s48{};
        for (u64 i = 0; i < n_s48; i++) {
            _elements_s48[i] = ref_ud[i / 24] * ar2[i / 12 % 2] * au1[i / 3 % 4] * au1_af1[i % 3];
        }
        return _elements_s48;
    }();

    constexpr std::array<u8, n_s48> inv_s48 = generate_table_inv<cube3, n_s48>(elements_s48);

    constexpr array_2d <u8, n_s48, n_s48> mul_s48 = generate_table_mul<cube3, n_s48>(elements_s48);

    constexpr u64 n_s16 = 16;

    constexpr std::array<cube3, n_s16> elements_s16 = []() -> std::array<cube3, n_s16> {
        std::array<cube3, 2> ref_ud{cube3::i(), REF_UD};
        std::array<cube3, 2> ar2{cube3::i(), AR1 * AR1};
        std::array<cube3, 4> au1{cube3::i(), AU1, AU1 * AU1, AU1.inv()};
        std::array<cube3, n_s16> _elements_s16{};
        for (u64 i = 0; i < n_s16; i++) {
            _elements_s16[i] = ref_ud[i / 8] * ar2[i / 4 % 2] * au1[i % 4];
        }
        return _elements_s16;
    }();

    constexpr std::array<u8, n_s16> inv_s16 = generate_table_inv<cube3, n_s16>(elements_s16);

    constexpr array_2d <u8, n_s16, n_s16> mul_s16 = generate_table_mul<cube3, n_s16>(elements_s16);
}

#endif
