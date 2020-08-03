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
    constexpr u64 n_cp = number_p<8>();
    constexpr u64 n_co = number_o<8, 3>();
    constexpr u64 n_ep = number_p<12>();
    constexpr u64 n_eo = number_o<12, 2>();

    struct cube3a {
        array_u8<8> cp;
        array_u8<8> co;
        array_u8<12> ep;
        array_u8<12> eo;

        static constexpr cube3a i() {
            return cube3a{i_p<8>(), i_o<8, 3>(), i_p<12>(), i_o<12, 2>()};
        }

        constexpr cube3a inv() const {
            const po<8, 3> &rc = inv_po<8, 3>(cp, co);
            const po<12, 2> &re = inv_po<12, 2>(ep, eo);
            return cube3a{std::get<0>(rc), std::get<1>(rc), std::get<0>(re), std::get<1>(re)};
        }

        constexpr cube3a operator*(const cube3a &other) const {
            const po<8, 3> &rc = mul_po<8, 3>(cp, co, other.cp, other.co);
            const po<12, 2> &re = mul_po<12, 2>(ep, eo, other.ep, other.eo);
            return cube3a{std::get<0>(rc), std::get<1>(rc), std::get<0>(re), std::get<1>(re)};
        }

        constexpr bool operator==(const cube3a &other) const {
            return (cp == other.cp) and (co == other.co) and (ep == other.ep) and (eo == other.eo);
        }
    };

    constexpr cube3a U1{{1, 3, 0, 2, 4, 5,  6, 7},
                        {0, 0, 0, 0, 0, 0,  0, 0},
                        {0, 1, 2, 3, 8, 11, 6, 7, 5, 9, 10, 4},
                        {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0}};
    constexpr cube3a D1{{0, 1, 2, 3, 6, 4, 7, 5},
                        {0, 0, 0, 0, 0, 0, 0, 0},
                        {0, 1, 2, 3, 4, 5, 9, 10, 8, 7, 6, 11},
                        {0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0}};
    constexpr cube3a R1{{4, 0, 2, 3, 5, 1, 6, 7},
                        {2, 1, 0, 0, 1, 2, 0, 0},
                        {9, 1, 2, 8, 4, 5, 6, 7, 0, 3, 10, 11},
                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0}};
    constexpr cube3a L1{{0, 1,  3,  7, 4, 5, 2, 6},
                        {0, 0,  1,  2, 0, 0, 2, 1},
                        {0, 11, 10, 3, 4, 5, 6, 7, 8, 9, 1, 2},
                        {0, 0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0}};
    constexpr cube3a F1{{2, 1, 6, 3, 0, 5, 4, 7},
                        {1, 0, 2, 0, 2, 0, 1, 0},
                        {4, 7, 2, 3, 1, 5, 6, 0, 8, 9, 10, 11},
                        {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0}};
    constexpr cube3a B1{{0, 5, 2, 1, 4, 7, 6, 3},
                        {0, 2, 0, 1, 0, 1, 0, 2},
                        {0, 1, 5, 6, 4, 3, 2, 7, 8, 9, 10, 11},
                        {0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,  0}};
    constexpr cube3a AU1{{1, 3, 0, 2, 5, 7,  4,  6},
                         {0, 0, 0, 0, 0, 0,  0,  0},
                         {3, 0, 1, 2, 8, 11, 10, 9, 5, 6, 7, 4},
                         {1, 1, 1, 1, 0, 0,  0,  0, 0, 0, 0, 0}};
    constexpr cube3a AR1{{4, 0,  6,  2, 5, 1, 7, 3},
                         {2, 1,  1,  2, 1, 2, 2, 1},
                         {9, 10, 11, 8, 7, 4, 5, 6, 0, 3, 2, 1},
                         {0, 0,  0,  0, 1, 1, 1, 1, 0, 0, 0, 0}};
    constexpr cube3a AF1{{2, 3, 6, 7, 0, 1, 4, 5},
                         {1, 2, 2, 1, 2, 1, 1, 2},
                         {4, 7, 6, 5, 1, 2, 3, 0, 11, 8, 9, 10},
                         {1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1}};
    constexpr cube3a REF_UD{{4, 5, 6, 7, 0, 1, 2, 3},
                            {3, 3, 3, 3, 3, 3, 3, 3},
                            {0, 1, 2, 3, 7, 6, 5, 4, 9, 8, 11, 10},
                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0}};

    constexpr cube3a U2 = U1 * U1;
    constexpr cube3a D2 = D1 * D1;
    constexpr cube3a R2 = R1 * R1;
    constexpr cube3a L2 = L1 * L1;
    constexpr cube3a F2 = F1 * F1;
    constexpr cube3a B2 = B1 * B1;
    constexpr cube3a AU2 = AU1 * AU1;
    constexpr cube3a AR2 = AR1 * AR1;
    constexpr cube3a AF2 = AF1 * AF1;

    constexpr cube3a U3 = U1.inv();
    constexpr cube3a D3 = D1.inv();
    constexpr cube3a R3 = R1.inv();
    constexpr cube3a L3 = L1.inv();
    constexpr cube3a F3 = F1.inv();
    constexpr cube3a B3 = B1.inv();
    constexpr cube3a AU3 = AU1.inv();
    constexpr cube3a AR3 = AR1.inv();
    constexpr cube3a AF3 = AF1.inv();

    constexpr std::array<cube3a, 18> cube3_base{
            U1, U2, U3, D1, D2, D3,
            R1, R2, R3, L1, L2, L3,
            F1, F2, F3, B1, B2, B3
    };

    const std::array<std::string, cube3_base.size()> cube3_base_name{
            "U", "U2", "U'", "D", "D2", "D'",
            "R", "R2", "R'", "L", "L2", "L'",
            "F", "F2", "F'", "B", "B2", "B'"
    };

    constexpr u64 n_s16 = 16;

    constexpr std::array<cube3a, n_s16> elements_s16 = []() -> std::array<cube3a, n_s16> {
        std::array<cube3a, 2> ref_ud{cube3a::i(), REF_UD};
        std::array<cube3a, 2> ar2{cube3a::i(), AR2};
        std::array<cube3a, 4> au1{cube3a::i(), AU1, AU2, AU3};
        std::array<cube3a, n_s16> _s16_to_cube3a{};
        for (u64 i = 0; i < n_s16; i++) {
            _s16_to_cube3a[i] = ref_ud[i / 8] * ar2[i / 4 % 2] * au1[i % 4];
        }
        return _s16_to_cube3a;
    }();

    constexpr std::array<u8, n_s16> inv_s16 = generate_table<u8, n_s16>(
            [](u64 i) -> u8 {
                return u8(find<cube3a, n_s16>(elements_s16, elements_s16[i].inv()));
            }
    );

    constexpr array_2d <u8, n_s16, n_s16> mul_s16 = generate_table<u8, n_s16, n_s16>(
            [](u64 i, u64 j) -> u8 {
                return u8(find<cube3a, n_s16>(elements_s16, elements_s16[i] * elements_s16[j]));
            }
    );
}

#endif
