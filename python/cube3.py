import functools
import random

from typing import NamedTuple, Tuple

from group import (
    T_P, T_O,
    p_number, int_to_p, p_to_int,
    gp_to_int,
    o_number, int_to_o, o_to_int,
    make_orbit, p_to_gp,
    i_p, inv_p, mul_p,
    inv_d2, mul_d2,
    inv_d3, mul_d3)

from utility import cache_data, generate_table_1d

"""
number of corner and edge blocks

      *3      5     *1
   11             8
*2      4     *0

       2             3

 1             0

      *7      6     *5
   10             9
*6      7     *4
"""

N_CP = p_number(8)
N_CO = o_number(8, 3)
N_EP = p_number(12)
N_EO = o_number(12, 2)


class Cube3a(NamedTuple):
    cp: T_P
    co: T_O
    ep: T_P
    eo: T_O

    def inv(self) -> 'Cube3a':
        cp_inv = inv_p(self.cp)
        ep_inv = inv_p(self.ep)
        return Cube3a(cp_inv,
                      tuple(map(inv_d3, mul_p(self.co, cp_inv))),
                      ep_inv,
                      tuple(map(inv_d2, mul_p(self.eo, ep_inv))))

    def __mul__(self, other: 'Cube3a') -> 'Cube3a':
        return Cube3a(mul_p(self.cp, other.cp),
                      tuple(map(mul_d3, mul_p(self.co, other.cp), other.co)),
                      mul_p(self.ep, other.ep),
                      tuple(map(mul_d2, mul_p(self.eo, other.ep), other.eo)))


class Cube3b(NamedTuple):
    cp: int
    co: int
    ep: int
    eo: int


def cube3a_to_b(a: Cube3a) -> Cube3b:
    return Cube3b(p_to_int(a.cp, 8),
                  o_to_int(a.co, 8, 3),
                  p_to_int(a.ep, 12),
                  o_to_int(a.eo, 12, 2))


def cube3b_to_a(b: Cube3b) -> Cube3a:
    return Cube3a(int_to_p(b.cp, 8),
                  int_to_o(b.co, 8, 3),
                  int_to_p(b.ep, 12),
                  int_to_o(b.eo, 12, 2))


def cube3b_to_int(b: Cube3b) -> int:
    return ((int(table_p8_ip[b.cp]) * N_CO + b.co) * N_EP + b.ep) * N_EO + b.eo


def i_cube3a() -> Cube3a:
    return Cube3a(i_p(8), (0,) * 8, i_p(12), (0,) * 12)


def random_cube3a(k: int = 40) -> Cube3a:
    return functools.reduce(Cube3a.__mul__, random.choices(base, k=k), IDENTITY)


def moves_name(moves: Tuple[int, ...]) -> str:
    return ' '.join(base_name[m] for m in moves)


IDENTITY = i_cube3a()
U1 = Cube3a((1, 3, 0, 2, 4, 5, 6, 7), (0, 0, 0, 0, 0, 0, 0, 0),
            (0, 1, 2, 3, 8, 11, 6, 7, 5, 9, 10, 4), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
D1 = Cube3a((0, 1, 2, 3, 6, 4, 7, 5), (0, 0, 0, 0, 0, 0, 0, 0),
            (0, 1, 2, 3, 4, 5, 9, 10, 8, 7, 6, 11), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
R1 = Cube3a((4, 0, 2, 3, 5, 1, 6, 7), (2, 1, 0, 0, 1, 2, 0, 0),
            (9, 1, 2, 8, 4, 5, 6, 7, 0, 3, 10, 11), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
L1 = Cube3a((0, 1, 3, 7, 4, 5, 2, 6), (0, 0, 1, 2, 0, 0, 2, 1),
            (0, 11, 10, 3, 4, 5, 6, 7, 8, 9, 1, 2), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
F1 = Cube3a((2, 1, 6, 3, 0, 5, 4, 7), (1, 0, 2, 0, 2, 0, 1, 0),
            (4, 7, 2, 3, 1, 5, 6, 0, 8, 9, 10, 11), (1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0))
B1 = Cube3a((0, 5, 2, 1, 4, 7, 6, 3), (0, 2, 0, 1, 0, 1, 0, 2),
            (0, 1, 5, 6, 4, 3, 2, 7, 8, 9, 10, 11), (0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0))
AU1 = Cube3a((1, 3, 0, 2, 5, 7, 4, 6), (0, 0, 0, 0, 0, 0, 0, 0),
             (3, 0, 1, 2, 8, 11, 10, 9, 5, 6, 7, 4), (1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0))
AR1 = Cube3a((4, 0, 6, 2, 5, 1, 7, 3), (2, 1, 1, 2, 1, 2, 2, 1),
             (9, 10, 11, 8, 7, 4, 5, 6, 0, 3, 2, 1), (0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0))
AF1 = Cube3a((2, 3, 6, 7, 0, 1, 4, 5), (1, 2, 2, 1, 2, 1, 1, 2),
             (4, 7, 6, 5, 1, 2, 3, 0, 11, 8, 9, 10), (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
REF_RL = Cube3a((2, 3, 0, 1, 6, 7, 4, 5), (3, 3, 3, 3, 3, 3, 3, 3),
                (1, 0, 3, 2, 4, 5, 6, 7, 11, 10, 9, 8), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

U2 = U1 * U1
D2 = D1 * D1
R2 = R1 * R1
L2 = L1 * L1
F2 = F1 * F1
B2 = B1 * B1
AU2 = AU1 * AU1
AR2 = AR1 * AR1
AF2 = AF1 * AF1

U3 = U1.inv()
D3 = D1.inv()
R3 = R1.inv()
L3 = L1.inv()
F3 = F1.inv()
B3 = B1.inv()
AU3 = AU1.inv()
AR3 = AR1.inv()
AF3 = AF1.inv()

base = (U1, U2, U3, D1, D2, D3, R1, R2, R3, L1, L2, L3, F1, F2, F3, B1, B2, B3)
base_name = ('U', 'U2', "U'", 'D', 'D2', "D'", 'R', 'R2', "R'", 'L', 'L2', "L'", 'F', 'F2', "F'", 'B', 'B2', "B'")

IDENTITY_B = cube3a_to_b(IDENTITY)

orbit_p8_ip = make_orbit(((0,), (1,), (2,), (3,), (4,), (5,), (6, 7)))

table_p8_ip = cache_data(
    'table_p8_ip',
    lambda: generate_table_1d(
        p_number(8),
        lambda i: gp_to_int(p_to_gp(int_to_p(i, 8), orbit_p8_ip), orbit_p8_ip.ns)))
