from typing import NamedTuple, Tuple

from group import (
    T_P, T_GP, T_O,
    p_number, int_to_p, p_to_int,
    gp_number, int_to_gp, gp_to_int,
    o_number, int_to_o, o_to_int,
    make_orbit, p_to_gp, p_to_pp,
    inv_p, mul_p,
    mul_d2,
    mul_d3)

from search import bfs_m3

from utility import cache_data, generate_table_2d

from cube3 import (
    Cube3a, moves_name,
    IDENTITY, U1, U2, U3, D1, D2, D3, R2, L2, F2, B2,
    base, base_name,
    table_p8_ip)

"""
G = <U, D, R, L, F, B>
H = <U, D, R2, L2, F2, B2>
I = <>

phase 0: G -> H
phase 1: H -> I
"""

orbit = make_orbit(((0, 1, 2, 3), (4, 5, 6, 7, 8, 9, 10, 11)))

##################################################
# phase 0
N_CO = o_number(8, 3)
N_EGP = gp_number(orbit.ns)
N_EO = o_number(12, 2)


class P0a(NamedTuple):
    co: T_O
    egp: T_GP
    eo: T_O

    def __mul__(self, other: 'Cube3a') -> 'P0a':
        return P0a(tuple(map(mul_d3, mul_p(self.co, other.cp), other.co)),
                   mul_p(self.egp, other.ep),
                   tuple(map(mul_d2, mul_p(self.eo, other.ep), other.eo)))


class P0b(NamedTuple):
    co: int
    egp: int
    eo: int

    def __mul__(self, base_index: int) -> 'P0b':
        return P0b(table_mul_p0_co[self.co, base_index],
                   table_mul_p0_egp[self.egp, base_index],
                   table_mul_p0_eo[self.eo, base_index])

    def adj(self) -> Tuple['P0b']:
        return tuple(self * i for i in range(len(p0_base)))


def cube3a_to_p0a(a: Cube3a) -> P0a:
    return P0a(a.co,
               p_to_gp(a.ep, orbit),
               a.eo)


def p0a_to_b(a: P0a) -> P0b:
    return P0b(o_to_int(a.co, 8, 3),
               gp_to_int(a.egp, orbit.ns),
               o_to_int(a.eo, 12, 2))


def p0b_to_a(b: P0b) -> P0a:
    return P0a(int_to_o(b.co, 8, 3),
               int_to_gp(b.egp, orbit.ns),
               int_to_o(b.eo, 12, 2))


def p0b_to_int(b: P0b) -> int:
    return (b.co * N_EGP + b.egp) * N_EO + b.eo


p0_moves_name = moves_name

p0_base = base
p0_base_name = base_name

IDENTITY_P0B = p0a_to_b(cube3a_to_p0a(IDENTITY))

table_mul_p0_co = cache_data(
    'table_mul_p0_co',
    lambda: generate_table_2d(
        N_CO, len(p0_base),
        lambda i, j: p0a_to_b(p0b_to_a(P0b(i, 0, 0)) * p0_base[j]).co))

table_mul_p0_egp = cache_data(
    'table_mul_p0_egp',
    lambda: generate_table_2d(
        N_EGP, len(p0_base),
        lambda i, j: p0a_to_b(p0b_to_a(P0b(0, i, 0)) * p0_base[j]).egp))

table_mul_p0_eo = cache_data(
    'table_mul_p0_eo',
    lambda: generate_table_2d(
        N_EO, len(p0_base),
        lambda i, j: p0a_to_b(p0b_to_a(P0b(0, 0, i)) * p0_base[j]).eo))

##################################################
# phase 1
N_CP = p_number(8)
N_E4P = p_number(4)
N_E8P = p_number(8)


class P1a(NamedTuple):
    cp: T_P
    e4p: T_P
    e8p: T_P

    def inv(self) -> 'P1a':
        return P1a(inv_p(self.cp),
                   inv_p(self.e4p),
                   inv_p(self.e8p))

    def __mul__(self, other: 'P1a') -> 'P1a':
        return P1a(mul_p(self.cp, other.cp),
                   mul_p(self.e4p, other.e4p),
                   mul_p(self.e8p, other.e8p))


class P1b(NamedTuple):
    cp: int
    e4p: int
    e8p: int

    def __mul__(self, base_index: int) -> 'P1b':
        return P1b(table_mul_p1_cp[self.cp, base_index],
                   table_mul_p1_e4p[self.e4p, base_index],
                   table_mul_p1_e8p[self.e8p, base_index])

    def adj(self) -> Tuple['P1b']:
        return tuple(self * i for i in range(len(p1_base)))


def cube3a_to_p1a(a: Cube3a) -> P1a:
    e4p, e8p = p_to_pp(a.ep, orbit)
    return P1a(a.cp, e4p, e8p)


def p1a_to_b(a: P1a) -> P1b:
    return P1b(p_to_int(a.cp, 8),
               p_to_int(a.e4p, 4),
               p_to_int(a.e8p, 8))


def p1b_to_a(b: P1b) -> P1a:
    return P1a(int_to_p(b.cp, 8),
               int_to_p(b.e4p, 4),
               int_to_p(b.e8p, 8))


def p1b_to_int(b: P1b) -> int:
    return (int(table_p8_ip[b.cp]) * N_E4P + b.e4p) * N_E8P + b.e8p


def p1_moves_name(moves: Tuple[int, ...]) -> str:
    return ' '.join(p1_base_name[m] for m in moves)


p1_base = tuple(map(cube3a_to_p1a, (U1, U2, U3, D1, D2, D3, R2, L2, F2, B2)))
p1_base_name = ('U', 'U2', "U'", 'D', 'D2', "D'", 'R2', 'L2', 'F2', 'B2')

IDENTITY_P1B = p1a_to_b(cube3a_to_p1a(IDENTITY))

table_mul_p1_cp = cache_data(
    'table_mul_p1_cp',
    lambda: generate_table_2d(
        N_CP, len(p1_base),
        lambda i, j: p1a_to_b(p1b_to_a(P1b(i, 0, 0)) * p1_base[j]).cp))

table_mul_p1_e4p = cache_data(
    'table_mul_p1_e4p',
    lambda: generate_table_2d(
        N_E4P, len(p1_base),
        lambda i, j: p1a_to_b(p1b_to_a(P1b(0, i, 0)) * p1_base[j]).e4p))

table_mul_p1_e8p = cache_data(
    'table_mul_p1_e8p',
    lambda: generate_table_2d(
        N_E8P, len(p1_base),
        lambda i, j: p1a_to_b(p1b_to_a(P1b(0, 0, i)) * p1_base[j]).e8p))

##################################################
table_dist_m3_p0 = cache_data(
    'table_dist_m3_p0',
    lambda: bfs_m3(IDENTITY_P0B, P0b.adj, N_CO * N_EGP * N_EO, p0b_to_int))

table_dist_m3_p1 = cache_data(
    'table_dist_m3_p1',
    lambda: bfs_m3(IDENTITY_P1B, P1b.adj, N_CP * N_E4P * N_E8P // 2, p1b_to_int))

##################################################
