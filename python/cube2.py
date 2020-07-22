import functools
import random

from typing import NamedTuple, Tuple

from group import (
    T_P, T_O,
    p_number, int_to_p, p_to_int,
    o_number, int_to_o, o_to_int,
    i_p, inv_p, mul_p,
    inv_d3, mul_d3)

from search import bfs, bfs_m3, ida_star, ida_star_m3

from utility import cache_data, generate_table_2d

"""
number of corner blocks

   3    1
2    0

   *    5
6    4
"""

N_CP = p_number(7)
N_CO = o_number(7, 3)


class Cube2a(NamedTuple):
    cp: T_P
    co: T_O

    def inv(self) -> 'Cube2a':
        cp_inv = inv_p(self.cp)
        return Cube2a(cp_inv,
                      tuple(map(inv_d3, mul_p(self.co, cp_inv))))

    def __mul__(self, other: 'Cube2a') -> 'Cube2a':
        return Cube2a(mul_p(self.cp, other.cp),
                      tuple(map(mul_d3, mul_p(self.co, other.cp), other.co)))


class Cube2b(NamedTuple):
    cp: int
    co: int

    def __mul__(self, base_index: int) -> 'Cube2b':
        return Cube2b(table_mul_cube2_cp[self.cp, base_index],
                      table_mul_cube2_co[self.co, base_index])

    def adj(self) -> Tuple['Cube2b']:
        return tuple(self * i for i in range(len(base)))


def cube2a_to_b(a: Cube2a) -> Cube2b:
    return Cube2b(p_to_int(a.cp, 7),
                  o_to_int(a.co, 7, 3))


def cube2b_to_a(b: Cube2b) -> Cube2a:
    return Cube2a(int_to_p(b.cp, 7),
                  int_to_o(b.co, 7, 3))


def cube2b_to_int(b: Cube2b) -> int:
    return b.cp * N_CO + b.co


def i_cube2a() -> Cube2a:
    return Cube2a(i_p(7), (0,) * 7)


def random_cube2a(k: int = 22) -> Cube2a:
    return functools.reduce(Cube2a.__mul__, random.choices(base, k=k), IDENTITY)


def moves_name(moves: Tuple[int, ...]) -> str:
    return ' '.join(base_name[m] for m in moves)


IDENTITY = i_cube2a()
U1 = Cube2a((1, 3, 0, 2, 4, 5, 6), (0, 0, 0, 0, 0, 0, 0))
R1 = Cube2a((4, 0, 2, 3, 5, 1, 6), (2, 1, 0, 0, 1, 2, 0))
F1 = Cube2a((2, 1, 6, 3, 0, 5, 4), (1, 0, 2, 0, 2, 0, 1))
U2 = U1 * U1
R2 = R1 * R1
F2 = F1 * F1
U3 = U1.inv()
R3 = R1.inv()
F3 = F1.inv()

base = (U1, U2, U3, R1, R2, R3, F1, F2, F3)
base_name = ('U', 'U2', "U'", 'R', 'R2', "R'", 'F', 'F2', "F'")

IDENTITY_B = cube2a_to_b(IDENTITY)

table_mul_cube2_cp = generate_table_2d(
    N_CP, len(base),
    lambda i, j: cube2a_to_b(cube2b_to_a(Cube2b(i, 0)) * base[j]).cp)

table_mul_cube2_co = generate_table_2d(
    N_CO, len(base),
    lambda i, j: cube2a_to_b(cube2b_to_a(Cube2b(0, i)) * base[j]).co)

table_dist_cube2 = cache_data(
    'table_dist_cube2',
    lambda: bfs(IDENTITY_B, Cube2b.adj, N_CP * N_CO, cube2b_to_int))

table_dist_m3_cube2 = cache_data(
    'table_dist_m3_cube2',
    lambda: bfs_m3(IDENTITY_B, Cube2b.adj, N_CP * N_CO, cube2b_to_int))


def solve(a: Cube2a, relax: int = 0) -> Tuple[Tuple[int, ...], ...]:
    result = []
    g = ida_star(cube2a_to_b(a),
                 IDENTITY_B,
                 Cube2b.adj,
                 lambda b: table_dist_cube2[cube2b_to_int(b)])
    for u in g:
        if isinstance(u, tuple):
            result.append(u)
        else:
            if len(result) > 0 and u > len(result[0]) + relax:
                break
    return tuple(result)


def solve_m3(a: Cube2a, relax: int = 0) -> Tuple[Tuple[int, ...], ...]:
    result = []
    g = ida_star_m3(cube2a_to_b(a),
                    IDENTITY_B,
                    Cube2b.adj,
                    lambda b: table_dist_m3_cube2[cube2b_to_int(b)])
    for u in g:
        if isinstance(u, tuple):
            result.append(u)
        else:
            if len(result) > 0 and u > len(result[0]) + relax:
                break
    return tuple(result)


def _test():
    print('test cube2 ...')

    a = random_cube2a()
    print(f'state={a}')

    result = solve(a, 2)
    result_m3 = solve_m3(a, 2)
    assert result == result_m3

    for moves in result[:10]:
        b = functools.reduce(Cube2a.__mul__, [base[m] for m in moves], IDENTITY)
        assert a * b == b * a == IDENTITY
        assert a.inv() == b and a == b.inv()
        print(f'{moves_name(moves)}, len={len(moves)}')

    print('test cube2 ok')


if __name__ == '__main__':
    random.seed(0)
    _test()
