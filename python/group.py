from typing import Callable, NamedTuple, Optional, Tuple

T_P = Tuple[int, ...]
T_GP = Tuple[int, ...]
T_PP = Tuple[T_P, ...]
T_O = Tuple[int, ...]


def p_number(m: int, n: Optional[int] = None) -> int:
    if n is None:
        n = m
    number = 1
    for i in range(n):
        number = number * (m - i)
    return number


def int_to_p(x: int, m: int, n: Optional[int] = None) -> T_P:
    if n is None:
        n = m
    x_c = x
    p = [-1] * n
    for i in range(n - 1, -1, -1):
        p[i] = x_c % (m - i)
        x_c = x_c // (m - i)
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            if p[j] >= p[i]:
                p[j] += 1
    return tuple(p)


def p_to_int(p: T_P, m: int, n: Optional[int] = None) -> int:
    if n is None:
        n = m
    p_c = list(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p_c[j] > p_c[i]:
                p_c[j] -= 1
    x = 0
    for i in range(n):
        x = x * (m - i) + p_c[i]
    return x


def gp_number(ns: Tuple[int, ...]) -> int:
    m = sum(ns)
    number = 1
    i = 0
    for n in ns:
        for j in range(n):
            number = number * (m - i) // (j + 1)
            i += 1
    return number


def int_to_gp(x: int, ns: Tuple[int, ...]) -> T_GP:
    m = sum(ns)
    number = gp_number(ns)
    ns_c = list(ns)
    x_c = x % number
    gp = []
    for i in range(m):
        for j in range(len(ns_c)):
            if ns_c[j] > 0:
                number_s = number * ns_c[j] // (m - i)
                if x_c < number_s:
                    number = number_s
                    ns_c[j] -= 1
                    gp.append(j)
                    break
                else:
                    x_c -= number_s
    return tuple(gp)


def gp_to_int(gp: T_GP, ns: Tuple[int, ...]) -> int:
    m = sum(ns)
    number = gp_number(ns)
    ns_c = list(ns)
    x = 0
    for i, j in enumerate(gp):
        for k in range(j):
            x += number * ns_c[k] // (m - i)
        number = number * ns_c[j] // (m - i)
        ns_c[j] -= 1
    return x


def o_number(m: int, n: int) -> int:
    return n ** (m - 1)


def int_to_o(x: int, m: int, n: int) -> T_O:
    x_c = x
    o = [-1] * m
    a = 0
    for i in range(m - 2, -1, -1):
        o[i] = x_c % n
        x_c = x_c // n
        a += o[i]
    o[-1] = -a % n
    return tuple(o)


def o_to_int(o: T_O, m: int, n: int) -> int:
    x = 0
    for i in range(m - 1):
        x = x * n + o[i]
    return x


class Orbit(NamedTuple):
    ns: Tuple[int, ...]
    os: Tuple[Tuple[int, ...], ...]
    x: Tuple[int, ...]
    y: Tuple[int, ...]


def make_orbit(os: Tuple[Tuple[int, ...], ...]) -> Orbit:
    ns = tuple(map(len, os))
    m = sum(ns)
    x = [-1] * m
    y = [-1] * m
    for i, o in enumerate(os):
        for j, k in enumerate(o):
            x[k] = i
            y[k] = j
    return Orbit(ns, os, tuple(x), tuple(y))


def p_to_gp(p: T_P, orbit: Orbit) -> T_GP:
    return tuple(orbit.x[i] for i in p)


def p_to_pp(p: T_P, orbit: Orbit) -> T_PP:
    return tuple(tuple(orbit.y[p[i]] for i in o)
                 for o in orbit.os)


def i_p(m: int) -> T_P:
    return tuple(range(m))


def inv_p(p: T_P) -> T_P:
    m = len(p)
    q = [-1] * m
    for i in range(m):
        q[p[i]] = i
    return tuple(q)


def mul_p(p: T_P, q: T_P) -> T_P:
    return tuple(p[i] for i in q)


def generate_pi_function(table_pi_to_p: Tuple[T_P, ...]) -> (
        Tuple[Callable[[int], T_P],
              Callable[[T_P], int],
              Callable[[int], int],
              Callable[[int, int], int]]):
    table_p_to_pi = {p: i for i, p in enumerate(table_pi_to_p)}
    table_inv = [table_p_to_pi[inv_p(p)] for p in table_pi_to_p]
    table_mul = [[table_p_to_pi[mul_p(p, q)] for q in table_pi_to_p]
                 for p in table_pi_to_p]

    def pi_to_p(x: int) -> T_P:
        return table_pi_to_p[x]

    def p_to_pi(p: T_P) -> int:
        return table_p_to_pi[p]

    def inv_pi(x: int) -> int:
        return table_inv[x]

    def mul_pi(x: int, y: int) -> int:
        return table_mul[x][y]

    return pi_to_p, p_to_pi, inv_pi, mul_pi


_, _, inv_d2, mul_d2 = generate_pi_function(
    ((0, 1),
     (1, 0)))

_, _, inv_d3, mul_d3 = generate_pi_function(
    ((0, 1, 2),
     (2, 0, 1),
     (1, 2, 0),
     (0, 2, 1),
     (2, 1, 0),
     (1, 0, 2)))
