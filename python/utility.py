import os.path
import pickle
import time

import numpy as np

from typing import Callable, TypeVar

A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')

_cache_dir = os.path.join(os.path.dirname(__file__), 'cache')


def cache_data(name: str, generate: Callable[[], A]) -> A:
    os.makedirs(_cache_dir, exist_ok=True)
    path = os.path.join(_cache_dir, name)
    if os.path.exists(path):
        print(f'load {name} ...')
        with open(path, 'rb') as f:
            a = pickle.load(f)
        print(f'load {name} ok')
    else:
        print(f'computer {name} ...')
        t0 = time.time()
        a = generate()
        t1 = time.time()
        print(f'computer {name} ok, time={t1 - t0:.3f}s')
        print(f'save {name} ...')
        with open(path, 'wb') as f:
            pickle.dump(a, f)
        print(f'save {name} ok')
    return a


def generate_table_mul(m: int, n: int,
                       int_to_a: Callable[[int], A],
                       int_to_b: Callable[[int], B],
                       c_to_int: Callable[[C], int],
                       mul: Callable[[A, B], C]
                       ) -> np.ndarray:
    table_mul = np.zeros((m, n), dtype=np.uint32)
    list_a = [int_to_a(i) for i in range(m)]
    list_b = [int_to_b(j) for j in range(n)]
    for i, a in enumerate(list_a):
        for j, b in enumerate(list_b):
            table_mul[i, j] = c_to_int(mul(a, b))
    return table_mul
