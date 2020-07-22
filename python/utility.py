import os.path
import pickle
import time

import numpy as np

from typing import Callable, TypeVar

A = TypeVar('A')

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


def generate_table_1d(m: int, f: Callable[[int], int]) -> np.ndarray:
    table = np.zeros(m, dtype=np.uint32)
    for i in range(m):
        table[i] = f(i)
    return table


def generate_table_2d(m: int, n: int, f: Callable[[int, int], int]) -> np.ndarray:
    table = np.zeros((m, n), dtype=np.uint32)
    for i in range(m):
        for j in range(n):
            table[i, j] = f(i, j)
    return table
