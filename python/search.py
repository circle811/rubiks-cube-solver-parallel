import numpy as np

from typing import Callable, Dict, Generator, List, Tuple, TypeVar, Union

A = TypeVar('A')
B = TypeVar('B')


class ArrayM4:
    __slots__ = ['size', 'a']

    def __init__(self, size: int, x: int = 0):
        assert 0 <= x < 4
        y = np.uint64(x)
        for i in range(1, 6):
            y = y | (y << np.uint64(2 ** i))
        self.size = size
        self.a = np.full((size + 31) // 32, y, dtype=np.uint64)

    def __getitem__(self, i: int) -> int:
        j = i // 32
        k = np.uint64(i % 32 * 2)
        return (self.a[j] >> k) & np.uint64(3)

    def __setitem__(self, i: int, x: int) -> None:
        assert 0 <= x < 4
        j = i // 32
        k = np.uint64(i % 32 * 2)
        self.a[j] = (self.a[j] & ~(np.uint64(3) << k)) | (np.uint64(x) << k)

    def to_array(self) -> np.ndarray:
        a = np.zeros(self.size, dtype=np.uint8)
        for i in range(self.size):
            a[i] = self[i]
        return a


def bfs_dict(start: A,
             adj: Callable[[A], Tuple[A, ...]],
             key: Callable[[A], B]
             ) -> Dict[B, int]:
    table_dist = {key(start): 0}
    layer = [start]
    depth = 0
    while len(layer) > 0:
        print(f'bfs_dict: depth={depth}, count={len(layer)}')
        depth += 1
        next_layer = []
        for a in layer:
            for b in adj(a):
                k = key(b)
                if k not in table_dist:
                    table_dist[k] = depth
                    next_layer.append(b)
        layer = next_layer
    return table_dist


def bfs(start: A,
        adj: Callable[[A], Tuple[A, ...]],
        size: int,
        index: Callable[[A], int]
        ) -> np.ndarray:
    table_dist = np.full(size, 255, dtype=np.uint8)
    table_dist[index(start)] = 0
    layer = [start]
    depth = 0
    while len(layer) > 0:
        print(f'bfs: depth={depth}, count={len(layer)}')
        depth += 1
        next_layer = []
        for a in layer:
            for b in adj(a):
                i = index(b)
                if table_dist[i] == 255:
                    table_dist[i] = depth
                    next_layer.append(b)
        layer = next_layer
    return table_dist


def bfs_m3(start: A,
           adj: Callable[[A], Tuple[A, ...]],
           size: int,
           index: Callable[[A], int]
           ) -> ArrayM4:
    table_dist = ArrayM4(size, 3)
    table_dist[index(start)] = 0
    layer = [start]
    depth = 0
    while len(layer) > 0:
        print(f'bfs_m3: depth={depth}, count={len(layer)}')
        depth += 1
        next_layer = []
        for a in layer:
            for b in adj(a):
                i = index(b)
                if table_dist[i] == 3:
                    table_dist[i] = depth % 3
                    next_layer.append(b)
        layer = next_layer
    return table_dist


def dfs(state: A,
        goal: A,
        adj: Callable[[A], Tuple[A, ...]],
        distance: Callable[[A], int],
        max_depth: int,
        moves: List[int],
        count: List[int]
        ) -> Generator[Union[Tuple[int], int], None, None]:
    count[0] += 1
    if len(moves) == max_depth:
        count[1] += 1
        if state == goal:
            count[2] += 1
            yield tuple(moves)
    elif len(moves) + distance(state) <= max_depth:
        for i, a in enumerate(adj(state)):
            moves.append(i)
            yield from dfs(a, goal, adj, distance, max_depth, moves, count)
            moves.pop()


def ida_star(start: A,
             goal: A,
             adj: Callable[[A], Tuple[A, ...]],
             distance: Callable[[A], int]
             ) -> Generator[Union[Tuple[int], int], None, None]:
    max_depth = distance(start)
    while True:
        yield max_depth
        print(f'ida_star: max_depth={max_depth}')
        moves = []
        count = [0, 0, 0]
        yield from dfs(start, goal, adj, distance, max_depth, moves, count)
        print(f'ida_star: max_depth={max_depth}, count={count}')
        max_depth += 1


def dfs_h(state: A,
          goal: A,
          adj: Callable[[A], Tuple[A, ...]],
          distance_m3: Callable[[A], int],
          depth: int
          ) -> int:
    if state == goal:
        return depth
    else:
        d_m3 = (distance_m3(state) - 1) % 3
        for a in adj(state):
            if distance_m3(a) == d_m3:
                return dfs_h(a, goal, adj, distance_m3, depth + 1)
        assert 0


def get_dist_z(start: A,
               goal: A,
               adj: Callable[[A], Tuple[A, ...]],
               distance_m3: Callable[[A], int]
               ) -> int:
    return dfs_h(start, goal, adj, distance_m3, 0)


def get_dist_h(dist_adj: int, dist_m3: int) -> int:
    return dist_adj + (dist_m3 - dist_adj + 1) % 3 - 1


def dfs_m3(state: A,
           goal: A,
           adj: Callable[[A], Tuple[A, ...]],
           distance_m3: Callable[[A], int],
           dist: int,
           max_depth: int,
           moves: List[int],
           count: List[int]
           ) -> Generator[Union[Tuple[int], int], None, None]:
    count[0] += 1
    if len(moves) == max_depth:
        count[1] += 1
        if state == goal:
            count[2] += 1
            yield tuple(moves)
    elif len(moves) + dist <= max_depth:
        for i, a in enumerate(adj(state)):
            moves.append(i)
            dist_a = get_dist_h(dist, distance_m3(a))
            yield from dfs_m3(a, goal, adj, distance_m3, dist_a, max_depth, moves, count)
            moves.pop()


def ida_star_m3(start: A,
                goal: A,
                adj: Callable[[A], Tuple[A, ...]],
                distance_m3: Callable[[A], int]
                ) -> Generator[Union[Tuple[int], int], None, None]:
    dist = get_dist_z(start, goal, adj, distance_m3)
    max_depth = dist
    while True:
        yield max_depth
        print(f'ida_star_m3: max_depth={max_depth}')
        moves = []
        count = [0, 0, 0]
        yield from dfs_m3(start, goal, adj, distance_m3, dist, max_depth, moves, count)
        print(f'ida_star_m3: max_depth={max_depth}, count={count}')
        max_depth += 1
