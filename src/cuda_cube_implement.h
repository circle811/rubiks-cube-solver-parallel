#ifndef _CUDA_CUBE_IMPLEMENT_H
#define _CUDA_CUBE_IMPLEMENT_H

#include "cuda_cube.h"

#ifdef CUDA_CUBE

#define CUDA_CHECK(EXP) cuda_check((EXP), __FILE__, __LINE__)

void cuda_check(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "%s:%d cuda error! code=%d, string=\"%s\"\n", file, line, code, cudaGetErrorString(code));
        abort();
    }
}

#endif

namespace cuda_cube {
#ifdef CUDA_CUBE

    void device_malloc(void **p, u64 size) {
        CUDA_CHECK(cudaMalloc(p, size));
    }

    void device_free(void *p) {
        CUDA_CHECK(cudaFree(p));
    }

    void device_set_zero(void *p, u64 size) {
        CUDA_CHECK(cudaMemset(p, 0, size));
    }

    void device_memcpy_h_to_d(void *dst, const void *src, u64 size) {
        CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
    }

    void device_memcpy_d_to_h(void *dst, const void *src, u64 size) {
        CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
    }

#else

    void device_malloc(void **p, u64 size) {
        *p = malloc(size);
    }

    void device_free(void *p) {
        free(p);
    }

    void device_set_zero(void *p, u64 size) {
        memset(p, 0, size);
    }

    void device_memcpy_h_to_d(void *dst, const void *src, u64 size) {
        memcpy(dst, src, size);
    }

    void device_memcpy_d_to_h(void *dst, const void *src, u64 size) {
        memcpy(dst, src, size);
    }

#endif

    template<typename _solver, u64 capacity>
    HOST_DEVICE
    void dfs_one(
            const _solver &s, const ida_star_node<_solver, capacity> &a, u64 n_moves,
            u64 &count, tuple<u64, t_moves<capacity>> &result, volatile bool &stop) {
        typedef ida_star_node<_solver, capacity> node;
        array<node, _solver::n_base * capacity> stack;
        stack[0] = a;
        u64 stack_size = 1;
        while (not stop and stack_size > 0) {
            stack_size--;
            node b = stack[stack_size];
            count++;
            if (b.moves.n == n_moves) {
                if (s.is_start(b.state)) {
                    result = {flag::solution | flag::optimum, b.moves};
                    stop = true;
                    break;
                }
            } else {
                u64 mask = b.moves.n == 0 ? u64(-1) : (*s.base_mask)[b.moves.a[b.moves.n - 1]];
                array<typename _solver::t_state, _solver::n_base> adj_b = s.adj(b.state);
                for (u64 i = _solver::n_base - 1; i < _solver::n_base; i--) {
                    if ((mask >> i) & u64(1)) {
                        typename _solver::t_state state_c = adj_b[i];
                        auto t = get_distance_hint<_solver>::call(s, state_c, b.hint);
                        auto dist_c = t.item0;
                        auto hint_c = t.item1;
                        if (b.moves.n + 1 + dist_c <= n_moves) {
                            node c{
                                    state_c,
                                    hint_c,
                                    t_moves<capacity>{u8(b.moves.n + 1), b.moves.a}
                            };
                            c.moves.a[b.moves.n] = i;
                            stack[stack_size] = c;
                            stack_size++;
                        }
                    }
                }
            }
        }
    }

#ifdef CUDA_CUBE

    template<typename _solver, u64 capacity>
    __global__
    void dfs_multi(
            const _solver &s, const ida_star_node<_solver, capacity> *nodes, u64 n_moves,
            const u64 *tasks, const u64 *split, u64 *count,
            tuple<u64, t_moves<capacity>> *result, volatile bool &stop) {
        u64 thread_id = blockIdx.x;
        u64 start = split[thread_id];
        u64 end = split[thread_id + 1];
        for (u64 i = start; not stop and i < end; i++) {
            u64 j = tasks[i];
            dfs_one(s, nodes[j], n_moves, count[j], result[thread_id], stop);
            u64 f = result[thread_id].item0;
            if (f & flag::solution) {
                break;
            }
        }
    }

    template<typename _solver, u64 capacity>
    void dfs_all(
            const _solver &s, const vector<ida_star_node<_solver, capacity>> &nodes, u64 n_thread, u64 n_moves,
            const vector<u64> &tasks, const vector<u64> &split, vector<u64> &count,
            vector<tuple<u64, t_moves<capacity>>> &result, volatile bool &stop) {
        dfs_multi<_solver, capacity><<<n_thread, 1>>>(
                s, &nodes[0], n_moves, &tasks[0], &split[0], &count[0], &result[0], stop);
    }

#else

    template<typename _solver, u64 capacity>
    void dfs_multi(
            const _solver &s, const ida_star_node<_solver, capacity> *nodes, u64 n_moves,
            const u64 *tasks, const u64 *split, u64 *count,
            tuple<u64, t_moves<capacity>> *result, volatile bool &stop, u64 thread_id) {
        u64 start = split[thread_id];
        u64 end = split[thread_id + 1];
        for (u64 i = start; not stop and i < end; i++) {
            u64 j = tasks[i];
            dfs_one(s, nodes[j], n_moves, count[j], result[thread_id], stop);
            u64 f = result[thread_id].item0;
            if (f & flag::solution) {
                break;
            }
        }
    }

    template<typename _solver, u64 capacity>
    void dfs_all(
            const _solver &s, const vector<ida_star_node<_solver, capacity>> &nodes, u64 n_thread, u64 n_moves,
            const vector<u64> &tasks, const vector<u64> &split, vector<u64> &count,
            vector<tuple<u64, t_moves<capacity>>> &result, volatile bool &stop) {
        for (u64 i = 0; i < n_thread; i++) {
            dfs_multi<_solver, capacity>(
                    s, &nodes[0], n_moves, &tasks[0], &split[0], &count[0], &result[0], stop, i);
        }
    }

#endif

    template
    void dfs_all<optx_solver, 20>(
            const optx_solver &s, const vector<ida_star_node<optx_solver, 20>> &nodes, u64 n_thread, u64 n_moves,
            const vector<u64> &tasks, const vector<u64> &split, vector<u64> &count,
            vector<tuple<u64, t_moves<20>>> &result, volatile bool &stop);

    template
    void dfs_all<opty_solver, 20>(
            const opty_solver &s, const vector<ida_star_node<opty_solver, 20>> &nodes, u64 n_thread, u64 n_moves,
            const vector<u64> &tasks, const vector<u64> &split, vector<u64> &count,
            vector<tuple<u64, t_moves<20>>> &result, volatile bool &stop);
}

#endif
