#ifndef _CUDA_CUBE_H
#define _CUDA_CUBE_H

#ifdef CUDA_CUBE
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace cuda_cube {
    // base
    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    // cuda
    void device_malloc(void **p, u64 size);

    void device_free(void *p);

    void device_set_zero(void *p, u64 size);

    void device_memcpy_h_to_d(void *dst, const void *src, u64 size);

    void device_memcpy_d_to_h(void *dst, const void *src, u64 size);

    template<typename T>
    void t_malloc(T *&p, u64 n = 1) {
        assert(p == nullptr);
        device_malloc(reinterpret_cast<void **>(&p), sizeof(T) * n);
    }

    template<typename T>
    void t_free(T *&p) {
        device_free(reinterpret_cast<void *>(p));
        p = nullptr;
    }

    template<typename T>
    void t_set_zero(T *p, u64 n = 1) {
        device_set_zero(reinterpret_cast<void *>(p), sizeof(T) * n);
    }

    template<typename T0, typename T1>
    void t_memcpy_h_to_d(T0 *dst, const T1 *src, u64 n = 1) {
        static_assert(sizeof(T0) == sizeof(T1), "");
        device_memcpy_h_to_d(reinterpret_cast<void *>(dst), reinterpret_cast<const void *>(src), sizeof(T0) * n);
    }

    template<typename T0, typename T1>
    void t_memcpy_d_to_h(T0 *dst, const T1 *src, u64 n = 1) {
        static_assert(sizeof(T0) == sizeof(T1), "");
        device_memcpy_d_to_h(reinterpret_cast<void *>(dst), reinterpret_cast<const void *>(src), sizeof(T0) * n);
    }

    // container
    template<typename T, u64 n>
    struct array {
        T _a[n];

        HOST_DEVICE
        constexpr T &operator[](u64 i) {
            return _a[i];
        }

        HOST_DEVICE
        constexpr const T &operator[](u64 i) const {
            return _a[i];
        }
    };

    template<typename T, u64 n0, u64 n1>
    using array_2d = array<array<T, n1>, n0>;

    template<typename ...TS>
    struct tuple {
    };

    template<typename T0, typename T1>
    struct tuple<T0, T1> {
        T0 item0;
        T1 item1;
    };

    template<typename T>
    struct device_ptr {
        T *_p;

        device_ptr() : _p(nullptr) {
            t_malloc(_p);
        }

        device_ptr(const device_ptr<T> &) = delete;

        device_ptr &operator=(const device_ptr<T> &) = delete;

        ~device_ptr() {
            t_free(_p);
        }

        HOST_DEVICE
        T *get() const {
            return _p;
        }

        HOST_DEVICE
        T &operator*() const {
            return *_p;
        }

        HOST_DEVICE
        T *operator->() const {
            return _p;
        }
    };

    template<typename T>
    struct vector {
        T *_p;
        u64 _n;

        vector(u64 __n) : _p(nullptr), _n(__n) {
            t_malloc(_p, __n);
        }

        vector(const vector<T> &) = delete;

        vector &operator=(const vector<T> &) = delete;

        ~vector() {
            t_free(_p);
        }

        HOST_DEVICE
        u64 size() const {
            return _n;
        }

        HOST_DEVICE
        constexpr T &operator[](u64 i) {
            return _p[i];
        }

        HOST_DEVICE
        constexpr const T &operator[](u64 i) const {
            return _p[i];
        }
    };

    // group
    template<typename U_G, typename u32, typename U_SS, u64 n_g, u64 n_sc, u64 n_sym, u64 n_base>
    struct table_conj_mul {
        array<u8, n_g> g_to_sym;
        array<u32, n_g> g_to_sc;
        array<U_G, n_sc> sc_to_g;
        array<U_SS, n_sc> sc_to_ss;
        array_2d<u8, n_sc, n_base> mul_sc_sym;
        array_2d<u32, n_sc, n_base> mul_sc_sc;
    };

    // search
    template<u64 _size>
    struct array_u2 {
        array<u64, (_size + 31) / 32> a;

        HOST_DEVICE
        u64 get(u64 i) const {
            u64 j = i / 32;
            u64 k = i % 32 * 2;
            return (a[j] >> k) & u64(3);
        }
    };

    HOST_DEVICE
    constexpr u64 computer_distance(u64 distance_m3, u64 distance_adj) {
        return distance_adj + (distance_m3 - distance_adj - 3) % 3 - 1;
    }

    template<typename _solver>
    HOST_DEVICE
    tuple<u64, typename _solver::t_hint> get_distance_hint(
            const _solver &s, const typename _solver::t_state &a, const typename _solver::t_hint &hint) {
        u64 d = computer_distance(s.distance_m3->get(s.state_to_int(a)), hint);
        return {d, d};
    }

    template<u64 capacity>
    struct t_moves {
        u8 n;
        array<u8, capacity> a;
    };

    struct flag {
        static constexpr u64 none = 0;
        static constexpr u64 solution = 1;
        static constexpr u64 optimum = 2;
        static constexpr u64 end = 4;
    };

    template<typename _solver, u64 capacity>
    struct ida_star_node {
        typename _solver::t_state state;
        typename _solver::t_hint hint;
        t_moves<capacity> moves;
    };

    template<typename _solver, u64 capacity>
    void dfs_all(
            const _solver &s, const vector<ida_star_node<_solver, capacity>> &nodes, u64 n_thread, u64 n_moves,
            const vector<u64> &tasks, const vector<u64> &split, vector<u64> &count,
            vector<tuple<u64, t_moves<capacity>>> &result, volatile bool &stop);

    // cube3
    constexpr u64 n_cube3_base = 18;
    constexpr u64 n_s16 = 16;
    constexpr u64 n_s3 = 3;

    struct p0es_solver {
        struct t_state {
            u8 sym;
            u16 co;
            u32 sc_egp_eo;
        };

        typedef u64 t_hint;

        static constexpr u64 n_co = 2187;
        static constexpr u64 n_egp = 11880;
        static constexpr u64 n_eo = 2048;
        static constexpr u64 n_sc_egp_eo = 1523864;
        static constexpr u64 n_state = n_co * n_sc_egp_eo;
        static constexpr u64 n_base = n_cube3_base;

        HOST_DEVICE
        static constexpr u64 state_to_int(const t_state &a) {
            return a.co * n_sc_egp_eo + a.sc_egp_eo;
        }

        HOST_DEVICE
        static constexpr bool is_start(const t_state &a) {
            return a.co == 0 and a.sc_egp_eo == 0;
        }

        device_ptr<array<u8, n_s16>> inv_s16;
        device_ptr<array_2d<u8, n_s16, n_s16> > mul_s16;
        device_ptr<array<u64, n_base>> base_mask;
        device_ptr<array_2d<u8, n_base, n_s16> > conj_base;
        device_ptr<array_2d<u16, n_co, n_s16>> conj_co;
        device_ptr<array_2d<u16, n_co, n_base>> mul_co;
        device_ptr<table_conj_mul<u32, u32, u16, n_egp * n_eo, n_sc_egp_eo, n_s16, n_base>> conj_mul_egp_eo;
        device_ptr<array_u2<n_state>> distance_m3;

        HOST_DEVICE
        array<t_state, n_base> adj(const t_state &a) const {
            const array<u16, n_base> &a_co = (*mul_co)[a.co];
            const array<u8, n_base> &a_sym = conj_mul_egp_eo->mul_sc_sym[a.sc_egp_eo];
            const array<u32, n_base> &a_sc = conj_mul_egp_eo->mul_sc_sc[a.sc_egp_eo];
            array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                u64 conj_i = (*conj_base)[i][(*inv_s16)[a.sym]];
                u64 sym = a_sym[conj_i];
                a_s[i] = t_state{
                        (*mul_s16)[sym][a.sym],
                        (*conj_co)[a_co[conj_i]][(*inv_s16)[sym]],
                        a_sc[conj_i]
                };
            }
            return a_s;
        }
    };

    struct c8_solver {
        typedef u64 t_state;
        typedef u64 t_hint;

        static constexpr u64 n_cp = 40320;
        static constexpr u64 n_co = 2187;
        static constexpr u64 n_state = n_cp * n_co;
        static constexpr u64 n_base = n_cube3_base;

        HOST_DEVICE
        static constexpr u64 state_to_int(const t_state &a) {
            return a;
        }

        HOST_DEVICE
        static constexpr bool is_start(const t_state &a) {
            return a == 0;
        }

        device_ptr<array<u64, n_base> > base_mask;
        device_ptr<array_2d<u16, n_cp, n_base>> mul_cp;
        device_ptr<array_2d<u16, n_co, n_base>> mul_co;
        device_ptr<array_u2<n_state>> distance_m3;

        HOST_DEVICE
        array<t_state, n_base> adj(const t_state &a) const {
            const array<u16, n_base> &a_cp = (*mul_cp)[a / n_co];
            const array<u16, n_base> &a_co = (*mul_co)[a % n_co];
            array<t_state, n_base> a_s{};
            for (u64 i = 0; i < n_base; i++) {
                a_s[i] = a_cp[i] * n_co + a_co[i];
            }
            return a_s;
        }
    };

    struct opt_solver {
        struct t_state {
            array<p0es_solver::t_state, n_s3> p0es;
            c8_solver::t_state c8;
        };

        typedef array<u8, n_s3 + 1> t_hint;

        static constexpr u64 n_base = n_cube3_base;

        HOST_DEVICE
        static constexpr bool is_start(const t_state &a) {
            for (u64 i = 0; i < n_s3; i++) {
                if (not p0es_solver::is_start(a.p0es[i])) {
                    return false;
                }
            }
            return c8_solver::is_start(a.c8);
        }

        device_ptr<array<u64, n_base>> base_mask;
        device_ptr<array_2d<u8, n_base, n_s3>> conj_base;
        p0es_solver p0es_s;
        c8_solver c8_s;

        HOST_DEVICE
        array<t_state, n_base> adj(const t_state &a) const {
            array<t_state, n_base> a_s{};
            for (u64 j = 0; j < n_s3; j++) {
                array<p0es_solver::t_state, n_base> adj_p = p0es_s.adj(a.p0es[j]);
                for (u64 i = 0; i < n_base; i++) {
                    a_s[i].p0es[j] = adj_p[(*conj_base)[i][j]];
                }
            }
            array<c8_solver::t_state, n_base> adj_c = c8_s.adj(a.c8);
            for (u64 i = 0; i < n_base; i++) {
                a_s[i].c8 = adj_c[i];
            }
            return a_s;
        }
    };

    template<>
    HOST_DEVICE
    tuple<u64, typename opt_solver::t_hint> get_distance_hint<opt_solver>(
            const opt_solver &s, const typename opt_solver::t_state &a, const typename opt_solver::t_hint &hint);
}

#endif
