#ifndef _CUDA_CUBE_ADAPTOR_H
#define _CUDA_CUBE_ADAPTOR_H

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube3_opt.h"
#include "cuda_cube.h"

namespace cuda_cube {
    template<typename _h_opt_solver, typename _d_opt_solver>
    struct g_opt_solver_manager {
        _d_opt_solver hd;
        device_ptr<_d_opt_solver> dd;

        g_opt_solver_manager(const _h_opt_solver &s) {
            t_memcpy_h_to_d(dd.get(), &hd);

            t_memcpy_h_to_d(hd.base_mask.get(), &_h_opt_solver::base_mask);
            t_memcpy_h_to_d(hd.conj_base.get(), &_h_opt_solver::conj_base);

            t_memcpy_h_to_d(hd.p0s_s.inv_s16.get(), &cube::_3::inv_s16);
            t_memcpy_h_to_d(hd.p0s_s.mul_s16.get(), &cube::_3::mul_s16);
            t_memcpy_h_to_d(hd.p0s_s.base_mask.get(), &_h_opt_solver::p0s_solver::base_mask);
            t_memcpy_h_to_d(hd.p0s_s.conj_base.get(), &_h_opt_solver::p0s_solver::conj_base);
            t_memcpy_h_to_d(hd.p0s_s.conj_co.get(), s.p0s_s.conj_co.get());
            t_memcpy_h_to_d(hd.p0s_s.mul_co.get(), s.p0s_s.mul_co.get());
            t_memcpy_h_to_d(hd.p0s_s.conj_mul_egp_eo.get(), s.p0s_s.conj_mul_egp_eo.get());
            t_memcpy_h_to_d(hd.p0s_s.distance_m3.get(), s.p0s_s.distance_m3.get());

            t_memcpy_h_to_d(hd.c8_s.base_mask.get(), &cube::_3::opt::c8_solver::base_mask);
            t_memcpy_h_to_d(hd.c8_s.mul_cp.get(), s.c8_s.mul_cp.get());
            t_memcpy_h_to_d(hd.c8_s.mul_co.get(), s.c8_s.mul_co.get());
            t_memcpy_h_to_d(hd.c8_s.distance_m3.get(), s.c8_s.distance_m3.get());
        }

        const _d_opt_solver &get() const {
            return *dd;
        }
    };

    template<typename _h_solver, typename _d_solver, u64 capacity>
    struct cuda_dfs {
        static constexpr char name[] = "cuda";

        typedef _d_solver solver;

        typedef cube::ida_star_node<_h_solver, capacity> h_node;
        typedef ida_star_node<_d_solver, capacity> node;

        static_assert(sizeof(h_node::state) == sizeof(node::state), "");
        static_assert(sizeof(h_node::hint) == sizeof(node::hint), "");
        static_assert(sizeof(h_node::moves) == sizeof(node::moves), "");
        static_assert(sizeof(h_node) == sizeof(node), "");

        const _d_solver &s;
        const u64 n_nodes;
        vector<node> nodes;
        const u64 n_thread;
        vector<u64> tasks;
        vector<u64> split;
        vector<u64> count;
        vector<tuple<u64, t_moves<capacity>>> result;
        std::vector<tuple<u64, t_moves<capacity>>> h_result;
        device_ptr<bool> stop;

        cuda_dfs(const _d_solver &_s, const std::vector<h_node> &h_nodes, u64 _n_thread) :
                s(_s), n_nodes(h_nodes.size()), nodes(n_nodes), n_thread(_n_thread),
                tasks(n_nodes), split(n_thread + 1), count(n_nodes),
                result(n_thread), h_result(n_thread), stop() {
            t_memcpy_h_to_d(&nodes[0], &h_nodes[0], n_nodes);
        }

        std::tuple<u64, cube::t_moves<capacity>> run(
                u64 n_moves,
                const std::vector<u64> &h_tasks, const std::vector<u64> &h_split, std::vector<u64> &h_count) {
            t_memcpy_h_to_d(&tasks[0], &h_tasks[0], n_nodes);
            t_memcpy_h_to_d(&split[0], &h_split[0], n_thread + 1);
            t_set_zero(&count[0], n_nodes);
            t_set_zero(&result[0], n_thread);
            t_set_zero(stop.get(), 1);
            dfs_all<_d_solver, capacity>(s, nodes, n_thread, n_moves, tasks, split, count, result, *stop);
            t_memcpy_d_to_h(&h_count[0], &count[0], n_nodes);
            t_memcpy_d_to_h(&h_result[0], &result[0], n_thread);
            for (u64 i = 0; i < n_thread; i++) {
                u64 f = h_result[i].item0;
                t_moves<capacity> moves = h_result[i].item1;
                if (f & flag::solution) {
                    return {f, reinterpret_cast<cube::t_moves<capacity> &>(moves)};
                }
            }
            return {cube::flag::none, cube::t_moves<capacity>{u8(0), {}}};
        }
    };
}

#endif
