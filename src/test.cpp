#include "cxxopts.hpp"

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube2.h"
#include "cube3.h"
#include "cube3_2p.h"
#include "cube3_opt.h"
#include "cube3_e12.h"

using namespace cube;
using namespace cube::_2;
using namespace cube::_3;
using namespace cube::_3::_2p;
using namespace cube::_3::opt;
using namespace cube::_3::e12;

template<typename _solver>
struct solved_check {
    static void call(const _solver &s, const typename _solver::t_cube &a) {
        assert(a == _solver::t_cube::i());
    }
};

template<typename _solver>
struct partial_check {
    static void call(const _solver &s, const typename _solver::t_cube &a) {
        assert(s.is_start(s.cube_to_state(a)));
    }
};

template<typename _solver, u64 capacity, typename check>
void test_one(u64 n_thread, u64 seed, u64 n_cube, u64 rand_n_moves, u64 max_n_moves, u64 max_n_solution) {
    std::cout << "##################################################" << std::endl;
    _solver s(n_thread);
    random_moves<capacity> rand(_solver::n_base, seed);
    for (u64 i = 0; i < n_cube; i++) {
        t_moves<capacity> moves_g{0, {}};
        if (i > 0) {
            moves_g = rand(rand_n_moves);
        }
        typename _solver::t_cube a = moves_to_cube<_solver, capacity>(moves_g);
        std::cout << "generation: " << moves_to_string<_solver, capacity>(moves_g) << std::endl;
        auto it = s.template solve<capacity>(a, max_n_moves);
        u64 n_solution = 0;
        while (n_solution < max_n_solution) {
            auto[f, moves] = it();
            if (f & flag::solution) {
                typename _solver::t_cube b = a * moves_to_cube<_solver, capacity>(moves);
                std::cout << f << " " << moves_to_string<_solver, capacity>(moves) << std::endl;
                check::call(s, b);
                n_solution++;
            } else if (f & flag::end) {
                break;
            }
        }
        std::cout << std::endl;
    }
}

void test(u64 n_thread, u64 seed, bool full) {
    constexpr u64 capacity = 20;
    u64 n_cube = 3;
    u64 rand_n_moves = 15;
    u64 max_n_moves = 14;
    u64 max_n_solution = 2;

    test_one<cube2_solver, capacity, solved_check<cube2_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<p0s_solver, capacity, partial_check<p0s_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<p1s_solver, capacity, solved_check<p1s_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<_2ps_solver, capacity, solved_check<_2ps_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<p0es_solver, capacity, partial_check<p0es_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<c8_solver, capacity, partial_check<c8_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    test_one<opt_solver, capacity, solved_check<opt_solver>>(
            n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

    if (full) {
        test_one<p0_solver, capacity, partial_check<p0_solver>>(
                n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

        test_one<p1_solver, capacity, solved_check<p1_solver>>(
                n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

        test_one<_2p_solver, capacity, solved_check<_2p_solver>>(
                n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);

        test_one<e12s_solver, capacity, partial_check<e12s_solver>>(
                n_thread, seed, n_cube, rand_n_moves, max_n_moves, max_n_solution);
    }

    test_one<cube2_solver, capacity, solved_check<cube2_solver>>(
            n_thread, seed, n_cube, rand_n_moves, 7, 100);

    test_one<_2ps_solver, capacity, solved_check<_2ps_solver>>(
            n_thread, seed, n_cube, rand_n_moves, 7, 100);
}

std::string usage = R"(Usage

Test:
    %s  --n_thread 4  --seed 0

Test Full
    %s  --n_thread 4  --seed 0  --full
)";

void parse_arg(int argc, char **argv, u64 &n_thread, u64 &seed, bool &full) {
    std::tuple<u64, u64> n_thread_t = {1, 256};

    cxxopts::Options option(argv[0], "Test");
    option.add_options()
            ("n_thread", "1~256", cxxopts::value<u64>(n_thread)->default_value("4"))
            ("seed", "0~max", cxxopts::value<u64>(seed)->default_value("0"))
            ("full", "bool", cxxopts::value<bool>(full)->default_value("false"))
            ("usage", "show usage")
            ("help", "show help");

    try {
        cxxopts::ParseResult result = option.parse(argc, argv);
        if (result.count("usage")) {
            u64 n = usage.size() + strlen(argv[0]) * 2;
            std::vector<char> buf(n, '\0');
            snprintf(&buf[0], n, usage.c_str(), argv[0], argv[0]);
            std::cout << &buf[0] << std::endl;
            exit(0);
        } else if (result.count("help")) {
            std::cout << option.help() << std::endl;
            exit(0);
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "error option: " << e.what() << std::endl;
        exit(1);
    }

    if (n_thread < std::get<0>(n_thread_t) or n_thread > std::get<1>(n_thread_t)) {
        std::cout << "error n_thread: " << n_thread << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv) {
    u64 n_thread;
    u64 seed;
    bool full;

    parse_arg(argc, argv, n_thread, seed, full);

    test(n_thread, seed, full);

    return 0;
}
