#include "cxxopts.hpp"

#include "base.h"
#include "group.h"
#include "search.h"
#include "cube3.h"
#include "cube3_2p.h"
#include "cube3_opt.h"
#include "cuda_cube_adaptor.h"
#include "interface.h"

using namespace cube;
using namespace cube::_3;
using namespace cube::_3::_2p;
using namespace cube::_3::opt;

constexpr u64 _2p_capacity = 29;
constexpr u64 opt_capacity = 20;

std::unique_ptr<_2ps_solver> _2p_s = nullptr;
std::unique_ptr<opt_solver> opt_s = nullptr;
std::unique_ptr<cuda_cube::opt_solver_manager> opt_s_m = nullptr;

struct solver_d {
    virtual ~solver_d() = default;

    virtual std::vector<std::vector<u8>> solve(const cube3 &a) = 0;
};

struct _2p_solver_d : solver_d {
    u64 _2p_n_moves;

    _2p_solver_d(u64 n_thread, u64 __2p_n_moves) : _2p_n_moves(__2p_n_moves) {
        if (_2p_s == nullptr) {
            _2p_s = std::make_unique<_2ps_solver>(n_thread);
        }
    }

    ~_2p_solver_d() override = default;

    std::vector<std::vector<u8>> solve(const cube3 &a) override {
        auto it = _2p_s->solve<_2p_capacity>(a);
        while (true) {
            auto[f, moves] = it();
            if (f & flag::solution) {
                std::cout << moves_to_string<_2ps_solver, _2p_capacity>(moves) << std::endl;
                if (f & flag::optimum or moves.n <= _2p_n_moves) {
                    return {std::vector<u8>(moves.a.begin(), moves.a.begin() + moves.n)};
                }
            } else if (f & flag::end) {
                return {};
            }
        }
    }
};

struct opt_solver_d : solver_d {
    u64 sym_mask_n_moves;
    u64 n_solution;

    opt_solver_d(u64 n_thread, u64 _sym_mask_n_moves, u64 _n_solution) :
            sym_mask_n_moves(_sym_mask_n_moves), n_solution(_n_solution) {
        if (opt_s == nullptr) {
            opt_s = std::make_unique<opt_solver>(n_thread);
        }
    }

    ~opt_solver_d() override = default;

    std::vector<std::vector<u8>> solve(const cube3 &a) override {
        std::vector<std::vector<u8>> solutions{};
        auto it = opt_s->solve<opt_capacity>(a, opt_capacity, sym_mask_n_moves);
        while (solutions.size() < n_solution) {
            auto[f, moves] = it();
            if (f & flag::solution) {
                std::cout << moves_to_string<opt_solver, opt_capacity>(moves) << std::endl;
                solutions.emplace_back(moves.a.begin(), moves.a.begin() + moves.n);
            } else if (f & flag::end) {
                break;
            }
        }
        return solutions;
    }
};

struct thread_opt_solver_d : solver_d {
    typedef thread_dfs<opt_solver, opt_capacity> parallel_dfs;

    const typename parallel_dfs::solver *pd_opt_s;
    std::string schedule;
    u64 n_pd_thread;
    u64 bfs_count;

    thread_opt_solver_d(const std::string &_schedule, u64 n_thread, u64 _bfs_count) :
            pd_opt_s(nullptr), schedule(_schedule), n_pd_thread(n_thread), bfs_count(_bfs_count) {
        if (opt_s == nullptr) {
            opt_s = std::make_unique<opt_solver>(n_thread);
        }
        pd_opt_s = opt_s.get();
    }

    ~thread_opt_solver_d() override = default;

    std::vector<std::vector<u8>> solve(const cube3 &a) override {
        u64 f;
        t_moves<opt_capacity> moves;
        if (schedule == "simple") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, simple_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else if (schedule == "linear") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, linear_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else if (schedule == "best") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, best_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else {
            assert(0);
        }
        if (f & flag::solution) {
            std::cout << moves_to_string<opt_solver, opt_capacity>(moves) << std::endl;
            return {std::vector<u8>(moves.a.begin(), moves.a.begin() + moves.n)};
        } else {
            return {};
        }
    }
};

struct cuda_opt_solver_d : solver_d {
    typedef cuda_cube::cuda_dfs<opt_solver, cuda_cube::opt_solver, opt_capacity> parallel_dfs;

    const typename parallel_dfs::solver *pd_opt_s;
    std::string schedule;
    u64 n_pd_thread;
    u64 bfs_count;

    cuda_opt_solver_d(const std::string &_schedule, u64 n_thread, u64 n_cuda_thread, u64 _bfs_count) :
            pd_opt_s(nullptr), schedule(_schedule), n_pd_thread(n_cuda_thread), bfs_count(_bfs_count) {
        if (opt_s == nullptr) {
            opt_s = std::make_unique<opt_solver>(n_thread);
        }
        if (opt_s_m == nullptr) {
            opt_s_m = std::make_unique<cuda_cube::opt_solver_manager>(*opt_s);
        }
        pd_opt_s = &opt_s_m->get();
    }

    ~cuda_opt_solver_d() override = default;

    std::vector<std::vector<u8>> solve(const cube3 &a) override {
        u64 f;
        t_moves<opt_capacity> moves;
        if (schedule == "simple") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, simple_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else if (schedule == "linear") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, linear_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else if (schedule == "best") {
            std::tie(f, moves) = parallel_ida_star<opt_solver, opt_capacity, parallel_dfs, best_schedule>::run(
                    *opt_s, *pd_opt_s, a, n_pd_thread, opt_capacity, bfs_count);
        } else {
            assert(0);
        }
        if (f & flag::solution) {
            std::cout << moves_to_string<opt_solver, opt_capacity>(moves) << std::endl;
            return {std::vector<u8>(moves.a.begin(), moves.a.begin() + moves.n)};
        } else {
            return {};
        }
    }
};

std::tuple<bool, std::string> read_cube(std::istream &in) {
    std::string s = "";
    while (true) {
        int c = in.get();
        if (c == std::char_traits<char>::eof()) {
            return {true, ""};
        } else if (c == ';') {
            return {false, s};
        } else {
            s += char(c);
        }
    }
}

std::string moves_to_string(const std::vector<u8> &moves) {
    std::string s = "(";
    u64 n = moves.size();
    for (u64 i = 0; i < n; i++) {
        s += cube3_base_name[moves[i]];
        if (i < n - 1) {
            s += " ";
        }
    }
    s += ")";
    return s;
}

std::string usage = R"(Usage

Two Phase Solver:
    %s  --algorithm two_phase  --n_thread 4  --2p_n_moves 24  --input example.txt  --output result.txt

Optimum Solver:
    %s  --algorithm optimum  --n_thread 4  --sym_n_moves 6  --n_solution=1  --input example.txt  --output result.txt

Thread Optimum Solver:
    %s  --algorithm thread_optimum  --schedule simple  --n_thread 4  --bfs_count=100000  --input example.txt  --output result.txt

CUDA Optimum Solver:
    %s  --algorithm cuda_optimum  --schedule simple  --n_thread 4  --n_cuda_thread 4096  --bfs_count=100000  --input example.txt  --output result.txt
)";

void parse_arg(
        int argc, char **argv,
        std::string &algorithm, std::string &schedule, u64 &n_thread, u64 &n_cuda_thread,
        u64 &_2p_n_moves, u64 &sym_n_moves, u64 &n_solution, u64 &bfs_count,
        std::string &input, std::string &output) {
    std::set<std::string> algorithm_set = {"two_phase", "optimum", "thread_optimum", "cuda_optimum"};
    std::set<std::string> schedule_set = {"simple", "linear", "best"};
    std::tuple<u64, u64> n_thread_t = {1, 256};
    std::tuple<u64, u64> n_cuda_thread_t = {1, 65536};
    std::tuple<u64, u64> _2p_n_moves_t = {0, 29};
    std::tuple<u64, u64> sym_n_moves_t = {0, 20};
    std::tuple<u64, u64> n_solution_t = {1, u64(-1)};
    std::tuple<u64, u64> bfs_count_t = {1, u64(-1)};

    cxxopts::Options option(argv[0], "Parallel Rubik's Cube Solver");
    option.add_options()
            ("algorithm", "[two_phase optimum thread_optimum cuda_optimum]",
             cxxopts::value<std::string>(algorithm))
            ("schedule", "[simple linear best]", cxxopts::value<std::string>(schedule)->default_value("simple"))
            ("n_thread", "1~256", cxxopts::value<u64>(n_thread)->default_value("4"))
            ("n_cuda_thread", "1~65536", cxxopts::value<u64>(n_cuda_thread)->default_value("4096"))
            ("2p_n_moves", "0~29", cxxopts::value<u64>(_2p_n_moves)->default_value("24"))
            ("sym_n_moves", "0~20", cxxopts::value<u64>(sym_n_moves)->default_value("6"))
            ("n_solution", "1~max", cxxopts::value<u64>(n_solution)->default_value("1"))
            ("bfs_count", "1~max", cxxopts::value<u64>(bfs_count)->default_value("100000"))
            ("input", "input file name", cxxopts::value<std::string>(input)->default_value(""))
            ("output", "output file name", cxxopts::value<std::string>(output)->default_value(""))
            ("usage", "show usage")
            ("help", "show help");

    try {
        cxxopts::ParseResult result = option.parse(argc, argv);
        if (result.count("usage")) {
            u64 n = usage.size() + strlen(argv[0]) * 4;
            std::vector<char> buf(n, '\0');
            snprintf(&buf[0], n, usage.c_str(), argv[0], argv[0], argv[0], argv[0]);
            std::cout << &buf[0] << std::endl;
            exit(0);
        } else if (result.count("help")) {
            std::cout << option.help() << std::endl;
            exit(0);
        } else if (result.count("algorithm") == 0) {
            std::cout << option.help() << std::endl;
            exit(0);
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "error option: " << e.what() << std::endl;
        exit(1);
    }

    if (algorithm_set.find(algorithm) == algorithm_set.end()) {
        std::cout << "error algorithm: " << algorithm << std::endl;
        exit(1);
    }

    if (schedule_set.find(schedule) == schedule_set.end()) {
        std::cout << "error schedule: " << schedule << std::endl;
        exit(1);
    }

    if (n_thread < std::get<0>(n_thread_t) or n_thread > std::get<1>(n_thread_t)) {
        std::cout << "error n_thread: " << n_thread << std::endl;
        exit(1);
    }

    if (n_cuda_thread < std::get<0>(n_cuda_thread_t) or n_cuda_thread > std::get<1>(n_cuda_thread_t)) {
        std::cout << "error n_cuda_thread: " << n_cuda_thread << std::endl;
        exit(1);
    }

    if (_2p_n_moves < std::get<0>(_2p_n_moves_t) or _2p_n_moves > std::get<1>(_2p_n_moves_t)) {
        std::cout << "error 2p_n_moves: " << _2p_n_moves << std::endl;
        exit(1);
    }

    if (sym_n_moves < std::get<0>(sym_n_moves_t) or sym_n_moves > std::get<1>(sym_n_moves_t)) {
        std::cout << "error sym_n_moves: " << sym_n_moves << std::endl;
        exit(1);
    }

    if (n_solution < std::get<0>(n_solution_t) or n_solution > std::get<1>(n_solution_t)) {
        std::cout << "error n_solution: " << n_solution << std::endl;
        exit(1);
    }

    if (bfs_count < std::get<0>(bfs_count_t) or bfs_count > std::get<1>(bfs_count_t)) {
        std::cout << "error bfs_count: " << bfs_count << std::endl;
        exit(1);
    }

    if (not std::filesystem::exists(input)) {
        std::cout << "file not exists: " << input << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv) {
    std::string algorithm;
    std::string schedule;
    u64 n_thread;
    u64 n_cuda_thread;
    u64 _2p_n_moves;
    u64 sym_n_moves;
    u64 n_solution;
    u64 bfs_count;
    std::string input;
    std::string output;

    parse_arg(argc, argv,
              algorithm, schedule, n_thread, n_cuda_thread,
              _2p_n_moves, sym_n_moves, n_solution, bfs_count,
              input, output);

    std::unique_ptr<std::ifstream> in_f = input.empty() ? nullptr :
                                          std::make_unique<std::ifstream>(input, std::ios::binary);
    std::unique_ptr<std::ofstream> out_f = output.empty() ? nullptr :
                                           std::make_unique<std::ofstream>(output, std::ios::binary);
    std::istream &in = input.empty() ? std::cin : *in_f;
    std::ostream &out = output.empty() ? std::cout : *out_f;

    std::unique_ptr<solver_d> sd = nullptr;
    if (algorithm == "two_phase") {
        sd = std::make_unique<_2p_solver_d>(n_thread, _2p_n_moves);
    } else if (algorithm == "optimum") {
        sd = std::make_unique<opt_solver_d>(n_thread, sym_n_moves, n_solution);
    } else if (algorithm == "thread_optimum") {
        sd = std::make_unique<thread_opt_solver_d>(schedule, n_thread, bfs_count);
    } else if (algorithm == "cuda_optimum") {
        sd = std::make_unique<cuda_opt_solver_d>(schedule, n_thread, n_cuda_thread, bfs_count);
    } else {
        assert(0);
    }

    cube3_interface cube3_i{};
    while (true) {
        if (input.empty()) {
            std::cout << "input:" << std::endl;
        }
        auto[eof, s] = read_cube(in);
        if (eof) {
            break;
        }
        auto[ok, error, scheme, a] = cube3_i.parse(s);
        if (ok) {
            std::vector<std::vector<u8>> solutions = sd->solve(a);
            out << "input:" << std::endl;
            out << cube3_interface::cube3_to_string(a, scheme) << std::endl;
            out << "solution: " << std::endl;
            for (const std::vector<u8> &moves: solutions) {
                out << moves_to_string(moves) << " " << moves.size() << "f" << std::endl;
            }
            out << std::endl;
        } else {
            out << "input:" << std::endl;
            out << s << std::endl;
            out << "error: " << error << std::endl;
            out << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
