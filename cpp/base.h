#ifndef _BASE_H
#define _BASE_H

#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace cube {
    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    template<typename T, u64 n0, u64 n1>
    using array_2d = std::array<std::array<T, n1>, n0>;

    template<typename T, u64 n>
    constexpr u64 array_find(const std::array<T, n> &a, const T &x) {
        for (u64 i = 0; i < n; i++) {
            if (a[i] == x) {
                return i;
            }
        }
        return u64(-1);
    }

    template<typename T, u64 n>
    constexpr bool array_eq(const std::array<T, n> &a, const std::array<T, n> &b) {
        for (u64 i = 0; i < n; i++) {
            if (a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }

    template<typename T, u64 n, u64 n_sub>
    constexpr std::array<T, n_sub> array_sub(const std::array<T, n> &a, const std::array<u64, n_sub> &index) {
        std::array<T, n_sub> b{};
        for (u64 i = 0; i < n_sub; i++) {
            b[i] = a[index[i]];
        }
        return b;
    }

    template<typename U, u64 n>
    constexpr U array_sum(const std::array<U, n> &a) {
        U s = 0;
        for (u64 i = 0; i < n; i++) {
            s += a[i];
        }
        return s;
    }

    template<typename T>
    std::string vector_to_string(const std::vector<T> &a) {
        std::string s = "(";
        u64 n = a.size();
        for (u64 i = 0; i < n; i++) {
            s += std::to_string(a[i]);
            if (i < n - 1) {
                s += " ";
            }
        }
        s += ")";
        return s;
    }

    template<typename U>
    constexpr U vector_sum(const std::vector<U> &a) {
        U s = 0;
        u64 n = a.size();
        for (u64 i = 0; i < n; i++) {
            s += a[i];
        }
        return s;
    }

    template<typename T>
    std::unique_ptr<T> cache_data(const std::string &name, const std::function<void(T &)> &init) {
        std::string dir = "cache/";
        std::string path = dir + name;
        std::unique_ptr<T> p = std::make_unique<T>();
        if (std::filesystem::exists(path)) {
            std::cout << "cache_data: load " << name << " ..." << std::endl;
            std::ifstream f{path, std::ios::binary};
            f.read(reinterpret_cast<char *>(p.get()), sizeof(T));
            std::cout << "cache_data: load " << name << " ok" << std::endl;
        } else {
            std::cout << "cache_data: compute " << name << " ..." << std::endl;
            auto t0 = std::chrono::steady_clock::now();
            init(*p);
            auto t1 = std::chrono::steady_clock::now();
            std::chrono::duration<double> d = t1 - t0;
            std::cout << "cache_data: compute " << name << " ok, time=" << d.count() << "s" << std::endl;
            std::cout << "cache_data: save " << name << " ..." << std::endl;
            std::filesystem::create_directory(dir);
            std::ofstream f{path, std::ios::binary};
            f.write(reinterpret_cast<const char *>(p.get()), sizeof(T));
            std::cout << "cache_data: save " << name << " ok" << std::endl;
        }
        return p;
    }
}

#endif
