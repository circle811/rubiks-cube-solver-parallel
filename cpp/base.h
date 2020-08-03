#ifndef _BASE_H
#define _BASE_H

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace cube {
    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    template<typename T, u64 _size_0, u64 _size_1>
    using array_2d = std::array<std::array<T, _size_1>, _size_0>;

    template<u64 _size>
    struct array_u8 : std::array<u8, _size> {
        constexpr bool operator==(const array_u8 &other) const {
            for (u64 i = 0; i < _size; i++) {
                if ((*this)[i] != other[i]) {
                    return false;
                }
            }
            return true;
        }
    };

    template<typename T, u64 _size>
    constexpr u64 find(const std::array<T, _size> &a, const T &x) {
        for (u64 i = 0; i < _size; i++) {
            if (a[i] == x) {
                return i;
            }
        }
        return u64(-1);
    }

    template<typename U, u64 _size>
    constexpr std::array<U, _size> generate_table(U (*f)(u64)) {
        std::array<U, _size> t{};
        for (U i = 0; i < _size; i++) {
            t[i] = f(i);
        }
        return t;
    }

    template<typename U, u64 _size_0, u64 _size_1>
    constexpr array_2d<U, _size_0, _size_1> generate_table(U (*f)(u64, u64)) {
        array_2d<U, _size_0, _size_1> t{};
        for (U i = 0; i < _size_0; i++) {
            for (U j = 0; j < _size_1; j++) {
                t[i][j] = f(i, j);
            }
        }
        return t;
    }

    template<typename T>
    std::unique_ptr<T> cache_data(const std::string &name, std::function<void(T &)> init) {
        std::string dir = "cache/";
        std::string path = dir + name;
        T *p = reinterpret_cast<T *>(operator new(sizeof(T)));
        if (std::filesystem::exists(path)) {
            std::cout << "cache_data: load " << name << " ..." << std::endl;
            std::ifstream f{path, std::ios::binary};
            f.read(reinterpret_cast<char *>(p), sizeof(T));
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
            f.write(reinterpret_cast<const char *>(p), sizeof(T));
            std::cout << "cache_data: save " << name << " ok" << std::endl;
        }
        return std::unique_ptr<T>(p);
    }
}

#endif
