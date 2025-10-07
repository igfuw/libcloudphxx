#pragma once
#include <thrust/device_vector.h>
#include <vector>
#include <cassert>


namespace libcloudphxx
{
  namespace lgrngn
  {
    template<typename vec_t>
    class tmp_vector_pool {
        struct entry {
            vec_t vec;
            bool in_use = false;
            entry(size_t n) : vec(n) {}
        };
        std::vector<entry> pool;
    public:
        tmp_vector_pool(size_t pool_size = 1) {
            for (size_t i = 0; i < pool_size; ++i)
                pool.emplace_back();
        }

        // Add a new vector to the pool
        // void add_vector(size_t vec_size) {
        //     pool.emplace_back(vec_size);
        // }
        void add_vector() {
            pool.emplace_back();
        }

        void resize(size_t n) {
            for (size_t i = 0; i < pool.size(); ++i) {
                pool[i].vec.resize(n);
            }
        }

        void reserve(size_t n) {
            for (size_t i = 0; i < pool.size(); ++i) {
                pool[i].vec.reserve(n);
            }
        }

        // Acquire an available vector, returns its index
        size_t acquire() {
            for (size_t i = 0; i < pool.size(); ++i) {
                if (!pool[i].in_use) {
                    pool[i].in_use = true;
                    return i;
                }
            }
            assert(false && "No available temporary vectors in pool!");
            return size_t(-1);
        }

        // Release a vector by index
        void release(size_t idx) {
            assert(idx < pool.size() && pool[idx].in_use);
            pool[idx].in_use = false;
        }

        // Access vector by index
        vec_t& get(size_t idx) {
            assert(idx < pool.size() && pool[idx].in_use);
            return pool[idx].vec;
        }

        std::pair<size_t, vec_t&> get() {
            const size_t idx = acquire();
            assert(idx < pool.size() && pool[idx].in_use);
            return std::make_pair(idx, pool[idx].vec);
        }

        // RAII guard
        class guard {
            tmp_vector_pool<vec_t>& pool;
            size_t idx;
            bool valid;
        public:
            guard(tmp_vector_pool<vec_t>& pool_)
                : pool(pool_), idx(pool_.acquire()), valid(true) {}
            ~guard() { if (valid) pool.release(idx); }
            guard(const guard&) = delete;
            guard& operator=(const guard&) = delete;
            guard(guard&& other) noexcept : pool(other.pool), idx(other.idx), valid(other.valid) {
                other.valid = false;
            }
            guard& operator=(guard&& other) noexcept {
                if (this != &other) {
                    if (valid) pool.release(idx);
                    pool = other.pool;
                    idx = other.idx;
                    valid = other.valid;
                    other.valid = false;
                }
                return *this;
            }
            // void release() {
            //   pool = nullptr;
            //   idx = 4444444;
            //   valid = false;
            // }
            vec_t& get() { return pool.get(idx); }
            // vec_t* operator->() { return &pool.get(idx); }
            vec_t& operator*() { return pool.get(idx); }
        };

        guard get_guard() {
            return guard(*this);
        }
    };
  };
};


// Usage example 1:
{
    tmp_vector_pool<float>::guard tmp(pool);
    thrust::device_vector<float>& tmp = guard.get();
    // Use *tmp or tmp.get() as a thrust::device_vector<float>&
} // Automatically released here

// Usage example 2:
size_t idx = pool.acquire();
thrust::device_vector<float>& tmp = pool.get(idx);
// Use tmp...
pool.release(idx); // Mark as available again
