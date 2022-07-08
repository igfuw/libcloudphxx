// helper classes to store arrays associated with grid refinement

#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    template<class T>
    class ref_common
    {
      private:
      const bool ref_flag; // true if refinement is actually done
      thrust_device::vector<T> arr, arr_ref; // actual data, arr_ref initialized only if refinement is actually done

      protected:
      // ctor
      ref_common(const bool ref_flag):
        ref_flag(ref_flag)
      {}

      auto begin()
      {
        return arr.begin();
      }
      auto begin_ref()
      {
        return ref_flag ? arr_ref.begin() : arr.begin();
      }
      auto end()
      {
        return arr.end();
      }
      auto end_ref()
      {
        return ref_flag ? arr_ref.end() : arr.end();
      }
      void resize(thrust_size_t n, thrust_size_t n_ref)
      {
        arr.resize(n);
        if(ref_flag) arr_ref.resize(n_ref);
      }
      void reserve(thrust_size_t n, thrust_size_t n_ref)
      {
        arr.reserve(n);
        if(ref_flag) arr_ref.reserve(n_ref);
      }
    }

    // for arrays of the size of the grid
    template <class T>
    class grid_ref : ref_common<T>
    {
      using parent_t = ref_common<T>;

      public:
      // ctor
      grid_ref(const int &n_cell, const int &n_cell_ref):
        parent_t(n_cell != n_cell_ref)
      {
        parent_t::resize(n_cell, n_cell_ref);
      }
    };

    // for arrays of the size of the number of particles
    template <class T>
    class part_ref
    {
      using parent_t = ref_common<T>;

      public:
      // ctor
      part_ref(const int &n_ref):
        parent_t(n_ref > 1)
      {}

      void resize(int n)
      {
        parent_t::resize(n, n);
      }
      void reserve(int n)
      {
        parent_t::reserve(n, n);
      }
    };
  };
};
