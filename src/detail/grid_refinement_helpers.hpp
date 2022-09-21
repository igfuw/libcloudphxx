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

      protected:
      // ctor
      ref_common(const bool ref_flag):
        ref_flag(ref_flag)
      {}

      ref_common() = delete;
    };

    
    // for single values of type T
    template<class T>
    class ref_val : ref_common<T>
    {
      private:
      using parent_t = ref_common<T>;
      T val, val_ref;

      public:
      ref_val(T val, T val_ref):
        val(val),
        val_ref(val_ref),
        parent_t(true)
       {}
      ref_val() = delete;

      T& get()
      {
        return val;
      }
      /*
      const T& get()
      {
        return val;
      }
      */
      T& get_ref()
      {
        //return ref_flag ? val_ref : val;
        return val_ref;
      }
      /*
      const T& get_ref()
      {
        //return ref_flag ? val_ref : val;
        return val_ref;
      }
      */
    };

    template<class T>
    class ref_arr_common : ref_common<T>
    {
      private:
      using parent_t = ref_common<T>;
      using vec_t = thrust_device::vector<T>;

      protected:
      using parent_t::parent_t;
      vec_t arr, arr_ref; // actual data, arr_ref initialized only if refinement is actually done

      void resize(thrust_size_t size, thrust_size_t size_ref)
      {
        arr.resize(size);
        if(this->ref_flag) arr_ref.resize(size_ref);
      }
      void resize(ref_val<thrust_size_t> size)
      {
        arr.resize(size.get());
        if(this->ref_flag) arr_ref.resize(size.get_ref());
      }
      void reserve(thrust_size_t size, thrust_size_t size_ref)
      {
        arr.reserve(size);
        if(this->ref_flag) arr_ref.reserve(size_ref);
      }
      void reserve(ref_val<thrust_size_t> size)
      {
        arr.reserve(size.get());
        if(this->ref_flag) arr_ref.reserve(size.get_ref());
      }

      public:

      auto begin()
      {
        return arr.begin();
      }

      auto begin_ref()
      {
        return this->ref_flag ? arr_ref.begin() : arr.begin();
      }

      auto end()
      {
        return arr.end();
      }

      auto end_ref()
      {
        return this->ref_flag ? arr_ref.end() : arr.end();
      }

      auto ptr()
      {
        return &arr;
      }

      auto ptr_ref()
      {
        return &arr_ref;
      }

      vec_t& get()
      {
        return arr;
      }

      vec_t& get_ref()
      {
        return arr_ref;
      }
    };

    // for arrays of the size of the grid that need values
    // in the normal grid and in the refined grid (or in only one of these)
    template <class T>
    class ref_grid : public ref_arr_common<T>
    {
      using parent_t = ref_arr_common<T>;

      public:
      // ctor
      ref_grid(const int &n_cell, const int &n_cell_ref):
        parent_t(n_cell != n_cell_ref)
      {
        parent_t::resize(n_cell, n_cell_ref);
      }

      ref_grid(ref_val<thrust_size_t> n):
        ref_grid(n.get(), n.get_ref())
        {}

      ref_grid() = delete;

      // grid size cant be 0, ensure that correct arrays are pointed to
      auto begin()
      {
	assert(this->arr.size() > 0);
        return parent_t::begin();
      }
      auto begin_ref()
      {
	assert((this->ref_flag && this->arr_ref.size() > 0) || (!this->ref_flag && this->arr.size() > 0));
        return parent_t::begin_ref();
      }
      auto end()
      {
	assert(this->arr.size() > 0);
        return parent_t::end();
      }
      auto end_ref()
      {
	assert((this->ref_flag && this->arr_ref.size() > 0) || (!this->ref_flag && this->arr.size() > 0));
        return parent_t::end_ref();
      }
    };

    // for arrays of the size of the number of particles that need
    // two independent values: for the normal and for the refined grid (e.g. ijk)
    template <class T>
    class ref_part : public ref_arr_common<T>
    {
      using parent_t = ref_arr_common<T>;

      public:
      // ctor
      ref_part(const int &n_ref):
        parent_t(n_ref > 1)
      {assert(n_ref > 0);}

      ref_part() = delete;

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
