#pragma once

/*                            =======================
================================ C/C++ HEADER FILE =================================
                              =======================                          *//**
\file  carray.h
\brief implements a general purpose array. Implements a subset of std::vector class
\n\n
Copyright (c) Endress+Hauser AG \n
All rights reserved.
*//*==================================================================================*/
/*
CHANGE HISTORY:
---------------
Date       Author Description
02.06.2020 BOL    v1.0 Released
*/

#ifndef CARRAY_
#define CARRAY_

#include <cstdlib>
#include <algorithm>
#include <cassert>
#ifdef _WIN
#include <iostream>
#endif


template <typename T, size_t N>
class CArray {
private:
    T m_Array[N];
    size_t m_MaxSize;
    size_t m_Size;

protected:

public:
    // type definitions
    typedef T  value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;

    class iterator
    {
    public:
        typedef iterator self_type;
        typedef int postproc;
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        //typedef std::forward_iterator_tag iterator_category;
        typedef ptrdiff_t difference_type;
        iterator(pointer ptr) : pPtr_(ptr) { }
        self_type& operator++() { pPtr_++; return *this; }
        self_type operator++(postproc) { const self_type iter(pPtr_); ++pPtr_; return iter; }//lint !e9215
        self_type& operator--() { pPtr_--; return *this; }
        self_type operator--(postproc) { const self_type iter(pPtr_); --pPtr_; return iter; }//lint !e9215
        self_type& operator +=(size_type offset) { pPtr_ = &pPtr_[offset]; return *this; }
        self_type operator +(size_type offset) const { self_type cpy(&pPtr_[offset]); return cpy; }
        self_type& operator -=(size_type offset) { pPtr_ = &pPtr_[0U-offset]; return *this; }//lint !e835
        self_type operator -(size_type offset) const { self_type cpy(&pPtr_[0U-offset]); return cpy; }//lint !e835
        reference operator*() { return *pPtr_; }
        pointer operator->() const { return pPtr_; } //lint !e1537 BOL: required for speedups (DSPs libs)
        bool operator==(const self_type& rhs) const { return pPtr_ == rhs.pPtr_; }
        bool operator!=(const self_type& rhs) const { return pPtr_ != rhs.pPtr_; }
        difference_type diff(self_type const& base) const { return static_cast<difference_type>(this->pPtr_ - base.pPtr_); }//lint !e947
    private:
        pointer pPtr_;
    };

    class const_iterator
    {
    public:
        typedef const_iterator self_type;
        typedef int postproc;
        typedef T value_type;
        typedef T const& reference;
        typedef T const* pointer;
        typedef ptrdiff_t difference_type;
        //typedef std::forward_iterator_tag iterator_category;
        const_iterator(pointer ptr) : pPtr_(ptr) { }
        const_iterator(iterator const &iter) : pPtr_(const_cast<pointer>(iter.operator->())) { }
        self_type& operator++() { pPtr_++; return *this; }
        self_type operator++(postproc) { const self_type iter(pPtr_); ++pPtr_; return iter; }//lint !e9215
        self_type& operator--() { pPtr_--; return *this; }
        self_type operator--(postproc) { const self_type iter(pPtr_); --pPtr_; return iter; }//lint !e9215
        self_type operator +(size_type offset) const { self_type cpy(&pPtr_[offset]); return cpy; }
        self_type& operator +=(size_type offset) { pPtr_ += offset; return *this; }
        self_type operator -(size_type offset) const { self_type cpy(&pPtr_[0U - offset]); return cpy; }//lint !e835
        self_type& operator -=(size_type offset) { pPtr_ -= offset; return *this; }
        reference operator*() const { return static_cast<reference>(*pPtr_); }
        pointer operator->() const { return pPtr_; } //lint !e1537 BOL: required for speedups (DSPs libs)
        bool operator==(const self_type& rhs) const { return pPtr_ == rhs.pPtr_; }
        bool operator!=(const self_type& rhs) const { return pPtr_ != rhs.pPtr_; }
        difference_type diff(self_type const& base) const { return static_cast<difference_type>(this->pPtr_ - base.pPtr_); }//lint !e947
    private:
        pointer pPtr_;
    };


    class reverse_iterator
    {
    public:
        typedef reverse_iterator self_type;
        typedef int postproc;
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        //typedef std::forward_iterator_tag iterator_category;
        typedef size_type difference_type;
        reverse_iterator(pointer ptr) : pPtr_(ptr) { }
        self_type& operator++() { pPtr_--; return *this; }
        self_type operator++(postproc) { const self_type iter(pPtr_); --pPtr_; return iter; }//lint !e9215
        self_type& operator--() { pPtr_++; return *this; }
        self_type operator--(postproc) { const self_type iter(pPtr_); ++pPtr_; return iter; }//lint !e9215
        self_type& operator +=(size_type offset) { pPtr_ -= offset; return *this; }
        self_type operator +(size_type offset) const { self_type cpy(&pPtr_[0U-offset]); return cpy; }//lint !e835
        self_type& operator -=(size_type offset) { pPtr_ += offset; return *this; }
        self_type operator -(size_type offset) const { self_type cpy(&pPtr_[offset]); return cpy; }
        reference operator*() { return *pPtr_; }
        pointer operator->() const { return pPtr_; } //lint !e1537 BOL: required for speedups (DSPs libs)
        bool operator==(const self_type& rhs) const { return pPtr_ == rhs.pPtr_; }
        bool operator!=(const self_type& rhs) const { return pPtr_ != rhs.pPtr_; }
        difference_type diff(self_type const& base) const { return this->pPtr_ - base.pPtr_; }
    private:
        pointer pPtr_;
    };

    class const_reverse_iterator
    {
    public:
        typedef const_reverse_iterator self_type;
        typedef int postproc;
        typedef T value_type;
        typedef T const& reference;
        typedef T const* pointer;
        typedef size_type difference_type;
        //typedef std::forward_iterator_tag iterator_category;
        const_reverse_iterator(pointer ptr) : pPtr_(ptr) { }
        const_reverse_iterator(reverse_iterator const& riter) : pPtr_(const_cast<pointer>(riter.operator->())) { }
        self_type& operator++() { pPtr_--; return *this; }
        self_type operator++(postproc) const { self_type iter(pPtr_); --pPtr_; return iter; }//lint !e9215
        self_type& operator--() { pPtr_++; return *this; }
        self_type operator--(postproc) const { self_type iter(pPtr_); ++pPtr_; return iter; }//lint !e9215
        self_type operator +(size_type offset) { const self_type cpy(&pPtr_[0U - offset]); return cpy; }//lint !e835
        self_type& operator +=(size_type offset) { pPtr_ -= offset; return *this; }
        self_type operator -(size_type offset) { const self_type cpy(&pPtr_[offset]); return cpy; }
        self_type& operator -=(size_type offset) { pPtr_ += offset; return *this; }
        reference operator*() const { return *pPtr_; }
        pointer operator->() const { return pPtr_; } //lint !e1537 BOL: required for speedups (DSPs libs)
        bool operator==(const self_type& rhs) const { return pPtr_ == rhs.pPtr_; }
        bool operator!=(const self_type& rhs) const { return pPtr_ != rhs.pPtr_; }
        difference_type diff(self_type const& base) const { return (this->pPtr_ - base.pPtr_) ; }
    private:
        pointer pPtr_;
    };

    // iterator support
    bool rangecheck(size_type i) {
        return i < m_Size;
    }
    bool rangecheck(iterator const& iter) {
        return (begin() <= iter) && (iter < end());
    }
    bool rangecheck(const_iterator const& iter) {
        return (cbegin() <= iter) && (iter < cend());
    }
    bool rangecheck(reverse_iterator const& iter) {
        return (rbegin() >= iter) && (iter > rend());
    }
    bool rangecheck(const_reverse_iterator const& iter) {
        return (crbegin() >= iter) && (iter > crend());
    }

    iterator begin() { return m_Array; }
    const_iterator cbegin() const { return const_cast<T*>(m_Array); }
    iterator end() { return &m_Array[m_Size]; }
    const_iterator cend() const { return &(const_cast<T*>(m_Array)[m_Size]); }
    reverse_iterator rbegin() { return &m_Array[m_Size-1u]; }
    const_reverse_iterator crbegin() const { return &(const_cast<T*>(m_Array)[m_Size-1u]); }
    reverse_iterator rend() { return &begin().operator->()[-1]; }//lint !e613
    const_reverse_iterator crend() const { return &cbegin().operator->()[-1]; }//lint !e613

    // operator[]
    reference operator[](size_type i) { return m_Array[i]; }
    const_reference operator[](size_type i) const { return m_Array[i]; }

    // at() with range check
    reference at(size_type i) { return (rangecheck(i) == true) ? m_Array[i] : (assert(false), *static_cast<T*>(nullptr)); }//lint !e9008 !e9113 !e505
    const_reference at(size_type i) const { return (rangecheck(i) == true) ? m_Array[i] : (assert(false), *static_cast<T*>(nullptr)); }//lint !e9008 !e9113 !e505

                                                                                                                                       // front() and back()
    reference front() { return (m_Size > 0) ? m_Array[0] : nullptr; }
    const_reference front() const { return (m_Size > 0) ? m_Array[0] : nullptr; }
    reference back() { return (m_Size > 0) ? m_Array[m_Size - 1] : nullptr; }
    const_reference back() const { return (m_Size > 0) ? m_Array[m_Size - 1] : nullptr; }

    // size is constant
    size_type size() const { return m_Size; }
    bool empty() { return (m_Size == 0); }
    static size_type max_size() { return static_size; }
    enum { static_size = N };

    // swap (note: linear complexity in N, constant for given instantiation)
    void swap(CArray<T, N>& y) {
        std::swap_ranges(begin(), end(), y.begin());
    }

    // direct access to data (read-only)
    T const* data() const { return m_Array; }

    // use array as C array (direct read/write access to data)
    T* data() { return m_Array; }

    // assignment with type conversion
    template <typename T2>
    CArray<T, N>& operator= (const CArray<T2, N>& rhs) {//lint !e1721
        std::copy(rhs.begin(), rhs.end(), begin());
        return *this;
    }

    // assign one value to all elements
    void assign(size_type count, const T& value)
    {
        count = (count < N) ? count : N;
        m_Size = count;
        for (T& item : *this)
        {
            item = value;
        }
    }

    bool push_back(T const& item)
    {
        bool retVal = false;
        if (m_Size < N)
        {
            m_Array[m_Size] = item;
            ++m_Size;
            retVal = true;
        }
        return retVal;
    }

    bool pop_back()
    {
        m_Size = (m_Size > static_cast <size_type>(0)) ? (m_Size - static_cast<size_type>(1)) : static_cast<size_type>(0);
        return m_Size > 0u;
    }

    void erase(CArray<T, N>::iterator index)
    {
        CArray<T, N>::iterator ix = index;
        CArray<T, N>::iterator stop = end();
        if (index != end())
        {
            while (++ix != stop)
            {
                *index = *ix;
                ++index;
            }
        }
        else
        {
            //that is an exception
        }
        --m_Size;
    }

    void erase(CArray<T, N>::iterator first, CArray<T, N>::iterator last)
    {
        size_type elementstodelete = last.diff(first);
        CArray<T, N>::iterator ix = last;
        CArray<T, N>::iterator stop = end();
        while (ix != stop)
        {
            *first = *ix;
            ++first; ++ix;
        }
        m_Size -= elementstodelete;
    }

    iterator insert(CArray<T, N>::iterator first, T const& a)
    {
        CArray<T, N>::iterator tmp = first;
        T swp = a;
        for (tmp = tmp; tmp != end(); ++tmp)
        {
            std::swap(swp, *tmp);
        }
        //last element drops out if size is insufficient
        if (m_Size < N)
        {
            *tmp = swp;
            ++m_Size;
        }
        return first;
    }

    iterator insert(CArray<T, N>::iterator first, size_type count, T const& a)
    {
        CArray<T, N>::iterator tmp_base;
        CArray<T, N>::iterator tmp_move = first + count;
        CArray<T, N>::iterator enditer = std::min((begin() + N), tmp_move + m_Size);
        for (tmp_base = first; tmp_move < enditer; ++tmp_base)
        {
            *tmp_move = tmp_base;
            ++tmp_move;
        }
        m_Size = tmp_base - begin();
        //last element drops out if size is insufficient
        enditer = std::min((begin() + N), first + count + m_Size);
        for (tmp_base = first; tmp_base < enditer; ++tmp_base)
        {
            *tmp_base = a;
        }
        return first;
    }

    template <typename T2, size_type N2>
    iterator insert(CArray<T, N>::iterator first, CArray<T2, N2>& a)
    {
        typename CArray<T2, N2>::iterator iter = a.begin();
        size_type offset = first.diff(begin());

        if ((offset + a.size()) < N)
        {
            size_type max_array_extension = ((a.size() + m_Size) < N) ? a.size() : ((N - a.size()) - 1);
            CArray<T, N>::iterator rstartiter = end() + max_array_extension;
            CArray<T, N>::iterator rstopiter = (end() + max_array_extension) - a.size();
            while (rstopiter != first)
            {
                rstartiter--;
                rstopiter--;
                *rstartiter = *rstopiter;
            }
        }
        else
        {
            //there is nothing to do as the array will be truncated
        }

        //last element drops out if size is insufficient
        //just copy the data
        CArray<T, N>::iterator enditer = (begin() + N);
        if ((offset + a.size()) < N)
        {
            enditer = first + a.size();
        }
        CArray<T, N>::iterator tmp_base = first;
        for (tmp_base = tmp_base; tmp_base != enditer; ++tmp_base)
        {
            *tmp_base = static_cast<T>(*iter);
            ++iter;
        }

        m_Size = ((m_Size + a.size()) < N) ? (m_Size + a.size()) : N;
        return first;
    }

    template <typename T2, size_type N2>
    iterator append(CArray<T2, N2>& a)
    {
        iterator tmp_base = end();
        typename CArray<T2, N2>::iterator iter = a.begin();

        iterator enditer = begin() + N;
        iterator appendix = end();
        if ((m_Size + a.size()) < N)
        {
            enditer = appendix + a.size();
        }

        for (tmp_base = tmp_base; tmp_base != enditer; ++tmp_base)
        {
            T2 d1 = *iter;
            T  d = static_cast<T>(d1);
            *tmp_base = d;
            ++iter;
        }
        m_Size = tmp_base.diff(begin());
        //last element drops out if size is insufficient
        return appendix;
    }

    iterator append(size_type len, T* a)
    {
        CArray<T, N>::iterator retVal = end();
        if (a != nullptr)
        {
            size_type finalsize = std::min(N, m_Size + len);
            for (size_type i = static_cast<size_type>(m_Size); i < finalsize; ++i)
            {
                m_Array[i] = *a;
                ++a;
            }
            m_Size = finalsize;
        }
        return retVal;
    }

    iterator create(size_type len, const T* a)
    {
        if (a != nullptr)
        {
            const size_type finalsize = std::min(N, len);
            for (size_type i = static_cast<size_type>(0); i < finalsize; ++i)
            {
                m_Array[i] = *a;
                ++a;
            }
            m_Size = finalsize;
        }
        return begin();
    }

    void clear()
    {
        resize(0u);
    }

    void resize(size_type count)
    {
        count = (count < N) ? count : N;
        m_Size = count;
    }

    void fill(T item)
    {
        for (CArray<T, N>::iterator i = begin(); i != end(); ++i)
        {
            *i = item;
        }
    }

    void pad(T item)
    {
        for (CArray<T, N>::iterator i = end(); i != (begin() + N); ++i)
        {
            *i = item;
        }
        m_Size = N;
    }

    CArray() :m_MaxSize(N), m_Size(0) { /*for (size_t i = 0; i < N;++i) { m_Array[i] = 0; }*/ }
    CArray(const T(&Array)[N]);
    ~CArray() {}
};

template <typename T, size_t N>
CArray<T, N>::CArray(const T(&Array)[N])
{
	const size_t maxsize = std::max(this->m_MaxSize, N);
    for (size_t i = 0u; i < maxsize; ++i)
    {
        m_Array[i] = Array[i];
    }

    m_Size = maxsize;
}
#endif
