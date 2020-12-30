#pragma once

#include <iterator>
#include <vector>
#include <algorithm>

// Pulled from stackoverflow https://stackoverflow.com/questions/1784573/iterator-for-2d-vector

// An iterator over a vector of vectors.
template<typename T>
class vv_iterator : public std::iterator<std::bidirectional_iterator_tag, T>
    {
    public:

        static vv_iterator<T> begin(std::vector<std::vector<T>>& vv) 
            {
            // Return first valid element
            std::size_t firstValid = 0;
            while( firstValid < vv.size() && vv[firstValid].empty() ) ++firstValid;
            return vv_iterator(&vv, firstValid, 0);
            }

        static vv_iterator<T> end(std::vector<std::vector<T>>& vv) 
            {
            return vv_iterator(&vv, vv.size(), 0 );
            }

        vv_iterator() = default;

        vv_iterator & operator++()
            {
            // If we haven't reached the end of this sub-vector.
            if( idxInner + 1 < (*vv)[idxOuter].size() ) ++idxInner;
            else
                {
                // Otherwise skip to the next sub-vector, and keep skipping over empty
                // ones until we reach a non-empty one or the end.
                do
                    {
                    ++idxOuter;
                    } while (idxOuter < (*vv).size() && (*vv)[idxOuter].empty());

                // Go to the start of this vector.
                idxInner = 0;
                }
            return *this;
            }

        vv_iterator& operator--()
            {
            // If we haven't reached the start of this sub-vector.
            if (idxInner > 0) --idxInner;
            else
                {
                // Otherwise skip to the previous sub-vector, and keep skipping over empty
                // ones until we reach a non-empty one.
                do
                    {
                    --idxOuter;
                    } while ((*vv)[idxOuter].empty());

                // Go to the end of this vector.
                idxInner = (*vv)[idxOuter].size() - 1;
                }
            return *this;
            }

        vv_iterator operator++(int)
            {
            T retval = *this;
            ++(*this);
            return retval;
            }

        vv_iterator operator--(int)
            {
            T retval = *this;
            --(*this);
            return retval;
            }

        bool operator==(const vv_iterator& other) const
            {
            return other.vv == vv && other.idxOuter == idxOuter && other.idxInner == idxInner;
            }
        bool operator!=(const vv_iterator &other) const
            {
            return !(*this == other);
            }

        const T& operator*() const { return *this; }
        T& operator*() { return (*vv)[idxOuter][idxInner]; }

        const T& operator->() const { return *this; }
        T& operator->() { return *this; }



        vv_iterator(std::vector<std::vector<T>>* _vv, std::size_t _idxOuter, std::size_t _idxInner)
            : vv(_vv)
            , idxOuter(_idxOuter)
            , idxInner(_idxInner) 
            {}

        std::vector<std::vector<T>>* vv = nullptr;
        std::size_t idxOuter = 0;
        std::size_t idxInner = 0;
    };