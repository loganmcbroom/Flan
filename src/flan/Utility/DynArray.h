#pragma once

namespace flan {

/*
The stl does not have a sequential data structure that can hold immovable types and decide its size at runtime, so we have to make one.
Counter for number of times I've tried to make this work: 4

As far as I can tell, you cannot design this by wrapping an array of const T. 
*/

// template<typename T>
// class dynarray {
// public:
//     using MutT = std::remove_const_t<T>;

//     dynarray() 
//         : vSize( 0 ) 
//         , v() 
//         {}

//     dynarray( size_t size ) 
//         : vSize( size )
//         , v( new MutT[size] )
//         {}

//     //dynarray( std::vector<MutT> && vec ) : v( vec ) {}
//     //dynarray( std::initializer_list<T> init ) : v() {}

//     const T & operator[]( size_t i ) const { return v[i]; }
//     size_t size() const { return vSize; }
//     bool empty() const { return size() == 0; }
//     const T * begin() const { return v; }
//     const T * end() const { return v + size(); }

// private:
//     const size_t vSize;
//     std::unique_ptr<MutT[]> v;
// };

}