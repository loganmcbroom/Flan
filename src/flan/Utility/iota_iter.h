#pragma once 

#include <iterator>

// This is a tool for allowing std algorithms to operate in parallel on an integer range
// For some reason you can't do that with iota_view in the ranges library
struct iota_iter
{
	using iterator_category = std::random_access_iterator_tag;
	using difference_type 	= int;
	using value_type 		= int;
	using pointer 			= int*;
	using reference 		= int&;

	iota_iter() : index( 0 ) {}
	iota_iter( int i ) : index( i ) {}
	operator int() const { return index; }
	int operator*() const { return index; }
	int operator[]( int i ) const { return i; }
	iota_iter operator++() { index++; return *this; }

	int index;
};