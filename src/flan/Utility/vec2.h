#pragma once

#include <array>
#include <complex>

#include "flan/defines.h"

namespace flan {

/** Simple two element math vector. Arithmetic overloads are defined component-wise.
 */
template<typename T>
struct vec2Base
	{
	T _x,_y;

	/** Default constructor, 0 initializes elements. */
	vec2Base() : _x( 0 ), _y( 0 ) {}

	/** Copy constructor */
	vec2Base( const vec2Base & _v ) : _x( _v.x() ), _y( _v.y() ) {}

	/** Secondary copy constructor for easier syntax */
	vec2Base( T x, T y ) : _x( x ), _y( y ) {}

	/** Construct from different T */
	template<typename O> 
	requires std::convertible_to<T, O>
	vec2Base( const vec2Base<O> & o ) : _x( o.x() ), _y( o.y() ) {}
	/** Convert to other T */
	template<typename O> 
	requires std::convertible_to<T, O>
	operator vec2Base<O>() const { return vec2Base<O>( x(), y() ); }

	/** Construct from std::complex */
	vec2Base( std::complex<T> c ) : _x( c.real() ), _y( c.imag() ) {}
	/** Convert to std::complex */
	operator std::complex<T>() const { return std::complex<T>( x(), y() ); }

	/** Construct from std::array */
	vec2Base( std::array<T,2> c ) : _x( c[0] ), _y( c[1] ) {}
	/** Convert to std::array */
	operator std::array<T,2>() const { return std::array<T,2>{ x(), y() }; }

	/** Construct from flan TF */
	vec2Base( TF c ) : _x( c.t ), _y( c.f ) {}
	/** Convert to flan TF */
	operator TF() const { return TF{ x(), y() }; }

	/** Construct from flan MF */
	vec2Base( MF c ) : _x( c.m ), _y( c.f ) {}
	/** Convert to flan MF */
	operator MF() const { return MF{ x(), y() }; }

	/** Construct from any pair */
	template<typename A, typename B>
	vec2Base( std::pair<A, B> c ) : _x( c.first ), _y( c.second ) {}
	/** Convert to any pair */
	template<typename A, typename B>
	operator std::pair<A, B>() const { return std::make_pair<A, B>( x(), y() ); }
 
	T & x() { return _x; }
	T & y() { return _y; }
	const T & x() const { return _x; }
	const T & y() const { return _y; }
	/** Identical to x() */
	T & t() { return _x; }
	/** Identical to y() */
	T & f() { return _y; }
	/** Identical to x() */
	const T & t() const { return _x; }
	/** Identical to y() */
	const T & f() const { return _y; }
	
	/** Scalar multiplication */
	vec2Base operator*( const T s ) const  { return vec2Base{ x() * s, y() * s }; }
	/** Component-wise multiplication */
	vec2Base operator*( const vec2Base & r ) const  { return vec2Base{ x() * r.x(), y() * r.y() }; }
	/** Scalar division */
	vec2Base operator/( const T s ) const  { return vec2Base{ x() / s, y() / s }; }
	/** Component-wise division */
	vec2Base operator/( const vec2Base & r ) const  { return vec2Base{ x() / r.x(), y() / r.y() }; }
	/** Component-wise modulus */
	vec2Base operator%( const vec2Base & r ) const  { return vec2Base{ fmod( x(), r.x() ), fmod( y(), r.y() ) }; }
	/** Component-wise addition */
	vec2Base operator+( const vec2Base & r ) const { return vec2Base{ x() + r.x(), y() + r.y() }; }
	/** Component-wise subtraction */
	vec2Base operator-( const vec2Base & r ) const { return vec2Base{ x() - r.x(), y() - r.y() }; }
	/** Component-wise negation */
	vec2Base operator-() const { return vec2Base{ -x(), -y() }; }

	T mag() const { return std::sqrt( x() * x() + y() * y() ); }
	};

using vec2 = vec2Base<float>;
using vec2d = vec2Base<double>;

};