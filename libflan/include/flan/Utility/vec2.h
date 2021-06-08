#pragma once

#include <array>
#include <complex>

namespace flan {

/** Simple two element math vector. Arithmetic overloads are defined component-wise.
 */
struct vec2
	{
	/** Default constructor, 0 initializes elements. */
	vec2() : v({0,0}) {}

	/** Copy constructor */
	vec2( const vec2 & _v ) : v( _v.v ) {}

	/** Secondary copy constructor for easier syntax */
	vec2( float x, float y ) : v( { x, y } ) {}

	/** Construct from std::complex */
	vec2( std::complex<float> c ) : v({ c.real(), c.imag() }) {}
	/** Convert to std::complex */
	operator std::complex<float>() const { return std::complex<float>( x(), y() ); }

	/** Construct from std::array */
	vec2( std::array<float,2> c ) : v({ c[0], c[1] }) {}
	/** Convert to std::array */
	operator std::array<float,2>() const { return std::array<float,2>{ x(), y() }; }
 
	float & x() { return v[0]; }
	float & y() { return v[1]; }
	const float & x() const { return v[0]; }
	const float & y() const { return v[1]; }
	/** Identical to x() */
	float & t() { return v[0]; }
	/** Identical to y() */
	float & f() { return v[1]; }
	/** Identical to x() */
	const float & t() const { return v[0]; }
	/** Identical to y() */
	const float & f() const { return v[1]; }
	
	/** Scalar multiplication */
	vec2 operator*( const float s ) const  { return vec2{ x() * s, y() * s }; }
	/** Component-wise multiplication */
	vec2 operator*( const vec2 & r ) const  { return vec2{ x() * r.x(), y() * r.y() }; }
	/** Component-wise division */
	vec2 operator/( const vec2 & r ) const  { return vec2{ x() / r.x(), y() / r.y() }; }
	/** Component-wise addition */
	vec2 operator+( const vec2 & r ) const { return vec2{ x() + r.x(), y() + r.y() }; }
	/** Component-wise subtraction */
	vec2 operator-( const vec2 & r ) const { return vec2{ x() - r.x(), y() - r.y() }; }
	/** Component-wise negation */
	vec2 operator-() const { return vec2{ -x(), -y() }; }

	float mag() const { return std::sqrt( x() * x() + y() * y() ); }

	std::array<float,2> v;
	};

};