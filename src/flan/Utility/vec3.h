#pragma once

#include <array>

#include "flan/defines.h"

namespace flan {

/** Simple three element math vector. Arithmetic overloads are defined component-wise.
 */
struct vec3
	{
	float _x,_y,_z;

	/** Default constructor, 0 initializes elements. */
	vec3() : _x( 0 ), _y( 0 ), _z( 0 ) {}

	/** Copy constructor */
	vec3( const vec3 & _v ) : _x( _v.x() ), _y( _v.y() ), _z( _v.z() ) {}

	/** Secondary copy constructor for easier syntax */
	vec3( float x, float y, float z ) : _x( x ), _y( y ), _z( z ) {}

	/** Construct from std::array */
	vec3( std::array<float,3> c ) : _x( c[0] ), _y( c[1] ), _z( c[2] ) {}
	/** Convert to std::array */
	operator std::array<float,3>() const { return std::array<float,3>{ x(), y(), z() }; }
 
	float & x() { return _x; }
	float & y() { return _y; }
	float & z() { return _z; }
	const float & x() const { return _x; }
	const float & y() const { return _y; }
	const float & z() const { return _z; }
	
	/** Scalar multiplication */
	vec3 operator*( const float s ) const  { return vec3{ x() * s, y() * s, z() * s }; }
	/** Component-wise multiplication */
	vec3 operator*( const vec3 & r ) const  { return vec3{ x() * r.x(), y() * r.y(), z() * r.z() }; }
	/** Component-wise division */
	vec3 operator/( const vec3 & r ) const  { return vec3{ x() / r.x(), y() / r.y(), z() / r.z() }; }
	/** Component-wise modulus */
	vec3 operator%( const vec3 & r ) const  { return vec3{ fmod( x(), r.x() ), fmod( y(), r.y() ), fmod( z(), r.z() ) }; }
	/** Component-wise addition */
	vec3 operator+( const vec3 & r ) const { return vec3{ x() + r.x(), y() + r.y(), z() + r.z() }; }
	/** Component-wise subtraction */
	vec3 operator-( const vec3 & r ) const { return vec3{ x() - r.x(), y() - r.y(), z() - r.z() }; }
	/** Component-wise negation */
	vec3 operator-() const { return vec3{ -x(), -y(), -z() }; }

	float mag() const { return std::sqrt( x() * x() + y() * y() + z() * z() ); }
	};

};