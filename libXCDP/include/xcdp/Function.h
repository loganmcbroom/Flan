#pragma once

#include <functional>
#include <complex>

#include "Utility.h"

namespace xcdp {

/** This is used as a base for all the vector (math vector) classes used throughout xcdp.
 *	These are primarily created by passing either lambda functions or constants to Audio or PVOC methods.
 */
template< typename I, typename O >
struct Function
{
	/** Callable constructor. This allows anything that can be called with the correct input and output type to be cast to xcdp::Function. */
	template< typename Callable >
	Function( Callable f_ ) : f( f_ ) {}
	Function( float t0	  ) : f( [t0]( I ){ return O( t0 ); } ) {}
	Function( double t0	  ) : f( [t0]( I ){ return O( t0 ); } ) {}
	Function( int t0	  ) : f( [t0]( I ){ return O( t0 ); } ) {}
	Function(			  ) : f( [  ]( I ){ return O( 0  ); } ) {}

	/** Function application */
	O operator()( I t ) const { return f(t); }

	/** Function composition */
	template< typename P >
	Function<P,O> operator()( const Function<P,I> g ) const { return [*this, g]( T t ){ return operator()( g(t) ); }; }

	/** Function multiplication */
    Function<I,O> operator*( const Function<I,O> r ) const { return [*this, r](I t){ return operator()(t) * r(t); }; }

	/** Function division */
	Function<I,O> operator/( const Function<I,O> r ) const { return [*this, r](I t){ return operator()(t) / r(t); }; }

	/** Function addition */
	Function<I,O> operator+( const Function<I,O> r ) const { return [*this, r](I t){ return operator()(t) + r(t); }; }

	/** Function subtraction */
	Function<I,O> operator-( const Function<I,O> r ) const { return [*this, r](I t){ return operator()(t) - r(t); }; }

	/** Function negation */
	Function<I,O> operator-() const { return [*this](I t){ return -operator()(t); }; }

	std::function< O ( I ) > f;
};

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

	std::array<float,2> v;
};

/** Real function of one real variable type. */
struct Func1x1 : public Function<float, float>
{
	template< typename Constructable >
	Func1x1( Constructable f ) : Function( f ) {}

	/** Exponentiate and evaluate 
	 *	\param t
	 */
	float exp( float t ) const { return std::pow( 2.0f, f(t) ); }

	/** Exponentiate as a function */
    Func1x1 exp() const { return [*this]( float t ){ return exp( t ); }; }

	/** Create and save a bmp of the function. 
	 * \param filename File path at which to save the created bmp.
	 * \param left Left side of graph window.
	 * \param right Right side of graph window.
	 * \param bottom Bottom side of graph window.
	 * \param top Top side of graph window.
	 * \param resolution Pixels per graph unit.
	 */
	void graph( const std::string & filename,
		float left = -5, float right = 5, float bottom = -5, float top = 5, size_t resolution = 128 ) const;

	/** Generate an ADSR envelope. Generated envelopes always range from 0 to 1.
	 *	The Exp parameters dictate how the envelope curves between fixed points. Using an Exp of 1 gives a linear movement.
	 *	Values between 0 and 1 give curves that move rapidly toward their destination before slowing down.
	 *	Values larger than 1 do the opposite. Other values will give unusual, but well defined, behaviour, so no clamping is applied.
	 *
	 * \param a Attack length in time.
	 * \param d Decay length in time.
	 * \param s Sustain length in time.
	 * \param r Release length in time.
	 * \param sLvl Sustain level.
	 * \param aExp Exponent of the attack curve.
	 * \param dExp Exponent of the decay curve.
	 * \param rExp Exponent of the release curve. 
	 */
	static Func1x1 ADSR( float a, float d, float s, float r, float sLvl,
		float aExp = 1, float dExp = 1, float rExp = 1 );

	/** This generates a periodic Function from the input function. 
	 *
	 *	\param wave The function to periodize. This will only be accessed on the domain [0,period).
	 *	\param min Bottom of the generated function.
	 *	\param max Top of the generated function.
	 *	\param period The period of the output.
	 */
	Func1x1 periodize( Func1x1 period = 1.0f );

	/** Generate a function that passes through a given set of points.
	 *
	 * \param points Points that the generated function must pass through.
	 * \param interp How the generated function should move between points.
	 */
	static Func1x1 interpolatePoints( const std::vector< std::pair< float, float > > points, 
		Interpolator interp = Interpolators::linear );
};

/** Real function of vec2 type. */
struct Func2x1 : public Function<vec2, float>
{
	template< typename Constructable >
	Func2x1( Constructable f ) : Function( f ) {}

	//Func2x1( std::function< float ( float, float ) > f ) 
	//	: Function( [f]( vec2 v ) { return f( v.x(), v.y() ); } ) {}

	/** Func2x1 can be constructed from Func1x1 by ignoring the second argument to the constructed Func2x1.  */
	Func2x1( Func1x1 f_ ) : Function( [f_]( vec2 v ){ return f_( v.x() ); } ) {}

	/** Helper for calling without converting parameters to vec2 */
	float operator()( float t, float f ) const { return Function::operator()( vec2{ t, f } ); }
};


/** vec2 function of vec2 type. */
struct Func2x2 : public Function<vec2, vec2>
{
	template< typename Constructable >
	Func2x2( Constructable f ) : Function( f ) {}

	/** Helper for calling without converting parameters to vec2 */
	vec2 operator()( float x, float y ) { return Function::operator()( { x, y } ); }
};

} // End namespace xcdp
