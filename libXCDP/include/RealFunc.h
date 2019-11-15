#pragma once

#include <functional>

namespace xcdp {

/*
 * Utility function abstraction for allowing doubles to be implicitly callable
 */
struct RealFunc 
{
	template< typename Callable >
	RealFunc( Callable f_ )
		: f( f_ )
		{}
	RealFunc( double t0 ) 
		: f( [t0]( double t ){ return t0; } ) 
		{}
	double operator()( double t ) const { return f(t); }

	void graph( const std::string & filename = std::string("tempRealFuncGraph.tga"),
		double left = -1, double right = 10, double bottom = -1, double top = 3.0, size_t resolution = 128 ) const;

private:
	std::function< double ( double ) > f;
};

RealFunc ADSR( double a, double d, double s, double r, double sLvl,
		double aExp = 1, double dExp = 1, double rExp = 1 );

} // End namespace xcdp