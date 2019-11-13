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
	std::function< double ( double ) > f;
};

} // End namespace xcdp