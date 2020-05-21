#include "xcdp/PVOC.h"

#include "xcdp/Function.h"

namespace xcdp {

//========================================================================
// Blur
//========================================================================

//PVOC PVOC::blur_blur( RealFunc blur ) const
//	{ 
//	return desample( blur, Interpolators::linear ); 
//	}

PVOC PVOC::blur_chorus( Func1x1 fspread ) const
	{
	return perturb( 0, fspread, Distributions::identity,  Distributions::normal );
	}

//========================================================================
// Combine
//========================================================================

PVOC PVOC::combine_cross( const PVOC & ampSource, Func1x1 amount ) const
	{
	return replaceAmplitudes( ampSource, amount );
	}

}