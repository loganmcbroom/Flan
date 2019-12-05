#pragma once

#include <functional>

namespace xcdp {

typedef std::function< double ( double mix, double a, double b ) > Interpolator;
namespace Interpolators 
	{
	extern const Interpolator 
		constant, 
		linear, 
		sine;
	}

typedef std::function< double ( double, double ) > Perturber;
namespace Perturbers {
extern const Perturber
	identity, 
	normalDist, 
	normalDistUp, 
	normalDistDown;
}
}