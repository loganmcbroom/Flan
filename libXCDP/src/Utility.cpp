#include "Utility.h"

#include <random>
#include <ctime>

static double pi = acos(-1);

namespace xcdp {

const Interpolator Interpolators::constant = []( double mix, double a, double b ) 
	{ 
	return a;
	};

const Interpolator Interpolators::linear = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - mix ) * a + mix * b;
	};

const Interpolator Interpolators::sine = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - cos( pi * mix ) ) / 2.0 * ( b - a ) + a; 
	};


///Code duplication?
const Perturber Perturbers::normalDist = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x + dist(rng) * amount;
	};

const Perturber Perturbers::normalDistUp = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x + abs(dist(rng) * amount);
	};

const Perturber Perturbers::normalDistDown = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x - abs(dist(rng) * amount);
	};

const Perturber Perturbers::identity = []( double x, double amount ){ return x; };

}