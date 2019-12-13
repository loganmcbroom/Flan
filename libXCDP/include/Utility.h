#pragma once

//Any small misc code can live here 

#include <functional>
#include <vector>
#include <string>
#include <array>

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

std::array<uint8_t,3> HSVtoRGB( int H, double S, double V );
bool writeBMP( const std::string & filename, const std::vector<std::vector<std::array<uint8_t,3>>> & data );

}