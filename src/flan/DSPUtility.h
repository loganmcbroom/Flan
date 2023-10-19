#pragma once

#include <vector>
#include <memory>
#include <functional>

#include "flan/defines.h"

namespace flan {

struct vec2;

// Returns the (x,y) coordinates of the peak
std::pair<float,float> parabolic_interpolation( float y0, float y1, float y2, int x1 );
std::pair<float,float> parabolic_interpolation( const std::vector<float> & d, int tau );
std::pair<float,float> parabolic_interpolation( std::function< float ( int ) > f, int tau );

std::vector<vec2> find_peaks( std::function< float ( int ) > data, int size, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> find_peaks( const std::vector<float> & data, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> find_valleys( std::function< float ( int ) > data, int size, int maxValleys = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> find_valleys( const std::vector<float> & data, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );

float mean( const std::vector<float> & data );
float mean( std::function< float ( int ) > data, int n );
vec2 mean_and_sd( const std::vector<float> & data );
vec2 mean_and_sd( std::function< float ( int ) > data, int n );

}