#pragma once

#include <vector>
#include <memory>
#include <functional>

#include "flan/defines.h"

/** This file contains generic dsp utility functions. 
 */

namespace flan {

struct vec2;

size_t autocorrelationFFTSize( size_t windowSize );
//std::vector<float> autocorrelation( const float * signal, Frame n, std::shared_ptr<FFTHelper> = nullptr );

// Returns the (x,y) coordinates of the peak
std::pair<float,float> parabolicInterpolation( float y0, float y1, float y2, int x1 );
std::pair<float,float> parabolicInterpolation( const std::vector<float> & d, int tau );
std::pair<float,float> parabolicInterpolation( std::function< float ( int ) > f, int tau );

std::vector<vec2> findPeaks( std::function< float ( int ) > data, int size, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> findPeaks( const std::vector<float> & data, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> findValleys( std::function< float ( int ) > data, int size, int maxValleys = -1, bool ampOrder = false, bool interpolate = true );
std::vector<vec2> findValleys( const std::vector<float> & data, int maxPeaks = -1, bool ampOrder = false, bool interpolate = true );

float mean( const std::vector<float> & data );
float mean( std::function< float ( int ) > data, int n );
vec2 meanAndStandardDeviation( const std::vector<float> & data );
vec2 meanAndStandardDeviation( std::function< float ( int ) > data, int n );

}