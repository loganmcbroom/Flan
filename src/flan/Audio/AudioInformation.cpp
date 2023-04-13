#include "flan/Audio/Audio.h"

#include <iostream>
#include <set>
#include <limits>
#include <numeric>

#include "flan/WindowFunctions.h"
#include "flan/FFTHelper.h"
#include "flan/DSPUtility.h"

#include "flan/Graph.h"

using namespace flan;

//============================================================================================================================================================
// Helpers
//============================================================================================================================================================

// Computes sum of squared differences for set number of lag values
static std::vector<float> compute_d( const float * in, Frame n )
	{
    const size_t fftSize = n; //std::pow( 2, 1 + (int) std::ceil( std::log2( n ) ) );
	auto fft = std::make_shared<FFTHelper>( fftSize, true, true, false );

    // Calculate power terms
	std::vector<float> powerTerms( n / 2 );
    powerTerms[0] = 0.0;
    for( int j = 0; j < n / 2; ++j ) 
        powerTerms[0] += in[j] * in[j];
    for( int tau = 1; tau < powerTerms.size(); ++tau ) 
        powerTerms[tau] = powerTerms[tau-1] - in[tau-1] * in[tau-1] + in[tau-1 + n / 2] * in[tau-1 + n / 2];  
		
    // Get modified autocorrelation. Note this needs two ffts, one of a full window and one of a zero padded half window
    
    // Get half fft
    std::vector<std::complex<float>> halfFFT( fft->complexBufferSize() );
    std::copy( in, in + n/2, fft->realBegin() );
    std::fill( fft->realBegin() + n / 2, fft->realEnd(), 0 );
    fft->r2cExecute();
    std::copy( fft->complexBegin(), fft->complexEnd(), halfFFT.begin() );

    // Get full fft
    std::copy( in, in + n, fft->realBegin() );
    fft->r2cExecute();

    // Multiply
    for( int i = 0; i < fft->complexBufferSize(); ++i ) 
		fft->getComplexBuffer()[i] *= std::conj( halfFFT[i] );

    // IFFT
    fft->c2rExecute();
    
    // Find difference function d.
	std::vector<float> d( n / 2 );
    for( int j = 0; j < n / 2; ++j ) 
        d[j] = powerTerms[0] + powerTerms[j] - 2 * fft->getRealBuffer()[j] / n;

	return d;
	}

static std::vector<float> compute_dPrime( const float * in, Frame n )
	{
	// Get d
	std::vector<float> d = compute_d( in, n );

	// Transform d into d' via cumulative mean normalization
    d[0] = 1;
    double sum = 0;
    for( int tau = 1; tau < d.size(); ++tau ) 
		{
        sum += d[tau];
        if( sum == 0 ) d[tau] = 1;
        else d[tau] *= tau / sum;
		} 

	return d;
	}

//static Frame detectAutocorrelationPeak( float * data, Frame start, Frame end, float confidenceRation = 0.5 )
//	{
//	const float confidenceRequired = data[0] * confidenceRation;
//
//	// Discard initial peak. If ac data is never negative, we don't really care about it
//	Frame searchFrame = 0;
//	while( searchFrame < end && data[searchFrame] > 0 )  ++searchFrame; 
//	if( searchFrame < start ) searchFrame = start;
//
//	// Find potential peak frames using local maxima
//	std::set<Frame> peakFrames;
//	const Frame windowSize = 128;
//	const Frame hopSize = windowSize / 2;
//	for( ; searchFrame < end; searchFrame += hopSize )
//		{
//		// If local maxima is close to the center, it's a fine guess
//		auto localMaxFrameIter = std::max_element( data + searchFrame, data + searchFrame + windowSize );
//		if( std::abs( std::distance( data + searchFrame + windowSize / 2, localMaxFrameIter ) ) <= windowSize / 4 )
//			{
//			const Frame localMaxFrame = std::distance( data, localMaxFrameIter );
//			peakFrames.insert( localMaxFrame );
//			}
//		}
//
//	// Find max peak among scaled peaks
//	auto distribution = [end]( Frame x )
//		{
//		return 1.0f - float( x ) / float( end );
//		};
//	Frame bestFrame = -1;
//	float peakAmp = 0;
//	for( auto f : peakFrames )
//		{
//		const float scaledAmp_c = data[f] * distribution( f );
//		if( scaledAmp_c > peakAmp && scaledAmp_c > confidenceRequired )
//			{
//			peakAmp = scaledAmp_c;
//			bestFrame = f;
//			}
//		}
//
//	return bestFrame;
//	}

std::vector<float> Audio::getTotalEnergy() const
	{
	std::vector<float> energies( getNumChannels() );
	std::transform( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumChannels() ), energies.begin(), [&]( Channel channel ) {
		return std::accumulate( channelBegin( channel ), channelEnd( channel ), 0.0f, []( float x, Sample s ){ return x + s * s; } );
		} );
	return energies;
	}

std::vector<float> Audio::getEnergyDifference( const Audio & other ) const
	{
	const Audio inv = other.invertPhase();
	return Audio::mix( std::vector<const Audio *>{ this, &inv } ).getTotalEnergy();
	}

//============================================================================================================================================================
// Wavelengths
//============================================================================================================================================================

float Audio::getLocalWavelength( Channel channel, Frame start, Frame windowSize ) const
	{
	const float absoluteCutoff = 0.2f;

	flan_PROCESS_START( 0 );

	// Get d' - Cumulative mean normalized difference function (search for YIN algorithm for details)
	const std::vector<float> dPrime = compute_dPrime( getSamplePointer( channel, start ), windowSize );
	
	// Get local minima of dPrime, sorted from largest to smallest
	const std::vector<vec2> minima = findValleys( dPrime );
	
	// For wavelength finding, we want absolute minimum, as long as it's low enough
	vec2 bestMinima = {0,-1};
	for( const vec2 & m : minima )
		if( m.y() < absoluteCutoff ) 
			bestMinima = m;
    
    return bestMinima.x();
	}

std::vector<float> Audio::getLocalWavelengths( Channel channel, Frame start, Frame end, Frame windowSize, Frame hop, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( std::vector<float>() );

    if( end == -1 ) end = getNumFrames();

    std::vector<Frequency> out;
    for( Frame frame = start; frame + windowSize < end; frame += hop )
		{
		flan_CANCEL_POINT( std::vector<float>() );
        out.push_back( getLocalWavelength( channel, frame, windowSize ) );
		}
        
    return out;
	}

float Audio::getAverageWavelength( Channel channel, float minActiveRatio, float maxLengthSigma, 
	Frame start, Frame end, Frame windowSize, Frame hop, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( 0 );
	return getAverageWavelength( getLocalWavelengths( channel, start, end, windowSize, hop, canceller ), minActiveRatio, maxLengthSigma );
	}

float Audio::getAverageWavelength( const std::vector<float> & locals, float minActiveRatio, float maxLengthSigma ) const
	{
	flan_PROCESS_START( 0 );

	// Copy valid wavelengths
	std::vector<float> validLengths;

	const int numValids = locals.size() - std::count( locals.begin(), locals.end(), -1 );
	if( numValids <= minActiveRatio * locals.size() ) return -1; // Check that a sufficient amount of data had a detectable wavelength

	std::copy_if( locals.begin(), locals.end(), std::back_inserter( validLengths ), []( float x ){ return x != -1; } );

	const vec2 msd = meanAndStandardDeviation( validLengths );
	if( maxLengthSigma != -1 && msd.y() > maxLengthSigma ) return -1; // Check that the lengths didn't vary too much

	return msd.x();
	}

//============================================================================================================================================================
// Frequencies
//============================================================================================================================================================

Frequency Audio::getLocalFrequency( Channel channel, Frame start, Frame windowSize, flan_CANCEL_ARG_CPP ) const
	{
	const float minimumProximityCutoff = 0.2f;
	const float absoluteCutoff = 0.15f;

    flan_PROCESS_START( 0 );

	// Get d' - Cumulative mean normalized difference function (search for YIN algorithm for details)
	const std::vector<float> dPrime = compute_dPrime( getSamplePointer( channel, start ), windowSize );
	
	// Get local minima of dPrime, sorted from largest to smallest
	const std::vector<vec2> minima = findValleys( dPrime, -1, true );

	// Filter out anything that isn't at least within 20% of the minimum minima.
	// Also filter anything that isn't below a global cutoff ( dPrime is normalized to some degree ).
	const float cutoff = std::min( minima.back().y() * ( 1.0f + minimumProximityCutoff ), absoluteCutoff );
	
	// For frequency finding, we want small wavelength to minimize octave errors
	vec2 bestMinima = { float( windowSize ),-1 };
	for( const vec2 & m : minima )
		if( m.y() < cutoff && m.x() < bestMinima.x() ) 
			bestMinima = m;
    
    return float( getSampleRate() ) / bestMinima.x();
    }

std::vector<Frequency> Audio::getLocalFrequencies( Channel channel, Frame start, Frame end, Frame windowSize, Frame hop, flan_CANCEL_ARG_CPP ) const
    {
    flan_PROCESS_START( std::vector<Frequency>() );

    if( end == -1 ) end = getNumFrames();

    std::vector<Frequency> out;
    for( Frame frame = start; frame + windowSize < end; frame += hop )
        out.push_back( getLocalFrequency( channel, frame, windowSize, canceller ) );
        
    return out;
    }




	