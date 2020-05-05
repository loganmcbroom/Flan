#pragma once

#include <functional>

#include "AudioBuffer.h"

namespace xcdp {

class PVOC;
class Spectrum;
struct RealFunc;

class Audio : public AudioBuffer
{
public:
	typedef const std::vector<Audio> & Vec;
	typedef std::function< Audio ( const Audio &, size_t ) > Mod;
	
	Audio();
	Audio( const AudioBuffer::Format & other );
	Audio( const std::string & filename );

	//======================================================
	//	Overloads
	//======================================================

	Audio operator+( const Audio & other ) const; //Mix
	Audio operator-() const; //Invert phase
	Audio operator-( const Audio other ) const; //Add negative

	//========================================================
	// Information
	//========================================================

	float getTotalEnergy() const;

	Audio graph( const std::string & filename, size_t width = 2048, size_t height = 512 ) const;

	//======================================================
	//	Conversions
	//======================================================

	PVOC convertToPVOC( const size_t frameSize = 2048, const size_t overlaps = 16 ) const;
	//PVOC convertToPVOC_old( const size_t frameSize = 2048, const size_t overlaps = 16 ) const;
	Spectrum convertToSpectrum() const;

	Audio convertToMidSide() const;
	Audio convertToLeftRight() const;

	Audio convertToStereo() const;
	Audio convertToMono() const;

	//===========================================================================================
	//	Procs
	//===========================================================================================

	Audio invertPhase() const;
	Audio modifyVolume( RealFunc volumeLevel ) const;
	Audio setVolume( RealFunc level ) const;
	Audio waveshape( RealFunc shaper ) const;
	Audio pan( RealFunc panAmount ) const;
	Audio widen( RealFunc widenAmount ) const;
	Audio iterate( size_t n, Audio::Mod mod = nullptr, bool fbIterate = false ) const;
	Audio reverse() const;
	Audio cut( float startTime, float endTime ) const;
	Audio repitch( RealFunc factor, size_t granul = 64, size_t qual = 0 ) const;
	Audio convolve( const std::vector<RealFunc> & ) const;
	Audio delay( float delayTime, size_t numDelays, float decayAmount = .5, Audio::Mod mod = nullptr, bool fbIterate = true ) const;
	Audio fades( float fadeTime = .05 ) const;
	Audio lowPass( RealFunc cutoff, size_t taps = 64 ) const;

	// Unfinished
	//Audio freeze( real_t freezeTime, real_t freezeLength, real_t delay, real_t rand = 0, Audio::Mod mod = nullptr ) const;
	// drunk walk / scramble
	// extend loop
	// extend repititions?

	//========================================================
	// Multi-In Procs
	//========================================================

	static Audio mix( Audio::Vec ins, 
		std::vector<RealFunc> balances = std::vector<RealFunc>(),
		std::vector<float> startTimes = std::vector<float>() );
	static Audio join( Audio::Vec ins );
};

} // End namespace xcdp
