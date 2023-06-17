#pragma once

#include <atomic>
#include <iostream>
#include <cmath>

namespace flan 
{

using Index = int;
using Second = float;
using Channel = int32_t;
using Frame = int32_t;
using Bin = int32_t;
using fFrame = float;
using fBin = float;
using Sample = float;
using Frequency = float;
using Magnitude = float;
using Pixel = int32_t;
using Harmonic = int32_t;
using FrameRate = float;
using Radian = float;
using Cycle = float;
using Amplitude = float;
using Decibel = float;

Amplitude decibelToAmplitude( Decibel d );
Decibel amplitudeToDecibel( Amplitude g );

const Radian pi = std::acos( -1.0f );
const Radian pi2 = pi * 2.0f;

};

/** These macros are used for logging function calls in debug mode, and checking for input nullity in all modes.
 */
#ifdef flan_LOG_FUNCTIONS
#define flan_FUNCTION_LOG std::cout << __FUNCTION__ << std::endl;
#else
#define flan_FUNCTION_LOG
#endif

#define flan_PROCESS_START( T ) flan_FUNCTION_LOG if( isNull() ) \
	{ std::cout << "Null input" << std::endl; return T; }

/** These macros are used to inject voluntary cancellation points into flan algorithms
 *	Cancellation points can be used to halt algorithms from other threads.
 */
#define flan_CANCELLABLE
#ifdef flan_CANCELLABLE
	static std::atomic<bool> default_canceller( false ); // Global flag used for default argument, doesn't change
	#define flan_CANCEL_POINT( T ) { if( canceller ) return T; }
	#define flan_CANCEL_ARG std::atomic<bool> & = default_canceller
	#define flan_CANCEL_ARG_CPP std::atomic<bool> & canceller
#else
	#define flan_CANCEL_POINT( T )
	#define flan_CANCEL_ARG void* = nullptr
	#define flan_CANCEL_ARG_CPP void*
#endif

