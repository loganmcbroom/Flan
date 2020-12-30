#pragma once

#include <atomic>

namespace xcdp 
{

using Time = float;
using Channel = int32_t;
using Frame = int32_t;
using Bin = int32_t;
using Sample = float;
using Frequency = float;
using Magnitude = float;
using Pixel = int32_t;

};

/** These macros are used for logging function calls in debug mode, and checking for input nullity in all modes.
 */
#ifndef NDEBUG
#define XCDP_FUNCTION_LOG std::cout << __FUNCTION__ << std::endl;;
#else
#define XCDP_FUNCTION_LOG
#endif

#define XCDP_PROCESS_START( T ) XCDP_FUNCTION_LOG if( isNull() ) \
	{ std::cout << "Null input" << std::endl; return T; }

/** These macros are used to inject voluntary cancellation points into xcdp algorithms
 *	Cancellation points can be used to halt algorithms from other threads.
 */
#define XCDP_CANCELLABLE
#ifdef XCDP_CANCELLABLE
	static std::atomic<bool> default_canceller( false ); // Global flag used for default argument, doesn't change
	#define XCDP_CANCEL_POINT( T ) { if( canceller ) return T; }
	#define XCDP_CANCEL_ARG std::atomic<bool> & = default_canceller
	#define XCDP_CANCEL_ARG_CPP std::atomic<bool> & canceller
#else
	#define XCDP_CANCEL_POINT( T )
	#define XCDP_CANCEL_ARG void* = nullptr
	#define XCDP_CANCEL_ARG_CPP void*
#endif

