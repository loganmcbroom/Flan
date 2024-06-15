#pragma once

#include <functional>
#include <complex>
#include <type_traits>
#include <random>
#include <ctime>
#include <algorithm>
#include <concepts>
#include <ranges>
#include <variant>

#include "flan/defines.h"
#include "flan/Utility/vec2.h"
#include "flan/Utility/Rect.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/execution.h"
#include "flan/Utility/buffer_access.h"
#include "flan/Utility/function_traits.h"
#include "bmp/bitmap_image.hpp"
#include "flan/Graph.h"
#include "flan/FunctionSample.h"

#undef min
#undef max
#undef clamp

namespace flan {

class Graph;

/** This is a simple wrapper around std::function that comes with some utilities and optimizations.
 *	These are primarily created by passing either lambda functions or constants to Audio or PV methods.
 */
template< typename I, typename O >
struct Function
{
	using StdFuncType = std::function< O ( I )>;
	using ReturnType = O;
	using ArgType = I;

	Function( const Function & ) 				= delete;
	Function & operator=( const Function & ) 	= delete;
	Function( Function && ) 					= default;
	Function & operator=( Function && ) 		= default;
	~Function() 								= default;

	template<typename T>
	requires std::convertible_to<T, O>
	Function( T t0 ) 
		: f( static_cast<O>( t0 ) )
		, execution_policy( ExecutionPolicy::Parallel_Unsequenced )
		{}

	template<typename T>
#ifndef __APPLE__
	// Christ apple clang, get your shit together
	requires std::convertible_to<T, StdFuncType>
#endif
	Function( T && f_, ExecutionPolicy policy = ExecutionPolicy::Parallel_Unsequenced ) 
		: f( std::move( f_ ) )
		, execution_policy( policy )
		{}

	Function copy() const 
		{
		if( is_constant() )
			return Function( std::get<O>( f ) );
		else
			return Function( std::get<StdFuncType>( f ), execution_policy );
		}

	bool is_constant() const 
		{
		return std::holds_alternative<O>( f );
		}

	ExecutionPolicy get_execution_policy() const 
		{ 
		return execution_policy; 
		}

	// Function() 
	// 	: f( []( I ){ return O( 0 ); } ) 
	// 	, execution_policy( ExecutionPolicy::Parallel_Unsequenced )
	// 	{}

	/** Function application */
	O operator()( I t ) const 
		{ 
		if constexpr ( std::is_copy_constructible_v<O> )
			{
			if( is_constant() )
				return std::get<O>( f );
			else 
				return std::get<StdFuncType>( f )( t ); 
			}
		else // If O isn't copy constructable, there shouldn't be a need for constants that I know of
			// Another way to do this is test if O has a copy method and use that for constants in this else
			{
			return std::get<StdFuncType>( f )( t ); 
			}
		}

	// static Function<I,O> uniformDistribution( Function<I,O> lowerBound, Function<I,O> upperBound )
	// 	{
	// 	static std::default_random_engine rng( static_cast<unsigned int>( std::time( nullptr ) ) );
	// 	return Function<I,O>( [lowerBound, upperBound]( I x ) -> O
	// 		{ 
	// 		return std::uniform_real_distribution<float>( lowerBound( x ), upperBound( x ) )( rng ); 
	// 		}, lowest_execution( lowerBound.get_execution_policy(), upperBound.get_execution_policy() ) );
	// 	}

	// static Function<I,O> normalDistribution( const Function<I,O> & mean, const Function<I,O> & sigma )
	// 	{
	// 	static std::default_random_engine rng( static_cast<unsigned int>( std::time( nullptr ) ) );
	// 	return Function<I,O>( [mean, sigma]( I x ) -> O
	// 		{ 
	// 		const O m = mean( x );
	// 		const O s = sigma( x );
	// 		if( s <= 0 ) return m;
	// 		return std::normal_distribution<float>( m, s )( rng ); 
	// 		}, lowest_execution( mean.get_execution_policy(), sigma.get_execution_policy() ) );
	// 	}

	template<typename Input = I>
	requires std::convertible_to<Input, float>
	Function periodize( float period ) const
		{
		if( is_constant() ) return copy();

		auto p_f = std::make_shared<Function>( copy() );
		return Function( [p_f, p = period]( float t )
			{
			return (*p_f)( std::fmod( t, p ) + ( t < 0.0f ? p : 0.0f ) );
			}, get_execution_policy() );
		}

	template<typename Input = I>
	requires std::convertible_to<Input, float>
	FunctionSample<O> sample( int start, int end, float scale ) const
		{
		if( is_constant() )
			return FunctionSample<O>( operator()(0), end - start );
		std::vector<O> out( end - start );
		runtime_execution_policy_handler( get_execution_policy(), [&]( auto policy ){
			std::for_each( FLAN_POLICY iota_iter( start ), iota_iter( end ), [&]( int x )
				{ 
				out[x-start] = operator()( x * scale ); 
				} ); 
			} );
		return FunctionSample<O>( std::move( out ) );
		}

	template<typename Input = I>
	requires std::convertible_to<Input, vec2>
	FunctionSample<O> sample( float x_start, float x_end, float x_scale, float y_start, float y_end, float y_scale ) const
		{
		const int x_size = std::ceil( x_end - x_start );
		const int y_size = std::ceil( y_end - y_start );
		if( is_constant() )
			return FunctionSample<O>( operator()({0,0}), x_size * y_size );
		std::vector<O> out( x_size * y_size );
		runtime_execution_policy_handler( get_execution_policy(), [&]( auto policy ){
			std::for_each( FLAN_POLICY iota_iter( x_start ), iota_iter( x_end ), [&]( int x ){ 
				for( int y = y_start; y < y_end; ++y )
					out[ buffer_access( (y-y_start), (x-x_start), y_size ) ] = operator()( vec2( x * x_scale, y * y_scale ) ); 
				} ); 
			} );
		return out;
		}

	template<typename Input = I>
	requires std::convertible_to<Input, vec2>
	std::vector<O> sample( int x_start, int x_end, float xScale, const std::vector<float> & yPositions ) const
		{
		const int xSize = x_end - x_start;
		const int ySize = yPositions.size();
		std::vector<O> out( xSize * ySize );
		runtime_execution_policy_handler( get_execution_policy(), [&]( auto policy ){
			std::for_each( FLAN_POLICY iota_iter( x_start ), iota_iter( x_end ), [&]( int x ){ 
				for( int y = 0; y < ySize; ++y )
					out[ buffer_access( y, (x-x_start), ySize ) ] = operator()( vec2( x * xScale, yPositions[y] ) ); 
				} ); 
			} );
		return out;
		}

	/** Create a graph of the function. 
	 * \param view The area of the cartesian plane that should be viewed.
	 * \param domain The interval on which the function should be graphed.
	 * \param width The bmp width.
	 * \param height The bmp height.
	 */
	template<typename Input = I, typename Output = O>
	requires std::convertible_to<Input, float> && std::convertible_to<Output, float>
	Graph convert_to_graph( 
		Rect view = { -5, -5, 5, 5 }, 
		Interval domain = Interval::R, 
		Pixel width = -1, 
		Pixel height = -1 
		) const
		{
		Graph g( width, height );
		g.set_view( view );
		g.fill_image( Color::White );
		g.draw_linear_grid( 1, 1, 0, Color( 200, 200, 200 ) );
		g.draw_axes( 0, Color::Black );
		g.draw_function( *this, domain, -1, Color::Black );
		return g;
		}

	/** Create and save a graph of the function. 
	 * \param filename The location to save the bmp.
	 * \param view The area of the cartesian plane that should be viewed.
	 * \param domain The interval on which the function should be graphed.
	 * \param width The bmp width.
	 * \param height The bmp height.
	 */
	template<typename Input = I, typename Output = O>
	requires std::convertible_to<Input, float> && std::convertible_to<Output, float>
	void save_to_bmp( 
		const std::string & filename, 
		Rect view = { -5, -5, 5, 5 }, 
		Interval domain = Interval::R, 
		Pixel width = -1, 
		Pixel height = -1 
		) const
		{
		auto bmp = convert_to_graph( view, domain, width, height );
		bmp.save_image( filename );
		}

	/** Constructor when argument is convertable to std::function< float ( float, float ) >.
	 */
	// template<typename T>
	// requires std::convertible_to< T, std::function< float ( float, float ) > >
	// Function( T f ) 
	// 	: Function( [f = std::move(f)]( vec2 v ) { return f( v.x(), v.y() ); } ) 
	// 	{}

	/** Constructor when argument is convertable to a std::function taking a pair of float convertable types
	 */
	// template<typename T, typename I1, typename I2>
	// requires std::convertible_to< T, std::function< float ( std::pair<I1, I2> ) > >
	// 	&& std::convertible_to< I1, float >
	// 	&& std::convertible_to< I2, float >
	// Function( T f ) 
	// 	: Function( [f = std::move(f)]( vec2 v ) 
	// 		{ 
	// 		return f( std::make_pair( 
	// 			static_cast<I1>( v.x() ), 
	// 			static_cast<I2>( v.y() ) ) ); 
	// 		} ) 
	// 	{}

	/* Call op taking two floats when vec2 would normally be needed. */
	template<typename Input = I>
	requires std::convertible_to<Input, vec2>
	float operator()( float x, float y ) const { return Function::operator()( vec2{ x, y } ); }

protected:
	std::variant<O, StdFuncType> f;
	const ExecutionPolicy execution_policy;
};

/** Generate an ADSR envelope. Generated envelopes always range from 0 to 1.
 *	The Exp parameters dictate how the envelope curves between fixed points. Using an Exp of 1 gives a linear movement.
 *	Values between 0 and 1 give curves that move rapidly toward their destination before slowing down.
 *	Values larger than 1 do the opposite. Other values will give unusual, but well defined, behaviour, so no clamping is applied.
 *
 * \param a Attack length in time.
 * \param d Decay length in time.
 * \param s Sustain length in time.
 * \param r Release length in time.
 * \param sLvl Sustain level.
 * \param aExp Exponent of the attack curve.
 * \param dExp Exponent of the decay curve.
 * \param rExp Exponent of the release curve. 
 */
Function<Second, Amplitude> ADSR( 
	Second attack_time, 
	Second decay_time, 
	Second sustain_time, 
	Second release_time, 
	Amplitude sustain_level,
	float attack_curve_exponent = 1, 
	float decay_curve_exponent = 1, 
	float release_curve_exponent = 1 
	);

// Each waveform has a period and amplitude of one.
namespace waveforms {

extern const std::function<Amplitude (Second)> sine; 	
extern const std::function<Amplitude (Second)> square; 	
extern const std::function<Amplitude (Second)> saw; 		
extern const std::function<Amplitude (Second)> triangle; 

}

} // End namespace flan
