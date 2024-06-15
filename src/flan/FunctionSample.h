#pragma once 

#include <concepts>
#include <numeric>

namespace flan {

template<typename T>
concept VectorSpace = requires( T a, T b, size_t c )
	{ 
	{a + b} -> std::same_as<T>; 
	{a * c} -> std::same_as<T>;
	};

template< typename I, typename O >
struct Function;

template<typename O>
struct FunctionSample
{
	std::variant<O, std::vector<O>> value;
	const size_t vec_size;

	FunctionSample( std::vector<O> && v )
		: value( v )
		, vec_size( v.size() )
		{}

	FunctionSample( O v, size_t n )
		: value( v )
		, vec_size( n )
		{}

	template<typename F>
	void for_each( const F & transformer )
		{
		if( std::holds_alternative<O>( value ) )
			{
			transformer( std::get<O>( value ) );
			}
		else
			{
			auto & v = std::get<std::vector<O>>( value );
			flan::for_each_i( v.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int i )
				{
				transformer( v[i] );
				} );
			}
		}

	template<typename F, typename T = liph::function_return_type<F>>
	FunctionSample<T> transform( const F & transformer ) const
		{
		if( std::holds_alternative<O>( value ) )
			{
			return FunctionSample<T>( transformer( std::get<O>( value ) ), vec_size );
			}
		else
			{
			auto & v = std::get<std::vector<O>>( value );
			std::vector<T> out( v.size() );
			flan::for_each_i( v.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int i )
				{
				out[i] = transformer( v[i] );
				} );
			return FunctionSample<T>( std::move( out ) );
			}
		}

	bool is_constant() const
		{
		return std::holds_alternative<O>( value );
		}

	O & get_constant()
		{
		return std::get<O>( value );
		}

	const O & get_constant() const
		{
		return std::get<O>( value );
		}

	std::vector<O> & get_vector()
		{
		return std::get<std::vector<O>>( value );
		}

	const std::vector<O> & get_vector() const
		{
		return std::get<std::vector<O>>( value );
		}

	template<typename C>
	//requires std::convertible_to<O, C>
	FunctionSample<C> convert_to() const
		{
		if( is_constant() )
			{
			return FunctionSample<C>( C( get_constant() ), size() );
			}
		else
			{
			std::vector<C> v;
			v.reserve( size() );
			for( auto & x : get_vector() )
				v.push_back( x );
			return std::move( v );
			}
		}

	size_t size() const
		{
		return vec_size;
		}

	O & operator[]( int n )
		{
		if( is_constant() ) return get_constant();
		else return get_vector()[n];
		}

	O operator[]( int n ) const
		{
		if( is_constant() ) return get_constant();
		else return get_vector()[n];
		}

	Function<Second, O> to_time_function( FrameRate sr )
		{
		return [sr, this]( Second t ){ return operator[]( t * sr ); };
		}
	
	// Function<TF, O> to_time_frequency_function( FrameRate frames_per_second, float bins_per_frequency, Bin num_bins )
	// 	{
	// 	return [=, this]( TF tf )
	// 		{ 
	// 		const Frame frame = tf.t * frames_per_second;
	// 		const Bin bin = tf.f * bins_per_frequency;
	// 		return operator[]( buffer_access( bin, frame, num_bins ) ); 
	// 		};
	// 	}

	template<VectorSpace T = O>
	T accumulate() const
		{
		if( is_constant() )
			return get_constant() * vec_size;
		else
			{
			auto & v = get_vector();
			if constexpr ( std::is_integral_v<T> )
				return std::accumulate( v.begin(), v.end(), T(0) );
			else
				return std::accumulate( v.begin(), v.end(), T() );
			}
		}

	template<VectorSpace T = O, class BinaryOperation>
	std::vector<T> exclusive_scan( T start, BinaryOperation binary_op ) const
		{
		std::vector<O> out( size() );
		const auto & vec = is_constant()? std::vector<O>( size(), get_constant() ) : get_vector();
		std::exclusive_scan( FLAN_PAR_UNSEQ vec.begin(), vec.end(), out.begin(), start, binary_op );
		return out;
		}

	O maximum() const
		{
		if( is_constant() ) return get_constant();
		return *std::max_element( get_vector().begin(), get_vector().end() );
		}

	template<class MapOp>
	float maximum( MapOp op ) const
		{
		if( is_constant() ) return op( get_constant() );
		return op( *std::ranges::max_element( get_vector(), std::ranges::less(), op ) );
		}
};

}