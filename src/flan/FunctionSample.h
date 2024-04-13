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

	template<VectorSpace T = O>
	std::vector<T> exclusive_scan( T start ) const
		{
		std::vector<O> out( size() );
		if( is_constant() )	flan::for_each_i( size(), ExecutionPolicy::Parallel_Unsequenced, [&out, val = get_constant()]( size_t i ){ out[i] = i * val; } );
		else std::exclusive_scan( FLAN_PAR_UNSEQ get_vector().begin(), get_vector().end(), out.begin(), start );
		return out;
		}
};

}