#pragma once

#include <vector>

#ifndef __APPLE__
#include <execution>
#define FLAN_PAR_UNSEQ std::execution::par_unseq,
#define FLAN_PAR_SEQ std::execution::par,
#define FLAN_POLICY policy,
#else
#define FLAN_PAR_UNSEQ
#define FLAN_PAR_SEQ
#define FLAN_POLICY
#endif

#include "iota_iter.h"

namespace flan {

enum class ExecutionPolicy 
	{
	Linear_Sequenced,
	Linear_Unsequenced,
	Parallel_Sequenced,
	Parallel_Unsequenced,
	};


// The stl uses template overloads on execution policy types to decide algorithm behaviour rather than taking an enum or something.
// Because of that, you can't just decide at runtime what execution policy to use without a wrapper like this.
// Flan functions carry a flan::ExecutionPolicy around and use this function to decide which stl policy overload should be used at runtime.
template<class F>
auto runtime_execution_policy_handler( ExecutionPolicy policy, F f )
	{
#ifdef __APPLE__
	return f( 0 );
#else
	switch( policy )
		{
		case ExecutionPolicy::Linear_Sequenced: 		return f( std::execution::seq		);
		case ExecutionPolicy::Linear_Unsequenced: 		return f( std::execution::unseq		);
		case ExecutionPolicy::Parallel_Sequenced: 		return f( std::execution::par		);
		case ExecutionPolicy::Parallel_Unsequenced: 	return f( std::execution::par_unseq	);

		default: 										return f( std::execution::seq 		);
		};
#endif
	}

ExecutionPolicy lowest_execution( const std::vector<ExecutionPolicy> & );

template<typename... Ts>
using AllExePolicies = typename std::enable_if_t<std::conjunction_v< std::is_convertible<Ts, ExecutionPolicy>... >>;
template<typename ... Ts, typename = AllExePolicies<Ts...> >
ExecutionPolicy lowest_execution( 
	Ts ... ins 
	)
	{
	std::vector<ExecutionPolicy> policies;
	( policies.push_back( ins ), ... ); // Fold expression
	return lowest_execution( policies );
	}

template<typename ... Fs >
ExecutionPolicy lowest_execution( 
	const Fs & ... ins 
	)
	{
	std::vector<ExecutionPolicy> policies;
	( policies.push_back( ins.get_execution_policy() ), ... ); // Fold expression
	return lowest_execution( policies );
	}

template<typename F>
void for_each_i( int end, ExecutionPolicy p, F f )
	{
	#ifdef __APPLE__
		std::for_each( iota_iter( 0 ), iota_iter( end ), f );
	#else 
		runtime_execution_policy_handler( p, [&]( auto policy ) 
			{
			std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( end ), f );
			} );
	#endif
	}

}