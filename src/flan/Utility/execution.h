#pragma once

#include <execution>

namespace flan 
{

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
auto runtimeExecutionPolicyHandler( ExecutionPolicy policy, F f )
	{
	switch( policy )
		{
		case ExecutionPolicy::Linear_Sequenced: 		return f( std::execution::seq		);
		case ExecutionPolicy::Linear_Unsequenced: 		return f( std::execution::unseq		);
		case ExecutionPolicy::Parallel_Sequenced: 		return f( std::execution::par		);
		case ExecutionPolicy::Parallel_Unsequenced: 	return f( std::execution::par_unseq	);

		default: 										return f( std::execution::seq 		);
		};
	}

bool isLinear( ExecutionPolicy a );
bool isParallel( ExecutionPolicy a );
bool isSequenced( ExecutionPolicy a );
bool isUnsequenced( ExecutionPolicy a );

ExecutionPolicy lowestExecution( ExecutionPolicy a, ExecutionPolicy b );

}