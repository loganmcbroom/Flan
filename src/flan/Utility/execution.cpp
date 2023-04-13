#include "execution.h"

namespace flan {

bool isLinear( ExecutionPolicy a )
	{
	return a == ExecutionPolicy::Linear_Sequenced || a == ExecutionPolicy::Linear_Unsequenced;
	}

bool isParallel( ExecutionPolicy a )
	{
	return ! isLinear( a );
	}

bool isSequenced( ExecutionPolicy a )
	{
	return a == ExecutionPolicy::Linear_Sequenced || a == ExecutionPolicy::Parallel_Sequenced;
	}

bool isUnsequenced( ExecutionPolicy a )
	{
	return ! isSequenced( a ); 
	}

ExecutionPolicy lowestExecution( ExecutionPolicy a, ExecutionPolicy b )
	{
	const bool linear = isLinear( a ) || isLinear( b );
	const bool sequenced = isSequenced( a ) || isSequenced( b );
	if( linear )
		{
		if( sequenced ) return ExecutionPolicy::Linear_Sequenced;
		else return ExecutionPolicy::Linear_Unsequenced;
		}
	else
		{
		if( sequenced ) return ExecutionPolicy::Parallel_Sequenced;
		else return ExecutionPolicy::Parallel_Unsequenced;
		}
	}

}