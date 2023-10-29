#include "execution.h"

namespace flan {

bool is_linear( ExecutionPolicy a )
	{
	return a == ExecutionPolicy::Linear_Sequenced || a == ExecutionPolicy::Linear_Unsequenced;
	}

bool is_parallel( ExecutionPolicy a )
	{
	return ! is_linear( a );
	}

bool is_sequenced( ExecutionPolicy a )
	{
	return a == ExecutionPolicy::Linear_Sequenced || a == ExecutionPolicy::Parallel_Sequenced;
	}

bool is_unsequenced( ExecutionPolicy a )
	{
	return ! is_sequenced( a ); 
	}

ExecutionPolicy lowest_execution( const std::vector<ExecutionPolicy> & ps )
	{
	#ifdef __APPLE__
		return ExecutionPolicy::Linear_Sequenced;
	#else
		bool linear = false;
		for( auto p : ps )
			if( is_linear( p ) )
				linear = true;

		bool sequenced = false;
		for( auto p : ps )
			if( is_sequenced( p ) )
				linear = true;

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
	#endif
	}

}