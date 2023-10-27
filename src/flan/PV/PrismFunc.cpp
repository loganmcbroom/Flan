#include "PrismFunc.h"

#include "flan/PV/PVBuffer.h"

using namespace flan;

PrismFunc::PrismFunc()
	: func( nullptr )
	, policy( ExecutionPolicy::Parallel_Unsequenced )
	, null( true )
	{
	}

MF PrismFunc::operator()( 
	int n, 
	Second t, 
	int h, 
	Frequency f, 
	const std::vector<Magnitude> & hM 
	) const 
	{ 
	if( is_null() ) return MF{ 0, 0 };
	return func(n, t, h, f, hM); 
	}

bool PrismFunc::is_null() const 
	{ 
	return null; 
	}

ExecutionPolicy PrismFunc::get_execution_policy() const 
	{ 
	return policy;
	}