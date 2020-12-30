#include "xcdp/Utility/View.h"

namespace xcdp {

void View::forceAspectOnX( float aspect )
	{
	const float expansion = aspect * V.aspect() / U.aspect();
	const float c = U.d.midpoint();
	U.d.x1 = U.x1() * expansion + c * ( 1.0f - expansion );
	U.d.x2 = U.x2() * expansion + c * ( 1.0f - expansion );
	}

void View::forceAspectOnY( float aspect )
	{
	const float expansion = U.aspect() / V.aspect() / aspect;
	const float c = U.r.midpoint();
	U.r.x1 = U.y1() * expansion + c * ( 1.0f - expansion );
	U.r.x2 = U.y2() * expansion + c * ( 1.0f - expansion );
	}

void View::growToAspect( float aspect )
	{
	if( aspect * V.aspect() / U.aspect() > 1.0f )
		forceAspectOnX( aspect );
	else
		forceAspectOnY( aspect );
	}

void View::shrinkToAspect( float aspect )
	{
	if( aspect * V.aspect() / U.aspect() <= 1.0f )
		forceAspectOnX( aspect );
	else
		forceAspectOnY( aspect );
	}

}