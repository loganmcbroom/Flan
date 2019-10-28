#pragma once

#include "Types.h"

#include "SpectralBuffer.h"

namespace xcdp::spec {

SpectralBuffer average( SpecInput in, int factor );

SpectralBuffer gate( SpecInput in, RealFunc cutoff );
SpectralBuffer gate( SpecInput in, double cutoff );

SpectralBuffer invert( SpecInput in, int lowerBound, int upperBound );

} // End namespace xcdp