#pragma once

#include "Types.h"

#include "SpectralBuffer.h"

namespace xcdp {

SpectralBuffer blur_average( SpecInput in, int factor );

} // End namespace xcdp