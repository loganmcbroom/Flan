#pragma once

//x should sit between 0 and 1

namespace xcdp{

struct RealFunc;

namespace Windows {

float Hann( float x );

/** Origin centered Hann DFT * 2 */
float HannDFT2( float K );

}}