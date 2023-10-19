#pragma once

/* 
Windows are functions with domain [0,1]. They can be used as analysis windows, or as envelopes, among other things.
*/

namespace flan{

namespace Windows {

float hann( float x );

}

}