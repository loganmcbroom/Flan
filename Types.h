#pragma once

#include <functional>

namespace xcdp {

class AudioBuffer;
class SpectralBuffer;

typedef std::function< double ( double time ) > RealFunc;
typedef std::vector<AudioBuffer> AudioBufferVec;
typedef const AudioBuffer & ProcInput;
typedef const AudioBufferVec & ProcInputVec;

typedef std::vector<SpectralBuffer> SpectralBufferVec;
typedef const SpectralBuffer & SpecInput;
typedef const SpectralBufferVec & SpecInputVec;

} // End namespace xcdp