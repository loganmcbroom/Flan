//#include "xcdp/Spectrum.h"

//#include <fftw3.h>
//
//#include "xcdp/Audio.h"
//
//namespace xcdp {
//
//Audio Spectrum::convertToAudio()
//	{
//	const size_t outNumSamples = ( getNumBins() - 1 ) * 2;
//
//	Audio::Format format;
//	format.numChannels = getNumChannels();
//	format.numFrames = outNumSamples;
//	format.sampleRate = getSampleRate();
//	Audio out( format );
//
//	//Allocate fftw buffers and plan
//	std::complex<double> * fftIn = ( std::complex<double> * ) fftw_alloc_complex( getNumBins() );
//	double *fftOut = fftw_alloc_real( outNumSamples );
//	fftw_plan plan = fftw_plan_dft_c2r_1d( int( outNumSamples ), (fftw_complex*) fftIn, fftOut, FFTW_ESTIMATE );
//
//	for( size_t channel = 0; channel < getNumChannels(); ++channel )
//		{
//		//Read buffer data into fftIn
//		for( size_t bin = 0; bin < getNumBins(); ++bin )
//			fftIn[bin] = getSpectra( channel, bin );
//
//		//ifft
//		fftw_execute( plan );
//
//		//Read fftOut into out audio
//		for( size_t sample = 0; sample < outNumSamples; ++sample )
//			out.setSample( channel, sample, fftOut[sample] );
//		}
//
//	fftw_destroy_plan( plan );
//	fftw_free( fftOut );
//	fftw_free( fftIn );
//
//	return out;
//	}



//}
