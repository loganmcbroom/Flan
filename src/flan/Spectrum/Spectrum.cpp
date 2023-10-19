//#include "flan/Spectrum.h"

//#include <fftw3.h>
//
//#include "flan/Audio/Audio.h"
//
//namespace flan {
//
//Audio Spectrum::convert_to_audio()
//	{
//	const size_t outNumSamples = ( get_num_bins() - 1 ) * 2;
//
//	Audio::Format format;
//	format.num_channels = get_num_channels();
//	format.num_frames = outNumSamples;
//	format.sample_rate = get_sample_rate();
//	Audio out( format );
//
//	//Allocate fftw buffers and plan
//	std::complex<double> * fftIn = ( std::complex<double> * ) fftw_alloc_complex( get_num_bins() );
//	double *fftOut = fftw_alloc_real( outNumSamples );
//	fftw_plan plan = fftw_plan_dft_c2r_1d( int( outNumSamples ), (fftw_complex*) fftIn, fftOut, FFTW_ESTIMATE );
//
//	for( size_t channel = 0; channel < get_num_channels(); ++channel )
//		{
//		//Read buffer data into fftIn
//		for( size_t bin = 0; bin < get_num_bins(); ++bin )
//			fftIn[bin] = getSpectra( channel, bin );
//
//		//ifft
//		fftw_execute( plan );
//
//		//Read fftOut into out audio
//		for( size_t sample = 0; sample < outNumSamples; ++sample )
//			out.set_sample( channel, sample, fftOut[sample] );
//		}
//
//	fftw_destroy_plan( plan );
//	fftw_free( fftOut );
//	fftw_free( fftIn );
//
//	return out;
//	}



//}
