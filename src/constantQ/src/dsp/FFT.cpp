/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Constant-Q library
    Copyright (c) 2013-2014 Queen Mary, University of London

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#include "FFT.h"

#include "MathUtilities.h"

#include <fftw3.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <execution>

class FFT::D
{
public:
    D(int n) : m_n(n) {
        m_kin  = (std::complex<float>*) fftwf_alloc_complex( m_n );
        m_kout = (std::complex<float>*) fftwf_alloc_complex( m_n );
        m_plan_forward = fftwf_plan_dft_1d( m_n, (fftwf_complex*) m_kin, (fftwf_complex*) m_kout, FFTW_FORWARD, FFTW_ESTIMATE );
        m_plan_inverse = fftwf_plan_dft_1d( m_n, (fftwf_complex*) m_kin, (fftwf_complex*) m_kout, FFTW_BACKWARD, FFTW_ESTIMATE );
    }

    ~D() {
        if( m_plan_forward ) fftwf_destroy_plan( m_plan_forward );
        if( m_plan_inverse ) fftwf_destroy_plan( m_plan_inverse );
        fftwf_free( m_kin );
        fftwf_free( m_kout );
    }

    void process(bool inverse, const std::vector<std::complex<float>> & in, std::vector<std::complex<float>> & out ) {
        std::copy( std::execution::par_unseq, in.begin(), in.end(), m_kin );
        fftwf_execute( inverse? m_plan_inverse : m_plan_forward );
        std::copy( std::execution::par_unseq, m_kout, m_kout + m_n, out.begin() );
        if( inverse ) 
            std::for_each( std::execution::par_unseq, out.begin(), out.end(), [scale = 1.0f / m_n]( std::complex<float> & z ){ z *= scale; } );
    }
    
private:
    const int m_n;
    std::complex<float> * m_kin;
    std::complex<float> * m_kout;
    fftwf_plan_s * m_plan_forward;
    fftwf_plan_s * m_plan_inverse;
};        

FFT::FFT(int n) :
    m_d( std::make_unique<D>( n ) )
{
}

FFT::~FFT()
{
}

void
FFT::process( bool inverse, const std::vector<std::complex<float>> & in, std::vector<std::complex<float>> & out )
{
    m_d->process( inverse, in, out );
}
    
class FFTReal::D
{
public:
    D(int n) : m_n(n) {
        if (n % 2) {
            throw std::invalid_argument("nsamples must be even in FFTReal constructor");
        }
        m_real = fftwf_alloc_real( m_n );
	    m_complex = (std::complex<float>*) fftwf_alloc_complex( m_n / 2 + 1 );
        m_plan_forward = fftwf_plan_dft_r2c_1d( m_n, m_real, (fftwf_complex*) m_complex, FFTW_ESTIMATE );
        m_plan_inverse = fftwf_plan_dft_c2r_1d( m_n, (fftwf_complex*) m_complex, m_real, FFTW_ESTIMATE );
    }

    ~D() {
        if( m_plan_forward ) fftwf_destroy_plan( m_plan_forward );
        if( m_plan_inverse ) fftwf_destroy_plan( m_plan_inverse );
        fftwf_free( m_real );
        fftwf_free( m_complex );
    }

    void forward( const std::vector<float> & real_in, std::vector<std::complex<float>> & out ) {
        for( int i = 0; i < m_n; ++i ) m_real[i] = real_in[i];
        fftwf_execute( m_plan_forward );
        for (int i = 0; i < m_n/2 + 1; ++i) out[i] = m_complex[i]; 
        for (int i = 1; i < m_n/2; ++i) out[m_n - i] = std::conj( out[i] );
    }

    // void forwardMagnitude(const float *ri, float *mo) {

    //     float *io = new float[m_n];

    //     forward(ri, mo, io);

    //     for (int i = 0; i < m_n; ++i) {
    //         mo[i] = sqrt(mo[i] * mo[i] + io[i] * io[i]);
    //     }

    //     delete[] io;
    // }

    void inverse( const std::vector<std::complex<float>> & in, std::vector<float> & real_out ) {
        for (int i = 0; i < m_n/2 + 1; ++i) m_complex[i] = in[i];
        fftwf_execute( m_plan_inverse );
        const float scale = 1.0 / m_n;
        for (int i = 0; i < m_n; ++i) real_out[i] = m_real[i] * scale;
    }

private:
    int m_n;
    std::complex<float> * m_complex;
    float * m_real;
    fftwf_plan_s * m_plan_forward;
    fftwf_plan_s * m_plan_inverse;
};

FFTReal::FFTReal(int n) :
    m_d( std::make_unique<D>( n ) ) 
{
}

FFTReal::~FFTReal()
{
}

void
FFTReal::forward( const std::vector<float> & real_in, std::vector<std::complex<float>> & out )
{
    m_d->forward( real_in, out );
}

//void
// FFTReal::forwardMagnitude(const float *ri, float *mo)
// {
//     m_d->forwardMagnitude(ri, mo);
// }

void
FFTReal::inverse( const std::vector<std::complex<float>> & in, std::vector<float> & realOut )
{
    m_d->inverse( in, realOut );
}


    
