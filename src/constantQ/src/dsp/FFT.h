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

#ifndef FFT_H
#define FFT_H

#include <memory>
#include <vector> 
#include <complex>

class FFT  
{
public:
    /**
     * Construct an FFT object to carry out complex-to-complex
     * transforms of size nsamples. nsamples does not have to be a
     * power of two.
     */
    FFT(int nsamples);
    ~FFT();

    /**
     * Carry out a forward or inverse transform (depending on the
     * value of inverse) of size nsamples, where nsamples is the value
     * provided to the constructor above.
     *
     * realIn and (where present) imagIn should contain nsamples each,
     * and realOut and imagOut should point to enough space to receive
     * nsamples each.
     *
     * imagIn may be NULL if the signal is real, but the other
     * pointers must be valid.
     *
     * The inverse transform is scaled by 1/nsamples.
     */
    void process( bool inverse, const std::vector<std::complex<float>> & in, std::vector<std::complex<float>> & out );
    
private:
    class D;
    std::unique_ptr<D> m_d;
};

class FFTReal
{
public:
    /**
     * Construct an FFT object to carry out real-to-complex transforms
     * of size nsamples. nsamples does not have to be a power of two,
     * but it does have to be even. (Use the complex-complex FFT above
     * if you need an odd FFT size. This constructor will throw
     * std::invalid_argument if nsamples is odd.)
     */
    FFTReal(int nsamples);
    ~FFTReal();

    /**
     * Carry out a forward real-to-complex transform of size nsamples,
     * where nsamples is the value provided to the constructor above.
     *
     * realIn, realOut, and imagOut must point to (enough space for)
     * nsamples values. For consistency with the FFT class above, and
     * compatibility with existing code, the conjugate half of the
     * output is returned even though it is redundant.
     */
    void forward( const std::vector<float> & real_in, std::vector<std::complex<float>> & out );

    /**
     * Carry out a forward real-to-complex transform of size nsamples,
     * where nsamples is the value provided to the constructor
     * above. Return only the magnitudes of the complex output values.
     *
     * realIn and magOut must point to (enough space for) nsamples
     * values. For consistency with the FFT class above, and
     * compatibility with existing code, the conjugate half of the
     * output is returned even though it is redundant.
     */
    //void forwardMagnitude(const float *realIn, float *magOut);

    /**
     * Carry out an inverse real transform (i.e. complex-to-real) of
     * size nsamples, where nsamples is the value provided to the
     * constructor above.
     *
     * realIn and imagIn should point to at least nsamples/2+1 values;
     * if more are provided, only the first nsamples/2+1 values of
     * each will be used (the conjugate half will always be deduced
     * from the first nsamples/2+1 rather than being read from the
     * input data).  realOut should point to enough space to receive
     * nsamples values.
     *
     * The inverse transform is scaled by 1/nsamples.
     */
    void inverse( const std::vector<std::complex<float>> & in, std::vector<float> & real_out );

private:
    class D;
    std::unique_ptr<D> m_d;
};    

#endif
