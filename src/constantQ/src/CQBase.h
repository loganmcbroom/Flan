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

#ifndef CQBASE_H
#define CQBASE_H

#include <vector>
#include <complex>

/// A single complex-valued sample.
typedef std::complex<float> Complex;

/// A series of real-valued samples ordered in time.
typedef std::vector<float> RealSequence;

/// A series of real-valued samples ordered by bin (frequency or similar).
typedef std::vector<float> RealColumn;

/// A series of complex-valued samples ordered in time.
typedef std::vector<Complex> ComplexSequence;

/// A series of complex-valued samples ordered by bin (frequency or similar).
typedef std::vector<Complex> ComplexColumn;

/// A matrix of real-valued samples, indexed by time then bin number.
typedef std::vector<RealColumn> RealBlock;

/// A matrix of complex-valued samples, indexed by time then bin number.
typedef std::vector<ComplexColumn> ComplexBlock;

#endif
