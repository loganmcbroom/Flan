#pragma once

/*
Utility functions for accessing single buffers storing multidimensional data.
*/

namespace flan {

int buffer_access( int smallPos, int bigPos, int smallSize );

int buffer_access( int smallPos, int mediumPos, int bigPos, int smallSize, int mediumSize );

}