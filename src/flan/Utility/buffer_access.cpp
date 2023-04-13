#include "buffer_access.h"

namespace flan {

int buffer_access( int smallPos, int bigPos, int smallSize ) {
	return bigPos 		* smallSize 
		 + smallPos;
	}

int buffer_access( int smallPos, int mediumPos, int bigPos, int smallSize, int mediumSize ) {
	return bigPos 		* smallSize * mediumSize 
		 + mediumPos 	* smallSize 
		 + smallPos;
	} 

}