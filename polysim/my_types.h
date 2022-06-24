#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include "machine.h"
#include <string>



#ifdef BIT_64

typedef  int64_t			fast_int;
typedef uint64_t			fast_uint;
typedef 	uint32_t		uint;

#else

typedef  int32_t			fast_int;
typedef uint32_t			fast_uint;

typedef 	uint32_t		uint;

#endif

typedef fast_uint     fuint;
typedef fast_int      fint;

typedef const fint&		cfint;
typedef const fuint&	cfuint;

typedef const int& 	 	cint;
typedef const uint& 	cuint;

typedef unsigned char uchar;

typedef const std::string& str_ref;


#endif
