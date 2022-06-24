#ifndef MACHINE_H
#define MACHINE_H

#ifdef _WIN32

#ifdef _WIN64
#define BIT_64
#endif


#endif

#ifdef __linux__

#ifdef _LP64
#define BIT_64
#endif

#endif



#endif
