#ifndef R250_H
#define R250_H

#include "my_types.h"

#define RANDOM_PREFETCH			256
#define RAND_MAX_INT			2147483647
#define RAND_MAX_DOUBLE			2147483647.0
#define RAND_NORMALIZE			4.65661287307739e-10

// template < class T >
class R250
{
private:
  uint	*dice,*arrayEnd,*pos,*other147,*other250,
   			array[RANDOM_PREFETCH];

 	static const std::string filename;
public:


	R250(cuint _seed = 1013904223);
// 			   ~R250();

//   unsigned int* GetPrefetchTable();
  inline uint	Rand();
	void	Refresh();
  void	Initialize(cuint _seed = 1013904223);
	void  Print(cuint limit = (RANDOM_PREFETCH > 500 ? 500 : RANDOM_PREFETCH));

	double Uniform();
	uint operator () () {return Rand();}

	void Save() const;
	void Load() ;

};


inline uint
R250::Rand()
{
	if(dice - array == RANDOM_PREFETCH) Refresh();
	return *(dice++);
}

inline double R250::Uniform()
{
	return ((double)Rand())*RAND_NORMALIZE;
}

/*
Number 	..	my type of random numbers
stack	..	list of 250 (ideally random) numbers
max	..	highest possible random number (binary: 11111...11 )
norm	..	shall be the double-version of the max-value
pos	..	actual position in the array (of actual random number)
other	..	other position for the central XOR-operation:
		stack[pos] = stack[pos-147] ^ stack[pos-250],
		where pos - 147 is equal pos + 103 =: other
		and pos - 250 is pos itself because the array is of length 250
		and we walk through periodically
bits	..	number of bits of the given "number" type
*/

#endif
