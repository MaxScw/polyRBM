
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
// #include "random.h"
#include "r250.h"
// #include "tools.h"

using namespace std;

const string R250::filename = "R250.state";

R250::R250(cuint _seed):
		   arrayEnd(&array[RANDOM_PREFETCH]),
		   pos(&array[250]),
		   other147(pos-147),
		   other250(pos-250)
{
	Initialize(_seed);

	Load();
// 	ifstream in("/dev/random");
// 	in.read((char*)(array),sizeof(uint)*RANDOM_PREFETCH);
}

// R250::~R250()
// {
// // 	delete []array;
// }

void
R250::Initialize(cuint _seed)
{
	srand(_seed);
	for (int i = 0; i < 250; i++)
	{
		array[i]     = rand();
	}
	Refresh();

// 	Print();
}

void
R250::Refresh()
{
	do
	{
	  *(pos++)		= *(other147++)^(*(other250++)); 	//the meaning of it all
	}
	while(pos != arrayEnd);
	pos			= array;

	do
	{
	  *(pos++)		= *(other147++)^(*(other250++)); 	//the meaning of it all
	}
	while(other147!=arrayEnd);

	other147			= array;
	do
	{
	  *(pos++)		= *(other147++)^(*(other250++)); 	//the meaning of it all
	}
	while(other250!=arrayEnd);

// 		 *pos		= (*other147)^(*other250); 	//the meaning of it all
// 		pos			++;
// 		other147++;
// 		other250++;
	dice 					= array;
	other250			= array;

}

// unsigned int*
// R250::GetPrefetchTable()
// {
// 	return &array[0];
// }


void
R250::Save() const
{
//  	ofstream file(filename.c_str());
// 	file << dice - array << endl;
// 	for (int i = 0 ; i < RANDOM_PREFETCH; ++i) file << array[i] << endl;
// 	file.close();
// 	cout << "r250 state written.\n";
}

void
R250::Load()
{
	if (0)
// 	string fname = filename+".backup";
//  	if(FileExists(filename))
	{
 		string fname = "R250.state";
	  
		std::ifstream file(fname.c_str());
		
		//filename = fn;

		unsigned offset;
		file >> offset; file.get(); dice = array + offset;
// 		cout << offset << endl;
		for (int i = 0 ; i < RANDOM_PREFETCH; ++i)
		{
			file >> array[i];
			file.get();
// 			cout << array[i] << endl;
		}

		file.close();
		cout << "loaded r250 state from " << fname << " .\n";
	}
	else
	{
// 		cout << "size of rng = " << sizeof(array) << endl;
//	RANDOM_PREFETCH

// 		cout << "loading r250 state from /dev/urandom...";
		
		ifstream in("/dev/urandom");
		in.read((char*)array,sizeof(array));
		
		for (uint i = 0; i < RANDOM_PREFETCH; ++i) array[i] &= RAND_MAX_INT;
		
// 		cout << "ready.\n";
	}
}

void
R250::Print(cuint limit)
{
	for (uint i = 0; i<limit; i++)
	{
		std::cout << "\n prefetch["
			 << std::setw(4)  << std::right << i << "] = "
			 << std::setw(14) << std::right << array[i] << ", ";
// 		printint(array[i]);
	}
}
