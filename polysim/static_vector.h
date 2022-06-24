#ifndef STATIC_VECTOR_H
#define STATIC_VECTOR_H

#include <iostream>
#include "stdlib.h"

template < class  T, unsigned _size > class StatVec;

template < class  T, unsigned _size >
std::ostream& operator << (std::ostream & stream, const StatVec<T,_size> & v);

/**
 * @brief A vector class, where the size is known at compile-time.
 * 
 * This can have benefits in adress calculation and memory layout.
 * The usage of the class should be analogous to the STL vector < T > class.
 * 
 * 
 * @param T the type of the stored data.
 * @param _size The number of elements of type T, which have to be stored.
 * 
 * 
 * 
 * **/
template < class  T, unsigned _size >
class StatVec
{
protected:
  
    /// @note The data array is not created dynamically using "new"
    T data[_size];
    
    template < class Container > 
    void checkSameLength( const Container& c ) const
    {
      if (c.size() != this->size())
      {
// 	std::ostringstream strm; 
//	strm << "Tried to combine StaticArray object with Container of different size.";
//	strm << "this->size() : " << this->size() << ", container size: " << c.size() << std::endl;
//	strm << "first elements: " << this->operator[](0) << ", "<<c[0]<<std::endl;
//	std::cout << strm; 
	exit(0);
// 	handleError(strm.str());
      }
    }

public:
  
  typedef T 	   value_type;
  typedef       T* iterator;
  typedef const T* const_iterator;

	StatVec(){ for (unsigned i = 0; i < _size; ++i) data[i] = T(0);  }

	StatVec(const T& value)
	{
		for (unsigned i = 0; i < _size; ++i) data[i] = value;
	}
	
 	template < class Container > StatVec(const Container& c)
	{
	   checkSameLength(c);
	   for (unsigned i = 0; i < this->size(); ++i) data[i] = c[i];
	}
	
	StatVec(const T* value)
	{
		for (unsigned i = 0; i < _size; ++i) { data[i] = value[i]; //std::cout << value[i] << std::endl;
		}
	}

	StatVec& operator = (const StatVec& rhs)
	{
		for (unsigned i = 0; i < _size; ++i)
			data[i] = rhs[i];
		return *this;
	}

	StatVec& operator = (const T& val)
	{
		for (unsigned i = 0; i < _size; ++i) data[i] = val;
		return *this;
	}


	StatVec& operator *= (const T& val)
	{
		for (unsigned i = 0; i < _size; ++i) data[i] *= val;
		return *this;
	}

	StatVec& operator += (const StatVec& rhs)
	{
		for (unsigned i = 0; i < _size; ++i) data[i] += rhs.data[i];
		return *this;
	}
	
	StatVec& operator -= (const StatVec& rhs)
	{
		for (unsigned i = 0; i < _size; ++i) data[i] -= rhs.data[i];
		return *this;
	}
	
	bool operator == (const StatVec& rhs) const
	{
		if (size() != rhs.size()) return false;
		else return std::equal(begin(),end(),rhs.begin());
		
// 		bool b(true);
// 		for (int i = 0; i != size(); ++i) b &= bool(data[i] == rhs[i]);
// 		return b;
	}
	
	bool operator != (const StatVec& rhs) const
	{
		return !( *this == rhs );
	}

// 	virtual ~StatVec() { }

	inline 			  T& operator[](int i)		    {return data[i];}
	inline const 	T& operator[](int i) const 	{return data[i];}

	unsigned capacity() const {return _size;}
	unsigned size()		  const {return _size;}


  T*       begin()        {return data;}
  const T* begin() const  {return data;}

  T*       end()          {return &data[_size];}
  const T* end()   const  {return &data[_size];}

	template < class T2, unsigned _S >
	friend std::ostream & operator << (std::ostream &, const StatVec<T2,_S>& );

};

template < class T, unsigned  _size>
std::ostream& operator << (std::ostream & stream, const StatVec<T,_size> & v)
{
//  for (unsigned i = 0; i < _size; i++) stream << v.data[i] << std::endl;
  for (unsigned i = 0; i < _size; i++) stream << '(' << v.data[i] << ") ";
	return stream;
}

#endif

