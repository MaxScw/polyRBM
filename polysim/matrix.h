#ifndef MATRIX_H
#define MATRIX_H

// #include <list>
#include <iostream>
#include "pair.h"
#include "tools.h"
#include "static_vector.h"
// #include <radixsort.cuh>

using namespace std;

// typedef uint64_t size_t;

// namespace matrix
// {
// typedef pair < size_t, size_t > index_pair;
// typedef pair < size_t, size_t > index_pair;
// template <
// template < class LinearStorage > FastMatrix
// {
// private:
// 	
// 	typedef LinearStorage::value_type T;
// 	
// 	LinearStorage data;
// 	
// public:
// 		
// 	typedef T value_type;
// 	

namespace new_matrix
{
  
  
template < class Indexer, class Storage > class Map
{
private:
	typedef typename Storage::value_type value_type;
	typedef typename Indexer::key_type key_type;

public:
  
	Map(const key_type& size, const value_type& val = 0)
	:idx(size)
	,data(idx.size(),val)
	{}
  
	const value_type& operator[] (const key_type& key) const 
	{
		return data[idx.getIndex(key)];
	}
	
	value_type& operator[] (const key_type& key) 
	{
		return const_cast< value_type& > ( this->operator[](key) );
	}
	
private: 
  
	Storage data;
	Indexer idx;
};

template < uint stop, uint start=0, uint step=1 > struct recursive_generic_loop
{
  template < class Func > inline static void exec(Func & f) {
	f(start); recursive_generic_loop <stop,start+step,step> ::exec(f);
  }
};

template < uint stop, uint step > struct recursive_generic_loop <stop,stop,step>
{
  template < class Func > inline static void exec(Func & f) {}
};

// template < uint stop, uint start=0, uint step=1 > struct generic_for_each
// {
//   for (uint i = 0
//   template < class Func > static void exec(Func & f) {
// 	f(start); recursive_generic_loop <stop,start+step,step> ::exec(f);
//   }
// };

// template < uint stop, uint start=0, uint step=1 > struct generic_for_each
// {
//   template < class Func > static void exec(Func & f) {
// 	f(start); generic_loop <stop,start+step,step> ::exec(f);
//   }
// };



// template < uint stop, uint step > struct generic_loop <stop,stop,step>
// {
//   template < class Func > static void exec(Func & f) {}
// };


// template < class Func, uint stop, uint start=0, uint step=1 > struct generic_for_each
// {
//   template < class T >
//   static void exec(T& data) const{
// 	Func f; (data[start]);
// 	generic_for_each<Func,stop,start+step,step>::operator()(data);
//   }
// };
// 
// template < class Func, uint stop, uint step > struct generic_for_each<Func,stop,stop,step>
// {
//   template < class T > void operator()(T& data) const{}
// };

// template < uint dim >
// struct CalcShifts
// {
//   uint total;
// //   const StatVec
//   
// }

/// calc idx ////////////////////////////////////////////////////////////////// 

template < uint dim, uint i > struct calc_grid_idx
{
  template < class T > inline size_t operator () (const T& shifts, const T& key) const
  {
	return (size_t(key[i])<<shifts[i]) + calc_grid_idx<dim,i+1>()(shifts,key);
  }
};

template < uint dim > struct calc_grid_idx < dim, dim >
{
  template < class T > inline size_t operator () (const T& shifts, const T& key) const
  {
	return 0;
  }
};

template < uint dim > struct calc_grid_idx < dim, 0 >
{
  template < class T > inline size_t operator () (const T& shifts, const T& key) const
  {
	return size_t(key[0]) + calc_grid_idx<dim,1>()(shifts,key);
  }
};

/// calc shifts ///////////////////////////////////////////////////////////////

template < uint max, uint i > struct calc_shifts
{
  template < class T > void operator () (T& shifts, const T& key) const
  {
	shifts[i] = shifts[i-1] + my_log2(key[i-1]);
	calc_shifts<max,i+1>()(shifts,key);
  }
};

template < uint max > struct calc_shifts < max, max >
{
  template < class T > void operator () (T& shifts, const T& key) const
  {}
};

template < uint max > struct calc_shifts < max, 0 >
{
  template < class T > void operator () (T& shifts, const T& key) const
  {
	shifts[0]=0;
	calc_shifts<max,1>()(shifts,key);
  }
};

/// grid_idx //////////////////////////////////////////////////////////////////

template < uint dim, bool periodic = false > class GridIndex;

template < uint D > class GridIndex < D, false >
{
public:
   typedef StatVec<uint,D> key_type;
   typedef key_type KT;
   
   enum{dim=D};
   
   GridIndex(const KT& max)
   {
 	 calc_shifts<D,0>()(shifts,max);
	 cout << "huhu" << shifts << endl;
   }
   
   size_t getIndex(const KT & key) const
   {
	 return calc_grid_idx<D,0>()(shifts,key);
   }
   
private:
  
   key_type shifts;

};
// template < class T, uint dim > class FastGrid
// {
//   
//   T* 	data;
//   Index shifts;
//   
//   
//   public:
//   
//   typedef StatVec < uint, dim > Index;
//   
//   FastGrid(const Index& size):
//   shifts(0)
//   {
// 	uint shift = 0;
// 	for(uint axis = 0; axis < dim; ++axis)
// 	{
// 	   shifts[axis] = shift;
// 	   shift+=log2(size[axis]);// <
// 	}
// 	
// 	
//   }
// };

class FastMatrixShape
{
	public:
		FastMatrixShape(const size_t& _d1, const size_t& _d2)
		:d1(_d1)
		,d2(_d2)
		,d2_shift(my_log2(d1-1l)+1l)
// 		,d2_total(size_t( 1 << d2_shift ) - 1 );
		{
// 					cout << "dim1 = " << dim1() << " t_dim2() = " << dim2_total();
		}
// 		my_log2(d0-1)+1

		size_t LinearIndex(const size_t& idx1, const size_t& idx2) const
		{
			return idx1 + (idx2 << d2_shift);
		}
		
		size_t size() const {return d1*d2;};
	
		const size_t& dim1() const { return d1; }
		const size_t& dim2() const { return d2; }
		
		virtual ~FastMatrixShape() {};
		
	protected:
		
		size_t 	total_size() const { return dim2_total()*dim1(); }
		
	private:
		
		size_t  dim2_total() const { return size_t( ( 1l << d2_shift )  ); }		
		
		size_t d1;
		size_t d2;
		
// 		size_t d2_total;
		size_t d2_shift;
		
};

///*****************************************************************************

template < class LinearStorage , class Shape = FastMatrixShape >
class SmartMatrix : public Shape
{
private: 	typedef typename LinearStorage::value_type T;
public:		typedef T	value_type;

	class ReturnHandle
	{
	private:
		SmartMatrix& m;
		const size_t& idx1;
		
		public:
		ReturnHandle(SmartMatrix& _m, const size_t& i):m(_m),idx1(i){}
			
		typename LinearStorage::iterator::reference operator[] (const size_t& idx2)
		{
		 return m(idx1,idx2);
		}
	};
	
	SmartMatrix(const size_t& _dim1, const size_t& _dim2, const T& val = T())
	:Shape(_dim1,_dim2)
	,data(Shape::total_size(),val)
	{
// 		*this = val;
	}
	
	virtual ~SmartMatrix() {}
	
	inline ReturnHandle operator [] (const size_t& i) 
	{
		return 	ReturnHandle(*this,i);
	}

	template < class T2 >
	inline typename LinearStorage::iterator::reference operator()	(const T2& x, const T2& y)
	{
		return *(data.begin() + Shape::LinearIndex(x,y));
	}
	
	void inline Set(const size_t & idx1, const size_t & idx2, const T& val)
	{
		*(data.begin() + Shape::LinearIndex(idx1,idx2)) = val;
// 		data.Set(Shape::LinearIndex(idx1,idx2),val);
	}

// 	template < class T2 >
// 	inline const T& operator []	(const pair < T2, T2 >& _c)	const
// 	{
// 		return *(data 	+ _c.second + _c.first*dim2());
// 	}

// 	template < class T2 >
// 	inline const T& operator []	(const Pair < T2 >& _c) const
// 	{
// 		return *(data 	+ _c.second + _c.first*dim2());
// 	}
// 
// 	inline const size_t&				dim1() 		const 								{return dimension.first;}
// 	inline const size_t&				dim2() 		const 								{return dimension.second;}
// 	inline const index_pair& 	dim() 		const 									{return dimension;}
// 	inline 			 size_t&				dim1()													{return dimension.first;}
// 	inline			 size_t&				dim2()													{return dimension.second;}
// 	inline			 index_pair& 	dim()															{return dimension;}
// 
// 	inline const size_t& size() const {return total;}
// 	
// 	void resize(uint nw, uint nh, const T& val = T())
// 	{
// 		delete[] data;
// 		dimension = index_pair(nw,nh);
// 		total			= dimension.first*dimension.second;
// 		data 			= new T[total];
// 		*this 		= val;
// 	}
// 
// 	template < class TR >
// 	SmartMatrix& operator = (const SmartMatrix < TR > & rhs)
// 	{
// 		if (size() >= rhs.size())
// 		{
// 			T* 	end 		= data + total;
// 			TR* rhs_ptr = rhs.data;
// 			for (T* lhs_ptr = data ; lhs_ptr != end; ++lhs_ptr)
// 			{
// 				*lhs_ptr = *rhs_ptr; ++rhs_ptr;
// 			}
// 			dimension = rhs.dim();
// 		}
// 	}
// 
// 	template < class T >
	SmartMatrix& operator = (const T& rhs)
	{
// 		T* end = data + total;
  	typename LinearStorage::iterator it(data.begin()),e(data.end());
		
		while(it != e ) { *it = rhs; it+=1;}
		
// 		for (T* lhs_ptr = data; lhs_ptr != end; ++lhs_ptr) *lhs_ptr = rhs;
// 			for (uint i = 0; i < this->dim1(); ++i)
// 				for (uint j = 0; j < this->dim2(); ++j) Set(i,j,rhs);
		return *this;
	}

private:
	LinearStorage data;
};

} //namespace new_matrix

///*****************************************************************************


template < class T > class MatrixDynamic
{
public:

	typedef T value_type;

	MatrixDynamic(const index_pair& _dim = index_pair(1,1), const T& val = T())
	:dimension(_dim)
	,total(dimension.first*dimension.second)
	,data(new T[total])
	{*this = val;}

	MatrixDynamic(const size_t& _dim1, const size_t& _dim2, const T& val = T())
	:dimension(_dim1,_dim2)
	,total(dimension.first*dimension.second)
	,data(new T[total])
	{
		*this = val;
	}
	
	template < class Iter >
	MatrixDynamic(const size_t& _dim1, const size_t& _dim2, Iter vi, Iter ve ) 
	:dimension(_dim1,_dim2)
	,total(dimension.first*dimension.second)
	,data(new T[total])
	{
		uint i = 0;
		while (vi != ve) {data[i] = *vi; ++vi; ++i;}
	}
	

	MatrixDynamic(const MatrixDynamic& rhs)
	:dimension(rhs.dimension)
	,total(rhs.total)
	,data(new T[total])
	{
		*this = rhs;
	}

	virtual ~MatrixDynamic()
	{
		delete[] data;
	}


// 	template < class T2 >
	inline T* operator [] (const uint& i)
	{
		return 	(data + i*dim2());
	}

	template < class T2 >
	inline T& operator()	(const T2& x, const T2& y)
	{
		return *(data + y + x*dim2());
	}

// 	template < class T2 >
// 	inline T& operator []	(const pair < T2, T2 >& _c)
// 	{
// 		return *(data 	+ _c.second + _c.first*dim2());
// 	}

	template < class T2 >
	inline T& operator []	(const Pair < T2 >& _c)
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

// 	template < class T2 >
	inline const T* operator [] (const uint& i)								const
	{
		return 	(data + i*dim2());
	}

	template < class T2 >
	inline const T& operator()	(const T2& x, const T2& y)	const
	{
		return *(data + y + x*dim2());
	}

// 	template < class T2 >
// 	inline const T& operator []	(const pair < T2, T2 >& _c)	const
// 	{
// 		return *(data 	+ _c.second + _c.first*dim2());
// 	}

	template < class T2 >
	inline const T& operator []	(const Pair < T2 >& _c)			const
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

	inline const size_t&				dim1() 		const 								{return dimension.first;}
	inline const size_t&				dim2() 		const 								{return dimension.second;}
	inline const index_pair& 	dim() 		const 									{return dimension;}
	inline 			 size_t&				dim1()													{return dimension.first;}
	inline			 size_t&				dim2()													{return dimension.second;}
	inline			 index_pair& 	dim()															{return dimension;}

	inline const size_t& size() const {return total;}
	
	void resize(uint nw, uint nh, const T& val = T())
	{
		delete[] data;
		dimension = index_pair(nw,nh);
		total			= dimension.first*dimension.second;
		data 			= new T[total];
		*this 		= val;
	}

	template < class TR >
	MatrixDynamic& operator = (const MatrixDynamic < TR > & rhs)
	{
		if (size() >= rhs.size())
		{
			T* 	end 		= data + total;
			TR* rhs_ptr = rhs.data;
			for (T* lhs_ptr = data ; lhs_ptr != end; ++lhs_ptr)
			{
				*lhs_ptr = *rhs_ptr; ++rhs_ptr;
			}
			dimension = rhs.dim();
		}
	}

	template < class TR >
	MatrixDynamic& operator = (const TR& rhs)
	{
		T* end = data + total;
		for (T* lhs_ptr = data; lhs_ptr != end; ++lhs_ptr) *lhs_ptr = rhs;
	}

protected:

	index_pair 	dimension;
	size_t			total;
	T* 					data;

};


///*****************************************************************************

template < class T, size_t X, size_t Y > class MatrixStatic
{
public:
	typedef T value_type;

	MatrixStatic(const T& val = T(0))
// 	:data(new T[X*Y])
	{
		*this = val;
	}

	template < class TR, size_t RX, size_t RY >
	MatrixStatic(const MatrixStatic<T,RX,RY>& rhs)
// 	:data(new T[X*Y])
	{
		*this = rhs;
// 		if (size() >= rhs.size())
// 		{
// 			TR* rhs_ptr 	= rhs.data;
// 			for (size_t i = 0;  i != X*Y; ++i)	{ data[i] = *rhs_ptr; ++rhs_ptr; }
// 		}
	}

	~MatrixStatic()
	{
// 		delete[] data;
	}

	template < class TR, size_t RX, size_t RY >
	MatrixStatic& operator = (const MatrixStatic < TR, RX, RY > & rhs)
	{
		if (size() >= rhs.size())
		{
// 			TR* rhs_ptr 	= rhs.data;
			for (size_t i = 0;  i != X*Y; ++i)	{ data[i] = rhs[i]; }
		}
	}

	MatrixStatic& operator = (const T& rhs)
	{
		for (size_t i = 0;  i!= X*Y; ++i)	data[i] = rhs;
        return *this;
	}

	template < class T2 >
	inline T& operator [] (const T2& i)
	{
		return 	*(data + i);
	}

	template < class T2 >
	inline T& operator()	(const T2& x, const T2& y)
	{
		return *(data + y + x*dim2());
	}

	template < class T2 >
	inline T& operator []	(const pair < T2, T2 >& _c)
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

	template < class T2 >
	inline T& operator []	(const Pair < T2 >& _c)
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

	template < class T2 >
	inline const T& operator [] (const T2& i)	const
	{
		return 	*(data + i);
	}

	template < class T2 >
	inline const T& operator()	(const T2& x, const T2& y)	const
	{
// 		cout << "returning " << *(data + y + x*dim2()) << endl;
		return *(data + y + x*dim2());
	}

	template < class T2 >
	inline const T& operator []	(const pair < T2, T2 >& _c)	const
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

	template < class T2 >
	inline const T& operator []	(const Pair < T2 >& _c)			const
	{
		return *(data 	+ _c.second + _c.first*dim2());
	}

	inline 			 size_t				dim1()	const {return X;}
	inline			 size_t				dim2()	const {return Y;}
	inline			 index_pair		dim()		const {return index_pair (X,Y);}

	inline 			 size_t				size() 	const {return X*Y;}

// 	T* begin () {return data;}
// 	T* end () {return data + X*Y;}

private:
// 	T* 					data;
		T data[X*Y];
};

template < class M >
inline bool limits(const int& x,const int& y, const M& m)
{
	if ( (x >= 0) && ( x < m.dim1() ) && ( y >= 0 ) && ( x < m.dim2() ))
		return true;
	else
		return false;
}

template < class M , class T >
inline bool limits(const pair < T , T > & ip, const M& m)
{
	return limits (ip.first, ip.second, m);
}

template < class T > ostream& operator <<
(ostream& str, const MatrixDynamic <T>& m)
{
	cout << endl;
	for (int i = 0; i != m.dim1(); ++i)
	{
		for (int j = 0; j != m.dim2(); ++j)
			cout << ' ' << m(i,j);
		cout << endl;
	}
	return str;
}

template < class T, size_t X, size_t Y > ostream& operator <<
(ostream& str, const MatrixStatic <T,X,Y> & m)
{
	cout << endl;
	for (int i = 0; i != m.dim1(); ++i)
	{
		for (int j = 0; j != m.dim2(); ++j)
			cout << setw(13) << right << m(i,j);
		cout << endl;
	}
	return str;
}


template < class Storage, class Shape > ostream& operator <<
(ostream& str,new_matrix::SmartMatrix <Storage,Shape> & m)
{
	cout << endl;
	for (int i = 0; i != m.dim1(); ++i)
	{
		for (int j = 0; j != m.dim2(); ++j)
		{
			cout << setw(2) << right ;
			if (m(i,j)) cout << m(i,j);
			else cout << ' ';
		}
		cout << endl;
	}
	return str;
}


// template < class T > matrix_to_ppm
// {
// 	
// 	
// 	
// 	
// };





// //
// }

#endif
