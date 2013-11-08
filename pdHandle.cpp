#include "pdHandle.h"

template<class T>
pdHandle<T>::pdHandle(T* p = 0) : ptr(p), use(new size_t(1)) {};

template<class T>
pdHandle<T>::~pdHandle()
{
	destroy();
}

template<class T>
inline T& pdHandle<T>::operator *()
{
	if (ptr)
	{
		return *ptr;
	}
	else
	{
		throw std::runtime_error("pdHandle: dereference of unbound handle");
	}
}

template<class T>
inline T* pdHandle<T>::operator ->()
{
	if (ptr)
	{
		return ptr;
	}
	else
	{
		throw std::runtime_error("pdHandle: access through unbound handle");
	}
}


template<class T>
const T& pdHandle<T>::operator *() const
{
	if (ptr)
	{
		return *ptr;
	}
	else
	{
		throw std::runtime_error("pdHandle: dereference of unbound handle");
	}
}

template<class T>
const T* pdHandle<T>::operator ->() const
{
	if (ptr)
	{
		return ptr;
	}
	else
	{
		throw std::runtime_error("pdHandle: access through unbound handle");
	}
}


//copy constructor
template<class T>
pdHandle<T>::pdHandle(const pdHandle& h): ptr(h.ptr), use(h.use)
{
	++*use;
}

//assignment operator
template<class T>
inline pdHandle<T>& pdHandle<T>::operator =(const pdHandle& rhs)
{
	++*rhs.use;
	destroy();
	ptr = rhs.ptr;
	use = rhs.use;
	return *this;
}

template<class T>
void pdHandle<T>::destroy()
{ 
	if (--*use == 0)
	{
		delete ptr;
		delete use;
	}
}

template<class T>
T* pdHandle<T>::GetPtr() const
{	
	if (ptr)
	{
		return ptr;
	}
	else
	{
		throw std::runtime_error("pdHandle: access through unbound handle");
	}
}
