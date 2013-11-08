/* 

This generic handle class provides pointer-like behavior. 
Access through an unbound Handle is checked and throws a runtime_error exception.
The object to which the pdHandle points is deleted when the last pdHandle goes away.
 Users should allocate new objects of type T and bind them to a pdHandle. Once an object is bound to a pdHandle,, the user must not delete that object.

*/

#ifndef PDHANDLE_H
#define PDHANDLE_H

#include <iostream>
using std::size_t;

template<class T> 
class pdHandle
{
public:
	// Constructor
	pdHandle(T* p);

	// Destructor
	virtual ~pdHandle();

	// Overloaded operators to support pointer behavior
	T& operator*();
	T* operator->();
	const T& operator*() const;
	const T* operator->() const;

	// Copy control: normal pointer behavior, but last pdHandle deletes the object
	pdHandle(const pdHandle&);
	pdHandle& operator=(const pdHandle&);

	T* GetPtr() const;

private:
	// Data
	T* ptr; // shared object
	size_t* use; // count of how many pdHandle point to *ptr

	// Function
	void destroy(); 
};

#endif