#include "pdBond.h"
#include "pdNode.h"

//// initial static member
//new_handler pdBond::currentHandler;

pdBond::pdBond() {}

pdBond::pdBond(const pdBond& bond) {}

//pdBond::pdBond(int id, pdNode* i, pdNode* j) : _id(id), _node_i(i), _node_j(j) {}

pdBond::pdBond(pdNode* i, pdNode* j) : _node_i(i), _node_j(j) {}

pdBond::~pdBond() 
{
	//_planes.clear();
}

//int pdBond::GetID() const
//{
//	return _id;
//}

pdNode* pdBond::GetNodeI() const
{
	return _node_i;
}

pdNode* pdBond::GetNodeJ() const
{
	return _node_j;
}

void pdBond::SetBreak(bool b)
{
	_break = b;
}

bool pdBond::IsBreak() const
{
	return _break;
}

//void pdBond::SetStretch(double str)
//{
//	_stretch = str;
//}

//double pdBond::GetStretch() const
//{
//	return _stretch;
//}

void pdBond::SetMinDistance(double d)
{
	_minDistance = d;
}

void pdBond::SetMaxDistance(double d)
{
	_maxDistance = d;
}

double pdBond::GetMinDistance() const
{
	return _minDistance;
}

double pdBond::GetMaxDistance() const
{
	return _maxDistance;
}

//void pdBond::AddPlane(pdSpacePlane* p)
//{
//	_planes.push_back(p);
//}
//
//int pdBond::GetNumberPlane()
//{
//	 // return the total number of planes this bond passes through
//	return int(_planes.size());
//}
//
//pdSpacePlane* pdBond::GetPlane(int i)
//{
//	return _planes[i];
//}

void pdBond::Print(ostream& os) const
{
	//os << "Bond (id = " << _id << "):" 
	os << "Bond:" 
		<< "\t" << "node_i_i = " << _node_i->GetID() << "\t" << "node_i_j = " << _node_j->GetID()
		//<< "\t" << "stretch = " << _stretch 
		<< "\t" << "break=" <<  _break;
}

ostream& operator << (ostream& os, const pdBond& b)
{
	b.Print(os);
	return os;
}

//new_handler pdBond::set_new_handler(new_handler p)
//{
//	new_handler oldHandler = currentHandler;
//	currentHandler = p;
//	return oldHandler;
//}
//
//void* pdBond::operator new(size_t size)
//{
//	// return a legitimate pointer when 0 bytes are requested
//	if (size == 0)
//	{
//		size = 1;
//	}
//
//	// if operator new is not defined in the derived class and it calls operator new in the base class
//	// let standard opertor new handle the request
//	if (size != sizeof(pdBond))
//	{
//		return ::operator new(size);
//	}
//
//	while (1)
//	{
//		// attempt to allocate size bytes
//
//		if 
//		{
//			// the allocation is successful
//			return ??;
//		}
//
//		// allocation was unsuccessful
//		// find out the current error-handling function
//		new_handler globalHandler = set_new_handler(0);
//		set_new_handler(globalHandler);
//
//		if (globalHandler)
//		{
//			(*globalHandler)();
//		}
//		else
//		{
//			throw std::bad_alloc;
//		}
//
//		void* memory;
//
//		try
//		{
//			memory = ::operator new(size);
//		}
//		catch (std::bad_alloc&)
//		{
//			std::set_new_handler(globalHandler);
//			throw;
//		}
//
//		std::set_new_handler(globalHandler);
//
//		return memory;
//	}
//}
//
//void pdBond::operator delete(void* rawMemory, size_t size)
//{
//	if (rawMemory == 0)
//	{
//		return;
//	}
//
//	if (size_t != sizeof(pdBond))
//	{
//		::operator delete(rawMemory);
//		return;
//	}
//
//	//deallocate the memory pointed to by rawMemeory
//
//	
//}
//
//void pdBond::outOfMemory()
//{
//	cerr << "pdBond:: unable to allocate requested memery" << endl;
//	abort();
//}