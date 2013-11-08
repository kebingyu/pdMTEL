#include "pdSlit.h"

pdSlit::pdSlit() {}

pdSlit::pdSlit(const pdSlit &slit) {}

pdSlit::pdSlit(int id, int type) : pdSpace(id, type) {}

pdSlit::~pdSlit() {}

void pdSlit::SetNorm(int norm)
{
	_norm = norm;
}

void pdSlit::SetHeadID(int i)
{
	_hea_id = i;
}

void pdSlit::SetTailID(int i)
{
	_tail_id = i;
}

void pdSlit::SetExt_head(double e)
{
	_ext_head = e;
}

void pdSlit::SetExt_tail(double e)
{
	_ext_tail = e;
}

int pdSlit::GetHeadID() const
{
	return _hea_id;
}

int pdSlit::GetTailID() const
{
	return _tail_id;
}

int pdSlit::GetNorm() const
{
	return _norm;
}

double pdSlit::GetExt_head() const
{
	return _ext_head;
}

double pdSlit::GetExt_tail() const
{
	return _ext_tail;
}

void pdSlit::Print(ostream& os) const
{
	pdSpace::Print(os);
	os << endl << "  Type: slit, norm= " << _norm 
		<< ", hea_id= " << "\t" << _hea_id << ", tail_id= " << "\t" << _tail_id 
		<< ", hea_extension= " << "\t" << _ext_head << ", tail_extension= " << "\t" << _ext_head << endl;
}