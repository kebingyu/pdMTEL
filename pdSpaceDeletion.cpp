#include "pdSpaceDeletion.h"

pdSpaceDeletion::pdSpaceDeletion() {}

pdSpaceDeletion::pdSpaceDeletion(const pdSpaceDeletion &bspace) {}

pdSpaceDeletion::pdSpaceDeletion(int id, int type) : pdSpace(id, type) {}

pdSpaceDeletion::~pdSpaceDeletion() {}

void pdSpaceDeletion::Print(ostream& os) const
{
	pdSpace::Print(os);
}