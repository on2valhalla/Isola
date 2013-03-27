

struct by_min
{
	bool operator() (const Node& lhs, const Node& rhs)
	{
		return lhs.evaluate() < rhs.evaluate();
	}
	
};

struct by_closeness
{
	bool operator() (const Node& lhs, const Node& rhs)
	{
		char Lx = abs(GETX(lhs.myIdx) - GETX(lhs.opIdx));
		char Ly = abs(GETX(lhs.myIdx) - GETX(lhs.opIdx));
		
		char Rx = abs(GETX(rhs.myIdx) - GETX(rhs.opIdx));
		char Ry = abs(GETX(rhs.myIdx) - GETX(rhs.opIdx));
		return (Lx + Ly) < (Rx + Ry);
	}
	
};