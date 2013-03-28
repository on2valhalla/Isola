

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



	
//	FROM VALID MOVE
//	SAVE THIS DAMN CODE JUST IN CASE FOR THE BLOCKED DIAGONAL
//	if( abs(moveIncr) == 7 || abs(moveIncr) == 9 )
//	{ // diagonal move
//		while (curIdx != *moveIdx)
//		{
//			// check for diagonal blocked path
//			if( board[ curIdx + TOIDX( SIGN( moveVec[0] ), 0) ] &&
//			   board[ curIdx + TOIDX( 0, SIGN( moveVec[1] ))] )
//				return false;
//			
//			// move to the next position
//			curIdx += moveIncr;
//			
//			// check if the position is visited
//			if( board[curIdx] )
//				return false;
//		}
//	}
//	else




char find_moves(const Node &node,
				vector<Node> &children,
				const char *dir,
				const char *player )
{
	// cout << "fastm:" <<endl;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx;
	char curIdx = fromIdx + *dir, n = 0;
	
	while ( node.board[curIdx] )
	{
		
		if( (curIdx % ROWSIZE == 7 &&
			 (*dir == NWEST || *dir == WEST || *dir == SWEST))
		   || (curIdx % ROWSIZE == 0 &&
			   (*dir == NEAST || *dir == EAST || *dir == SEAST)))
			break; //edge of table reached
		
		Node child;
		if(*player == MY)
			child = Node(node.board, &curIdx, &node.opIdx);
		else
			child = Node(node.board, &node.myIdx, &curIdx);
		
//		child.evaluate();
		children.push_back(child);
		n++;
		curIdx += *dir;
	}
	
	return n;
}

char get_children( const Node &node,
				  vector<Node> &children,
				  const char *player )
{
	// cout << "getc:" <<endl;
	// number of children found
	char n = 0;
	
	
	n += find_moves(node, children, &NORTH, player);
	n += find_moves(node, children, &SOUTH, player);
	n += find_moves(node, children, &EAST, player);
	n += find_moves(node, children, &WEST, player);
	n += find_moves(node, children, &NWEST, player);
	n += find_moves(node, children, &NEAST, player);
	n += find_moves(node, children, &SWEST, player);
	n += find_moves(node, children, &SEAST, player);
	
	return n;
}