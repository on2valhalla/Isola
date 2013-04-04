// ISOLA GAME
// ----------

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))


int main()
{
// Establish initial board state:
	bool board[8][8];
	char opPos[2];
	char myPos[2];

	memset(board, true);

// Take input to decide if you are first or second player:
	if( firstPlayer )
		make_move( myPos, opPos );
	else
		take_move( );
}


take_move()
{
	// take keyboard input for opponents move
	char opToPos[2] = std_in; 

	// test opponents move for validity
	if( !valid_move( opPos, opToPos ) )
		fail("Opponent made invalid move");

	make_move()
}

bool valid_move( char curPos[2], char toPos[2] )
{
	// check if given position is off the board or visited
	if ( off_board( toPos ) || board[toPos[0]][toPos[1]] ) )
		return false; // to-position is invalid

	// get vectors representing the direction and entire move
	char moveVec[2] = { toPos[0] - curPos [0], toPos[1] - curPos[1] };
	char moveDir[2] = { (moveVec[0] > 0) - (moveVec[1] < 0), 
						(moveVec[0] > 0) - (moveVec[1] < 0) };

	if (moveVec[0] != 0 && moveVec[1] !=0)
	{ // diagonal move

		// test for a valid diagonal
		if( abs(moveVec[0]) != abs(moveVec[1]) )
			return false;

		for( ; curPos[0] == toPos[0], curPos[1] == toPos[1]; )
		{
			// check for diagonal blocked path
			if( board[curPos[0] + moveDir[0]][curPos[1]] &&
				board[curPos[0]][curPos[1] + moveDir[1]])
				return false;

			// move to the next position
			curPos[0] += moveDir[0]; curPos[1] += moveDir[1];

			// check if the position is visited
			if( board[curPos[0]][curPos[1]] )
				return false;
		}
	}
	else
	{ // horizontal or vertical move
		for( ; curPos[0] == toPos[0], curPos[1] == toPos[1]; )
		{
			// move to the next position
			curPos[0] += moveDir[0]; curPos[1] += moveDir[1];

			// check if the position is visited
			if( board[curPos[0]][curPos[1]] )
				return false;
		}
	}

	// valid move
	return true;
}


// Get Children
get_children( &parentPos, &posHeap)
{
	&newPos, &lastPos;

	moveVectors = { N, NE, NW, S, SE, SW };

	for (each moveVec in moveVectors)
	{
		newPos = parentPos + moveVec;
		lastPos = parentPos;

		while( valid_move(lastPos, newPos) )
		{
			posHeap.add(newPos);
			newPos += moveVec;
		}
	}
}



char[2] make_move(double maxTime)
{
	// find universal cutoff time from current system and max

	Node v = starting board, myPos, opPos;

	for( int d = 0; d < MAX_DEPTH; d++)
	{
		// if the current time is within 3% of cutoff time, break

		v = negamax_a_b(v, alpha, beta, d);
	}

	
	return v;
}















