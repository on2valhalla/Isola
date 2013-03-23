/*
 * Created by Jason Carlisle Mann (jcm2207@columbia.edu)
 * Isola board game playing program
 *
 * Uses the Negamax variant of the minimax search with alpha/beta pruning.
 * Also implements a history heuristic and transposition table.
 *
 */


#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atoi */

#include <iostream>		/* cin/cout */
#include <bitset>		/* bitset */
#include <string>		/* string */



using namespace std;

static const size_t ROWSIZE = 8;
static const size_t BOARDSIZE = ROWSIZE * ROWSIZE;

int usage()
{
	cout << "USAGE: $ isola PLAYER\n" <<
			"Where PLAYER is [1/2] representing "<<
			"the player this instance is"<< endl;
	return 1;
}

int draw_board(bitset<BOARDSIZE> &board, char myPos[2], char opPos[2])
{
	cout << "X represents me, 0 represents you" <<endl;

	for (size_t i = 0; i < board.size() ; i += ROWSIZE)
	{

	}

	return 0;
}






int main(int argc, char *argv[])
{
	/* USAGE
	*
	* $ isola PLAYER
	*
	*	where SIZE is a single int board size (normally should be 8)
	*	and player is [1/2] representing the player this instance is
	*
	*/

	if( argc != 2 )
		return usage();

	char playerNum;
	try
	{ // take in command line argument
		playerNum = argv[1][0] - '0';
	} 
	catch (int e)
	{
		return usage();
	}

	// setup the board: init to zeros, add starting positions
	bitset<BOARDSIZE> initialBoard;
	initialBoard.reset();
	initialBoard.set(0, 1);
	initialBoard.set(BOARDSIZE-1, 1);

	cout << initialBoard.to_string() << endl;


	// completed sucessfully
	return 0;
}










