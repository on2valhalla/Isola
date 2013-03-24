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
#include <assert.h>     /* assert */

#include <iostream>		/* cin/cout */
#include <bitset>		/* bitset */
#include <string>		/* string */
#include <queue>		/* priority queue */


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define TOIDX(x,y) ((x) + ((y) * 8))
#define SIGN(x) (((x) > 0) - ((x) < 0))
#define GETX(idx) ((idx) % 8)
#define GETY(idx) ((idx) / 8)




using namespace std;

static const char ROWSIZE = 8;
static const char BOARDSIZE = ROWSIZE * ROWSIZE;
static const char	NEAST = -7,
					NORTH = -8,
					NWEST = -9,
					WEST = -1,
					EAST = 1,
					SWEST = 7,
					SOUTH = 8,
					SEAST = 9;

typedef bitset<BOARDSIZE> BitBoard;

class Node 
{
public:
	BitBoard board;
	char opIdx;
	char myIdx;
	double heuristic;

	Node(const BitBoard & board, const char *myIdx, const char *opIdx)
	{
		this->board = board;
		this->opIdx = *opIdx;
		this->myIdx = *myIdx;
		this->heuristic = this->evaluate();
	}

	Node( const Node &other )
	{
		this->board = other.board;
		this->opIdx = other.opIdx;
		this->myIdx = other.myIdx;
		this->heuristic = other.heuristic;
	}

	Node()
	{
		opIdx = -1;
		myIdx = -1;
	}

	~Node()
	{	}

	double evaluate()
	{
		return myIdx;
	}

	bool operator() (const Node& lhs, const Node& rhs)
	{
		return lhs.heuristic > rhs.heuristic;
	}

	friend ostream& operator<< (ostream &o, const Node &n)
	{
		o << "myP: " << (int) n.myIdx << "\topP: " 
				<< (int) n.opIdx << "\th: " << n.heuristic;
		return o;
	}


};

int usage()
{
	cout << "USAGE: $ isola PLAYER\n" <<
			"Where PLAYER is [1/2] representing "<<
			"the player this instance is"<< endl;
	return 1;
}

bool fail(string message)
{
	cout << "\n\nERROR!!!!\n" << message << endl;
	return false;
}

void move(BitBoard &board, char *myIdx, char *opIdx, const Node &child)
{
	*myIdx = child.myIdx;
	*opIdx = child.opIdx;
	board.set(*myIdx, 1);
}


bool valid_move( const BitBoard &board, const char *oldIdx, 
					const char *moveIdx )
{
	char curIdx = *oldIdx;

	// check if given position is off the board or visited
	if ( *moveIdx >= BOARDSIZE || *moveIdx < 0 || board[*moveIdx] )
		return false; // to-position is invalid

	// get vectors representing the entire move and direction
	char moveVec[2] = { GETX(*moveIdx) - GETX(curIdx), 
						GETY(*moveIdx) - GETY(curIdx) };
	char moveIncr = TOIDX( SIGN(moveVec[0]), SIGN(moveVec[1]) );


	if( abs(moveIncr) == 7 || abs(moveIncr) == 9 )
	{ // diagonal move
		while (curIdx != *moveIdx)
		{
			// check for diagonal blocked path
			if( board[ curIdx + TOIDX( SIGN( moveVec[0] ), 0) ] &&
				board[ curIdx + TOIDX( 0, SIGN( moveVec[1] ))] )
				return false;

			// move to the next position
			curIdx += moveIncr;

			// check if the position is visited
			if( board[curIdx] )
				return false;
		}
	}
	else if( abs(moveIncr) == 1 || abs(moveIncr) == 8 )
	{ // horizontal or vertical move
		while (curIdx != *moveIdx)
		{
			// move to the next position
			curIdx += moveIncr;

			// check if the position is visited
			if( board[curIdx] )
				return false;
		}
	}
	else
		return false; //invalid move

	// valid move
	return true;
}

char find_moves(const Node &node, 
				priority_queue<Node, vector<Node>, Node> &children,
				const char *dir )
{
	char curIdx = node.myIdx + *dir;
	char n = 0;

	while ( valid_move( node.board, &node.myIdx, &curIdx ) )
	{

		if( (curIdx % ROWSIZE == 7 && 
				(*dir == NWEST || *dir == WEST || *dir == SWEST))
			|| (curIdx % ROWSIZE == 0 && 
				(*dir == NEAST || *dir == EAST || *dir == SEAST)))
			break; //edge of table reached

		BitBoard childBoard (node.board);
		childBoard.set(curIdx, 1);
		Node child (childBoard, &curIdx, &node.opIdx);
		children.push(child);
		cout << "d: " << (int)*dir << "\t" << child << endl;
		n++;
		curIdx += *dir;
	}

	return n;
}

char get_children( const Node &node, 
				priority_queue<Node, vector<Node>, Node> &children )
{
	// number of children found
	char n = 0;

	n += find_moves(node, children, &NORTH);
	n += find_moves(node, children, &SOUTH);
	n += find_moves(node, children, &EAST);
	n += find_moves(node, children, &WEST);
	n += find_moves(node, children, &NWEST);
	n += find_moves(node, children, &NEAST);
	n += find_moves(node, children, &SWEST);
	n += find_moves(node, children, &SEAST);

	return n;
}

// Draws a Isola board represented by a bitset, when supplied
// with the bitset and two current positions as indexes in the bitset
int draw_board(const Node &node)
{
	// cout << "A represents me, B represents you, X a visited space\n\n";

	char i, j, idx;
	cout << "   |";
	for ( i = 1; i <= ROWSIZE; i++ )
		cout << " " << (int)i << " |";
	cout<< endl <<"___";
	for ( i = 0; i < ROWSIZE; i++ )
		cout << "____";

	cout << "_" << endl;

	for ( i = 0; i < ROWSIZE; i++ )
	{
		cout << " " << (int)i+1 << " |";

		for ( j = 0; j < ROWSIZE; j++ )
		{
			// get idxs for the bitset *MACRO* (col, row)
			idx = TOIDX(j, i);

			if( idx == node.myIdx )			cout << "*A*";
			else if( idx == node.opIdx )	cout << "*B*";
			else if( node.board[idx] )		cout << " X ";
			else 						cout << "   ";
			// else if(idx < 10)
			// 	cout << " " << (int)idx << " ";
			// else
			// 	cout << " " << (int)idx;

			cout << "|";
		}

		if (i < ROWSIZE -1)
		{
			cout << endl << "---|";
		
			for ( j = 0; j < ROWSIZE; j++ )
				cout << "———";

			cout << "———————|";
		}

		cout << endl;
	}

	cout << "---|";

	for ( i = 0; i < ROWSIZE; i++ )
		cout << "———";

	cout << "———————|" << endl;

	return 0;
}

bool take_move( BitBoard &board, char *oldIdx )
{
	cout << "take move:" <<endl;
	priority_queue<Node, vector<Node>, Node> children;
	char def = 0;
	char n = get_children(Node(board, oldIdx, &def), children);
	if(n == 0)
		return false;

	string move;
	cout << "Please enter a move (X Y):  ";
	getline(cin, move);

	size_t x = move.find_first_of("12345678");
	size_t y = move.find_first_of("12345678", x+1);
	if (x == string::npos || y == string::npos)
		return fail( "Incorrect move format. " + move);


	char newIdx = TOIDX(move[x] - '0' - 1, move[y] - '0' - 1);

	if( !valid_move(board, oldIdx, &newIdx) )
		return fail( "Invalid move was attempted.");

	*oldIdx = newIdx;
	board.set(newIdx, 1);
	return true;
}


bool play( BitBoard &board, char *myIdx, char *opIdx )
{
	Node curNode (board, myIdx, opIdx);
	priority_queue< Node, vector<Node>, Node> children;


	// cout << endl << (int)size << endl;
	// while ( !children.empty() )
	// {
	// 	cout  << "\t" << children.top() <<endl;
	// 	children.pop();
	// }
	// cout << endl;

	while( get_children(curNode, children) > 0)
	{
		curNode = children.top();
		while( !children.empty() )
			children.pop();

		move(curNode.board, &curNode.myIdx, &curNode.opIdx, curNode);
		draw_board(curNode);

		cout.flush();
		if( !take_move(curNode.board, &curNode.opIdx))
			return true;

		draw_board(curNode);
	}

	return false;
}




int main(int argc, char *argv[])
{
	/* USAGE
	*
	* $ isola PLAYER
	*
	*	where SIZE is a single int board size (normally should be 8)
	*	and player is [1 or 2] representing the player this instance is
	*
	*/

	char playerNum;
	try
	{ // take in command line argument
		assert(argc == 2);
		playerNum = argv[1][0] - '0';
	} 
	catch (int e)
	{
		return usage();
	}

	// setup the board: init to zeros, add starting positions
	BitBoard board;
	board.reset();

	char myIdx, opIdx;
	if(playerNum == 1)
	{
		myIdx = TOIDX(0, 0);
		opIdx = TOIDX(ROWSIZE-1, ROWSIZE-1);
	}
	else
	{
		myIdx = TOIDX(0, 0);
		opIdx = TOIDX(ROWSIZE-1, ROWSIZE-1);
	}

	board.set(myIdx, 1);
	board.set(opIdx, 1);

	Node initNode (board, &myIdx, &opIdx);

	draw_board(initNode);


	if(playerNum == 2)
	{
		take_move(board, &opIdx);
	}

	bool win = play(board, &myIdx, &opIdx);

	if ( win )
		cout << "I win" << endl;
	else
		cout << "You win" << endl;


	//---------------------------------------------
	// TESTING CODE
	//---------------------------------------------
	// cout << "\n\nValid Move Test:" << endl;

	// myIdx = 29;
	// board.set(myIdx, 1);
	// char valid[] = {2, 5,
	// 				11, 13, 15,
	// 				20, 21, 22,
	// 				26, 28, 30, 31,
	// 				36, 38,
	// 				43, 45, 47,
	// 				50, 53,
	// 				57, 61 };
	// char invalid[] = {1, 3, 4, 12, 23, 29, 35, 
	// 				42, 44, 46, 51, 52, 54, 55,
	// 				58, 59, 60, 62, 63, 64, 65, 67, -1};
	// draw_board(board, &myIdx, &opIdx);

	// for(size_t i = 0; i < sizeof(valid); i++)
	// 	cout << (int) (valid[i]) << ": " 
	// 		<< valid_move(board, &myIdx, &valid[i]) << endl;



	// cout << "\n\ninvalid Move Test:" << endl;
	// for(size_t i = 0; i < sizeof(invalid); i++)
	// 	cout << (int) (invalid[i]) << ": " 
	// 		<< valid_move(board, &myIdx, &invalid[i]) << endl;



	// cout << "\n\nblocked Move Test:" << endl;
	// board.set(28, 1);
	// board.set(21, 1);
	// board.set(37, 1);
	// board.set(30, 1);
	// for(size_t i = 0; i < sizeof(valid); i++)
	// 	cout << (int) (valid[i]) << ": " 
	// 		<< valid_move(board, &myIdx, &valid[i]) << endl;

	//---------------------------------------------


	// completed sucessfully
	return 0;
}










