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
#include <sys/time.h>

#include <iostream>		/* cin/cout */
#include <bitset>		/* bitset */
#include <string>		/* string */
#include <queue>		/* priority queue */
#include <limits>		/* min and max for types */
#include <ctime>		/* for system time */
#include <sstream>		/* stringstream */


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define TOIDX(x,y) ((x) + ((y) * 8))
#define SIGN(x) (((x) > 0) - ((x) < 0))
#define GETX(idx) ((idx) % 8 + 1)
#define GETY(idx) ((idx) / 8 + 1)




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
static const char MY = 1;
static const char OP = -1;
//Keep Even?
static const char MAXDEPTH = 7;
static const char MAXSECS = 58;
static const int MININT = numeric_limits<int>::min() + 10;
static const int MAXINT = numeric_limits<int>::max() - 10;



typedef bitset<BOARDSIZE> BitBoard;

class Node;
struct by_max;
struct by_min;

int usage();
bool fail(string);
int draw_board(const Node&);


bool valid_move( const BitBoard&, const char*,
				const char*);


char get_children( const Node&,	vector<Node>&, const char* player);
char find_moves(const Node&, vector<Node>&, const char*, const char* );


char fast_moves(const Node&, const char*, const char*);
char fast_children(const Node&, const char*);


bool take_move( BitBoard&, char* );
bool play( BitBoard&, char*, char* );


Node search_root(Node &, int&);
int alpha_beta(Node&, int, int, const char*, const timeval&, int);

bool game_over(Node&);
bool lookup(const BitBoard &board, int *hashValue);
bool store(const BitBoard &board);


string display_time( const timeval &tv );
timeval split_time( const timeval &tv, char split);
timeval add_time( const timeval &tv, const char *seconds);
void add_time( timeval &to, const timeval &add);
bool past_time(const timeval &cur, const timeval &end);
timeval diff(const timeval &start, const timeval &end);


class Node
{
public:
	BitBoard board;
	char opIdx;
	char myIdx;
	int heuristic;
	bool eval;
	
	Node(const BitBoard & board, const char *myIdx, const char *opIdx)
	{
		this->board = board;
		this->board.set(*myIdx, 1);
		this->board.set(*opIdx, 1);
		this->opIdx = *opIdx;
		this->myIdx = *myIdx;
		this->eval = false;
	}
	
	Node( const Node &other )
	{
		// cout << "node copied" << endl;
		this->board = other.board;
		this->opIdx = other.opIdx;
		this->myIdx = other.myIdx;
		this->heuristic = other.heuristic;
		this->eval = other.eval;
	}
	
	Node()
	{
		this->eval = false;
	}
	
	~Node()
	{	}
	
	int evaluate()
	{
		if(!eval)
		{
			int opMoves = fast_children(*this, &OP);
			int myMoves = fast_children(*this, &MY);
			heuristic = (opMoves == 0) ? numeric_limits<int>::max() - 10:
							(myMoves == 0) ? numeric_limits<int>::min() + 10:
							(myMoves - opMoves*2 + 50) * (int)board.count();
			eval = true;
		}
		return heuristic;
	}
	
	int evaluate() const
	{
		if (!eval)
		{
			cout << "Node not evaluated: " << *this << endl;
			throw "Node not evaluated";
		}
		return heuristic;
	}
	
	friend ostream& operator<< (ostream &o, const Node &n)
	{
		o << "myP: " << (int) n.myIdx << "\topP: "
		<< (int) n.opIdx << "\th: " << n.heuristic;
		return o;
	}
};







struct by_max
{
	bool operator() (const Node& lhs, const Node& rhs)
	{
		return lhs.evaluate() > rhs.evaluate();
	}
	
};

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


bool valid_move( const BitBoard &board, const char *oldIdx,
				const char *moveIdx )
{
	// cout << "valid:" <<endl;
	char curIdx = *oldIdx;
	
	// check if given position is off the board or visited
	if ( *moveIdx >= BOARDSIZE || *moveIdx < 0 || board[*moveIdx] )
		return false; // to-position is invalid
	
	// get vectors representing the entire move and direction
	char x = GETX(*moveIdx) - GETX(curIdx);
	char y = GETY(*moveIdx) - GETY(curIdx);
	char moveVec[2] = { x, y };
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

bool game_over(Node &node)
{
	if(node.evaluate() >= numeric_limits<int>::max()-20
	   || node.evaluate() <= numeric_limits<int>::min() + 20)
		return true;
	else
		return false;
}

// Draws a Isola board represented by a bitset, when supplied
// with the bitset and two current positions as indexes in the bitset
int draw_board(const Node &node)
{
	// cout << "A represents me, B represents you, X a visited space\n\n";
	// cout << "draw:" <<endl;
	
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




char find_moves(const Node &node,
				vector<Node> &children,
				const char *dir,
				const char *player )
{
	// cout << "fastm:" <<endl;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx;
	char curIdx = fromIdx + *dir, n = 0;
	
	while ( valid_move( node.board, &fromIdx, &curIdx ) )
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
		
		child.evaluate();
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


char fast_moves(const Node &node, const char *dir, const char *player)
{
	// cout << "fastm:" <<endl;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx;
	char curIdx = fromIdx + *dir, n = 0;
	
	while ( valid_move( node.board, &fromIdx, &curIdx ) )
	{
		
		if( (curIdx % ROWSIZE == 7 &&
			 (*dir == NWEST || *dir == WEST || *dir == SWEST))
		   || (curIdx % ROWSIZE == 0 &&
			   (*dir == NEAST || *dir == EAST || *dir == SEAST)))
			break; //edge of table reached
		
		n++;
		curIdx += *dir;
	}
	
	return n;
}

char fast_children(const Node &node, const char *player)
{
	// cout << "fastc:" <<endl;
	// number of children found
	char n = 0;
	
	n += fast_moves(node, &NORTH, player);
	n += fast_moves(node, &SOUTH, player);
	n += fast_moves(node, &EAST, player);
	n += fast_moves(node, &WEST, player);
	n += fast_moves(node, &NWEST, player);
	n += fast_moves(node, &NEAST, player);
	n += fast_moves(node, &SWEST, player);
	n += fast_moves(node, &SEAST, player);
	
	return n;
}








bool take_move( BitBoard &board, char *oldIdx )
{
	// cout << "take:" <<endl;
	char def = 0, newIdx = 0,
	n = fast_children(Node(board, oldIdx, &def), &MY);
	if(n == 0)
		return false;
	
	cout << "numKids Op: " << (int) n << endl;
	while( true )
	{
		string move;
		cout << "Please enter a move (X Y):  ";
		getline(cin, move);
		
		size_t x = move.find_first_of("12345678");
		size_t y = move.find_first_of("12345678", x+1);
		if (x == string::npos || y == string::npos)
		{
			cout << "Incorrect move format, try again." <<endl;
			continue;
		}
		
		newIdx = TOIDX(move[x] - '0' - 1, move[y] - '0' - 1);
		
		if( !valid_move(board, oldIdx, &newIdx) )
		{
			cout << "Invalid move, try again." <<endl;
			continue;
		}
		
		break;
	}
	
	// success
	*oldIdx = newIdx;
	board.set(newIdx, 1);
	return true;
}

bool play( const Node &initNode )
{
	// cout << "play:" <<endl;
	Node curNode = initNode;
	curNode.evaluate();
	int alpha = 251;
	
	while( fast_children(curNode, &MY) > 0)
	{
		// cout << "playloop:" <<endl;
		// cout <<"before "<< curNode <<endl;
		curNode = search_root(curNode, alpha);
		// cerr <<"\n\n\n\n\n--------------------------------------"<< curNode <<endl;
		
		cout << "\n\n";
		draw_board(curNode);
		
		if( !take_move(curNode.board, &curNode.opIdx))
			return true;
		
		cout << "\n\n";
		draw_board(curNode);
	}
	
	return false;
}




timeval diff(const timeval &start, const timeval &end)
{
	timeval temp;
	if ((end.tv_usec-start.tv_usec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_usec = 1000000+end.tv_usec-start.tv_usec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_usec = end.tv_usec-start.tv_usec;
	}
	return temp;
}

bool past_time(const timeval &cur, const timeval &end)
{
	timeval tmp = diff(cur, end);
	long t = tmp.tv_sec*1000000 + tmp.tv_usec;
	// cerr << t << "\tcur: " << display_time(cur) << "\tend: "
	// 				<< display_time(end) <<"\ttmp: " << display_time(tmp) <<endl;
	return t <= 0;
}

timeval add_time( const timeval &tv, const char *seconds)
{
	timeval temp;
	temp.tv_sec = tv.tv_sec + *seconds;
	temp.tv_usec = tv.tv_usec;
	return temp;
}

void add_time( timeval &to, const timeval &add)
{
	to.tv_sec += add.tv_sec;
	to.tv_usec += add.tv_usec;
	if( to.tv_usec >= 1000000)
	{
		to.tv_usec -= 1000000;
		to.tv_sec++;
	}
}

timeval split_time( const timeval &tv, char split)
{
	timeval tmp;
	long usecs = ((tv.tv_sec*1000000) + tv.tv_usec) / split;
	tmp.tv_sec = usecs / 1000000;
	tmp.tv_usec = usecs % 1000000;
	return tmp;
}

string display_time( const timeval &tv )
{
	std::string number;
	std::stringstream strstream;
	strstream << tv.tv_sec << "." << tv.tv_usec << " seconds";
	strstream >> number;
	return number;
}




Node search_root(Node &initNode, int &alpha)
{
	// set time limit
	timeval begin, end, tv;
	gettimeofday(&begin, NULL);
	end = add_time(begin, &MAXSECS);
	
	
	vector<Node> children;
	int numKids = get_children(initNode, children, &MY);
	sort(children.begin(), children.end(), by_max());
	// random_shuffle(children.begin(),children.end());
	
	int i = 0, d = 0, value, beta = MAXINT, bestIdx = 0, lastIdx = 0;
	for(; ; d++)
	{
		try{
			for (; i < numKids; i++)
			{
				// gettimeofday(&tv, NULL);
				// tmp = split_time( diff(tv, end), numKids - i);
				// add_time(tv, tmp );
				// cerr << "split: " << display_time( tmp ) << endl;
				value = -1 * alpha_beta(children[i],
										-1*beta,
										-1*alpha,
										&OP,
										end,
										d);
				
				cerr << "last: " << value
				<< "  pos: " << GETX((int)children[i].myIdx) << ","
				<< GETY((int)children[i].myIdx);
				
				if(value == MAXINT)
					return children[i];
				
				if (value > alpha)
				{
					alpha = value;
					bestIdx = i;
				}
				
				cerr <<"\t\tbest: " << alpha
				<< "  pos: " << GETX((int)children[bestIdx].myIdx) << ","
				<< GETY((int)children[bestIdx].myIdx)
				<< endl;
				
								
			}
		}catch (...)
		{
			cerr << "ended because of time cutoff" << endl;
			break;
		}
		
		lastIdx = bestIdx;
		i = 0;
//		alpha -= 100;
		
		gettimeofday(&tv, NULL);
		tv = diff(begin, tv);
		cerr << "time: " << display_time( tv )
			<< "\tdepth:" << d << endl <<endl;
		tv = diff(tv, end);
		if( tv.tv_sec < (int)MAXSECS / 2)
			break;
	}
	
	gettimeofday(&tv, NULL);
	cerr << "completed in: " << display_time( diff(begin, tv))
	<< endl<< endl<< endl;
	return children[lastIdx];
	
}

int alpha_beta(Node &node, int alpha, int beta, const char *player,
			   const timeval &end, int depth)
{
	// for(int i = 1; i <= depth; i += 3)
	// 	cerr << " ";
	// cerr << "d:" << depth << "\t" << node << "\t\n\tend: " << display_time(end) << "\n";
	
	
	struct timeval tv;
	gettimeofday(&tv, NULL);
	if( past_time(tv, end) )
		throw "Ran out of time";
	
	if(game_over(node) || depth == 0)
	{
		// cerr << "\t\tleafNode:\t" << node<< endl;
		return node.evaluate() * *player;
	}
	
	// cerr << endl;
	
	vector<Node> children;
	char numKids = get_children(node, children, player);
//	sort(children.begin(), children.end(), by_max());
	// random_shuffle(children.begin(),children.end());
	
	int i = 0, value = MININT;
	for ( ; i < numKids; i++)
	{
		// gettimeofday(&tv, NULL);
		// add_time(tv, split_time( diff(tv, end), numKids - i));
		
		value = MAX(value,
					-1*(alpha_beta(children[i],
								   -1 * beta,
								   -1 * alpha,
								   (*player == MY) ? &OP : &MY,
								   end,
								   depth - 1)));
		// cerr << "\t\tresult: " << result <<endl;
		
		if( value >= beta)
		{
			return beta;
		}
		
		if( value > alpha)
			alpha = value;
	}
	// cerr << "alpha:\t" << alpha << endl;
	return alpha;
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
	
	
	
	
	if(playerNum == 2)
	{
		take_move(board, &opIdx);
	}
	
	Node initNode (board, &myIdx, &opIdx);
	draw_board(initNode);
	
	bool win = play(initNode);
	
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
	
	
	
	
	// vector<Node> children, grandchildren, great;
	// int n = get_children(initNode, children, &MY);
	// sort(children.begin(), children.end(), by_max());
	// Node min(initNode), max(initNode), max2(initNode);
	// max.heuristic = numeric_limits<int>::min()+10;
	// for (int i = 0; i < n; i++)
	// {
	// 	min.heuristic = numeric_limits<int>::max()-10;
	// 	cout << children[i] <<endl;
	// 	int m = get_children(children[i], grandchildren, &OP);
	// 	sort(grandchildren.begin(), grandchildren.end(), by_min());
	// 	for (int j =0; j < m; j++)
	// 	{
	// 		max2.heuristic = numeric_limits<int>::min()+10;
	// 		cout << "\t" << grandchildren[j] <<endl;
	// 		int o = get_children(grandchildren[i], great, &MY);
	// 		sort(great.begin(), great.end(), by_max());
	// 		for (int k =0; k < o; k++)
	// 		{
	// 			cout << "\t\t" << great[k] <<endl;
	// 			max2 = (great[k].heuristic > max2.heuristic)
	// 						? great[k] : max2;
	// 		}
	// 		cout << "\t\tmax2 " << max2 <<endl;
	
	// 		min = (max2.heuristic < min.heuristic)
	// 					? max2 : min;
	// 		great.clear();
	// 	}
	// 	cout << "\tmin " << min <<endl;
	
	// 	max = (min.heuristic > max.heuristic)
	// 				? min : max;
	// 	cout << "max " << max <<endl <<endl;
	
	// 	// for (int j =0; j < m; j++)
	// 	// 	cout << "\t" << grandchildren[j] <<endl;
	// 	grandchildren.clear();
	// }
	
	// cout << "min " << min << "\tmax " << max << endl;
	
	
	
	
	
	
	
	
	
	
	
	//---------------------------------------------
	
	
	// completed sucessfully
	return 0;
}










