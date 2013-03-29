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
#include <unordered_map> /* hash set */
#include <algorithm>	/* heap functions */


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define TOIDX(y,x) ((x) + ((y) * ROWSIZE))
#define SIGN(x) (((x) > 0) - ((x) < 0))
#define GETX(idx) ((idx) % ROWSIZE + 1)
#define GETY(idx) ((idx) / ROWSIZE + 1)




using namespace std;



static const char ROWSIZE = 8;
static const char BOARDSIZE = ROWSIZE * ROWSIZE;
static const char	NEAST = -ROWSIZE + 1,
					NORTH = -ROWSIZE,
					NWEST = -ROWSIZE - 1,
					WEST = -1,
					EAST = 1,
					SWEST = ROWSIZE - 1,
					SOUTH = ROWSIZE,
					SEAST = ROWSIZE + 1;

static const char DIRECTIONS[] = { NWEST,
									WEST,
									SWEST,
									NEAST,
									EAST,
									SEAST,
									NORTH,
									SOUTH};

static const char MY = 1;
static const char OP = -1;
static const char MAXDEPTH = 40;
static const int MAXSECS = 60;
static const int MININT = numeric_limits<int>::min() + 40;
static const int MAXINT = numeric_limits<int>::max() - 40;

static const char CONNCOMPLIMIT = ROWSIZE * ROWSIZE / 3;
static const int CLOSEDCOMPONENT = MAXINT / 64;


typedef bitset<BOARDSIZE> BitBoard;
class Node;
struct HashEntry;
struct by_max;
struct by_closeness;


string me, op, block = "|||";
char connCompSplitDepth = numeric_limits<char>::max();
vector<Node> moves;





int usage();
bool game_over(const Node&);
bool lookup(const Node &node);
bool store(const Node &node, const char &scoreType,
		   const int &score, const char &bestIdx);
int draw_board(const Node&);

string display_time( const timeval &tv );
timeval add_time( const timeval &tv, const int *seconds);
void add_time( timeval &to, const timeval &add);
bool past_time(const timeval &cur, const timeval &end);
timeval diff(const timeval &start, const timeval &end);



int evaluate_node( const Node &node );
bool valid_move( const BitBoard&, const char*, const char*);

char get_children( const Node&,	vector<Node>&, const char* player);
char fast_children(const Node&, const char*);


bool take_move( Node& );
bool play( Node& );


Node search_root(Node &, int&);
int alpha_beta(Node&, int, int, const char*, const timeval&, int);
int alpha_beta_transpo(Node&, int, int, const char*, const timeval&, int);

Node isolated_move(const Node &node);





class Node
{
public:
	BitBoard board;
	char opIdx;
	char myIdx;
	int heuristic;
	bool eval;
	
	Node(const BitBoard & board, const char &myIdx, const char &opIdx)
		:eval(false), myIdx(myIdx), opIdx(opIdx), board(board)
	{
		this->board.set(myIdx, 1);
		this->board.set(opIdx, 1);
	}
	
	Node( const Node &other )
		:eval(other.eval), myIdx(other.myIdx),
		opIdx(other.opIdx), board(other.board),
		heuristic(other.heuristic)
	{
		// cout << "node copied" << endl;
	}
	
	Node()
		:eval(false)
	{
	}
	
	friend ostream& operator<< (ostream &o, const Node &n)
	{
		o << "myP: " << (int) n.myIdx << "\topP: "
		<< (int) n.opIdx << "\th: " << n.heuristic;
		return o;
	}
};


struct HashEntry
{
	char scoreType;
	int score;
	// Relies on static ordering of returned moves
	// so it should be checked before sorting/shuffling
	char bestMoveIdx;
};






struct by_max
{
	bool operator() (const Node& lhs, const Node& rhs)
	{
		return evaluate_node(lhs) > evaluate_node(rhs);
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
	bool operator() (const char& lhs, const char& rhs)
	{
		char Lx = abs(GETX(lhs) - GETX(lhs));
		char Ly = abs(GETX(lhs) - GETX(lhs));
		
		char Rx = abs(GETX(rhs) - GETX(rhs));
		char Ry = abs(GETX(rhs) - GETX(rhs));
		return (Lx + Ly) < (Rx + Ry);
	}
	
};





int connect_comp( const Node &node, bool &opWithMe, const char &player)
{
	char pos, newPos,
	myIdx = (player == MY) ? node.myIdx : node.opIdx;
	int n = 0, i = 0;
	vector<char> h;
	unordered_map<char, bool> m;
	
	h.push_back(myIdx);
	make_heap(h.begin(), h.end(), by_closeness());
	
	while (!h.empty())
	{
		pos = h.front();
		pop_heap(h.begin(), h.end(), by_closeness()); h.pop_back();
		for (i = 0; i < sizeof(DIRECTIONS); i++)
		{
			newPos = pos + DIRECTIONS[i];
			if (m[newPos] || newPos < 0 || newPos >= BOARDSIZE)
				continue;
			
			m[newPos] = true;
			
			if( !node.board[newPos])
			{
				h.push_back(newPos); push_heap(h.begin(),h.end());
				n++;
			}
		}
		if (n == BOARDSIZE - node.board.count())
			break;
	}
	
	
	return n;
}


bool op_same_comp( const Node &node)
{
	char pos, newPos;
	int n = 0, i;
	vector<char> h;
	unordered_map<char, bool> m;
	
	h.push_back(node.myIdx);
	make_heap(h.begin(), h.end(), by_closeness());
	
	while (!h.empty())
	{
		pos = h.front();
		pop_heap(h.begin(), h.end(), by_closeness()); h.pop_back();
		for (i = 0; i < sizeof(DIRECTIONS); i++)
		{
			newPos = pos + DIRECTIONS[i];
			if (m[newPos] || newPos < 0 || newPos >= BOARDSIZE)
				continue;
			
			m[newPos] = true;
			
			if (newPos == node.opIdx)
				return true;
			if( !node.board[newPos])
			{
				h.push_back(newPos); push_heap(h.begin(),h.end());
				n++;
			}
		}
		if (n == BOARDSIZE - node.board.count())
			break;
	}
	
	
	return false;
}

int evaluate_node( const Node &node )
{
	if (!node.eval)
	{
		const_cast<Node&>(node).eval = true;
		
		if (node.board.count() > CONNCOMPLIMIT)
		{
			bool opWithMe = false;
			if (node.board.count() < connCompSplitDepth)
			{
				opWithMe = op_same_comp(node);
			}
			if (!opWithMe)
			{
				int mySpace = connect_comp(node, opWithMe, MY);
				int opSpace = connect_comp(node, opWithMe, OP);
				
				return const_cast<Node &>(node).heuristic =
				(opSpace == 0) ? MAXINT + mySpace:
				(mySpace == 0) ? MININT + mySpace:
				(mySpace > opSpace) ?
				CLOSEDCOMPONENT * (int)node.board.count():
				-1 * CLOSEDCOMPONENT * (int)node.board.count();
			}
		}
		int opMoves = fast_children(node, &OP);
		int myMoves = fast_children(node, &MY);
		return const_cast<Node &>(node).heuristic =
		(opMoves == 0) ? MAXINT + myMoves:
		(myMoves == 0) ? MININT + myMoves:
		(myMoves - opMoves*3 + 60) * (int)node.board.count();
	}
	return node.heuristic;
}

int usage()
{
	cout << "USAGE: $ isola PLAYER\n" <<
	"Where PLAYER is [1/2] representing "<<
	"the player this instance is"<< endl;
	return 1;
}

bool game_over(const Node &node)
{
	int v = evaluate_node(node);
	if(v == MAXINT || v == MININT)
		return true;
	else
		return false;
}

Node undo_last()
{
	Node prevNode = moves.back();
	moves.pop_back();
	return prevNode;
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

timeval add_time( const timeval &tv, const int *seconds)
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

string display_time( const timeval &tv )
{
	std::string number;
	std::stringstream strstream;
	strstream << tv.tv_sec << "." << tv.tv_usec << " seconds";
	strstream >> number;
	return number;
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
	char moveVec[2] = { y, x };
	char moveIncr = TOIDX( SIGN(moveVec[0]), SIGN(moveVec[1]) );
	
	if( abs(moveIncr) == 1 || abs(moveIncr) == ROWSIZE
	   || abs(moveIncr) == ROWSIZE - 1 || abs(moveIncr) == ROWSIZE + 1)
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


// Draws a Isola board represented by a bitset, when supplied
// with the bitset and two current positions as indexes in the bitset
int draw_board(const Node &node)
{
	// cout << "A represents me, B represents you, X a visited space\n\n";
	// cout << "draw:" <<endl;
	
	char i, j, idx;
	cout << "\n\n   |";
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
			idx = TOIDX(i, j);
			
			if( idx == node.myIdx )			cout << me;
			else if( idx == node.opIdx )	cout << op;
			else if( node.board[idx] )		cout << block;
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
	
	cout << "———————|\n\nMove: " << node.board.count() -2 << endl;
	
	return 0;
}




char get_children( const Node &node,
				  vector<Node> &children,
				  const char *player )
{
	// cout << "getc:" <<endl;
	// number of children found
	char n = 0, i = 0, curIdx=0;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx;
	
	for ( ; i < sizeof(DIRECTIONS); i++) {
		// cout << "fastm:" <<endl;
		curIdx = fromIdx + DIRECTIONS[i];
		
		while ( curIdx >= 0 && curIdx < BOARDSIZE && !node.board[curIdx] )
		{
			
			if(	(curIdx % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdx % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			Node child;
			// cerr << "node: " << node << endl << (int)curIdx<<endl;
			if(*player == MY)
				child = Node(node.board, curIdx, node.opIdx);
			else
				child = Node(node.board, node.myIdx, curIdx);
			
			//		child.evaluate();
			children.push_back(child);
			n++;
			curIdx += DIRECTIONS[i];
		}
	}
	
	return n;
}

char fast_children(const Node &node, const char *player)
{
	// cout << "getc:" <<endl;
	// number of children found
	unsigned char n = 0, i = 0;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx,
	curIdx;
	
	for ( ; i < sizeof(DIRECTIONS); i++) {
		// cout << "fastm:" <<endl;
		curIdx = fromIdx + DIRECTIONS[i];
		// cout << curIdx << endl;
		while ( curIdx >= 0 && curIdx < BOARDSIZE && !node.board[curIdx] )
		{
			
			if(	(curIdx % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdx % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			n++;
			curIdx += DIRECTIONS[i];
		}
	}
	
	return n;
}

char fast_children_both(const Node &node, const char *player, char *myMoves, char *opMoves)
{
	// cout << "getc:" <<endl;
	// number of children found
	char i = 0, curIdxMY, curIdxOP;
	char fromIdx = (*player == MY) ? node.myIdx : node.opIdx;
	
	for ( ; i < sizeof(DIRECTIONS); i++) {
		// cout << "fastm:" <<endl;
		curIdxMY = fromIdx + DIRECTIONS[i];
		curIdxOP = fromIdx + DIRECTIONS[i];
		
		while ( !node.board[curIdxMY] )
		{
			
			if(	curIdxMY < 0 || curIdxMY >= BOARDSIZE
			   || (curIdxMY % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdxMY % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			myMoves++;
			curIdxMY += DIRECTIONS[i];
		}
		while ( !node.board[curIdxOP] )
		{
			
			if(	curIdxOP < 0 || curIdxOP >= BOARDSIZE
			   || (curIdxOP % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdxOP % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			opMoves++;
			curIdxOP += DIRECTIONS[i];
		}
	}
	
	return myMoves - opMoves;
}








bool take_move( Node &node )
{
	// cout << "take:" <<endl;
	char newIdx = 0,
	n = fast_children(node, &OP);
	
	//	cout << "numKids Op: " << (int) n << endl;
	while( n != 0 )
	{
		string move;
		cout << "Please enter a move (row col):  ";
		getline(cin, move);
		
		size_t undo = move.find_first_of("$");
		if (undo != string::npos)
		{
			Node prevNode = undo_last();
			node.board = prevNode.board;
			node.myIdx = prevNode.myIdx;
			node.opIdx = prevNode.opIdx;
			node.heuristic = prevNode.heuristic;
			node.eval = node.eval;
			draw_board(node);
			continue;
		}
		
		size_t y = move.find_first_of("12345678");
		size_t x = move.find_first_of("12345678", y+1);
		if (x == string::npos || y == string::npos)
		{
			cout << "Incorrect move format, try again." <<endl;
			continue;
		}
		
		newIdx = TOIDX(move[y] - '0' - 1, move[x] - '0' - 1);
		
		if( !valid_move(node.board, &node.opIdx, &newIdx) )
		{
			cout << "Invalid move, try again." <<endl;
			continue;
		}
		
		moves.push_back(node);
		node.opIdx = newIdx;
		node.board.set(newIdx, 1);
		break;
	}
	
	// success
	n = fast_children(node, &OP);
	if(n == 0)
		return false;
	return true;
}




Node search_root(Node &initNode, int &alpha)
{
	// set time limit
	timeval begin, end, tv, tmp;
	gettimeofday(&begin, NULL);
	end = add_time(begin, &MAXSECS);
	
	
	vector<Node> children;
	int numKids = get_children(initNode, children, &MY);
	sort(children.begin(), children.end(), by_max());
	// random_shuffle(children.begin(),children.end());
	
	if (initNode.board.count() < connCompSplitDepth)
	{
		bool opWithMe = op_same_comp(initNode);
		if(!opWithMe)
		{
			connCompSplitDepth = initNode.board.count();
		}
	}
	
	int i = 0, d = 0, value, beta = MAXINT, bestIdx = 0, lastIdx = 0;
	for(; d < MAXDEPTH ; d += 2)
	{
		try{
			for (; i < numKids; i++)
			{
				value = -1 * alpha_beta(children[i],
										-1*beta,
										-1*alpha,
										&OP,
										end,
										d);
				
				cerr << "last: " << value
				<< "  pos: " << GETY((int)children[i].myIdx) << ","
				<< GETX((int)children[i].myIdx);
				
				
				if(value == MAXINT)
					return children[i];
				
				if (value > alpha)
				{
					alpha = value;
					bestIdx = i;
				}
				
				cerr <<"\t\tbest: " << alpha
				<< "  pos: " << GETY((int)children[bestIdx].myIdx) << ","
				<< GETX((int)children[bestIdx].myIdx)
				<< endl;
				
			}
			if(alpha == MININT)
				return children[bestIdx];
			
		}catch (...)
		{
			alpha *= .9;
//			alpha = MININT;
			cerr << "ended because of time cutoff" << endl;
			break;
		}
		cerr << "best: " << alpha
		<< "  pos: " << GETY((int)children[bestIdx].myIdx) << ","
		<< GETX((int)children[bestIdx].myIdx)<<endl;
		
		lastIdx = bestIdx;
		i = 0;
		alpha *= .9;
//		alpha = MININT;
		
		gettimeofday(&tmp, NULL);
		tv = diff(begin, tmp);
		cerr << "time: " << display_time( tv )
		<< "\tdepth:" << d << endl <<endl;
		tv = diff(tmp, end);
		if( tv.tv_sec < (3 * (long)MAXSECS / 4))
		{
			cerr << "early cuttoff:\t" << tv.tv_sec
					<< "\t" << (3*(long)MAXSECS/4) << endl;
			break;
		}
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
		return evaluate_node(node) * *player;
	}
	
	// cerr << endl;
	
	vector<Node> children;
	char numKids = get_children(node, children, player);
	//	sort(children.begin(), children.end(), by_max());
	random_shuffle(children.begin(),children.end());
	
	int i = 0, value = MININT;
	for ( ; i < numKids; i++)
	{
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

int alpha_beta_transpo(Node &node, int alpha, int beta, const char *player,
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
		return evaluate_node(node) * *player;
	}
	
	// cerr << endl;
	
//	HashEntry entry;
//	bool hash_hit = lookup(node);
	
	vector<Node> children;
	char numKids = get_children(node, children, player);
	//	sort(children.begin(), children.end(), by_max());
	random_shuffle(children.begin(),children.end());
	
	int i = 0, value = MININT;
	for ( ; i < numKids; i++)
	{
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










bool play( Node &curNode )
{
	// cout << "play:" <<endl;
	evaluate_node(curNode);
	int alpha = MININT;
	
	cout << "\n\n";
	draw_board(curNode);
	cout << "\n\n";
	
	while( fast_children(curNode, &MY) > 0)
	{
		// cout << "playloop:" <<endl;
		cout <<"before: "<< curNode <<endl << endl;
		
//		alpha = MININT;
		
		curNode = search_root(curNode, alpha);
		
		cout << "\n\n";
		draw_board(curNode);
		cout << "\n\n";
		
		
		if(evaluate_node(curNode) == MININT)
			return false;
		cout <<"after: "<< curNode <<endl << endl <<endl;
		
		
		if( !take_move(curNode))
			return true;
		
		cout << "\n\n";
		draw_board(curNode);
		cout << "\n\n";
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
		me = "*X*";
		op = "*O*";
	}
	else
	{
		myIdx = TOIDX(ROWSIZE-1, ROWSIZE-1);
		opIdx = TOIDX(0, 0);
		me = "*O*";
		op = "*X*";
	}
	


	//TESTING/////////

	myIdx = TOIDX(6, 6);
	opIdx = TOIDX(7, 6);
	me = "*X*";
	op = "*O*";

	board.set(myIdx, 1);
	board.set(opIdx, 1);
	board.set(TOIDX(0, 0), 1);
	board.set(TOIDX(ROWSIZE-1, ROWSIZE-1), 1);
	board.set(TOIDX(1, 4), 1);
	board.set(TOIDX(1, 6), 1);
	board.set(TOIDX(2, 2), 1);
	board.set(TOIDX(2, 3), 1);
	board.set(TOIDX(2, 4), 1);
	board.set(TOIDX(2, 5), 1);
	board.set(TOIDX(2, 6), 1);
	board.set(TOIDX(2, 7), 1);
	board.set(TOIDX(3, 1), 1);
	board.set(TOIDX(3, 5), 1);
	board.set(TOIDX(3, 6), 1);
	board.set(TOIDX(3, 7), 1);
	board.set(TOIDX(4, 2), 1);
	board.set(TOIDX(4, 3), 1);
	board.set(TOIDX(4, 6), 1);
	board.set(TOIDX(5, 4), 1);
	board.set(TOIDX(5, 5), 1);
	board.set(TOIDX(6, 7), 1);
	
	
	Node initNode (board, myIdx, opIdx);
	cerr << initNode << endl;
	moves.push_back(initNode);
	
	//	cerr << connect_comp(initNode, opWithMe, MY) << " " << opWithMe << endl;
	
	
	if(playerNum == 2)
	{
		draw_board(initNode);
		take_move(initNode);
	}
	
	
	bool win = play(initNode);
	
	cout << "\n\n";
	draw_board(initNode);
	cout << "\n\n";
	
	if ( win )
		cout << "!!!!!!!    I win    !!!!!!!!!" << endl;
	else
		cout << ":( :(  You win  :( :(" << endl;
	
	
	
	// completed sucessfully
	return 0;
}










