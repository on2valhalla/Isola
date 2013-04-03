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
#include <fstream>		/* file operations */
#include <bitset>		/* bitset */
#include <string>		/* string */
#include <queue>		/* priority queue */
#include <limits>		/* min and max for types */
#include <ctime>		/* for system time */
#include <sstream>		/* stringstream */
#include <unordered_map> /* hash set */
#include <algorithm>	/* heap functions */
//#include <boost/serialization/serialization.hpp>
//#include <boost/serialization/map.hpp>
//#include <boost/serialization/bitset.hpp>
//#include <boost/serialization/unordered_map.hpp>
//#include <boost/serialization/dense_hash_map.hpp>
//#include <boost/serialization/sparse_hash_map.hpp>
//#include "dense_hash_map.hpp"
//#include "sparse_hash_map.hpp"
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <google/dense_hash_map>
//#include <google/sparse_hash_map>
//#include "sparsehash/dense_hash_map"


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define TOIDX(y,x) ((x) + ((y) * ROWSIZE))
#define SIGN(x) (((x) > 0) - ((x) < 0))
#define GETX(idx) ((idx) % ROWSIZE + 1)
#define GETY(idx) ((idx) / ROWSIZE + 1)




using namespace std;

//-------------------GLOBAL CONSTANTS-----------------------------------------
//
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
static const int MININT = numeric_limits<int>::min() + 100;
static const int MAXINT = numeric_limits<int>::max() - 100;

static const char CONNCOMPLIMIT = ROWSIZE * ROWSIZE / 4;
static const int CLOSEDCOMPONENT = MAXINT / 64;

static const char LOWERBOUND = 2;
static const char UPPERBOUND = 3;
static const char EXACTSCORE = 4;
//
//------------------------------------------------------------------------------



//-----DECLARATIONS-------------------------------------------------------------
//


class Node;
struct HashEntry;
struct by_max;
struct by_closeness;


typedef bitset<BOARDSIZE> BitBoard;
typedef unordered_map<Node, HashEntry, hash<Node>, equal_to<Node> > NodeMap;
//typedef google::dense_hash_map<Node, HashEntry, hash<Node>, equal_to<Node> > NodeMap;
//typedef google::sparse_hash_map<Node, HashEntry, hash<Node>, equal_to<Node> > NodeMap;
//typedef unordered_map<size_t, HashEntry> NodeMap;
//typedef google::dense_hash_map<size_t, HashEntry> NodeMap;
//typedef google::sparse_hash_map<size_t, HashEntry> NodeMap;


int usage();
bool game_over(const Node&);
int draw_board(const Node&);

string display_time( const timeval &tv );
timeval add_time( const timeval &tv, const int *seconds);
bool past_time(const timeval &cur, const timeval &end);
timeval diff(const timeval &start, const timeval &end);



int evaluate_node( const Node &node );
bool valid_move( const BitBoard&, const char*, const char*);

char get_children( const Node&,	vector<Node>&, const char* player);
char fast_children_both(const Node &node, char *myMoves, char *opMoves);


bool take_move( Node& , int&);
bool play( Node& );


Node search_root(Node &, int&);
int alpha_beta(Node&, int, int, const char*, const timeval&, char);


bool lookup(const Node &node, HashEntry &entry);
void store(const Node &node, const char &scoreType,
		   const int &score, const char &depth);
int alpha_beta_transpo(Node&, int, int, const char*, const timeval&, char);
//
//-----------------------------------------------------------------------------




//-----CLASSES AND STRUCTS------------------------------------------------------
//
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
		o << "  myPos: " << GETY((int)n.myIdx) << ","
		<< GETX((int)n.myIdx) << "  opPos: " << GETY((int)n.opIdx) << ","
		<< GETX((int)n.opIdx) << "\th: " << evaluate_node(n);
		return o;
	}
	
	//
	//    friend class boost::serialization::access;
	//    // When the class Archive corresponds to an output archive, the
	//    // & operator is defined similar to <<.  Likewise, when the class Archive
	//    // is a type of input archive the & operator is defined similar to >>.
	//    template<class Archive>
	//    void serialize(Archive & ar, const unsigned int version)
	//    {
	//        ar & board;
	//        ar & myIdx;
	//        ar & opIdx;
	//    }
};


struct HashEntry
{
	char scoreType;
	int score;
	// Relies on static ordering of returned moves
	// so it should be checked before sorting/shuffling
	char depth;
	
	//
	//    friend class boost::serialization::access;
	//    // When the class Archive corresponds to an output archive, the
	//    // & operator is defined similar to <<.  Likewise, when the class Archive
	//    // is a type of input archive the & operator is defined similar to >>.
	//    template<class Archive>
	//    void serialize(Archive & ar, const unsigned int version)
	//    {
	//        ar & scoreType;
	//        ar & score;
	//        ar & depth;
	//    }
	//
};

namespace std {
	// Hash combination emulates from Boost library
	template<>
	class hash<Node> {
	public:
		size_t operator()(const Node &n) const
		{
			hash<bitset<BOARDSIZE> > bHash;
			hash<char > cHash;
			size_t hash = bHash(n.board);
			hash ^= cHash(n.myIdx)
			+ 0x9e3779b9 + (hash << 6) + (hash >> 2);
			hash ^= cHash(n.opIdx)
			+ 0x9e3779b9 + (hash << 6) + (hash >> 2);
			return hash;
		}
	};
	template<> class equal_to<Node>
	{
	public:
		bool operator() (const Node& lhs, const Node& rhs) const
		{
			return lhs.board == rhs.board
			&& lhs.myIdx == rhs.myIdx
			&& lhs.opIdx == rhs.opIdx;
		}
		
	};
}


// COMPARISON STRUCTS
struct by_max
{
	bool operator() (const Node& lhs, const Node& rhs)
	{
		return evaluate_node(lhs) > evaluate_node(rhs);
	}
	
};
//
//-----------------------------------------------------------------------------





//-----GLOBAL VARIABLES---------------------------------------------------------
//
string me, op, block = "|||";

char connCompSplitDepth = numeric_limits<char>::max();
bool splitBoards = false;

int maxSecs = 60;

vector<Node> moves;

NodeMap transpos(100000000);
//
//-----------------------------------------------------------------------------










int usage()
{
	cout << "USAGE: $ isola PLAYER TIME_LIMIT\n" <<
	"Where PLAYER is [1/2] representing the player this instance is"
	<< "\n and TIME_LIMIT is the maximum allotted time" endl;
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




int connect_comp( const Node &node, const char &player)
{
	char pos, newPos,
	myIdx = (player == MY) ? node.myIdx : node.opIdx;
	int n = 0, i = 0;
	vector<char> h;
	unordered_map<char, bool> m;
	
	h.push_back(myIdx);
	
	while (!h.empty())
	{
		pos = h.back();
		h.pop_back();
		for (i = 0; i < sizeof(DIRECTIONS); i++)
		{
			newPos = pos + DIRECTIONS[i];
			if (m[newPos] || newPos < 0 || newPos >= BOARDSIZE
				|| (newPos % ROWSIZE == ROWSIZE - 1 && i <= 2)
				|| (newPos % ROWSIZE == 0 && (i > 2 && i <= 5)))
				continue;
			
			m[newPos] = true;
			
			if( !node.board[newPos])
			{
				h.push_back(newPos);
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
	char pos, newPos, i;
	vector<char> h;
	unordered_map<char, bool> m;
	
	h.push_back(node.myIdx);
	
	while (!h.empty())
	{
		pos = h.back();
		h.pop_back();
		for (i = 0; i < sizeof(DIRECTIONS); i++)
		{
			newPos = pos + DIRECTIONS[i];
			if (m[newPos] || newPos < 0 || newPos >= BOARDSIZE
				|| (newPos % ROWSIZE == ROWSIZE - 1 && i <= 2)
				|| (newPos % ROWSIZE == 0 && (i > 2 && i <= 5)))
				continue;
			
			m[newPos] = true;
			
			if (newPos == node.opIdx)
				return true;
			if( !node.board[newPos])
			{
				h.push_back(newPos);
			}
		}
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
			if (!splitBoards)
			{
				if (!op_same_comp(node))
				{
					int mySpace = connect_comp(node, MY);
					int opSpace = connect_comp(node, OP);
					//					cerr << "\tmyspace: " << mySpace
					//					<< "\topspace: " << opSpace << "\t" << node <<endl;
					
					return const_cast<Node &>(node).heuristic =
					(opSpace == 0) ? MAXINT + mySpace:
					(mySpace == 0) ? MININT - opSpace + 40:
					(mySpace > opSpace) ?
					CLOSEDCOMPONENT * (int)node.board.count() - opSpace:
					-1 * CLOSEDCOMPONENT * (int)node.board.count() - opSpace;
				}
			}
		}
		char myMoves=0, opMoves=0;
		fast_children_both(node, &myMoves, &opMoves);
		
		return const_cast<Node &>(node).heuristic =
		(opMoves == 0) ? MAXINT + myMoves:
		(myMoves == 0) ? MININT - opMoves + 40:
		(myMoves - opMoves*3 + 60) * (int)node.board.count();
	}
	return node.heuristic;
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
			children.push_back(move(child));
			n++;
			curIdx += DIRECTIONS[i];
		}
	}
	
	return n;
}


char fast_children_both(const Node &node, char *myMoves, char *opMoves)
{
	// cout << "getc:" <<endl;
	// number of children found
	char i = 0, curIdxMY, curIdxOP;
	
	for ( ; i < sizeof(DIRECTIONS); i++) {
		// cout << "fastm:" <<endl;
		curIdxMY = node.myIdx + DIRECTIONS[i];
		curIdxOP = node.opIdx + DIRECTIONS[i];
		
		while ( !node.board[curIdxMY] )
		{
			
			if(	curIdxMY < 0 || curIdxMY >= BOARDSIZE
			   || (curIdxMY % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdxMY % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			(*myMoves)++;
			curIdxMY += DIRECTIONS[i];
		}
		while ( !node.board[curIdxOP] )
		{
			
			if(	curIdxOP < 0 || curIdxOP >= BOARDSIZE
			   || (curIdxOP % ROWSIZE == ROWSIZE - 1 && i <= 2)
			   || (curIdxOP % ROWSIZE == 0 && (i > 2 && i <= 5)))
				break; //edge of table reached
			
			(*opMoves)++;
			curIdxOP += DIRECTIONS[i];
		}
	}
	
	return (*myMoves) - (*opMoves);
}









Node search_root(Node &initNode, int &alpha)
{
	// set time limit
	timeval begin, end, tv, tmp;
	gettimeofday(&begin, NULL);
	end = add_time(begin, &maxSecs);
	
	
	vector<Node> children;
	int numKids = get_children(initNode, children, &MY);
	sort(children.begin(), children.end(), by_max());
	// random_shuffle(children.begin(),children.end());
	
	if (!splitBoards)
	{
		splitBoards = !op_same_comp(initNode);
		cerr << "boards are split: " << splitBoards << endl;
	}
	
	int i = 0, d = 0, value, beta = MAXINT, bestIdx = 0, lastIdx = 0;
	for(; d < MAXDEPTH ; d += 1)
	{
		try{
			for (; i < numKids; i++)
			{
				//				value = -1 * alpha_beta(children[i],
				//										-1*beta,
				//										-1*alpha,
				//										&OP,
				//										end,
				//										d);
				value = -1 * alpha_beta_transpo(children[i],
												-1*beta,
												-1*alpha,
												&OP,
												end,
												d);
				
				cerr << "last: " << value
				<< "  pos: " << GETY((int)children[i].myIdx) << ","
				<< GETX((int)children[i].myIdx);
				
				
				
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
			//			if(alpha == MININT)
			//				return children[bestIdx];
			if(value >= MAXINT)
				return children[bestIdx];
			
		}catch (...)
		{
			//			alpha *= .9;
			alpha = MININT;
			cerr << "ended because of time cutoff" << endl;
			break;
		}
		cerr << "best: " << alpha
		<< "  pos: " << GETY((int)children[bestIdx].myIdx) << ","
		<< GETX((int)children[bestIdx].myIdx)<<endl;
		
		lastIdx = bestIdx;
		i = 0;
		//	alpha *= .9;
		alpha = MININT;
		
		gettimeofday(&tmp, NULL);
		tv = diff(begin, tmp);
		cerr << "time: " << display_time( tv )
		<< "\tdepth:" << d << endl <<endl;
		//		tv = diff(tmp, end);
		//		if( tv.tv_sec < EARLYCUTOFF)
		//		{
		//			cerr << "early cuttoff:\t" << tv.tv_sec
		//			<< "\t" << (long)maxSecs/4 << endl;
		//			break;
		//		}
	}
	
	gettimeofday(&tv, NULL);
	cerr << "completed in: " << display_time( diff(begin, tv))
	<< endl<< endl<< endl;
	return children[lastIdx];
	
}


int alpha_beta_transpo(Node &node, int alpha, int beta, const char *player,
					   const timeval &end, char depth)
{
	// for(int i = 1; i <= depth; i += 3)
	// 	cerr << " ";
	// cerr << "d:" << depth << "\t" << node << "\t\n\tend: " << display_time(end) << "\n";
	
	
	struct timeval tv;
	gettimeofday(&tv, NULL);
	if( past_time(tv, end) )
		throw "Ran out of time";
	
	
	//	if(game_over(node) || depth == 0)
	if(depth == 0)
	{
		// cerr << "\t\tleafNode:\t" << node<< endl;
		return evaluate_node(node) * *player;
	}
	
	HashEntry entry;
	bool hash_hit = lookup(node, entry);
	if(hash_hit && entry.depth >= depth)
	{
		switch(entry.scoreType)
		{
			case EXACTSCORE:
				return entry.score;
				break;
			case LOWERBOUND:
				if(alpha < entry.score)
					alpha = entry.score;
				break;
			case UPPERBOUND:
				if(beta > entry.score)
					beta = entry.score;
				break;
		}
		if(alpha >= beta)
			return entry.score;
	}
	
	
	// cerr << endl;
	
	
	
	vector<Node> children;
	char numKids = get_children(node, children, player);
	//	sort(children.begin(), children.end(), by_max());
	random_shuffle(children.begin(),children.end());
	
	int i = 0, value, best = MININT;
	for ( ; i < numKids; i++)
	{
		value = -1*(alpha_beta_transpo(children[i],
									   -1 * beta,
									   -1 * alpha,
									   (*player == MY) ? &OP : &MY,
									   end,
									   depth - 1));
		// cerr << "\t\tresult: " << result <<endl;
		
		if( value > best)
			best = value;
		if( best > alpha)
			alpha = best;
		if( best >= beta)
			break;
	}
	// cerr << "alpha:\t" << alpha << endl;
	
	if(best <= alpha) // a lowerbound value
		store(node, LOWERBOUND, best, depth);
	else if(best >= beta) // an upperbound value
		store(node, UPPERBOUND, best, depth);
	else // a true minimax value
		store(node, EXACTSCORE, best, depth);
	return best;
}







bool lookup(const Node &node, HashEntry &entry)
{
	NodeMap::const_iterator got = transpos.find(node);
	//	NodeMap::const_iterator got = transpos.find(hash<Node>()(node));
	if( got == transpos.end())
		return false;
	else
	{
		entry = got->second;
		return true;
	}
}

void store(const Node &node, const char &scoreType,
		   const int &score, const char &depth)
{
	HashEntry entry;
	entry.scoreType = scoreType;
	entry.score = score;
	entry.depth = depth;
	transpos[node] = entry;
	//	transpos[hash<Node>()(node)] = entry;
}






bool take_move( Node &node, int &alpha)
{
	// cout << "take:" <<endl;
	char newIdx = 0, myMoves = 0, opMoves = 0;
	fast_children_both(node, &myMoves, &opMoves);
	
	//	cout << "numKids Op: " << (int) n << endl;
	while( opMoves != 0 )
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
			alpha = MININT;
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
	opMoves = 0;
	fast_children_both(node, &myMoves, &opMoves);
	if(opMoves == 0)
		return false;
	return true;
}

bool make_move( Node &node, char player)
{
	// cout << "take:" <<endl;
	char newIdx = 0, myMoves = 0, opMoves = 0, localPlayer = player;
	fast_children_both(node, &myMoves, &opMoves);
	
	//	cout << "numKids Op: " << (int) n << endl;
	while( true )
	{
		string move;
		cout << "Please enter a move (row col):  ";
		getline(cin, move);
		
		size_t begin = move.find_first_of("#");
		if (begin != string::npos)
		{
			break;
		}
		
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
		
		char *curIdx = (localPlayer == MY) ? &node.myIdx : &node.opIdx;
		
		size_t y = move.find_first_of("12345678");
		size_t x = move.find_first_of("12345678", y+1);
		if (x == string::npos || y == string::npos)
		{
			cout << "Incorrect move format, try again." <<endl;
			continue;
		}
		
		newIdx = TOIDX(move[y] - '0' - 1, move[x] - '0' - 1);
		
		if( !valid_move(node.board, curIdx, &newIdx) )
		{
			cout << "Invalid move, try again." <<endl;
			continue;
		}
		
		*curIdx = newIdx;
		node.board.set(newIdx, 1);
		moves.push_back(node);
		draw_board(node);
		localPlayer *= -1;
	}
	
	int fakealpha;
	if (player != localPlayer)
		take_move(node, fakealpha);
	
	// success
	opMoves = 0;
	fast_children_both(node, &myMoves, &opMoves);
	if(opMoves == 0)
		return false;
	return true;
}


bool play( Node &curNode )
{
	// cout << "play:" <<endl;
	evaluate_node(curNode);
	int alpha = MININT;
	
	cout << "\n\n";
	draw_board(curNode);
	cout << "\n\n";
	
	char myMoves = 0, opMoves = 0;
	fast_children_both(curNode, &myMoves, &opMoves);
	
	while( myMoves > 0)
	{
		//		cout << "playloop:" <<endl;
		//		cout <<"before: "<< curNode <<endl << endl;
		
		//		if ( alpha > CLOSEDCOMPONENT )
		alpha = MININT;
		
		curNode = search_root(curNode, alpha);
		
		cerr <<"\n\nafter: "<< curNode <<"\ttranspo: "<< transpos.size()
		<<"\tbuckets: "<< transpos.bucket_count()
		<<endl << endl <<endl;
		
		cout << "\n\n";
		draw_board(curNode);
		cout << "\n\n";
		
		
		cout << "I MOVE TO:  (" << GETY(curNode.myIdx)
		<< " " << GETX(curNode.myIdx) << ")\n" << endl;
		
		if(evaluate_node(curNode) == MININT)
			return false;
		
		
		if( !take_move(curNode, alpha))
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
	 * $ isola PLAYER TIME_LIMIT
	 *
	 *	Where PLAYER is [1 or 2] representing the player this instance is
	 *
	 */
	
	
	// only needed for dense hashmap
	//	transpos.set_empty_key(Node());
	
	char playerNum;
	//	string filename;
	try
	{ // take in command line argument
		assert(argc >= 3);
		playerNum = argv[1][0] - '0';
		maxSecs = atoi(argv[2]);
		//		makemove = argv[2][0];
		
		//		filename = argv[2];
		//
		//		string read;
		//		cout << "do you want to read from file (y/n):  ";
		//		getline(cin, read);
		//
		//		if (read == "y")
		//			read_table(filename);
		
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
	//
	//	myIdx = TOIDX(6, 6);
	//	opIdx = TOIDX(7, 6);
	//	me = "*X*";
	//	op = "*O*";
	//
	//	board.set(myIdx, 1);
	//	board.set(opIdx, 1);
	//	board.set(TOIDX(0, 0), 1);
	//	board.set(TOIDX(ROWSIZE-1, ROWSIZE-1), 1);
	//	board.set(TOIDX(1, 4), 1);
	//	board.set(TOIDX(1, 6), 1);
	//	board.set(TOIDX(2, 2), 1);
	//	board.set(TOIDX(2, 3), 1);
	//	board.set(TOIDX(2, 4), 1);
	//	board.set(TOIDX(2, 5), 1);
	//	board.set(TOIDX(2, 6), 1);
	//	board.set(TOIDX(2, 7), 1);
	//	board.set(TOIDX(3, 1), 1);
	//	board.set(TOIDX(3, 5), 1);
	//	board.set(TOIDX(3, 6), 1);
	//	board.set(TOIDX(3, 7), 1);
	//	board.set(TOIDX(4, 2), 1);
	//	board.set(TOIDX(4, 3), 1);
	//	board.set(TOIDX(4, 6), 1);
	//	board.set(TOIDX(5, 4), 1);
	//	board.set(TOIDX(5, 5), 1);
	//	board.set(TOIDX(6, 7), 1);
	
	
	Node initNode (board, myIdx, opIdx);
	//	cerr << initNode << endl;
	moves.push_back(initNode);
	
	string makeMove;
	cout << "Would you like to start manually?";
	getline(cin, makeMove);
	
	size_t begin = makeMove.find_first_of("y");
	if (begin != string::npos)
	{
		make_move(initNode, (playerNum == 1) ? MY : OP);
	}
	
	//	cerr << connect_comp(initNode, opWithMe, MY) << " " << opWithMe << endl;
	
	int alpha = MININT;
	if(playerNum == 2)
	{
		draw_board(initNode);
		take_move(initNode, alpha);
	}
	
	
	bool win = play(initNode);
	
	cout << "\n\n";
	draw_board(initNode);
	cout << "\n\n";
	
	if ( win )
		cout << "!!!!!!!    I win    !!!!!!!!!" << endl;
	else
		cout << ":( :(  You win  :( :(" << endl;
	
	//	serialize_table(filename);
	// completed sucessfully
	return 0;
}








