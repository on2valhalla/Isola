Jason Mann
jcm2207

USAGE: 

$ make all

$ ./isola PLAYER TIME_LIMIT
-Where PLAYER is (1 or 2) with 1 corresponding to the x player and 2-O
-Where TIME_LIMIT is number of seconds as a maximum limit


Algorithm Background:
My implementation of the minimax algorithm with alpha-beta pruning follows the Negamax/Negascout pattern with iterative deepening search.  I have also implemented a hash table lookup that stores whether a node found a lower/upper bound, or exact score, as well as the depth the node was explored to previously. I decided on this combination after testing different combinations of algorithms including the above and MTDf, killer heuristics, history heuristics, aspiration windows and a couple of others. The combination of iterative deepening, negamax, and a transposition table beat any other combination by at least 2 plys.

My heuristic is two part: one for when both players are in the same component, and one for rating boards in which the players are in different components. The first heuristic is of the form (myMoves - 3*opMoves) * filledSpaces. myMoves and opMoves are respectively the amount of possible valid movements from the position of the computer and opponent at any one state. filledSpaces is the number of spaces on the board that are filled. The constant 3 was determined through testing to be the optimal combination in early games. After a quarter of the board is filled, the heuristic includes a check to see whether or not both players are in the same connected component. If they are not, the evaluation function instead returns a judgement of the size of the computers space vs the size of the opponents space using a value that is the maximum or minimum integer divided by the amount of spaces on the board. This value is then also adjusted to give preference to depth.


Defense of evaluation function:
1)	This heuristic function provides a good approximation of goodness, while it not only rates the computers current position against its opponent, but it also considers whether a space is closed off or not. This allows it to virtually look ahead more plys than are computationally possible, predicting the winner or loser in practice up to 25 ply ahead.

2)	By having a bipartite heuristic, I have chosen a balance of speed versus accuracy, allowing for further look ahead early in the game, while preventing being closed in a smaller space in the mid to end game.

3)	I determined this combination was worthwhile initially by human v. human games that I played myself, while I was testing possible options. Later I tested this heuristic against many other variations of itself, as well as completely separate heuristics keeping everything else the same.

4)	Besides the heuristic efficiency of not having to search as deep by virtually looking ahead, and the implementations of more specific algorithms above, I have sought out efficiency in other aspects of the program. Firstly, the greatest efficiency of both time and space was to use std::bitset for the board, which allows storage and modification of independent bits. This was a difference of 54 BYTES PER STATE when compared to an array of chars (the next smallest storage)! Saving that much memory allocation, as well as bitset's speedy access, also contributed to time efficiency. Also if at all possible all arguments to functions are passed by const reference, and if not const then references. Also the char data type is used anywhere that low numbers such as indicies are needed.



Extensions:
I guess you can label std::bitset as an extension, I also used std::unordered_map for my transposition table, while on clic machines the c++0x standard unordered_map beat both google's dense_hash_map and boost's unordered_map in speed tests (though these results are not shown elsewhere). I borrowed boost library's hash combination algorithm, to combine the std:hash<bitset> and the hashes for the char type position indexes. One final speedup is that I initialize the transposition table to 150,000,000 buckets, which is the amount for 75% capacity of the average end-game # of nodes stored. While this takes up quite a bit of memory outright, it saves time that would be spent rehashing.