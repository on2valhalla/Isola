#include <stdio>      /* printf, fgets */
#include <stdlib>     /* atoi */







int usage()
{
	cout << "USAGE: $ isola SIZE PLAYER\n" <<
			"Where SIZE is a single int board" <<
			" size (normally should be 8)\n" <<
			"and player is [1/2] representing "<<
			"the player this instance is"<< endl;
	return 1;
}






int main(int argc, char *argv[])
{
	/* USAGE
	*
	* $ isola SIZE PLAYER
	*
	*	where SIZE is a single int board size (normally should be 8)
	*	and player is [1/2] representing the player this instance is
	*
	*/

	if( argc < 3 )
		return usage();


	const int BOARDSIZE = 



	// completed sucessfully
	return 0;
}