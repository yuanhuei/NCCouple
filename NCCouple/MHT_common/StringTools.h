
#ifndef _StringTools_
#define _StringTools_

#include <fstream>
#include <sstream>
#include <string>

//This function can be used to extract content between a couple of symbols, 
//such as '(' and ')' from a given file.
//function returns:
//	0 for succussful extraction, and
//	1 for failure.
int GetBySymbol
(
std::ifstream& infile,
std::stringstream& outstringstream,
char leftChar,
char rightChar
);


//This function is designed to extract content between a couple of symbols, 
//such as '(' and ')' from a given stringstream.
//function returns:
//	0 for succussful extraction, and
//	1 for failure.
int GetBySymbol
(
std::stringstream& instringstream,
std::stringstream& outstringstream,
char leftChar,
char rightChar
);

bool getline(std::stringstream& sstream, std::string& outstring, char symbol);

#endif
