/*---------------------------------------------------------------------------*\
File Name:
    StringTools.h

Description:
    Manipulations on files and strings

    Author:		Kong Ling
    Date:   2016-09-10
\*---------------------------------------------------------------------------*/

#include "../MHT_common/StringTools.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

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
)
{
	char onechar;
	outstringstream.str("");
	//find the position of the first left parathese
	while (infile.get(onechar))
	{
		if (onechar == leftChar) break;
	}

	if (leftChar == rightChar)
	{
		//For example you want to get the content between a pair of double quotes,
		//in which the left char (\") is the same as the right one
		while (infile.get(onechar))
		{
			if (onechar == rightChar)
			{
				return 0;
			}
			outstringstream << onechar;
		}
	}
	else
	{
		//For example you want to get the content between a pair of brackes {} or parentheses (),
		//where the left char (\") differs from the right one
		int bracketNum = 1;
		//read content until reaching the corresponding right parathese
		while (infile.get(onechar))
		{
			if (onechar == leftChar)
			{
				bracketNum++;
			}
			else if (onechar == rightChar)
			{
				bracketNum--;
			}
			if (bracketNum == 0)
			{
				return 0;
			}
			outstringstream << onechar;
		}
	}
	//empty output stringstream and return failure information if corresponding right was not found.
	outstringstream.str("");
	return 1;
}

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
)
{
	char onechar;
	outstringstream.str("");
	//find the position of the first left parathese
	while (instringstream.get(onechar))
	{
		if (onechar == leftChar) break;
	}

	if (leftChar == rightChar)
	{
		//For example you want to get the content between a pair of double quotes,
		//in which the left char (\") is the same as the right one
		while (instringstream.get(onechar))
		{
			if (onechar == rightChar)
			{
				return 0;
			}
			outstringstream << onechar;
		}
	}
	//For example you want to get the content between a pair of brackes {} or parentheses (),
	//where the left char (\") differs from the right one
	else
	{
		int bracketNum = 1;
		//read content until reaching the corresponding right parathese
		while (instringstream.get(onechar))
		{
			if (onechar == leftChar)
			{
				bracketNum++;
			}
			else if (onechar == rightChar)
			{
				bracketNum--;
			}
			if (bracketNum == 0)
			{
				return 0;
			}
			outstringstream << onechar;
		}
	}
	//empty output stringstream and return failure information if corresponding right was not found.
	outstringstream.str("");
	return 1;
}

bool getline(std::stringstream& sstream, std::string& outstring, char symbol)
{
	char onechar;
    std::string stTemp;
	while (sstream.get(onechar))
	{
		if (onechar == symbol)
		{
			outstring = stTemp;
			return true;
		}
		stTemp.push_back(onechar);
	}
	if (stTemp.size() > 0)
	{
		outstring = stTemp;
		return true;
	}
	else
	{
		return false;
	}
}
