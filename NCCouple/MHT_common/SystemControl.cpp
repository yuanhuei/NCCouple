
#include "../MHT_common/SystemControl.h"

//print error information and stop
void FatalError(const std::string info)
{
	std::cout << LogLevel(LogLevel::mlError, std::string("FATAL ERROR:\n") + info) << std::endl;
	system("pause");
	exit(1);
}

//print waring information and pause
void WarningPause(const std::string info)
{
	std::cout << LogLevel(LogLevel::mlWarning, std::string("WARNING:\n") + info) << std::endl;
	std::cout << "Do you want to continue?[y/n]";
	std::string yesOrNo;
	while (true)
	{
		std::cin >> yesOrNo;
		if ("y" == yesOrNo)
		{
			break;
		}
		else if ("n" == yesOrNo)
		{
			exit(1);
		}
		else
		{
			std::cout << "Please input \'y\' or \'n\': " << std::endl;
		}
	}
	return;
}

//print waring information but continue
void WarningContinue(const std::string info)
{
	std::cout << LogLevel(LogLevel::mlWarning, std::string("WARNING:\n") + info) << std::endl;
	return;
}