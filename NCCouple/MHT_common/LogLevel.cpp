
#include "../MHT_common/LogLevel.h"

#if defined(_BasePlatformWinddows_)
	LogLevel::SystemPlatform LogLevel::sp_sysPlatform = spWindows;
#else
	LogLevel::SystemPlatform LogLevel::sp_sysPlatform = spLinux;
#endif

LogLevel::LogLevel(const MsgLevel& lv, const std::string& stMsg)
	:
	m_lv(lv),
	m_stMsg(stMsg)
{

#if defined(_BasePlatformWinddows_)
	m_handle = GetStdHandle(STD_OUTPUT_HANDLE);

	switch (m_lv)
	{
	case mlError:
		SetConsoleTextAttribute(m_handle, FG_RED);
		break;
	case mlWarning:
		SetConsoleTextAttribute(m_handle, FG_YELLOW);
		break;
	case mlOK:
		SetConsoleTextAttribute(m_handle, FG_GREEN);
		break;
	default:
		SetConsoleTextAttribute(m_handle, FG_GRAY);
		break;
	}
#else
	switch (m_lv)
	{
	case mlError:
		printf("\033[01;31m");
		break;
	case mlWarning:
		printf("\033[01;33m");
		break;
	case mlOK:
		printf("\033[01;32m");
		break;
	default:
		printf("\033[0m");
		break;
	}
#endif
}

LogLevel::~LogLevel()
{
#if defined(_BasePlatformWinddows_)
	SetConsoleTextAttribute(m_handle, FG_GRAY);
#else
	printf("\033[0m");
#endif
}

