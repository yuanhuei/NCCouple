/*---------------------------------------------------------------------------*\
Class:
	LogLevel
File Name:
	LogLevel.h
Description:
	A basic class for change screen color.
	1.Change screen output information color;
	2.This class allow user to compile not only in Windows, but also in Linux
	3.Color : RED,GREEN,BLUE,GRAY,YELLOW

	Author:		ShuaiZhang
	Revisor:	ZhiHaoJia
	Modified Date: Long ago
\*---------------------------------------------------------------------------*/

#pragma   once

#ifndef _LogLevel_
#define _LogLevel_

#include <iostream>
#include <string>
#include "stdio.h"

#include "../MHT_common/Configuration.h"

#if defined(_BasePlatformWinddows_)
#include <Windows.h>

#define	FG_RED			FOREGROUND_INTENSITY | FOREGROUND_RED
#define	FG_GREEN		FOREGROUND_INTENSITY | FOREGROUND_GREEN
#define	FG_BLUE			FOREGROUND_INTENSITY | FOREGROUND_BLUE
#define	FG_GRAY			FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE
#define FG_YELLOW		FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN

#endif

class	LogLevel
{
public:

	enum	MsgLevel
	{
		mlError,
		mlWarning,
		mlOK,
		mlInfo	
	};

private:

	enum SystemPlatform
	{
		spWindows,
		spLinux,
		spNoType
	};

	static SystemPlatform sp_sysPlatform;
	
#if defined(_BasePlatformWinddows_)

	HANDLE		m_handle;

#endif

	MsgLevel	m_lv;

	std::string	m_stMsg;

public:
	LogLevel(const MsgLevel& lv, const std::string& stMsg);
		

	~LogLevel();


	// "<<" overload
	inline friend std::ostream& operator << (std::ostream& out, const LogLevel& rhs)
	{
		out << rhs.m_stMsg << "\t";
		return out;
	}
};



#endif
