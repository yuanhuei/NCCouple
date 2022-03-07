/*---------------------------------------------------------------------------*\
File Name:
	SystemControl.h

Description:
	Display warning and error informations and
	pause or end the program correspondingly.

	Author:		Kong Ling
	Date: 2016-11-23
\*---------------------------------------------------------------------------*/
#pragma once

#ifndef _SystemControl_
#define _SystemControl_

#include <iostream>
#include <string>
#include <stdlib.h>

#include "../MHT_common/LogLevel.h"


//print error information and stop
void FatalError(const std::string info);

//print waring information and pause
void WarningPause(const std::string info);

//print waring information but continue
void WarningContinue(const std::string info);

#endif
