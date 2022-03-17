#include "Logger.h"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/basic_file_sink.h>
extern int g_iMpiID;
extern int g_iNumProcs;
static auto my_logger = spdlog::basic_logger_mt("basic_logger","log.txt");

void Logger::LogInfo(std::string info,bool bOnlyInLogFile) {
	
	if (g_iMpiID == 0 && g_iNumProcs > 1)
		info = "O Process Log Informataion:" + info;
	if(bOnlyInLogFile==false)
		spdlog::info(info);
	my_logger->info(info);
	return;
}

void Logger::LogWarn(std::string warn) {
	if (g_iMpiID == 0 && g_iNumProcs > 1)
		warn = "O Process Log Informataion:" + warn;
	spdlog::warn(warn);

}

void Logger::LogError(std::string err) {
	if (g_iMpiID == 0 && g_iNumProcs > 1)
		err = "O Process Log Informataion:" + err;
	spdlog::error(err);

	exit(1);
	return;
}

void Logger::LogInfotoFile(std::string info)
{

	my_logger->info(info);
}