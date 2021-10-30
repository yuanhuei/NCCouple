#include "Logger.h"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/basic_file_sink.h>

static auto my_logger = spdlog::basic_logger_mt("basic_logger", "log.txt");
void Logger::LogInfo(std::string info) {
	spdlog::info(info);
	return;
}

void Logger::LogWarn(std::string warn) {
	spdlog::warn(warn);
	return;
}

void Logger::LogError(std::string err) {
	spdlog::error(err);
	exit(1);

	return;
}

void Logger::LogInfotoFile(std::string info)
{

	my_logger->info(info);
}