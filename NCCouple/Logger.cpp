#include "Logger.h"
#include "spdlog/spdlog.h"

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