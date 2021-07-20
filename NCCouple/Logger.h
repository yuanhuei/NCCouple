#ifndef LOGGER_HEADER
#define LOGGER_HEADER

#include <iostream>
#include <string>
#include <memory>

template <class... Args>
static std::string FormatStr(const std::string& fmtStr, Args... args) {
	auto size_buf = std::snprintf(nullptr, 0, fmtStr.c_str(), args ...) + 1;
	std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

	if (!buf)
		return std::string("");

	std::snprintf(buf.get(), size_buf, fmtStr.c_str(), args ...);
	return std::string(buf.get(), buf.get() + size_buf - 1);
}

class Logger
{
public:
	static void LogInfo(std::string info);
	static void LogWarn(std::string warn);
	static void LogError(std::string error);

private:
	Logger() = delete;
};

#endif