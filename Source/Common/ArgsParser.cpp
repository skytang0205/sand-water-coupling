#include "ArgsParser.h"

#include <algorithm>
#include <iostream>

#include <cstddef>
#include <cstdlib>
#include <cstring>

namespace PhysX {

std::string ArgsParser::generateUsage() const
{
	std::string usage;
	// Usage header.
	usage += fmt::format("Usage: {}", progName);
	for (const auto &arg : args)
		if (arg->isMandatory())
			usage += fmt::format(" {}", arg->getShortDesc());
	for (const auto &arg : args)
		if (!arg->isMandatory())
			usage += fmt::format(" [{}]", arg->getShortDesc());
	usage += "\nOptions:\n";
	// Usage body.
	std::size_t maxWidth = 0;
	for (const auto &arg : args) maxWidth = std::max(maxWidth, arg->getName().length());
	for (const auto &arg : args) {
		if (arg->getFlag()) usage += fmt::format("  -{}, ", arg->getFlag());
		else usage += fmt::format("{:^6}", "");
		usage += fmt::format("--{:<{}}{}\n", arg->getName(), maxWidth + 4, arg->getDesc());
	}
	// Return usage.
	return usage;
}

void ArgsParser::parse(const int argc, char *const argv[])
{
	if (argc == 0) reportError("Empty command line");

	addArgument<bool>("help", '?', "print this message", false);
	progName = argv[0];
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			ArgDataBase *arg = nullptr;
			if (argv[i][1] == '-') arg = findArgByName(argv[i] + 2);
			else arg = findArgByFlag(argv[i][1]);
			if (arg) {
				// Handle options without value, assuming default_value == false.
				if (arg->getType() == typeid(bool) && !arg->isMandatory()) arg->parseValue("1");
				else {
					if (i + 1 == argc) reportError(fmt::format("Missing value for option {}", argv[i]));
					else if (!arg->parseValue(argv[i + 1]))
						reportError(fmt::format("Invalid value {} for option {}", argv[i + 1], argv[i]));
					i++;
				}
			}
			else reportError(fmt::format("Invalid option {}", argv[i]));
		}
		else extraArgs.push_back(argv[i]);
	}
	if (std::any_cast<bool>(getValueByName("help"))) {
		std::cerr << generateUsage();
		std::exit(0);
	}
	for (const auto &arg : args)
		if (arg->isMandatory() && !arg->isSet()) reportError(fmt::format("Unassigned argument {}", arg->getName()));
}

void ArgsParser::parse(const char *cmdLine)
{
	const std::size_t len = std::strlen(cmdLine);
	char *buffer = new char[len + 1];
	std::memcpy(buffer, cmdLine, len + 1);

	int argc = 0;
	for (std::size_t i = 0; i < len; i++) {
		if (!std::isspace(buffer[i])) {
			argc++;
			while (i + 1 < len && !std::isspace(buffer[i + 1])) i++;
		}
		else buffer[i] = 0;
	}
	if (argc == 0) reportError("Empty command line");

	char **argv = new char *[argc];
	argc = 0;
	for (std::size_t i = 0; i < len; i++) {
		if (buffer[i]) {
			argv[argc++] = buffer + i;
			while (i + 1 < len && buffer[i + 1]) i++;
		}
	}

	parse(argc, argv);

	delete[] argv;
	delete[] buffer;
}

ArgDataBase *ArgsParser::findArgByName(const std::string &name) const
{
	for (const auto &arg : args)
		if (name == arg->getName()) return arg.get();
	return nullptr;
}

ArgDataBase *ArgsParser::findArgByFlag(const char flag) const
{
	for (const auto &arg : args)
		if (flag == arg->getFlag()) return arg.get();
	return nullptr;
}

void ArgsParser::reportError(const std::string &msg)
{
	std::cerr << fmt::format("Error: [ArgsParser] {}.\n{}", msg, generateUsage());
	std::exit(1);
}

} // namespace PhysX
