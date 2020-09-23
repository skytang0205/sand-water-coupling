#include "ArgsParser.h"

#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>

ArgDataBase *ArgsParser::Find_Arg_by_Name(const std::string &name) const
{
	for (const auto &arg : args)
		if (name == arg->Get_Name()) return arg.get();
	return nullptr;
}

ArgDataBase *ArgsParser::Find_Arg_by_Flag(char flag) const
{
	for (const auto &arg : args)
		if (flag == arg->Get_Flag()) return arg.get();
	return nullptr;
}

void ArgsParser::Throw(const std::string &msg)
{
	std::cerr << "Error: [ArgsParser] " << msg << "." << std::endl << Generate_Usage();
	std::exit(1);
}

std::string ArgsParser::Generate_Usage() const
{
	std::ostringstream oss;
	// Usage header.
	oss << "Usage: " << prog_name;
	for (const auto &arg : args)
		if (arg->Is_Mandatory())
			oss << " " << arg->Get_Short_Desc();
	for (const auto &arg : args)
		if (!arg->Is_Mandatory())
			oss << " [" << arg->Get_Short_Desc() << "]";
	oss << std::endl << "Options:" << std::endl;
	// Usage body.
	size_t maxWidth = 0;
	for (const auto &arg : args) maxWidth = std::max(maxWidth, arg->Get_Name().length());
	for (const auto &arg : args)
	{
		if (arg->Get_Flag()) oss << "  -" << arg->Get_Flag() << ", ";
		else oss << "      ";
		oss << "--" << arg->Get_Name();
		oss << std::string(maxWidth + 4 - arg->Get_Name().length(), ' ');
		oss << arg->Get_Desc() << std::endl;
	}
	// Return usage as string.
	return oss.str();
}

void ArgsParser::Parse(int argc, char *argv[])
{
	Add_Argument<bool>("help", '?', "print this message", false);
	prog_name = argv[0];
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] == '-')
		{
			ArgDataBase *arg = nullptr;
			if (argv[i][1] == '-') arg = Find_Arg_by_Name(argv[i] + 2);
			else arg = Find_Arg_by_Flag(argv[i][1]);
			if (arg)
			{
				// Handle options without value, assuming default_value == false.
				if (arg->Get_Type() == typeid(bool) && !arg->Is_Mandatory()) arg->Parse_Value(std::istringstream("1"));
				else
				{
					if (i + 1 == argc) Throw(std::string("Missing value for option ") + argv[i]);
					else if (!arg->Parse_Value(std::istringstream(argv[i + 1])))
						Throw(std::string("Invalid value ") + argv[i + 1] + " for option " + argv[i]);
					i++;
				}
			}
			else Throw(std::string("Invalid option ") + argv[i]);
		}
		else extra_args.push_back(argv[i]);
	}
	if (std::any_cast<bool>(Get_Value_by_Name("help"))) {
		std::cerr << Generate_Usage();
		std::exit(0);
	}
	for (const auto &arg : args)
		if (arg->Is_Mandatory() && !arg->Is_Set()) Throw(std::string("Unassigned argument ") + arg->Get_Name());
}
