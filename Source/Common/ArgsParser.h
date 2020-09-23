#pragma once

#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <memory>
#include <any>
#include <typeinfo>

class ArgDataBase
{
protected:

	const std::type_info &type;
	std::string name;
	char flag;
	std::string desc;
	bool mandatory;
	bool set = false;

public:

	ArgDataBase(const std::type_info &_type, const std::string &_name, char _flag, const std::string _desc, bool _mandatory) :
		type(_type),
		name(_name),
		flag(_flag),
		desc(_desc),
		mandatory(_mandatory)
	{
		if (desc == "") desc = "<no description available>";
	}

	~ArgDataBase() = default;

	const std::type_info &Get_Type() const { return type; }
	const std::string &Get_Name() const { return name; }
	char Get_Flag() const { return flag; }
	const std::string &Get_Desc() const { return desc; }
	bool Is_Mandatory() const { return mandatory; }
	bool Is_Set() const { return set; }
	virtual std::string Get_Short_Desc() const = 0;
	virtual bool Parse_Value(std::istringstream iss) = 0;
	virtual std::any Get_Value() const = 0;
};

template <typename ArgType>
class ArgData : public ArgDataBase
{
protected:

	ArgType set_value;
	ArgType default_value;

public:

	ArgData(const std::string &_name, char _flag, const std::string _desc) :
		ArgDataBase(typeid(ArgType), _name, _flag, _desc, true),
		set_value(ArgType()),
		default_value(ArgType())
	{ }

	ArgData(const std::string &_name, char _flag, const std::string _desc, const ArgType &_default_value) :
		ArgDataBase(typeid(ArgType), _name, _flag, _desc, false),
		set_value(ArgType()),
		default_value(_default_value)
	{ }

	~ArgData() = default;

	virtual std::string Get_Short_Desc() const { return "--" + name + (mandatory ? "" : "=" + std::to_string(default_value)); }
	virtual bool Parse_Value(std::istringstream iss) { return set = true, static_cast<bool>(iss >> set_value); }
	virtual std::any Get_Value() const { return set ? set_value : (mandatory ? std::any() : default_value); }
};

template <> std::string ArgData<std::string>::Get_Short_Desc() const { return "--" + name + (mandatory ? "" : "=" + default_value); }

class ArgsParser
{
protected:
	
	std::string prog_name;
	std::vector<std::unique_ptr<ArgDataBase>> args;
	std::vector<std::string> extra_args;

	ArgDataBase *Find_Arg_by_Name(const std::string &name) const;
	ArgDataBase *Find_Arg_by_Flag(char flag) const;
	void Throw(const std::string &msg);

public:

	ArgsParser() = default;
	ArgsParser(const ArgsParser &rhs) = delete;
	ArgsParser operator=(const ArgsParser &rhs) = delete;
	~ArgsParser() = default;

	template <typename ArgType>
	void Add_Argument(const std::string &name, char flag = 0, const std::string &desc = "", const std::any &default_value = std::any())
	{
		if (default_value.has_value())
			args.push_back(std::unique_ptr<ArgDataBase>(new ArgData<ArgType>(name, flag, desc, std::any_cast<ArgType>(default_value))));
		else
			args.push_back(std::unique_ptr<ArgDataBase>(new ArgData<ArgType>(name, flag, desc)));
	}

	std::string Generate_Usage() const;

	void Parse(int argc, char *argv[]);
	void Parse(char *cmd_line);

	std::any Get_Value_by_Name(const std::string &name) const { return Find_Arg_by_Name(name)->Get_Value(); }
	std::any Get_Value_by_Flag(char flag) const { return Find_Arg_by_Flag(flag)->Get_Value(); }
};
