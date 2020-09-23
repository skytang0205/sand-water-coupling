#include "Simulator.h"
#include "ArgsParser.h"

int main()
{
	ArgsParser args_parser;
	args_parser.Add_Argument<int>("search", 's');
	args_parser.Parse("main.exe -s 54");
	return 0;
}
