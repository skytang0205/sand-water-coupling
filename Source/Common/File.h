#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <windows.h>

namespace File
{
	// ***** Binary IO *****

	template <class Type>
	inline void Write_Binary(std::ostream &output, const Type &data) { output.write(reinterpret_cast<const char *>(&data), sizeof(Type)); }

	template <class Type>
	inline void Read_Binary(std::istream &input, Type &data) { input.read(reinterpret_cast<char *>(&data), sizeof(Type)); }

	template <class Type>
	inline void Write_Binary_Array(std::ostream &output, const Type *array, const std::size_t n) { if (n > 0) output.write(reinterpret_cast<const char *>(array), n * sizeof(Type)); }

	template <class Type>
	inline void Read_Binary_Array(std::istream &input, Type *array, const std::size_t n) { if (n > 0) input.read(reinterpret_cast<char *>(array), n * sizeof(Type)); }

	// ***** File operations *****

	inline bool Directory_Exists(const char *dir_name)
	{
		auto attr = GetFileAttributesA(dir_name);
		return ((attr != -1) && (attr & FILE_ATTRIBUTE_DIRECTORY));
	}

	inline void Create_Directory(const std::string &dir_name)
	{
		if (!Directory_Exists(dir_name.c_str())) {
			std::string::size_type pos = 0;
			do {
				pos = dir_name.find_first_of("\\/", pos + 1);
				if (!Directory_Exists(dir_name.substr(0, pos).c_str())) {
					if (!CreateDirectoryA(dir_name.substr(0, pos).c_str(), nullptr) && ERROR_ALREADY_EXISTS != GetLastError()) {
						std::cerr << "Error: [File] Create directory " << dir_name << " failed." << std::endl;
						std::exit(1);
					}
				}
			} while (pos != std::string::npos);
		}
	}

	inline bool File_Exists(const std::string &file_name)
	{
		return std::ifstream(file_name.c_str()).good();
	}

	inline std::string File_Extension_Name(const std::string &file_name)
	{
		auto p = file_name.rfind('.');
		if (p != std::string::npos) return file_name.substr(p + 1);
		else return "";
	}
};
