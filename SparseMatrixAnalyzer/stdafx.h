#pragma once

#include<vector>
#include<array>
#include<unordered_set>
#include<set>
#include<functional>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<sys/stat.h>
#include<string>
#include<algorithm>
#include<memory>
#include<unordered_map>
#include<unordered_map>
#include <locale>
#include <clocale>
#include <cfloat>
#include <cmath>

using namespace std;

#define RETCODE_OK 0
#define RETCODE_NOMEM 1
#define RETCODE_NOFILE 2
#define RETCODE_OUTOFRANGE 3
#define RETCODE_ERROR 10

#define PI 3.1415926535897932

inline ostream& operator < (ostream& file, const double& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}
inline ostream& operator < (ostream& file, const int& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}

inline istream& operator > (istream& file, double& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

inline istream& operator > (istream& file, int& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}
