/* 
 * File:   MyTemplates.h
 * Author: aleksandra
 *
 * Created on 2 czerwiec 2012, 11:59
 */

#ifndef MYTEMPLATES_H
#define	MYTEMPLATES_H
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

template <class T>

T string2num(std::string& s,
std::ios_base& (*f)(std::ios_base&)) {
	
	size_t position = s.find("-X");//usuwam -X, ktore moze sie pojawic w plikach

	if (position != std::string::npos) {
		s.erase(position, 2);
	}
	
	
	T t;
	std::istringstream iss(s);
	if (!(iss >> f >> t).fail())
		return t;
	//

	return 0;
	//throw new Exception("String2num function: Wrong input type!"+s);
}

template <class T>
std::string num2string( T &t)
{
    std::stringstream iss;
    iss<<t;
    std::string result;
    iss>>result;
    return result;
}
#endif	/* MYTEMPLATES_H */

