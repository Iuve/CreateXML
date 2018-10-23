/* 
 * File:   PeriodicTable.h
 * Author: aleksandra
 * //source: http://www.cerim.pl/index.php?page=content/poradnik/pierwiastki.php
 * Created on 2 czerwiec 2012, 23:52
 */

#ifndef PERIODICTABLE_H
#define	PERIODICTABLE_H
#include<string>

class PeriodicTable {
public:
	static int GetAtomicNumber(std::string name);
	static std::string GetAtomicNameCap(int atomicNumber);
	static void UpToLow(std::string &str);
private:
    PeriodicTable();
    PeriodicTable(const PeriodicTable& orig);
	virtual ~PeriodicTable();

	static const std::string elements[];
	static const std::string elementsCap[];
	static const int nrOfElements = 118;

};


#endif	/* PERIODICTABLE_H */

