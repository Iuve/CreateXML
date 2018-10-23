/* 
 * File:   PeriodicTable.cpp
 * Author: aleksandra
 * 
 * Created on 2 czerwiec 2012, 23:52
 */

//#include "Exception.hh"


#include "PeriodicTable.hh"
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>

const std::string PeriodicTable::elements[] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F","Ne", "Na", "Mg",
"Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
"As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
"I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", 
"Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", 
"Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cp", 
"Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"};

const std::string PeriodicTable::elementsCap[] = {"H", "HE", "LI", "BE", "B", "C", "N", "O", "F","NE", "NA", "MG",
"AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE",
"AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE",
"I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", 
"TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", 
"NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CP", 
"UUT", "UUQ", "UUP", "UUH", "UUS", "UUO"};

void PeriodicTable::UpToLow(std::string &str)
{
    int i=1;
    while(str[i]!='\0')
    {
	    if (str[i] >= 0x41 && str[i] <= 0x5A)
		str[i] = str[i] + 0x20;
    i++;
    }
}

int PeriodicTable::GetAtomicNumber(std::string atomicName)
{
	int atomicNumber(0);
	UpToLow(atomicName);
	std::string symbol;
	std::stringstream sym;
	sym << atomicName;
	sym >> symbol;
	for (int i = 0; i < nrOfElements; i++)
		if (elements[i] == symbol)
		{
			atomicNumber = i + 1;
		}

	return atomicNumber;
}

std::string PeriodicTable::GetAtomicNameCap(int atomicNumber)
{
	if (atomicNumber < 1 || atomicNumber > nrOfElements)
//		throw new Exception("PeriodicTable::getAtomicNameCap - invalid agrument");
		std::cout << "PeriodicTable::getAtomicNameCap - invalid agrument" << std::endl;
	return elementsCap[atomicNumber - 1];
}
