#include "pugixml.hpp"
#include "MyTemplates.hh"
#include "PeriodicTable.hh"
#include "FermiDistribution.hh"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

const double g_delayedNeutronPercentage(6.58 / 100.);
const double g_SnEnergy = 7053.;

using namespace std;

string toStringPrecision(double input,int n)
{
    stringstream stream;
    stream << fixed << setprecision(n) << input;
    return stream.str();
}

string RemoveExtension(string fileName)
{
	while(fileName.back() != '.')
		fileName.pop_back();
	fileName.pop_back();	
	
	return fileName;
}

void ZeroT12(double* T12)
{
	if(*T12 < pow(10., -50.))
		*T12 = 0.;
}

void CreateDefaultDecayXML(string defaultXmlFilename, string xmlMFilename, string xmlDFilename,
	int atomicMassM, int atomicNumberM, int atomicNumberD, int atomicNumberDN = 0, string xmlDNFilename = "_")
{
	ofstream outDefaultXML(defaultXmlFilename.c_str());
	if (!outDefaultXML.is_open())
//		throw IOException("Warning message: The file " + (string) defaultXmlFilename + " is not open!");
		cout << "Warning message: The file " + (string) defaultXmlFilename + " is not open!" << endl;
	
	if(atomicNumberDN == 0){
		outDefaultXML << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		outDefaultXML << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlMFilename + "\"/>" << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlDFilename + "\"/>" << endl;
		outDefaultXML << "<StartLevel AtomicNumber=\"" << atomicNumberM << "\" AtomicMass=\"" << atomicMassM << "\" Energy=\"0.0\"/>" << endl;
		outDefaultXML << "<StopLevel AtomicNumber=\"" << atomicNumberD << "\" AtomicMass=\"" << atomicMassM << "\" Energy=\"0.0\"/>" << endl;
		outDefaultXML << "<EventLength Value=\"0.5\" TimeUnit=\"US\"/>" << endl;
		outDefaultXML << "<CycleLength Value=\"15\" TimeUnit=\"M\"/>" ;
	}
	else{
		outDefaultXML << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		outDefaultXML << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlMFilename + "\"/>" << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlDFilename + "\"/>" << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlDNFilename + "\"/>" << endl;
		outDefaultXML << "<StartLevel AtomicNumber=\"" << atomicNumberM << "\" AtomicMass=\"" << atomicMassM << "\" Energy=\"0.0\"/>" << endl;
		outDefaultXML << "<StopLevel AtomicNumber=\"" << atomicNumberDN << "\" AtomicMass=\"" << atomicMassM - 1 << "\" Energy=\"0.0\"/>" << endl;
		outDefaultXML << "<EventLength Value=\"0.5\" TimeUnit=\"US\"/>" << endl;
		outDefaultXML << "<CycleLength Value=\"15\" TimeUnit=\"M\"/>" ;
	}
	
	outDefaultXML.close();
	cout << "Default file Decay.xml created." << endl;
}

int main(int argc,char** argv){
	
	int atomicNumber(0), atomicMass(0), whichBeta(0);
	double qVal_(0.);
	double maxBetaEnergy, energyLevel;
	double totalBetaIntensity(0.), betaIntensityUnderSn(0.), betaIntensityAboveSn(0.);
	double normalizedTotalBetaIntensity(0.);
	double averageGammaEnergy(0.), averageBetaEnergy(0.);
	double betaIntensity(0.);
	string atomicName;
	string ensDecayFilename, ensNeutronFilename;
	string xmlDNFilename;
	bool ifNeutronsInDecay = false;
	
	if(argc == 2)
	{
		ensDecayFilename = argv[1];
	}
	else if(argc == 3)
	{
		ensDecayFilename = argv[1];
		ensNeutronFilename = argv[2];
		ifNeutronsInDecay = true;
	}
	else
	{
		cout << "Wrong number of input files. Should be 1 for data without neutrons or 2 for data including neutron file." << endl;
		return 0;
	}
	
	pugi::xml_document docMother; // Mother nucleus
	pugi::xml_node nodeNuclideM = docMother.append_child("Nuclide");
	pugi::xml_document docDaughter; // Daughter nucleus
	pugi::xml_node nodeNuclideD = docDaughter.append_child("Nuclide");

	ifstream inDecay(ensDecayFilename.c_str());
	if (!inDecay.is_open())
//		throw IOException("Warning message: The file " + (string) ensDecayFilename + " is not open!");
		cout << "Warning message (first): The file " + (string) ensDecayFilename + " is not open!" << endl;
	
	string decayInfoSummaryFilename = RemoveExtension(ensDecayFilename) + "_summary.txt";
		
	ofstream outSummary(decayInfoSummaryFilename.c_str());
	if (!outSummary.is_open())
//		throw IOException("Warning message: The file " + (string) decayInfoSummaryFilename + " is not open!");
		cout << "Warning message: The file " + (string) decayInfoSummaryFilename + " is not open!" << endl;
	
	string additionalFilename = RemoveExtension(ensDecayFilename) + "_additional.txt";
	
	ofstream outAdditional(additionalFilename.c_str());
	if (!outAdditional.is_open())
//		throw IOException("Warning message: The file " + (string) additionalFilename + " is not open!");
		cout << "Warning message: The file " + (string) additionalFilename + " is not open!" << endl;
	
	string data, buff, stringParity;
	getline(inDecay, buff); //first line - some information about decay
	getline(inDecay, buff);

	//decay to this:
	data = buff.substr(0, 3);
	atomicMass = string2num <int> (data, std::dec);
	int atomicMassDaughter = atomicMass;

	data = buff.substr(3, 2);
	atomicName = data;
	atomicNumber = PeriodicTable::GetAtomicNumber(atomicName);
	string atomicNameDaughter = atomicName;
	
	nodeNuclideD.append_attribute("AtomicNumber") = atomicNumber;
	nodeNuclideD.append_attribute("AtomicMass") = atomicMass;
	nodeNuclideD.append_attribute("QBeta"); //nuclideD Qvalue?
	
	pugi::xml_node nodeLevelM;
	pugi::xml_node nodeLevelD;
	
	vector<double> gammaIntensitiesFromLevels;
	bool unlockGamma(false);
	
	while (!inDecay.eof())
	{
		getline(inDecay, buff);
		if (buff.size() < 80)
			break;
			
		if (buff[6]!='c')
		{	
			if(buff[7]=='B' && buff[5]==' ' && buff[6]==' ')
			{
				data = buff.substr(21, 8);
				betaIntensity = string2num <double>(data, std::dec);
				totalBetaIntensity += betaIntensity;
			}
			if(buff[7]=='E' && buff[5]==' ' && buff[6]==' ')
			{
				data = buff.substr(21, 8);
				betaIntensity = string2num <double>(data, std::dec);
				totalBetaIntensity += betaIntensity;
				data = buff.substr(31,8);
				double ECIntensity = string2num <double>(data, std::dec);
				totalBetaIntensity += ECIntensity;
			}
			if(buff[7]=='L' && buff[5]==' ' && buff[6]==' ') //each level
			{				
				unlockGamma = true;
				gammaIntensitiesFromLevels.push_back(0.);
			}
			if(buff[7]=='G' && buff[5]==' ' && buff[6]==' ' && unlockGamma)
			{
				data = buff.substr(21, 8);
				double gammaIntensity = string2num <double>(data, std::dec);
				if(gammaIntensity > 0)
					gammaIntensitiesFromLevels.back() += gammaIntensity;
				else
					gammaIntensitiesFromLevels.back() += 1;
			}
		}
	}
	inDecay.close();
//	cout << "First end of " + (string) ensDecayFilename + " file." << endl;
	inDecay.clear();
	inDecay.open(ensDecayFilename.c_str());
	if (!inDecay.is_open())
//		throw IOException("Warning message: The file " + (string) ensDecayFilename + " is not open!");
		cout << "Warning message: The file " + (string) ensDecayFilename + " is not open!" << endl;
		
	getline(inDecay, buff); //first line - some information about decay
	getline(inDecay, buff);
	
	int gammaIntensityPointer(-1);
	unlockGamma = false;
	
	outAdditional << "#Energy betaIntensity" << endl;
	
	while (!inDecay.eof())
	{
		getline(inDecay, buff);
		if (buff.size() < 80)
			break;
			
		if (buff[6]!='c')
		{
			if(buff[7]=='P')//parent level
			{
				data=buff.substr(9,10);
				double energy = string2num <double>(data, std::dec);
								
				//spin and parity shenanigans
				double spin;
				string parity;
				stringParity=buff.substr(21,18); 
				stringParity.erase(remove(stringParity.begin(), stringParity.end(), 
				' '), stringParity.end());
				if( (stringParity[0] != '(') && (stringParity.length() <= 4) )
				{
					if(stringParity[1] == '/')
					{
						double left = int(stringParity[0]) - '0';
						double right = int(stringParity[2]) - '0';
						spin = left / right;
						parity = stringParity[3];
					}
					else
					{
						spin = int(stringParity[0]) - '0';
						parity = stringParity[1];
					}	
				}
				
				
				data=buff.substr(39,10);
				istringstream iss(data);
				double T12;
				iss >> T12;
				ZeroT12(&T12);
				string timeUnit;
				iss >> timeUnit;
				
				//T12 uncertainty
				data=buff.substr(49,6);
				double dT12 = string2num <double>(data, std::dec);
				
				
				nodeLevelM = nodeNuclideM.append_child("Level");
				nodeLevelM.append_attribute("Energy").set_value(toStringPrecision(energy,2).c_str());;
				if( parity == "-" || parity == "+" )
				{
					nodeLevelM.append_attribute("Spin") = spin;
					nodeLevelM.append_attribute("Parity") = parity.c_str();
				}
				else
				{
					nodeLevelM.append_attribute("Spin");
					nodeLevelM.append_attribute("Parity");
				}
				if(!stringParity.empty())
					nodeLevelM.append_attribute("SpinParity") = stringParity.c_str();
				nodeLevelM.append_attribute("HalfLifeTime").set_value(toStringPrecision(T12,2).c_str());
				if(dT12 != 0.)
					nodeLevelM.append_attribute("d_T12") = dT12;
				nodeLevelM.append_attribute("TimeUnit") = timeUnit.c_str();
				nodeLevelM.append_attribute("Origin") = "Database";
				
				T12 = 0.; //apparently that variable can have the same address in further code
				// as other double T12 so it needs to be zeroed here
				
				data=buff.substr(64,10);
				qVal_ = string2num <double>(data, std::dec);
				if(qVal_<10)
					qVal_ = qVal_*1000;
				qVal_ += energy; //total energy for decay takes parent level en
				data=buff.substr(74,2);
				double dqVal_ = string2num <double>(data, std::dec);
				nodeNuclideM.append_attribute("QBeta") = qVal_;
				nodeNuclideM.append_attribute("d_QBeta") = dqVal_;
				outSummary << "Q val (parent level): " << qVal_ << endl;
			}
			
			if(buff[7]=='L' && buff[5]==' ' && buff[6]==' ') //each level
			{	
				gammaIntensityPointer++;
				unlockGamma = true;
				
				data = buff.substr(9, 10);
				energyLevel = string2num <double>(data, std::dec);
				
				//spin and parity shenanigans
				double spin;
				string parity;
				stringParity=buff.substr(21,18);
				stringParity.erase(remove(stringParity.begin(), stringParity.end(), 
				' '), stringParity.end()); 
				if( (stringParity[0] != '(') && (stringParity.length() <= 4) )
				{
					if(stringParity[1] == '/')
					{
						double left = int(stringParity[0]) - '0';
						double right = int(stringParity[2]) - '0';
						spin = left / right;
						parity = stringParity[3];
					}
					else
					{
						spin = int(stringParity[0]) - '0';
						parity = stringParity[1];
					}	
				}
	
				data=buff.substr(39,9);
				istringstream iss(data);
				double T12;
				iss >> T12;
				ZeroT12(&T12);
				string timeUnit;
				iss >> timeUnit;
				
				//T12 uncertainty
				data=buff.substr(49,6);
				double dT12 = string2num <double>(data, std::dec);
				
				maxBetaEnergy = qVal_ - energyLevel;
				
				nodeLevelD = nodeNuclideD.append_child("Level");
				nodeLevelD.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());;
				if( parity == "-" || parity == "+" )
				{
					nodeLevelD.append_attribute("Spin") = spin;
					nodeLevelD.append_attribute("Parity") = parity.c_str();
				}
				else
				{
					nodeLevelD.append_attribute("Spin");
					nodeLevelD.append_attribute("Parity");
				}
				if(!stringParity.empty())
					nodeLevelD.append_attribute("SpinParity") = stringParity.c_str();
				nodeLevelD.append_attribute("HalfLifeTime").set_value(toStringPrecision(T12,2).c_str());;
				nodeLevelD.append_attribute("TimeUnit") = timeUnit.c_str();	
				if(dT12 != 0.)
					nodeLevelD.append_attribute("d_T12") = dT12;
				nodeLevelD.append_attribute("Origin") = "Database";
				
				T12 = 0.; //apparently that variable can have the same address in further code
				// as other double T12 so it needs to be zeroed here	
			}
			
			if(buff[7]=='B' && buff[5]==' ' && buff[6]==' ')
			{
				whichBeta = -1;
				
				data = buff.substr(21, 8);
				betaIntensity = string2num <double>(data, std::dec);
				betaIntensity *= 100. / totalBetaIntensity;
				normalizedTotalBetaIntensity += betaIntensity; //just to check if it is 100 at the end
				
				if(ifNeutronsInDecay)
					betaIntensity *= (1. - g_delayedNeutronPercentage);
				
				data = buff.substr(29,2);
				double intensityError;
				intensityError = string2num <double>(data, std::dec);
				string logft;
				logft = buff.substr(43,6);
				
				if(energyLevel <= g_SnEnergy)
					betaIntensityUnderSn += betaIntensity;
				else
					betaIntensityAboveSn += betaIntensity;
				
				pugi::xml_node nodeTransition = nodeLevelM.append_child("Transition");
				nodeTransition.append_attribute("Type") = "B-";
				nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(maxBetaEnergy,2).c_str());
				nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(betaIntensity,6).c_str());
				if(intensityError != 0.)
					nodeTransition.append_attribute("d_Intensity") = intensityError;
				nodeTransition.append_attribute("Origin") = "Database";
				
				pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
				nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());
				nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
				nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
				
				averageGammaEnergy += energyLevel*betaIntensity/(100-g_delayedNeutronPercentage*100);
				
				FermiDistribution* fermiDist = new FermiDistribution(atomicNumber + whichBeta, maxBetaEnergy, whichBeta);
				double averageLvlBetaEnergy = fermiDist->GetAverageBetaEnergy();
				averageBetaEnergy += betaIntensity*averageLvlBetaEnergy/100.;
				
				//outAdditional << logft << endl;
				outAdditional << energyLevel << " " << betaIntensity << endl;
				//outAdditional << energyLevel << " " << betaIntensityUnderSn + betaIntensityAboveSn << endl;
			}
			
			if(buff[7]=='E'&& buff[5]==' ' && buff[6]==' ')
			{
				whichBeta = +1;
				
				data = buff.substr(21, 8);
				double betaPlusIntensity = string2num <double>(data, std::dec);
				data = buff.substr(29,2);
				double intensityError;
				intensityError = string2num <double>(data, std::dec);
				data=buff.substr(31,8);
				double ECIntensity = string2num <double>(data, std::dec);
				
				betaPlusIntensity *= 100. / totalBetaIntensity;
				ECIntensity *= 100. / totalBetaIntensity;
				
				if(betaPlusIntensity != 0)
					normalizedTotalBetaIntensity += betaPlusIntensity; //just to check if it is 100 at the end
				normalizedTotalBetaIntensity += ECIntensity; //just to check if it is 100 at the end
					
				data=buff.substr(64,10);
				//double totalIntensity = string2num <double>(data, std::dec);
				
				if(betaPlusIntensity != 0)
				{
					pugi::xml_node nodeTransition = nodeLevelM.append_child("Transition");
					nodeTransition.append_attribute("Type") = "B+";
					nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(maxBetaEnergy - 1022.,2).c_str());
					//nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(maxBetaEnergy,2).c_str());
					nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(betaPlusIntensity,6).c_str());
					if(intensityError != 0.)
						nodeTransition.append_attribute("d_Intensity") = intensityError;
					nodeTransition.append_attribute("Origin") = "Database";		
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
				}
				
				//FermiDistribution* fermiDist = new FermiDistribution(atomicNumber + whichBeta, maxBetaEnergy, whichBeta);
				//double averageLvlBetaEnergy = fermiDist->GetAverageBetaEnergy();
				//averageBetaEnergy += betaIntensity*averageLvlBetaEnergy/100.;
				
				if(ECIntensity != 0)
				{
					pugi::xml_node nodeTransition = nodeLevelM.append_child("Transition");
					nodeTransition.append_attribute("Type") = "EC";
					nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(maxBetaEnergy,2).c_str());
					nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(ECIntensity,6).c_str());
					nodeTransition.append_attribute("Origin") = "Database";	
					
					pugi::xml_node nodeConversion = nodeTransition.append_child("ElectronCapture");
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
				}
		
			}
			
			if(buff[7]=='E'&& (buff[5]=='S'||buff[5]=='2') && buff[6]==' ')
			{
				//data = buff.substr(9, 3);
				//if (data == "CK=")
				//{
				pugi::xml_node nodeTransition = nodeLevelM.last_child();
				if( !nodeTransition.child("ElectronCapture").empty() )
				{
					pugi::xml_node nodeCapture = nodeTransition.child("ElectronCapture");
					
					data = buff.substr(9, 74); //all line
					stringstream ss;
					ss << data;

					while (!ss.eof())
					{
						string singleCoef;
						getline(ss, singleCoef, '$');//"$" - separator
						stringstream sss; //data in file: CK=num or CK+=num
						sss << singleCoef;

						string type;
						double value;
						getline(sss, type, '=');
						sss >> value;
						
						if(type.back() == '+') type.pop_back();
						
						if(type == "EAV") continue;
						
						//double eCC = nodeTransition.attribute("ElectronConversionCoefficient").as_double();
						//ensdf format may be a reason to write some data 2 times, e.g. NC=
						pugi::xml_attribute checkDoubles = nodeCapture.last_attribute();
						if( checkDoubles.name() == type )
						{
							//eCC -= checkDoubles.as_double();
							checkDoubles.set_value(value);
						}
						else
							nodeCapture.append_attribute(type.c_str()) = value;
							
						/*
						if(eCC == 0. || eccAuxiliaryIndicator)
						{
							eccAuxiliaryIndicator = true;
							eCC += value;
							nodeTransition.attribute("ElectronConversionCoefficient").set_value(eCC);
						}
						*/
					}
				}
			}	
			
				if(buff[7]=='G' && buff[5]==' ' && buff[6]==' ' && unlockGamma)
				{
					
					data = buff.substr(9, 10);
					double gammaEnergy = string2num <double>(data, std::dec);
					data = buff.substr(21, 8);
					double gammaIntensity = string2num <double>(data, std::dec);
					data = buff.substr(29, 2);
					double intensityError;
						intensityError = string2num <double>(data, std::dec);
					data = buff.substr(55, 7);
					double electronConversionCoefficient = string2num <double>(data, std::dec);
					
					if(gammaIntensity > 0)
					{
						gammaIntensity *= 100. / gammaIntensitiesFromLevels.at(gammaIntensityPointer);
					}
					else
						gammaIntensity = 100. / gammaIntensitiesFromLevels.at(gammaIntensityPointer);

					pugi::xml_node nodeTransition = nodeLevelD.append_child("Transition");
					nodeTransition.append_attribute("Type") = "G";
					nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(gammaEnergy,2).c_str());
					nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(gammaIntensity,6).c_str());
					if(intensityError != 0.)
						nodeTransition.append_attribute("d_Intensity") = intensityError;
					nodeTransition.append_attribute("Origin") = "Database";
					if(electronConversionCoefficient != 0)
					{
						pugi::xml_node nodeConversion = nodeTransition.append_child("ElectronConversionCoefficient");
						nodeConversion.append_attribute("Total") = electronConversionCoefficient;
					}
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(energyLevel - gammaEnergy,2).c_str());
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
										
					//data = buff.substr(31,10);			
					//(*(gamma.end() - 1))->SetGammaMultiplicity(data);

				}
			
			
			if (buff[7] == 'G' && (buff[5] == 'S' || buff[5] == '2') && buff[6] == ' ') //electron conversion
			{
				data = buff.substr(9, 3);
				//cout << data << endl;
				//bool testKC = data == "KC=";
				//bool testNC = data == "NC=";
				if (data == "KC=" || data == "NC=" || data == "CC=")
				{
					pugi::xml_node nodeTransition = nodeLevelD.last_child();
					pugi::xml_node nodeConversion;
					
					if(data == "CC=")
						nodeConversion = nodeTransition.append_child("ElectronConversionCoefficient");
					else
						nodeConversion = nodeTransition.child("ElectronConversionCoefficient");
						
					//double eCC = nodeConversion.attribute("Total").as_double();
					
					data = buff.substr(9, 74); //all line
					stringstream ss;
					ss << data;

					while (!ss.eof())
					{
						string singleCoef;
						getline(ss, singleCoef, '$');//"$" - separator
						stringstream sss; //data in file: KC=num or KC+=num
						sss << singleCoef;

						string type;
						double value;
						getline(sss, type, '=');
						sss >> value;
						
						if(type.back() == '+') type.pop_back(); // type without '+'
						if(type == "CC") type = "Total";
						
						//ensdf format may be a reason to write some data 2 times, e.g. NC=
						pugi::xml_attribute checkDoubles = nodeConversion.last_attribute();
						if( checkDoubles.name() == type )
						{
							checkDoubles.set_value(value);
						}
						else
							nodeConversion.append_attribute(type.c_str()) = value;
					}
				}
			}			
		}
	}
	inDecay.close();
//	cout << "Second end of " + (string) ensDecayFilename + " file." << endl;
	outSummary << "normalizedTotalBetaIntensity: " << normalizedTotalBetaIntensity << endl;
	
	string atomicNameMother;
	nodeNuclideM.append_attribute("AtomicNumber") = atomicNumber + whichBeta;
	nodeNuclideM.append_attribute("AtomicMass") = atomicMass;
	atomicNameMother = PeriodicTable::GetAtomicNameCap(atomicNumber + whichBeta);
					
	if(ifNeutronsInDecay){
		int atomicNumberN(0), atomicMassN(0); //neutrons emmiter
		
		ifstream inNeutron(ensNeutronFilename.c_str());
		if (!inNeutron.is_open())
	//		throw IOException("Warning message: The file " + (string) ensNeutronFilename + " is not open!");
			cout << "Warning message: The file " + (string) ensNeutronFilename + " is not open!" << endl;
		
		getline(inNeutron, buff); //first line - some information about decay
		getline(inNeutron, buff);

		//decay to this:
		data = buff.substr(0, 3);
		atomicMassN = string2num <int> (data, std::dec);

		data = buff.substr(3, 2);
		atomicName = data;
		atomicNumberN = PeriodicTable::GetAtomicNumber(atomicName);
	
		pugi::xml_document docDaughterNeutron; // Daughter of the Daughter, neutrons emmiter
		pugi::xml_node nodeNuclideDN = docDaughterNeutron.append_child("Nuclide");
		nodeNuclideDN.append_attribute("AtomicNumber") = atomicNumberN;
		nodeNuclideDN.append_attribute("AtomicMass") = atomicMassN;
		nodeNuclideDN.append_attribute("QBeta"); //nuclideDN Qvalue?
		
		pugi::xml_node nodeLevelDN;
		
		double neutronQVal(0.);
		double totalNeutronIntensity(0.), normalizedTotalNeutronIntensity(0.);
		
		while (!inNeutron.eof()) 
		{
			getline(inNeutron, buff);
			if (buff.size() < 80)
				break;
			
			if (buff[6]!='c')
			{
				if(buff[7]=='D' && buff[8]=='N' && buff[5]==' ' && buff[6]==' '){
					data = buff.substr(21, 8);
					double neutronIntensity = string2num <double>(data, std::dec);
					totalNeutronIntensity += neutronIntensity;
								
				}
			}
		}
		inNeutron.close();
//		cout << "End of " + (string) ensNeutronFilename + " file - first." << endl;
		
		double intensityConverter = 100. * g_delayedNeutronPercentage / totalNeutronIntensity;	
		outSummary << "Total ENSDF beta intensity: " << totalBetaIntensity << ", total ENSDF neutron intensity: "
		<< totalNeutronIntensity << ", intensity converter: " << intensityConverter << endl;
		
		inNeutron.clear();
		inNeutron.open(ensNeutronFilename.c_str());

		if (!inNeutron.is_open())
	//		throw IOException("Warning message: The file " + (string) ensNeutronFilename + " is not open!");
			cout << "Warning message: The file " + (string) ensNeutronFilename + " is not open!" << endl;
		
		getline(inNeutron, buff);//first line - some information about decay
		getline(inNeutron, buff);
		
		while (!inNeutron.eof()) 
		{
			getline(inNeutron, buff);
			if (buff.size() < 80)
				break;
			
			if (buff[6]!='c')
			{
				if(buff[7]=='P')//parent level, different spins for 87Br ???
				{
					data=buff.substr(9,10);
					double energy = string2num <double>(data, std::dec);
					// create new level? check if it exists?
					data=buff.substr(64,10);
					neutronQVal = string2num <double>(data, std::dec);
				}
				
				if(buff[7]=='L' && buff[5]==' ' && buff[6]==' ') //each level
				{
					data = buff.substr(9, 10);
					energyLevel = string2num <double>(data, std::dec);
					
					//data=buff.substr(21,14);
					//(*(levels_.end() - 1))->SetLevelSpin(data);
		
					/*data=buff.substr(39,9); // no info, had to ?? use different method in general (function)
					istringstream iss(data);
					double T12;
					iss >> T12;
					string timeUnit;
					iss >> timeUnit;*/
					
					nodeLevelDN = nodeNuclideDN.append_child("Level");
					nodeLevelDN.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());
					nodeLevelDN.append_attribute("Spin"); // to do in the future
					nodeLevelDN.append_attribute("Parity"); // to do in the future
					nodeLevelDN.append_attribute("HalfLifeTime");
					nodeLevelDN.append_attribute("TimeUnit");
					nodeLevelDN.append_attribute("Origin") = "Database";
				}
				
				if(buff[7]=='D' && buff[8]=='N' && buff[5]==' ' && buff[6]==' '){
					data = buff.substr(9, 10);
					double neutronEnergy = string2num <double>(data, std::dec);
				
					data = buff.substr(21, 8);
					double neutronIntensity = string2num <double>(data, std::dec);
				
					data=buff.substr(31,8);
					double neutronEnergyLevel= string2num <double>(data, std::dec);
					
					maxBetaEnergy = neutronQVal - neutronEnergy;
					neutronIntensity *= intensityConverter;
					normalizedTotalNeutronIntensity += neutronIntensity;
					
					nodeLevelD = nodeNuclideD.append_child("Level");
					nodeLevelD.append_attribute("Energy").set_value(toStringPrecision(neutronEnergyLevel,2).c_str());
					nodeLevelD.append_attribute("Spin"); // to do in the future
					nodeLevelD.append_attribute("Parity"); // to do in the future
					nodeLevelD.append_attribute("HalfLifeTime");
					nodeLevelD.append_attribute("TimeUnit");
					nodeLevelD.append_attribute("Origin") = "Database";
					
					pugi::xml_node nodeTransition = nodeLevelD.append_child("Transition");
					nodeTransition.append_attribute("Type") = "N";
					nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(neutronEnergy,2).c_str());
					nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(neutronIntensity,6).c_str());
					nodeTransition.append_attribute("Origin") = "Database";
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(energyLevel,2).c_str());
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumberN;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMassN;
					
					nodeTransition = nodeLevelM.append_child("Transition");
					nodeTransition.append_attribute("Type") = "B-";
					nodeTransition.append_attribute("TransitionQValue").set_value(toStringPrecision(maxBetaEnergy,2).c_str());
					nodeTransition.append_attribute("Intensity").set_value(toStringPrecision(neutronIntensity,6).c_str());
					nodeTransition.append_attribute("Origin") = "Database";
					
					nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy").set_value(toStringPrecision(neutronEnergyLevel,2).c_str());
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
					
					FermiDistribution* fermiDist = new FermiDistribution(atomicNumber + whichBeta, maxBetaEnergy, whichBeta);
					double averageLvlBetaEnergy = fermiDist->GetAverageBetaEnergy();
					averageBetaEnergy += neutronIntensity*averageLvlBetaEnergy/100.;
					
					outAdditional << neutronEnergyLevel << " " << neutronIntensity << endl;
					//outAdditional << "/gun/energy " << neutronEnergy << " keV" << endl;
					//outAdditional << "/run/beamOn " << floor(neutronIntensity*1000000) << endl;
				}
			}
		}
		inNeutron.close();
//		cout << "End of " + (string) ensNeutronFilename + " file - second." << endl;
		
		string atomicNameDaughterNeutron = PeriodicTable::GetAtomicNameCap(atomicNumberN);
		xmlDNFilename = num2string(atomicMassN) + atomicNameDaughterNeutron + ".xml";
		std::cout << "Saving result: " << docDaughterNeutron.save_file(xmlDNFilename.c_str()) << std::endl;
		
		outSummary << "normalizedTotalNeutronIntensity: " << normalizedTotalNeutronIntensity << endl;
	}
	
	string xmlDFilename = num2string(atomicMassDaughter) + atomicNameDaughter + ".xml";
	std::cout << "Saving result: " << docDaughter.save_file(xmlDFilename.c_str()) << std::endl;
	string xmlMFilename = num2string(atomicMassDaughter) + atomicNameMother + ".xml";
	std::cout << "Saving result: " << docMother.save_file(xmlMFilename.c_str()) << std::endl;
	
	string defaultXmlFilename = "Decay.xml";
	if(ifNeutronsInDecay)
		CreateDefaultDecayXML(defaultXmlFilename, xmlMFilename, xmlDFilename, atomicMassDaughter, 
		atomicNumber + whichBeta, atomicNumber, atomicNumber, xmlDNFilename);
	else
		CreateDefaultDecayXML(defaultXmlFilename, xmlMFilename, xmlDFilename, atomicMassDaughter, 
		atomicNumber + whichBeta, atomicNumber);
	
	outSummary << endl << "beta-gamma intensity smaller than Sn: " << betaIntensityUnderSn << endl;
	outSummary << "beta-gamma intensity bigger than Sn: " << betaIntensityAboveSn << endl;
	outSummary << "summed: " << betaIntensityAboveSn + betaIntensityUnderSn << endl;
	outSummary << "ratio in percents (including delayed neutrons levels): " << betaIntensityAboveSn * (1. - g_delayedNeutronPercentage) << endl;
	outSummary << "Average gamma intensity: " << averageGammaEnergy << endl;
	outSummary << "Average beta intensity: " << averageBetaEnergy << endl;
	outSummary.close();
	
	outAdditional.close();

return 0;
}
