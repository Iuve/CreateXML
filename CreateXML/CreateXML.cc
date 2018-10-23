#include "pugixml.hpp"
#include "MyTemplates.hh"
#include "PeriodicTable.hh"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

const double g_delayedNeutronPercentage(2.6 / 100.);
const double g_SnEnergy = 5515.;

using namespace std;

void ZeroT12(double* T12)
{
	if(*T12 < pow(10., -50.))
		*T12 = 0.;
}

void CreateDefaultDecayXML(string defaultXmlFilename, string xmlMFilename, string xmlDFilename,
	int atomicMassM, int atomicNumberM, int atomicNumberD, int atomicNumberDN = 0, string xmlDNFilename = 0L)
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
		outDefaultXML << "<StopLevel AtomicNumber=\"" << atomicNumberD << "\" AtomicMass=\"" << atomicMassM << "\" Energy=\"0.0\"/>" ;
	}
	else{
		outDefaultXML << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		outDefaultXML << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlMFilename + "\"/>" << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlDFilename + "\"/>" << endl;
		outDefaultXML << "<NuclideFile FileName=\"" + xmlDNFilename + "\"/>" << endl;
		outDefaultXML << "<StartLevel AtomicNumber=\"" << atomicNumberM << "\" AtomicMass=\"" << atomicMassM << "\" Energy=\"0.0\"/>" << endl;
		outDefaultXML << "<StopLevel AtomicNumber=\"" << atomicNumberDN << "\" AtomicMass=\"" << atomicMassM - 1 << "\" Energy=\"0.0\"/>" ;
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
	string atomicName;
	string ensDecayFilename, ensNeutronFilename;
	//string additionalFilename = "someInformation.txt";
	string xmlDNFilename;
	
	if(argc == 2)
	{
		ensDecayFilename = argv[1];
	}
	else if(argc == 3)
	{
		ensDecayFilename = argv[1];
		ensNeutronFilename = argv[2];
		#define NEUTRON
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
		string decayInfoSummaryFilename = ensDecayFilename + "_summary";
		
	ofstream outSummary(decayInfoSummaryFilename.c_str());
	if (!outSummary.is_open())
//		throw IOException("Warning message: The file " + (string) decayInfoSummaryFilename + " is not open!");
		cout << "Warning message: The file " + (string) decayInfoSummaryFilename + " is not open!" << endl;
			
//	ofstream outAdditional(additionalFilename.c_str());
//	if (!outAdditional.is_open())
//		throw IOException("Warning message: The file " + (string) additionalFilename + " is not open!");
//		cout << "Warning message: The file " + (string) additionalFilename + " is not open!" << endl;
	
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
	
	while (!inDecay.eof())
	{
		getline(inDecay, buff);
		if (buff.size() < 80)
			break;
			
		if (buff[6]!='c')
		{	
			if(buff[7]=='B' && buff[5]==' ' && buff[6]==' ')
			{
				whichBeta = -1;
				
				data = buff.substr(21, 8);
				double betaIntensity = string2num <double>(data, std::dec);
				totalBetaIntensity += betaIntensity;
			}
		}
	}
	inDecay.close();
	//cout << "First end of " + (string) ensDecayFilename + " file." << endl;
	inDecay.clear();
	inDecay.open(ensDecayFilename.c_str());
	if (!inDecay.is_open())
//		throw IOException("Warning message: The file " + (string) ensDecayFilename + " is not open!");
		cout << "Warning message: The file " + (string) ensDecayFilename + " is not open!" << endl;
		
	getline(inDecay, buff); //first line - some information about decay
	getline(inDecay, buff);
	
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
								
				//data=buff.substr(21,18); 
				//parentLevel_->SetSpinAndParity(data);
				
				data=buff.substr(39,9);
				istringstream iss(data);
				double T12;
				iss >> T12;
				ZeroT12(&T12);
				string timeUnit;
				iss >> timeUnit;
				
				
				nodeLevelM = nodeNuclideM.append_child("Level");
				nodeLevelM.append_attribute("Energy") = energy;
				nodeLevelM.append_attribute("Spin"); // to do in the future
				nodeLevelM.append_attribute("Parity"); // to do in the future
				nodeLevelM.append_attribute("HalfLifeTime") = T12;
				nodeLevelM.append_attribute("TimeUnit") = timeUnit.c_str();
				
				data=buff.substr(64,10);
				qVal_ = string2num <double>(data, std::dec);
				if(qVal_<10)
					qVal_ = qVal_*1000;
				qVal_ += energy; //total energy for decay takes parent level en
				nodeNuclideM.append_attribute("QBeta") = qVal_;
				outSummary << "Q val (parent level): " << qVal_ << endl;
			}
			
			if(buff[7]=='L' && buff[5]==' ' && buff[6]==' ') //each level
			{				
				data = buff.substr(9, 10);
				energyLevel = string2num <double>(data, std::dec);
				
				//data=buff.substr(21,14);
				//(*(levels_.end() - 1))->SetLevelSpin(data);
				//energyError = string2num <double>(data, std::dec);
				
				//stringParity = buff.substr(21,14);
	
				data=buff.substr(39,9);
				istringstream iss(data);
				double T12;
				iss >> T12;
				ZeroT12(&T12);
				string timeUnit;
				iss >> timeUnit;
				
				maxBetaEnergy = qVal_ - energyLevel;
				
				nodeLevelD = nodeNuclideD.append_child("Level");
				nodeLevelD.append_attribute("Energy") = energyLevel;
				nodeLevelD.append_attribute("Spin"); // to do in the future
				nodeLevelD.append_attribute("Parity"); // to do in the future
				nodeLevelD.append_attribute("HalfLifeTime") = T12;
				nodeLevelD.append_attribute("TimeUnit") = timeUnit.c_str();		
			}
			
			if(buff[7]=='B' && buff[5]==' ' && buff[6]==' ')
			{
				whichBeta = -1;
				
				data = buff.substr(21, 8);
				double betaIntensity = string2num <double>(data, std::dec);
				betaIntensity *= 100. / totalBetaIntensity;
				normalizedTotalBetaIntensity += betaIntensity; //just to check if it is 100 at the end
				#ifdef NEUTRON
					betaIntensity *= (1. - g_delayedNeutronPercentage);
				#endif	
				
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
				nodeTransition.append_attribute("TransitionQValue") = maxBetaEnergy;
				nodeTransition.append_attribute("Intensity") = betaIntensity;
				
				pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
				nodeTargetLevel.append_attribute("Energy") = energyLevel;
				nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
				nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;

			}
			/*
			if(buff[7]=='E'&& buff[5]==' ' && buff[6]==' ' && levels_.size() > 0)
			{
				
				data = buff.substr(21, 8);
				double betaPlusIntensity = string2num <double>(data, std::dec);
				data=buff.substr(31,8);
				double ECIntensity= string2num <double>(data, std::dec);
				data=buff.substr(64,10);
				//double totalIntensity = string2num <double>(data, std::dec);
				(*(levels_.end() - 1))->SetBetaFeedingFunction(betaPlusIntensity+ECIntensity);
				(*(levels_.end() - 1))->SetBetaPlus(betaPlusIntensity, ECIntensity);
				//data=buff.substr(41,8);
				//double betaLogft = string2num <double>(data, std::dec);
				//(*(levels_.end() - 1))->SetBetaLogft(betaLogft);

		
			}
			
			if(buff[7]=='E'&& (buff[5]=='S'||buff[5]=='2') && buff[6]==' ')
			{
				data = buff.substr(9, 3);
				if (data == "CK=")
				{
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
						(*(levels_.end() - 1))->GetBeta()->SetECCoef(type, value);
					}

				}
			}	*/
			
				if(buff[7]=='G' && buff[5]==' ' && buff[6]==' ')
				{
					data = buff.substr(9, 10);
					double gammaEnergy = string2num <double>(data, std::dec);
					data = buff.substr(21, 8);
					double gammaIntensity = string2num <double>(data, std::dec);
					data = buff.substr(29, 2);
					data = buff.substr(55, 7);
					double electronConversionCoefficient = string2num <double>(data, std::dec);

					pugi::xml_node nodeTransition = nodeLevelD.append_child("Transition");
					nodeTransition.append_attribute("Type") = "G";
					nodeTransition.append_attribute("TransitionQValue") = gammaEnergy;
					nodeTransition.append_attribute("Intensity") = gammaIntensity;	
					nodeTransition.append_attribute("ElectronConversionCoefficient") = electronConversionCoefficient;	
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy") = energyLevel - gammaEnergy;
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
					
					//data = buff.substr(31,10);			
					//(*(gamma.end() - 1))->SetGammaMultiplicity(data);

				}
			
			
			if (buff[7] == 'G' && (buff[5] == 'S' || buff[5] == '2') && buff[6] == ' ') //electron conversion
			{
				data = buff.substr(9, 3);
				if (data == "KC=" || "NC=")
				{
					pugi::xml_node nodeTransition = nodeLevelD.last_child();
					
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
						
						string typ1, typ2;
						typ1 = type[0];
						typ2 = type[1];
						string typ = typ1 + typ2; // type without '+'
						
						nodeTransition.append_attribute(typ.c_str()) = value;
						//nodeTransition.append_attribute("shellElectronConvCoefType") = type;
						//(*(gamma.end() - 1))->SetShellElectronConvCoef(type, value);
					}

				}
			}			
		}
	}
	inDecay.close();
//	cout << "Second end of " + (string) ensDecayFilename + " file." << endl;
	outSummary << "normalizedTotalBetaIntensity: " << normalizedTotalBetaIntensity << endl;
	
	string atomicNameMother;
	if(whichBeta == -1){
		nodeNuclideM.append_attribute("AtomicNumber") = atomicNumber - 1;
		nodeNuclideM.append_attribute("AtomicMass") = atomicMass;
		atomicNameMother = PeriodicTable::GetAtomicNameCap(atomicNumber - 1);
	}
					
	#ifdef NEUTRON
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
					nodeLevelDN.append_attribute("Energy") = energyLevel;
					nodeLevelDN.append_attribute("Spin"); // to do in the future
					nodeLevelDN.append_attribute("Parity"); // to do in the future
					nodeLevelDN.append_attribute("HalfLifeTime");
					nodeLevelDN.append_attribute("TimeUnit");			
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
					nodeLevelD.append_attribute("Energy") = neutronEnergyLevel;
					nodeLevelD.append_attribute("Spin"); // to do in the future
					nodeLevelD.append_attribute("Parity"); // to do in the future
					nodeLevelD.append_attribute("HalfLifeTime");
					nodeLevelD.append_attribute("TimeUnit");
					
					pugi::xml_node nodeTransition = nodeLevelD.append_child("Transition");
					nodeTransition.append_attribute("Type") = "N";
					nodeTransition.append_attribute("TransitionQValue") = neutronEnergy;
					nodeTransition.append_attribute("Intensity") = neutronIntensity;
					
					pugi::xml_node nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy") = energyLevel;
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumberN;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMassN;
					
					nodeTransition = nodeLevelM.append_child("Transition");
					nodeTransition.append_attribute("Type") = "B-";
					nodeTransition.append_attribute("TransitionQValue") = maxBetaEnergy;
					nodeTransition.append_attribute("Intensity") = neutronIntensity;
					
					nodeTargetLevel = nodeTransition.append_child("TargetLevel");
					nodeTargetLevel.append_attribute("Energy") = neutronEnergyLevel;
					nodeTargetLevel.append_attribute("AtomicNumber") = atomicNumber;
					nodeTargetLevel.append_attribute("AtomicMass") = atomicMass;
				}
			}
		}
		inNeutron.close();
//		cout << "End of " + (string) ensNeutronFilename + " file - second." << endl;
		
		string atomicNameDaughterNeutron = PeriodicTable::GetAtomicNameCap(atomicNumberN);
		xmlDNFilename = num2string(atomicMassN) + atomicNameDaughterNeutron + ".xml";
		std::cout << "Saving result: " << docDaughterNeutron.save_file(xmlDNFilename.c_str()) << std::endl;
		
		outSummary << "normalizedTotalNeutronIntensity: " << normalizedTotalNeutronIntensity << endl;
	#endif
	
	string xmlDFilename = num2string(atomicMassDaughter) + atomicNameDaughter + ".xml";
	std::cout << "Saving result: " << docDaughter.save_file(xmlDFilename.c_str()) << std::endl;
	string xmlMFilename = num2string(atomicMassDaughter) + atomicNameMother + ".xml";
	std::cout << "Saving result: " << docMother.save_file(xmlMFilename.c_str()) << std::endl;
	
	string defaultXmlFilename = "Decay.xml";
	CreateDefaultDecayXML(defaultXmlFilename, xmlMFilename, xmlDFilename, atomicMassDaughter, 
	atomicNumber, atomicNumber + whichBeta, atomicNumber + whichBeta, xmlDNFilename);
	
	outSummary << endl << "beta-gamma intensity smaller than Sn: " << betaIntensityUnderSn << endl;
	outSummary << "beta-gamma intensity bigger than Sn: " << betaIntensityAboveSn << endl;
	outSummary << "summed: " << betaIntensityAboveSn + betaIntensityUnderSn << endl;
	outSummary << "ratio in percents (including delayed neutrons levels): " << betaIntensityAboveSn * (1. - g_delayedNeutronPercentage) << endl;
	outSummary.close();

return 0;
}
