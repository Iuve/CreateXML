#DESTDIR=$(MTAS_BIN)
DESTDIR=.
CINCLUDEDIRS = -Iinclude -std=c++11
CPPFLAGS += $(CINCLUDEDIRS)
SRC = src/
CreateXML: CreateXML.o PeriodicTable.o pugixml.o FermiDistribution.o
	g++ CreateXML.o PeriodicTable.o pugixml.o FermiDistribution.o -o $(DESTDIR)/CreateXML
	
PeriodicTable.o:
	g++ $(CPPFLAGS) -c $(SRC)PeriodicTable.cc

pugixml.o:
	g++ $(CPPFLAGS) -c $(SRC)pugixml.cpp

FermiDistribution.o:
	g++ $(CPPFLAGS) -c $(SRC)FermiDistribution.cc
					
CreateXML.o:
	g++ $(CPPFLAGS) -c CreateXML.cc
		
clean:
	rm -f *.o *~	


