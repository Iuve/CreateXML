#DESTDIR=$(MTAS_BIN)
DESTDIR=.
CINCLUDEDIRS = -Iinclude
CPPFLAGS += $(CINCLUDEDIRS)
SRC = src/
CreateXML: CreateXML.o PeriodicTable.o pugixml.o
	g++ CreateXML.o PeriodicTable.o pugixml.o -o $(DESTDIR)/CreateXML
	
PeriodicTable.o:
	g++ $(CPPFLAGS) -c $(SRC)PeriodicTable.cc

pugixml.o:
	g++ $(CPPFLAGS) -c $(SRC)pugixml.cpp
					
CreateXML.o:
	g++ $(CPPFLAGS) -c CreateXML.cc
		
clean:
	rm -f *.o *~	


