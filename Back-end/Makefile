# Makefile
all: final.exe

# Run Python script to generate CSV before compiling C++ code
maturity_days_and_rates.csv: cubic.py
	python cubic.py

main.o: main.cpp TimeContingentCashFlows.h TermStructureHoLee.h TermStructure.h date.h maturity_days_and_rates.csv
	g++ -std=c++98 -g -Wall -c main.cpp -o main.o

date.o: date.cpp date.h
	g++ -std=c++98 -g -Wall -c date.cpp -o date.o

TermStructure.o: TermStructure.cpp TermStructure.h
	g++ -std=c++98 -g -Wall -c TermStructure.cpp -o TermStructure.o

TermStructureHoLee.o: TermStructureHoLee.cpp TermStructure.h TermStructureHoLee.h
	g++ -std=c++98 -g -Wall -c TermStructureHoLee.cpp -o TermStructureHoLee.o

TimeContingentCashFlows.o: TimeContingentCashFlows.cpp TermStructure.h TermStructureHoLee.h TimeContingentCashFlows.h
	g++ -std=c++98 -g -Wall -c TimeContingentCashFlows.cpp -o TimeContingentCashFlows.o

final.exe: main.o TermStructure.o TermStructureHoLee.o TimeContingentCashFlows.o
	g++ -std=c++98 main.o TermStructure.o TermStructureHoLee.o TimeContingentCashFlows.o -o final.exe

clean:
	rm -f *.o *.exe



# all: final.exe

# main.o: main.cpp TimeContingentCashFlows.h TermStructureHoLee.h TermStructure.h 
# 	g++ -std=c++98 -g -Wall -c main.cpp -o main.o

# TermStructure.o: TermStructure.cpp TermStructure.h 
# 	g++ -std=c++98 -g -Wall -c TermStructure.cpp -o TermStructure.o

# TermStructureHoLee.o: TermStructureHoLee.cpp TermStructure.h TermStructureHoLee.h
# 	g++ -std=c++98 -g -Wall -c TermStructureHoLee.cpp -o TermStructureHoLee.o

# TimeContingentCashFlows.o: TimeContingentCashFlows.cpp TermStructure.h TermStructureHoLee.h TimeContingentCashFlows.h
# 	g++ -std=c++98 -g -Wall -c TimeContingentCashFlows.cpp -o TimeContingentCashFlows.o

# final.exe: main.o TermStructure.o TermStructureHoLee.o TimeContingentCashFlows.o
# 	g++ -std=c++98 main.o TermStructure.o TermStructureHoLee.o TimeContingentCashFlows.o -o final.exe

# clean:
# 	rm -f *.o *.exe
