CC = gcc
CFLAGS = -Wall -O3 -g

exe: main.o GlobalVariables.o StructsRelated.o TargetsComputation.o Feasibility.o Outputs.o GeneralUsageFunctions.o Aging.o EA.o MutationsCrossover.o
	$(CC) $(CFLAGS) -o exe main.o GlobalVariables.o StructsRelated.o TargetsComputation.o Feasibility.o Outputs.o GeneralUsageFunctions.o Aging.o EA.o MutationsCrossover.o -lm

main.o: main.c GlobalVariables.h StructsRelated.h TargetsComputation.h Feasibility.h Outputs.h GeneralUsageFunctions.h Aging.h EA.h Constants.h Structs.h
	$(CC) $(CFLAGS) -c main.c

EA.o: EA.c EA.h Feasibility.h MutationsCrossover.h StructsRelated.h Constants.h Structs.h
	$(CC) $(CFLAGS) -c EA.c

Outputs.o: Outputs.c Outputs.h Constants.h GlobalVariables.h Structs.h
	$(CC) $(CFLAGS) -c Outputs.c

TargetsComputation.o: TargetsComputation.c TargetsComputation.h
	$(CC) $(CFLAGS) -c TargetsComputation.c

MutationsCrossover.o: MutationsCrossover.c MutationsCrossover.h GlobalVariables.h Constants.h Structs.h
	$(CC) $(CFLAGS) -c MutationsCrossover.c

GeneralUsageFunctions.o: GeneralUsageFunctions.c GeneralUsageFunctions.h Constants.h
	$(CC) $(CFLAGS) -c GeneralUsageFunctions.c

Feasibility.o: Feasibility.c Feasibility.h Constants.h Structs.h
	$(CC) $(CFLAGS) -c Feasibility.c

Aging.o: Aging.c Aging.h GlobalVariables.h
	$(CC) $(CFLAGS) -c Aging.c

GlobalVariables.o: GlobalVariables.c GlobalVariables.h
	$(CC) $(CFLAGS) -c GlobalVariables.c

StructsRelated.o: StructsRelated.c StructsRelated.h Structs.h GlobalVariables.h Constants.h
	$(CC) $(CFLAGS) -c StructsRelated.c

clean: 
	rm -f exe *.o Makefile~
