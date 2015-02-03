all:fair

hello_mpi:
	rm -f -r *.o
	mpicc -c hello_mpi.c
	mpicc hello_mpi.o -o hello_mpi
	rm -f -r *.o

send_recv:
	rm -f -r *.o
	mpicc -std=gnu99 -c send_recv.c
	mpicc send_recv.o -o send_recv
	rm -f -r *.o

comm2: comm2.c
	rm -f -r *.o
	rm -f coom2
	mpicc -std=gnu99 -c comm2.c
	mpicc comm2.o -o comm2
	rm -f -r *.o

clean:
	rm -f -r *.o
	rm -f comm2

pi: pi.c
	rm -f -r *.o	
	rm -f pi
	mpicc -std=gnu99 -c pi.c
	mpicc pi.o -o pi
	rm -f -r *.o
	
benmark: benmark.c
	rm -f -r *.o	
	rm -f benmark
	mpicc -std=gnu99 -c benmark.c
	mpicc benmark.o -o benmark
	rm -f -r *.o	
	
fair: fair.c
	rm -f -r *.o	
	rm -f fair
	mpicc -std=gnu99 -c fair.c
	mpicc fair.o -o fair
	rm -f -r *.o	
