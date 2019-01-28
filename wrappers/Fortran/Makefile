default: mui_3df.cpp unit_test.f90
	@echo "Generating F-wrapper object file..."
	mpic++ -std=c++11 -c mui_3df.cpp -o mui_3df.o -g
	@echo "Compiling unit test code..."
	mpifort  -std=f2003 mui.f90 unit_test.f90 mui_3df.o  -o unit_test.x -g -lstdc++ -lmpi_cxx
	@echo "Launching unit test code..."
	#mpirun  -tag-output -np 1 ./unit_test.x mpi://domain1/interface :  -np 1 ./unit_test.x mpi://domain2/interface
	mpirun  -np 1 ./unit_test.x mpi://domain1/interface :  -np 1 ./unit_test.x mpi://domain2/interface
	@echo "Usage: "
	@echo "Step 1) Compile with a Fortran compiler & link"
clean:
	rm *.o *.mod unit_test.x
