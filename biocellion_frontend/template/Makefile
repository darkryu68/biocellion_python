include /users/baguilar/biocellion1.1/Makefile.common
include /users/baguilar/biocellion1.1/Makefile.model

INCLUDE += -I /users/baguilar/biocellion1.1/libmodel/include

CXXFLAG += -fPIC

libmodelroutine.a: model_routine_config.o model_routine_grid.o model_routine_agent.o model_routine_mech_intrct.o model_routine_output.o
	$(AR) rcs libmodelroutine.a model_routine_config.o model_routine_grid.o model_routine_agent.o model_routine_mech_intrct.o model_routine_output.o

model_routine_config.o: model_routine_config.cpp
	$(CXX) $(CXXFLAG) $(INCLUDE) -c model_routine_config.cpp -o model_routine_config.o

model_routine_grid.o: model_routine_grid.cpp
	$(CXX) $(CXXFLAG) $(INCLUDE) -c model_routine_grid.cpp -o model_routine_grid.o

model_routine_agent.o: model_routine_agent.cpp
	$(CXX) $(CXXFLAG) $(INCLUDE) -c model_routine_agent.cpp -o model_routine_agent.o

model_routine_mech_intrct.o: model_routine_mech_intrct.cpp
	$(CXX) $(CXXFLAG) $(INCLUDE) -c model_routine_mech_intrct.cpp -o model_routine_mech_intrct.o

model_routine_output.o: model_routine_output.cpp
	$(CXX) $(CXXFLAG) $(INCLUDE) -c model_routine_output.cpp -o model_routine_output.o

all:
	make

clean:
	$(RM) *.ilo
	$(RM) *.o
	$(RM) *.a

clean_all:
	$(RM) *.ilo
	$(RM) *.o
	$(RM) *.a
	
