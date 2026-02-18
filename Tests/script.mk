#
# Single-test makefile template for script usage
SHELL=/bin/bash

MPI=mpirun -np 2
SIESTA=../../../siesta

# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(SIESTA))),"mpirun")
MPI=
endif
ifeq ($(strip $(firstword $(SIESTA))),"mpiexec")
MPI=
endif

label=work

.PHONY: completed
completed: completed_$(label)

completed_$(label):
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -f script.sh ] ; then cp -f script.sh $(label) ; fi
	@echo "    ==> Running script with SIESTA as ${SIESTA}"
	@(cd $(label) ; sh script.sh "$(MPI) $(SIESTA)")  && touch completed
	@if [ -f completed ] ; then \
           echo "    ===> Script finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi

xmlcheck: completed
	@echo "---- xmllint check $(name).xml ..."
	xmllint $(name).xml > /dev/null

clean:
	@echo ">>>> Cleaning $(name) [label=$(label)] test..."
	rm -rf $(label) completed_$(label) $(name).out $(name).xml
