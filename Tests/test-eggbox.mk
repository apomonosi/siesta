#
# Single-test makefile template
#

$(label)-completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for ps in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file $$ps ..." ;\
          ln ../Pseudos/$$ps $(label)/$$ps ;\
	done 

	@echo "    ==> Running SIESTA as $(SIESTA)"
	@(cd $(label) ; ../../Scripts/eggbox_size.sh $(SIESTA) $(name) ) \
		&& touch $(label)-completed
	@if [ -f $(label)-completed ]; then \
	   echo "    ===> SIESTA finished";\
	   else \
	   echo " **** Test $(name) did not complete successfully";\
	fi

#${SIESTA} 2>&1 > $(name).out < ../$(name).fdf
