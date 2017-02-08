ll:
	( cd SRC          ; make )
	( cd TESTS_Gen    ; mkdir -p OUT  ; make all)
	( cd TESTS_Lap    ; mkdir -p OUT  ; make all)
	( cd TESTS_GenEv  ; mkdir -p OUT  ; make all)
clean:
	( cd SRC    ; make clean)
	( cd TESTS_Gen    ; make cleanall)
	( cd TESTS_GenEv  ; make cleanall)
	( cd TESTS_Lap    ; make cleanall)

docs:
	 ( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

