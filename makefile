# makefile of Seebeck Coefficient

##################################################################################

clean:
	rm -r *.out *.o *.dat 


seebeck:ZEROIN.F90 thermcd.for	SURFPH.F90 SPLINT.FOR SPLINE.FOR GaussLeg.f90 GAULEG.FOR \
	Bi.f90 CURRENT.FOR CURRENT_EP.FOR SCURRENT.FOR FERMIR.FOR FERMIL.FOR FZ0SQ.F90 \
	FZPNSQ.F90 main.for

	pgf90 -o seebeck.out ZEROIN.F90 thermcd.for SURFPH.F90 SPLINT.FOR SPLINE.FOR \
 	GaussLeg.f90 GAULEG.FOR Bi.f90 CURRENT.FOR CURRENT_EP.FOR SCURRENT.FOR   \
 	FERMIR.FOR FERMIL.FOR FZ0SQ.F90 FZPNSQ.F90 main.for



