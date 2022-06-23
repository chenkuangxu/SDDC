FCOMP    = mpif90
OPTS     = -c -O3 
#OPTS     = -c -O0 -g -traceback -fpe:0 -check all -fpstkchk

SWP	= swplist
LINKOPTS = -O3 -o
OBJS = setup2d.o main.o readinput.o readGRIDinfo.o connectivity.o \
	connect_bdry.o map_proc_match.o map_cyc_match_rem.o getIVCELL_proc.o \
	getIVCELL_cyc.o map_cyc_loc.o initsetup.o mapder.o calcjacob.o \
	setinitialcond.o tecplotter2dsetup.o tecplotter.o \
	compflux.o getrusanov.o interfaceflux_all.o procintflux_all.o \
	cycremflux_all.o cyclocflux_all.o compresid.o \
	iterations.o write_data.o read_data_restart.o bcflux.o AV.o
	

sd3dhexa:$(OBJS)
	$(FCOMP) $(LINKOPTS) ./rundir/sd2d_MHD $(OBJS)

clean:
	 rm  *.o *.mod
#	$(RM) $(EXEC) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .F .c .f .swp .f90

.F.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<

.f.o:
	$(FCOMP)  $(OPTS) $(DEFINES) $<

.f90.o:
	$(FCOMP)  $(OPTS) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<
