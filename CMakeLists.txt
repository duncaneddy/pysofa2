cmake_minimum_required(VERSION 3.12)
project(_sofa_c)


# Add any additional include directories here
include_directories(src)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

if (WIN32)
  if (BUILD_SHARED_LIBS)
	  if (MINGW OR MSYS OR CYGWIN)
	  # only for the microsoft compiler
	  elseif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
		MESSAGE(STATUS "Exporting all symbols on windows")
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
	  endif()
  endif()
endif()

if(WIN32)
  set(LIBM_LIBRARIES "")
else(WIN32)
  if(MINGW)
    set(LIBM_LIBRARIES "")
  else(MINGW)
    find_library(LIBM_LIBRARIES m)
  endif(MINGW)
endif(WIN32)


add_library(${PROJECT_NAME}	SHARED
	src/a2af.c
	src/a2tf.c
	src/ab.c
	src/ae2hd.c
	src/af2a.c
	src/anp.c
	src/anpm.c
	src/apcg.c
	src/apcg13.c
	src/apci.c
	src/apci13.c
	src/apco.c
	src/apco13.c
	src/apcs.c
	src/apcs13.c
	src/aper.c
	src/aper13.c
	src/apio.c
	src/apio13.c
	src/atcc13.c
	src/atccq.c
	src/atci13.c
	src/atciq.c
	src/atciqn.c
	src/atciqz.c
	src/atco13.c
	src/atic13.c
	src/aticq.c
	src/aticqn.c
	src/atio13.c
	src/atioq.c
	src/atoc13.c
	src/atoi13.c
	src/atoiq.c
	src/bi00.c
	src/bp00.c
	src/bp06.c
	src/bpn2xy.c
	src/c2i00a.c
	src/c2i00b.c
	src/c2i06a.c
	src/c2ibpn.c
	src/c2ixy.c
	src/c2ixys.c
	src/c2s.c
	src/c2t00a.c
	src/c2t00b.c
	src/c2t06a.c
	src/c2tcio.c
	src/c2teqx.c
	src/c2tpe.c
	src/c2txy.c
	src/cal2jd.c
	src/cp.c
	src/cpv.c
	src/cr.c
	src/d2dtf.c
	src/d2tf.c
	src/dat.c
	src/dtdb.c
	src/dtf2d.c
	src/eceq06.c
	src/ecm06.c
	src/ee00.c
	src/ee00a.c
	src/ee00b.c
	src/ee06a.c
	src/eect00.c
	src/eform.c
	src/eo06a.c
	src/eors.c
	src/epb.c
	src/epb2jd.c
	src/epj.c
	src/epj2jd.c
	src/epv00.c
	src/eqec06.c
	src/eqeq94.c
	src/era00.c
	src/fad03.c
	src/fae03.c
	src/faf03.c
	src/faju03.c
	src/fal03.c
	src/falp03.c
	src/fama03.c
	src/fame03.c
	src/fane03.c
	src/faom03.c
	src/fapa03.c
	src/fasa03.c
	src/faur03.c
	src/fave03.c
	src/fk5hip.c
	src/fk5hz.c
	src/fk45z.c
	src/fk52h.c
	src/fk54z.c
	src/fk425.c
	src/fk524.c
	src/fw2m.c
	src/fw2xy.c
	src/g2icrs.c
	src/gc2gd.c
	src/gc2gde.c
	src/gd2gc.c
	src/gd2gce.c
	src/gmst00.c
	src/gmst06.c
	src/gmst82.c
	src/gst00a.c
	src/gst00b.c
	src/gst06.c
	src/gst06a.c
	src/gst94.c
	src/h2fk5.c
	src/hd2ae.c
	src/hd2pa.c
	src/hfk5z.c
	src/icrs2g.c
	src/ir.c
	src/jd2cal.c
	src/jdcalf.c
	src/ld.c
	src/ldn.c
	src/ldsun.c
	src/lteceq.c
	src/ltecm.c
	src/lteqec.c
	src/ltp.c
	src/ltpb.c
	src/ltpecl.c
	src/ltpequ.c
	src/makefile
	src/moon98.c
	src/num00a.c
	src/num00b.c
	src/num06a.c
	src/numat.c
	src/nut00a.c
	src/nut00b.c
	src/nut06a.c
	src/nut80.c
	src/nutm80.c
	src/obl06.c
	src/obl80.c
	src/p06e.c
	src/p2pv.c
	src/p2s.c
	src/pap.c
	src/pas.c
	src/pb06.c
	src/pdp.c
	src/pfw06.c
	src/plan94.c
	src/pm.c
	src/pmat00.c
	src/pmat06.c
	src/pmat76.c
	src/pmp.c
	src/pmpx.c
	src/pmsafe.c
	src/pn.c
	src/pn00.c
	src/pn00a.c
	src/pn00b.c
	src/pn06.c
	src/pn06a.c
	src/pnm00a.c
	src/pnm00b.c
	src/pnm06a.c
	src/pnm80.c
	src/pom00.c
	src/ppp.c
	src/ppsp.c
	src/pr00.c
	src/prec76.c
	src/pv2p.c
	src/pv2s.c
	src/pvdpv.c
	src/pvm.c
	src/pvmpv.c
	src/pvppv.c
	src/pvstar.c
	src/pvtob.c
	src/pvu.c
	src/pvup.c
	src/pvxpv.c
	src/pxp.c
	src/refco.c
	src/rm2v.c
	src/rv2m.c
	src/rx.c
	src/rxp.c
	src/rxpv.c
	src/rxr.c
	src/ry.c
	src/rz.c
	src/s00.c
	src/s00a.c
	src/s00b.c
	src/s06.c
	src/s06a.c
	src/s2c.c
	src/s2p.c
	src/s2pv.c
	src/s2xpv.c
	src/sepp.c
	src/seps.c
	src/sp00.c
	src/starpm.c
	src/starpv.c
	src/sxp.c
	src/sxpv.c
	src/t_sofa_c.c
	src/taitt.c
	src/taiut1.c
	src/taiutc.c
	src/tcbtdb.c
	src/tcgtt.c
	src/tdbtcb.c
	src/tdbtt.c
	src/tf2a.c
	src/tf2d.c
	src/tpors.c
	src/tporv.c
	src/tpsts.c
	src/tpstv.c
	src/tpxes.c
	src/tpxev.c
	src/tr.c
	src/trxp.c
	src/trxpv.c
	src/tttai.c
	src/tttcg.c
	src/tttdb.c
	src/ttut1.c
	src/ut1tai.c
	src/ut1tt.c
	src/ut1utc.c
	src/utctai.c
	src/utcut1.c
	src/xy06.c
	src/xys00a.c
	src/xys00b.c
	src/xys06a.c
	src/zp.c
	src/zpv.c
	src/zr.c
)

# Link the libraries
target_link_libraries(${PROJECT_NAME} ${LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "src/sofa.h;src/sofam.h")
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN")
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

install(TARGETS ${PROJECT_NAME}
        EXCLUDE_FROM_ALL
		DESTINATION ${PY_BUILD_CMAKE_MODULE_NAME}
		COMPONENT python_module)