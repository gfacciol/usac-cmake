%
% Created by Gabriele Facciolo on 22/04/16.
%
mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/grammar.c   -outdir MEX
mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/libconfig.c -outdir MEX
mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/scanctx.c   -outdir MEX
mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/scanner.c   -outdir MEX
mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/strbuf.c    -outdir MEX

mex -c -DHAVE_XLOCALE_H -I../ext/libconfig-1.4.5/lib ../ext/libconfig-1.4.5/lib/libconfigcpp.cc -outdir MEX

mex -c -I../src/ -I../ext/libconfig-1.4.5/lib ../src/config/ConfigFileReader.cpp       -outdir MEX
mex -c -I../src/ -I../ext/libconfig-1.4.5/lib ../src/config/ConfigParamsHomog.cpp      -outdir MEX     
mex -c -I../src/ -I../ext/libconfig-1.4.5/lib ../src/config/ConfigParams.cpp           -outdir MEX     
mex -c -I../src/ -I../ext/libconfig-1.4.5/lib ../src/config/ConfigParamsFundmatrix.cpp -outdir MEX
mex -c -I../src/ ../src/utils/MathFunctions.cpp           -outdir MEX
mex -c -I../src/ ../src/utils/FundmatrixFunctions.cpp     -outdir MEX
mex -c -I../src/ ../src/utils/HomographyFunctions.cpp     -outdir MEX


mex -I../ext/libconfig-1.4.5/lib -I../src/  MEX/MEX_usac.cpp  MEX/*.o -llapack 
%delete MEX/*.o
