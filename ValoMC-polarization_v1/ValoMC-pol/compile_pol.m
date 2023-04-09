%% For Mie Theory mex file
mex cpp/mie/MTmex.cpp COMPFLAGS='$COMPFLAGS -O3'
% There may be some warnings while compiling this just ignore them.
%%For Monte Carlo
% For MC without paralelization
mex cpp/3d/MC3Dmex.cpp COMPFLAGS='$COMPFLAGS -O3'
% for MC with paralelization
%mex -DUSE_OMP cpp/3d/MC3Dmex.cpp CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"