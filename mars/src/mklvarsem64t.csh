#! /bin/csh

setenv MKLROOT "/opt/intel/Compiler/11.1/064/mkl"

if ($?INCLUDE) then
    setenv INCLUDE "${MKLROOT}/include:$INCLUDE"
else
    setenv INCLUDE "${MKLROOT}/include"
endif

if ($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH "${MKLROOT}/lib/em64t:$LD_LIBRARY_PATH"
else
    setenv LD_LIBRARY_PATH "${MKLROOT}/lib/em64t"
endif

if ($?MANPATH) then
    setenv MANPATH "${MKLROOT}/man/en_US:$MANPATH"
else
    setenv MANPATH "${MKLROOT}/man/en_US:`manpath`"
endif

if ($?LIBRARY_PATH) then
    setenv LIBRARY_PATH "${MKLROOT}/lib/em64t:$LIBRARY_PATH"
else
    setenv LIBRARY_PATH "${MKLROOT}/lib/em64t"
endif

if ($?CPATH) then
    setenv CPATH "${MKLROOT}/include:$CPATH"
else
    setenv CPATH "${MKLROOT}/include"
endif

if ($?FPATH) then
    setenv FPATH "${MKLROOT}/include:$FPATH"
else
    setenv FPATH "${MKLROOT}/include"
endif

if ($?NLSPATH) then
    setenv NLSPATH "${MKLROOT}/lib/em64t/locale/%l_%t/%N:$NLSPATH"
else
    setenv NLSPATH "${MKLROOT}/lib/em64t/locale/%l_%t/%N"
endif
