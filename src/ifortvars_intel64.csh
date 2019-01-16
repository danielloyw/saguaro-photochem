#
# Copyright  (C) 1985-2009 Intel Corporation. All rights reserved.
#
# The information and source code contained herein is the exclusive property
# of Intel Corporation and may not be disclosed, examined, or reproduced in
# whole or in part without explicit written authorization from the Company.
#

#! /bin/tcsh

if !($?PATH) then
    setenv PATH /opt/intel/Compiler/11.1/064/bin/intel64
else
    setenv PATH /opt/intel/Compiler/11.1/064/bin/intel64:${PATH}
endif
if ( "`uname`" == "Darwin" ) then
   # DYLD_LIBRARY_PATH is used on MAC OS*
   if !($?DYLD_LIBRARY_PATH) then
       setenv DYLD_LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib
   else
       setenv DYLD_LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib:${DYLD_LIBRARY_PATH}
   endif

   if !($?LIBRARY_PATH) then
       setenv LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib
   else
       setenv LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib:${LIBRARY_PATH}
   endif

   if !($?NLSPATH) then
       setenv NLSPATH /opt/intel/Compiler/11.1/064/lib/locale/en_US/%N
   else
       setenv NLSPATH /opt/intel/Compiler/11.1/064/lib/locale/en_US/%N:${NLSPATH}
   endif
   if !($?INTEL_LICENSE_FILE) then
       setenv INTEL_LICENSE_FILE "/opt/intel/Compiler/11.1/064/licenses:/opt/intel/licenses:${HOME}/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses"
   else
       setenv INTEL_LICENSE_FILE "/opt/intel/Compiler/11.1/064/licenses:/opt/intel/licenses:${HOME}/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses:${INTEL_LICENSE_FILE}"
   endif
   if !($?MANPATH) then
       setenv MANPATH /opt/intel/Compiler/11.1/064/man/en_US:`manpath`
   else
       setenv MANPATH /opt/intel/Compiler/11.1/064/man/en_US:${MANPATH}
   endif
else
   if !($?LD_LIBRARY_PATH) then
       setenv LD_LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib/intel64
   else
       setenv LD_LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib/intel64:${LD_LIBRARY_PATH}
   endif
   if !($?LIBRARY_PATH) then
       setenv LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib/intel64
   else
       setenv LIBRARY_PATH /opt/intel/Compiler/11.1/064/lib/intel64:${LIBRARY_PATH}
   endif
   if !($?NLSPATH) then
       setenv NLSPATH /opt/intel/Compiler/11.1/064/lib/intel64/locale/%l_%t/%N
   else
       setenv NLSPATH /opt/intel/Compiler/11.1/064/lib/intel64/locale/%l_%t/%N:${NLSPATH}
   endif
   if !($?INTEL_LICENSE_FILE) then
       setenv INTEL_LICENSE_FILE "/opt/intel/Compiler/11.1/064/licenses:/opt/intel/licenses:${HOME}/intel/licenses"
   else
       setenv INTEL_LICENSE_FILE "/opt/intel/Compiler/11.1/064/licenses:/opt/intel/licenses:${HOME}/intel/licenses:${INTEL_LICENSE_FILE}"
   endif
   if ($?LC_CTYPE) then
     if ("$LC_CTYPE" != "") then
       set LANGUAGE_TERRITORY=`echo $LC_CTYPE | sed s/\\..\*//`
     endif
   else if ($?LANG) then
     if ("$LANG" != "") then
       set LANGUAGE_TERRITORY=`echo $LANG | sed s/\\..\*//`
     endif
   endif
   endif
   if ("$LANGUAGE_TERRITORY" != "en_US") then
      if ("$LANGUAGE_TERRITORY" != "ja_JP") then
         set LANGUAGE_TERRITORY="en_US"
      endif
   endif
   if !($?MANPATH) then
        setenv MANPATH /opt/intel/Compiler/11.1/064/man/${LANGUAGE_TERRITORY}:`manpath`
   else
        setenv MANPATH /opt/intel/Compiler/11.1/064/man/${LANGUAGE_TERRITORY}:${MANPATH}
   endif
endif
