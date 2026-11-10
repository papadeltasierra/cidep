ftn = fl /G2 /4I2 /4Nt /Ox /FPa /c

mode1.obj: $*.for
    $(ftn) $*.for

screen.obj:   $*.for
    $(ftn) $*.for

scrnbox.obj:  $*.for
    $(ftn) $*.for

scrnplot.obj: $*.for
    $(ftn) $*.for

scrncurs.obj: $*.for
    $(ftn) $*.for

mode1.exe:    mode1.obj screen.obj \
              scrnplot.obj scrnbox.obj scrncurs.obj \
              c:\gph\asm\graph.lib
    link @mode1.lnk

