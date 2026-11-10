ftn = fl /G2 /4I2 /4Nt /Ox /FPa /c

..\mode1\scrnbox.obj:  ..\mode1\$*.for
    $(ftn) /Fo..\mode1\$*.obj ..\mode1\$*.for

..\mode1\scrnplot.obj: ..\mode1\$*.for
    $(ftn) /Fo..\mode1\$*.obj ..\mode1\$*.for

..\mode1\scrnstck.obj: ..\mode1\$*.for
    $(ftn) /Fo..\mode1\$*.obj ..\mode1\$*.for

..\mode1\scrncurs.obj: ..\mode1\$*.for
    $(ftn) /Fo..\mode1\$*.obj ..\mode1\$*.for

..\mode1\screen.obj:   ..\mode1\$*.for
    $(ftn) /Fo..\mode1\$*.obj ..\mode1\$*.for

cidep.obj:    $*.for
    $(ftn) $*.for

stm.obj:      $*.for
    $(ftn) $*.for

clr_rpm.obj:  $*.for
    $(ftn) $*.for

stick.obj:    $*.for
    $(ftn) $*.for

pol_fld.obj:  $*.for
    $(ftn) $*.for

tm_pol.obj:   $*.for
    $(ftn) $*.for

stm_pol.obj:  $*.for
    $(ftn) $*.for

lineshpe.obj: $*.for
    $(ftn) $*.for

plots.obj:    $*.for
    $(ftn) $*.for

exchline.obj: $*.for
    $(ftn) $*.for

files.obj:    $*.for
    $(ftn) $*.for

msg_inpt.obj: $*.for
    $(ftn) $*.for

calc_etc.obj: $*.for
    $(ftn) $*.for

cidep.exe:    cidep.obj    stm.obj      clr_rpm.obj  stick.obj     \
              pol_fld.obj  tm_pol.obj   stm_pol.obj  lineshpe.obj  \
              plots.obj    exchline.obj files.obj    msg_inpt.obj  \
              calc_etc.obj
    link @cidep.lnk

