module nml_mod

contains

  !============================================================
  subroutine make_namelist(filename, name, tmpfile)
    !============================================================

    implicit none
    
    integer, parameter :: NCH = 200, NML_LEN=100000

    character(len=*), intent(in)         :: name
    character(len=*), intent(in)         :: filename,tmpfile
    character(len=NCH)                   :: read_str
    character(len=NML_LEN)               :: nml_str
    character                            :: lead
    !============================================================

    OPEN(1, file=TRIM(filename))
    OPEN(2,file=TRIM(tmpfile))
    WRITE(2,*) '&'//TRIM(ADJUSTL(name))
    DO 
       READ(1,'(a)', end=100) read_str
       READ(read_str,'(a1)') lead
       IF(lead /= '#'.AND.lead /= ' ') THEN
          WRITE(2,*) ADJUSTL(read_str)
       ENDIF
    ENDDO
100 continue
    CLOSE(1)
    WRITE(2,*) '/'
    CLOSE(2)

    return
  end subroutine make_namelist

end module nml_mod
