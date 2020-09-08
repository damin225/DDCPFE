INCLUDE 'SetUpCrystalProps.f90'
!-----------------------Usage of UEXTERNALDB---------------------------
!	called once each at the beginning of the analysis, at the beginning of each increment, at the end of each increment, and at the end of the analysis (in addition, the user subroutine is also called once at the beginning of a restart analysis);

!	can be used to communicate between other software and user subroutines within Abaqus/Standard;

!	can be used to open external files needed for other user subroutines at the beginning of the analysis and to close those files at the end of the analysis;

!	can be used to calculate or read history information at the beginning of each increment. This information can be written to user-defined COMMON block variables or external files for use during the analysis by other user subroutines; and

!	can be used to write the current values of the user-calculated history information to external files.


!**********************************************************
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!**********************************************************
    use numtype
    use workdir
    use timing
    implicit none


    INTEGER :: LOP,LRESTART,KSTEP,KINC, length, lenjobname, lenoutdir
    REAL(KIND=rkind)::TIME(2),DTIME

!---------------Timing--------------
    CHARACTER(8) :: DATE
    CHARACTER(10) :: CURRENTTIME
    CHARACTER(5) :: ZONE
    INTEGER,DIMENSION(8) :: VALUES
!-----------------------------------
    logical :: deldexists,comp
    integer :: STAT

!-----------------------------------------------------------------------
    
    IF (LOP == 0 .or. LOP == 4) THEN !Do the following if it is at the beginning of a new analysis (LOP=0) or a restart analysis (LOP=4)
        
        if (debug==1) write(*,*) 'uexternaldb called with lop=', lop 
        FirstIncr=.true.

!		initialize the time counter of the dense solver and sparse solver        
        tss=0.d0
        ths=0.d0
        tsj=0.d0        
 
        call system_clock(count_rate=cr)
        call system_clock(count_max=cm)  
        call system_clock(c1) 
     
!		obtain the simulation director
        call getoutdir( outdir,lenoutdir )
        length = index(outdir,' ') - 1
        outdir = outdir(1:length)//'/'
        
        call getjobname( jobname,lenjobname )   
        call DATE_AND_TIME(DATE,CURRENTTIME,ZONE,VALUES)
        write(*,"(A30,A,A,I2,A,I2,A,I2,A,I3,A)")&
             'Start Analysis of job ', TRIM(ADJUSTL(jobname)), ': ',  VALUES(5),&
            ':', VALUES(6), ':', VALUES(7), ' on', VALUES(3), '^th'
    
        call SetUpCrystalProps( )            

    ELSE IF (LOP == 3) THEN !END OF ANALYSIS  
        
        call system_clock(c2)    
        write(*,*) 'tss= ', tss, 'ths= ', ths, 'tsj= ', tsj
        CALL DATE_AND_TIME(DATE,CURRENTTIME,ZONE,VALUES)
        write(*,"(A30,A,A,I2,A,I2,A,I2,A,I3,A)")&
             'End Analysis of job ', TRIM(ADJUSTL(jobname)), ': ',  VALUES(5),&
            ':', VALUES(6), ':', VALUES(7), ' on', VALUES(3), '^th'  
            Write(*,'(A30, I6, A, I6, A, f12.5, A)') 'Total walltime: ', floor((c2-c1)/REAL(cr)/3600),  &
                        ' hours', floor((c2-c1)/REAL(cr)/60-floor((c2-c1)/REAL(cr)/3600)*60), ' minutes',                       &
                    ((c2-c1)/REAL(cr)- floor((c2-c1)/REAL(cr)/60)*60), ' seconds.'                 
    ELSE                    
        FirstIncr=.false.               
    ENDIF
    
    RETURN

END SUBROUTINE UEXTERNALDB
