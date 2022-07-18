program compONTvBS
        !Programa para comparar marcas de metilaciÃƒÂ³n de ONT con BS.
        implicit none
        integer :: io,i,j,nlinesONT, nlinesBS
        integer, allocatable :: startONT(:),startBS(:),endBS(:),endONTcorrected(:),endONT(:), &
                                startONTcorrected(:), col1(:),col2(:),col3(:),covBSunmeth(:),covBSmeth(:),covONT(:),covBS(:), &
                                chrONT(:), chrBS(:)
        real, allocatable :: freqBScorrected(:),freqONT(:), freqBS(:)


write (*,*) "*********************************************************************************"

write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) "                                       /;    ;\"
write (*,*) "                                   __  \\____//"
write (*,*) "                                  /{_\_/   `'\____"
write (*,*) "                                  \___   (o)  (o  }"
write (*,*) "       _____________________________/          :--'  "
write (*,*) "   ,-,'`@@@@@@@@       @@@@@@         \_    `__\"
write (*,*) "  ;:(  @@@@@@@@@        @@@             \___(o'o)"
write (*,*) "  :: )  @@@@          @@@@@@        ,'@@(  `===='     "
write (*,*) "  :: : @@@@@:          @@@@         `@@@:"
write (*,*) "  :: \  @@@@@:       @@@@@@@)    (  '@@@'"
write (*,*) "  ;; /\      /`,    @@@@@@@@@\   :@@@@@)"
write (*,*) "  ::/  )    {_----------------:  :~`,~~;"
write (*,*) " ;;'`; :   )                  :  / `; ;"
write (*,*) ";;;; : :   ;                  :  ;  ; : "             
write (*,*) "`'`' / :  :                   :  :  : :"
write (*,*) "    )_ \__;      ;          :_ ;  \_\       `,',' "
write (*,*) "    :__\  \    * `,'*         \  \  :  \   *  8`;'*  *"
write (*,*) "        `^'     \ :/           `^'  `-^-'   \v/ :  \/ "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " "
write (*,*) " ______________________________________________________________ "
write (*,*) "|CUIDAO! Este programa asume archivos de entrada sin cabecera. |"
write (*,*) "|El cromosoma X ha sido renombrado a cromosoma 30.             |"
write (*,*) "|El cromosoma mitocondrial ha sido renombrado a cromosoma 31.  |"
write (*,*) "|______________________________________________________________|"
write (*,*) "********************************************************************************"

!Read files

    open (10,file="methylation_renamed.tsv",form='formatted',  status='old')
    open (20,file="Bisulfite.txt", form='formatted',status='old')
    open (30,file="ComparedSites.tsv",form="formatted",status="new")

!Count lines
    nlinesONT=0
        do
		    read(10,*,iostat=io)
		    if (io.ne.0) exit
		    nlinesONT=nlinesONT+1
	    end do
	rewind(10)

    nlinesBS=0
        do
		    read(20,*,iostat=io)
		    if (io.ne.0) exit
		    nlinesBS=nlinesBS+1
	    end do
	rewind(20)

!Allocate variables
    allocate(chrONT(nlinesONT),chrBS(nlinesBS),startONT(nlinesONT),startBS(nlinesBS),startONTcorrected(nlinesONT), &
             endBS(nlinesBS),endONT(nlinesONT),endONTcorrected(nlinesONT),col1(nlinesONT),col2(nlinesONT),col3(nlinesONT),&
             freqBS(nlinesBS),freqONT(nlinesONT),freqBScorrected(nlinesBS),covBS(nlinesBS),&
             covBSunmeth(nlinesBS),covBSmeth(nlinesBS),covONT(nlinesONT))

!Assign cols

    do i=1,nlinesONT
	    read (10,*) chrONT(i),startONT(i),endONT(i),col1(i),covONT(i),col3(i),freqONT(i)
    end do

    do i=1,nlinesBS
	    read (20,*) chrBS(i),startBS(i),covBSunmeth(i),covBSmeth(i),freqBS(i)
    end do

    write (*,*) "ONT file has ",nlinesONT," methylated sites."
    write (*,*) "BS file has ",nlinesBS," methylated sites."



    do i=1, nlinesONT
        startONTcorrected(i)=startONT(i)+1
        endONTcorrected(i)=endONT(i)+1
    enddo

    do i=1,nlinesBS
        endBS(i)=startBS(i)
        freqBScorrected(i)=freqBS(i)/100
        covBS(i)=covBSmeth(i)+covBSunmeth(i)
    enddo



write (*,*) " __________________________________________________________________________"   
write (*,*) "|Files loaded. Frequencies and positions transformed. Starting comparison. |"
write (*,*) "|__________________________________________________________________________|"
     i=1 !nlinesONT
     j=1 !nlinesBS
     12 continue
     if (i.gt.nlinesONT) then
        write(*,*) "Finalizado archivo de ONT en la iteración= ", i
        goto 15
     end if  
     if	(j.gt.nlinesBS) then
        write(*,*) "Fin del archivo de BS en la iteración= ", j
        goto 15
     end if

if (chrBS(j)==chrONT(i)) then
    if (startONTcorrected(i)<=startBS(j) .and. endONTcorrected(i)>=startBS(j)) then
        !write(30,'(i4,i15,2f12.6,2i6)') chrBS(j),startBS(j), freqBScorrected(j),freqONT(i),covBS(j),covONT(i)
        write(30,*) chrBS(j),startBS(j), freqBScorrected(j),freqONT(i),covBS(j),covONT(i)
        j=j+1
    else if (startONTcorrected(i)>startBS(j)) then
        j=j+1
    else if (endONTcorrected(i)<startBS(j)) then
        i=i+1
    end if    
else if (chrBS(j)<chrONT(i)) then
    j=j+1
else if (chrBS(j)>chrONT(i)) then
    i=i+1
end if 

goto 12

15 continue

write (*,*) 'All done! Go check your matches!'

end program compONTvBS

