!test program to compute components of the bispectrum covariance from lensing 

program bisvar 
  Use SepBispectrum
  implicit none
  !integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359

  character(80) :: Folder1, Folder2, Folder3
  character(120) :: Clfile, Cllfile

  real(dl), pointer :: Cl(:,:), Cll(:,:), invCll(:,:,:),Cllm(:,:,:)
  real(dl), pointer :: pClpp(:,:)
  integer :: lmax, lmin, l1, l2, l3, l1b, l2b, l3b
  integer :: min_l, max_l, min_lb, max_lb, l_min, l_max
  integer :: i,j,m1,p1,q1,m2,p2,q2
  !wigner 3j
  real(dl)  :: atj(0:20000),atj2(0:20000)
  real(dl), pointer :: a3j(:,:), a3joC(:,:)
  real(dl), pointer :: bispectrum(:,:,:,:,:),bispectrum_ISWlens(:,:,:,:,:), bis(:,:),bis_ISWlens(:,:)

  real(dl) :: CMB2COBEnorm = 7428350250000.d0
  real(dl) :: DB,detCovCV,tmpPrefac

  type(bfs) :: P
  type(bfs) :: P_ISWlens
  real(dl) :: DSNGauss, DSNonGauss, TotSumGauss, TotSumNGauss,TotSumNGaussBD4,DSNonGaussP
  real(dl) :: SumNGauss, SumGauss, TotNoise
  real(dl) :: sigsq, sigsqX, fnl,fnlb, fnlISW, fnlISWb
  real(dl) :: TotSumGauss_outer, TotSumNGauss_outer
  real(dl) :: Det,TempCovCV,TotSumCV,DetISWLens,TotSumCVISWLens,DetLensCross,TotSumCVLensCross
  real(dl) :: DetFishCV, DefnlMarCV
  real(dl) :: alpha, beta, tempfac,tempfacFcM(2,2)

  logical :: AWigner

  integer :: testl = 70
  integer :: ellar(512)
  real(dl):: dellar(512)  !multiples of 32
  integer :: intmax
  logical :: want_ISW_correction = .false.

  real(dl) ::stepsum, DBx
  real(dl), allocatable ::  tempB3(:,:,:,:)

  integer ::  shape,nfields,minfields
  character(120) :: alphabetafile, alphabetaPolfile,tensDir

  shape = 1
  minfields = 1
  nfields = 2

  !you can see the effect of removing ISW-lensing helps in reducing the extra covariance 
  want_ISW_correction = .true. 

  !various binning schemes
!!$  do i  = 1, 256
!!$     if (i .le. 19) then
!!$        ellar(i) = i+1
!!$        dellar(i) = 1.d0
!!$     elseif (i .le. 75 .and. i .ge. 20) then 
!!$        ellar(i) = ellar(i-1) + 4
!!$        dellar(i) = 4.d0
!!$     elseif (i .le. 257 .and. i .ge. 76) then
!!$        ellar(i)  = ellar(i-1) + 20
!!$        dellar(i) = 20.d0
!!$     endif
!!$     !write(*,*) 'ell:', ellar(i)
!!$  enddo
  do i  = 1, 256
     if (i .le. 47) then
        ellar(i) = i+1
        dellar(i) = 1.d0
     elseif (i .le. 85 .and. i .ge. 48) then 
        ellar(i) = ellar(i-1) + 4
        dellar(i) = 4.d0
     elseif (i .le. 110 .and. i .ge. 86) then
        ellar(i)  = ellar(i-1) + 12
        dellar(i) = 12.d0
     elseif (i .le. 175 .and. i .ge. 111) then
        ellar(i)  = ellar(i-1) + 24
        dellar(i) = 24.d0
     else
        ellar(i)  = ellar(i-1) + 50
        dellar(i) = 50.d0
     endif
     !write(*,*) 'ell:', ellar(i)
  enddo
  do i  = 1, 512
     if (i .le. 399) then
        ellar(i) = i + 1
        dellar(i) = 1.d0
     else
        ellar(i)  = ellar(i-1) + 2
        dellar(i) = 2.d0
     endif
     !write(*,*) 'ell:', ellar(i)
  enddo
  !stop
  !∆`=  1  for`≤50,  ∆`=  4  for50< `≤200,  ∆`=  12  for  200< `≤500,  ∆`=  24for 500< `≤2000,  and finally ∆`= 40 for` >2000       

  !this is for reading in the files
  lmax = 5000

  lmin = 2
  allocate(Cl(4,2:lmax))
  allocate(Cll(4,2:lmax))
  allocate(invCll(2,2,2:lmax))
  allocate(Cllm(2,2,2:lmax))
  allocate(pClpp(3,2:lmax))

  Folder1 = 'SOspectra/'
  Clfile = trim(Folder1)//trim('SOspectra_lenspotentialCls.dat')
  Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
  !from Alex:
  !Clfile = trim(Folder1)//trim('cosmo2017_10K_acc3_lenspotentialCls.dat')
  !Cllfile = trim(Folder1)//trim('cosmo2017_10K_acc3_lensedCls.dat')
  !open(unit=17,file = Clfile, status='old')
  !open(unit=18,file = Cllfile, status='old')



  call getenv('SCRATCHDIR',tensDir)
  !Will:
  tensDir = TRIM(tensDir)//'/Data/alphaBetaDir/'
  !Daan:
  !tensDir = '/mnt/raid-cita/meerburg/SO_forecasts/alphabeta'
  write(*,*) tensDir      
  !allocate and read bessel transforms
  !you should check the subroutine to see if your file is directed correctly
  !it is now set to my directory 
  !allocate and read bessel transforms
  !you should check the subroutine to see if your file is directed correctly
  !it is now set to my directory
  !Will
  alphabetafile = TRIM(tensDir)//'/l_r_alpha_beta.txt.MAX4000'
  alphabetaPolfile = TRIM(tensDir)//'/l_r_alpha_beta_Pol.txt.MAX4000'
  !Daan
  !alphabetafile = TRIM(tensDir)//'l_r_alpha_beta_new_Lmax5000.txt'
  !alphabetaPolfile = TRIM(tensDir)//'l_r_gamma_delta_new_Lmax5000.txt'
  !P%flagDoWigner=flagDoWigner
  if (shape .eq. 5) then
     P%flagDoWigner = 5
  endif

  P%flagDoWigner=0
  P_ISWlens%flagDoWigner = 5
  call allocate_besseltransforms(P,alphabetafile,alphabetaPolfile,cllFile,ClFile)
  call allocate_besseltransforms(P_ISWlens,alphabetafile,alphabetaPolfile,cllFile,ClFile)
  P_ISWlens%flagDoWigner = 0

  open(unit=17,file = Clfile, status='old')
  open(unit=18,file = Cllfile, status='old')
  Cllm(:,:,:) = 0 
  invCll(:,:,:) = 0
  do j = 1, lmax
     !#    L    TT             EE             BB             TE 
     !#    L    TT             EE             BB             TE             PP             TP             EP
     if (j .eq. 1) then
        read(17,*)
        read(18,*)
        cycle
     endif

     read(17,*) l1, Cl(1:4,j),pClpp(1:3,j)
     read(18,*) l1, Cll(1:4,j)
     Cll(1:4,j) = 2.*pi*Cll(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     Cl(1:4,j) = 2.*pi*Cl(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     pClpp(1,j) = 2.*pi*pClpp(1,j)/(l1*(l1+1.))**2
     pClpp(2:3,j) = 2.*pi*pClpp(2:3,j)/(l1*(l1+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

    detCovCV = (Cll(1,j)*Cll(2,j)-Cll(4,j)**2) 
    !inverse covariance in CV limit 
    invCll(1,1,j) = Cll(2,j)/detCovCV
    invCll(1,2,j) = -Cll(4,j)/detCovCV
    invCll(2,1,j) = -Cll(4,j)/detCovCV
    invCll(2,2,j) = Cll(1,j)/detCovCV


    Cllm(1,1,j) = Cll(1,j)!/detCovCV
    Cllm(1,2,j) = Cll(4,j)!/detCovCV
    Cllm(2,1,j) = Cll(4,j)!/detCovCV
    Cllm(2,2,j) = Cll(2,j)!/detCovCV


     !write(*,*) l1,Cll(1:4,j) 
  enddo
  close(17)
  close(18)

  intmax = 39
  !lmax = ellar(intmax)
  !lmax = 2000
  lmin = 2

  write(*,*) 'lmin:', lmin
  write(*,*) 'lmax:', lmax

  !setting everything to zero:




  open(unit=12,file='BispectrumCovariance_ISWmarg_1.0.dat', status = 'replace')
  !want to use analytical approximation of wigner3J (slightly faster)
  AWigner = .True.
  
  !if(want_ISW_correction) then   
  do j = 200, 1000, 400
     TotSumCV = 0.d0
     TotSumCVISWLens = 0.d0
     TotSumCVLensCross = 0.d0
     !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
     !$OMP PRIVATE(l1,l2,l3, min_l,max_l,m1,p1,q1,m2,p2,q2), &
     !$OMP PRIVATE(Det,TempCovCV,atj,bis,bis_ISWlens,tmpPrefac) &
     !$OMP PRIVATE(DetISWLens,DetLensCross,fnlISW,fnl,bispectrum,bispectrum_ISWlens) &
     !$OMP REDUCTION(+:TotSumCV,TotSumCVISWLens,TotSumCVLensCross) 
     do l1 = lmin, j
        allocate(bis(nfields**3,j))
        allocate(bis_ISWlens(nfields**3,j))

        allocate(bispectrum(nfields,nfields,nfields,1,j))
        allocate(bispectrum_ISWlens(nfields,nfields,nfields,1,j))       
        do l2 =  max(lmin,l1), j
           min_l = max(abs(l1-l2),l2)
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if
           max_l = min(j,l1+l2)

           bis = 0
           bis_ISWlens = 0
           bispectrum = 0
           bispectrum_ISWlens = 0
           call get_bispectrum_sss(P,l1,l2,2,j,shape,nfields,bis)
           call get_bispectrum_sss(P_ISWlens,l1,l2,2,j,5,nfields,bis_ISWlens)
           call GetThreeJs(atj(abs(l2-l1)),l1,l2,0,0)
           call reshapeBispectrum(bis,bispectrum,1,minfields,nfields)
           call reshapeBispectrum(bis_ISWlens,bispectrum_ISWlens,1,minfields,nfields)
           do l3=min_l,max_l, 2 !sum has to be even
              tmpPrefac = (atj(l3)*prefactor(l1,l2,l3)*.5)**2/tr(l1,l2,l3)
              do m2  = minfields,nfields !T,E (8 terms only)
                do p2 = minfields,nfields !T,E
                  do q2 = minfields,nfields !T,E
                    do m1  = minfields,nfields !T,E (8 terms only)
                      do p1 = minfields,nfields !T,E
                        do q1 = minfields,nfields !T,E
                          TotSumCVLensCross = TotSumCVLensCross+bispectrum_ISWlens(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                          invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum(m2,p2,q2,1,l3)*tmpPrefac
                          TotSumCVISWLens = TotSumCVISWLens+bispectrum_ISWlens(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                          invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum_ISWlens(m2,p2,q2,1,l3)*tmpPrefac
                          TotSumCV = TotSumCV +bispectrum(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                          invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum(m2,p2,q2,1,l3)*tmpPrefac
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              !signal squared (in SW limit) 
              !fnl = floc(l1,l2,l3)*atj(l3)*prefactor(l1,l2,l3)
              ! fnl = bis(1,l3)*atj(l3)*prefactor(l1,l2,l3)*.5 ! 1/2 to account for term in bis
              ! Det = fnl*fnl
              ! !auto primordial 
              ! TempCovCV = 1.d0/Cll(1,l1)/Cll(1,l2)/Cll(1,l3)

              ! !fnl auto 
              ! TotSumCV = TotSumCV + Det*TempCovCV/tr(l1,l2,l3)

              ! !fnlISW = fPhiISW(l1,l2,l3,pClpp(2,l2),Cll(1,l3),1) + fPhiISW(l1,l3,l2,pClpp(2,l3),Cll(1,l2),1) + fPhiISW(l2,l1,l3,pClpp(2,l2),Cll(1,l3),1) + &
              ! !     fPhiISW(l3,l2,l1,pClpp(2,l2),Cll(1,l1),1) + fPhiISW(l3,l1,l2,pClpp(2,l1),Cll(1,l2),1) + fPhiISW(l2,l3,l1,pClpp(2,l3),Cll(1,l1),1)
              ! fnlISW = bis_ISWlens(1,l3)*.5 ! 1/2 to account for term in bis
              ! DetISWLens = fnlISW*fnlISW*atj(l3)**2*prefactor(l1,l2,l3)**2!*.25
              ! !auto lensing
              ! TotSumCVISWLens = TotSumCVISWLens + DetISWLens*TempCovCV/tr(l1,l2,l3)

              ! !ISW-lensing x primordial and ISW-reinization x primordial (any shape)
              ! DetLensCross =fnlISW*fnl*atj(l3)*prefactor(l1,l2,l3)

              ! TotSumCVLensCross = TotSumCVLensCross + DetLensCross*TempCovCV/tr(l1,l2,l3)                             

           enddo !l3 loop
        enddo !l2 loop
        deallocate(bis)
        deallocate(bis_ISWlens)
        deallocate(bispectrum)
        deallocate(bispectrum_ISWlens)
     enddo !L1 loop
     !$OMP END PARAllEl DO


     DetFishCV = TotSumCVISWLens*TotSumCV -TotSumCVLensCross**2
     alpha = TotSumCVISWLens/DetFishCV !C/det
     beta = -TotSumCVLensCross/DetFishCV !A/det 
     if(.not. want_ISW_correction) then
        alpha=1./TotSumCV
        beta = 0.
     endif
     write(*,*) 'lensing-ISW-fnl_local correlation coefficient:', TotSumCVLensCross/TotSumCV**(1./2)/TotSumCVISWLens**(1./2)
     write(*,*) 'Fisher local error:', 1/TotSumCV**(1./2)
     write(*,*) 'Fisher ISW-lensing error:', 1/TotSumCVISWLens**(1./2)
     write(*,*) 'alpha', alpha, 'beta', beta
     !endif

!!!!!PDM jul 2019
!!!!!testing; with Fisher code I find fnl_template error of 16.658 for lmax = 500 and for 44.511 lmax = 250

     !call get_bispectrum_sss(P,l1,l2,2,lmax,shape,nfields,bis)

     DB = 0.d0
     DSNGauss = 0.d0
     DSNonGauss = 0.d0
     SumGauss = 0.d0
     SumNGauss = 0.d0
     TotNoise = 0.d0
     TotSumGauss_outer = 0.d0
     TotSumNGauss_outer = 0.d0
     do i = 1, intmax !multiples of 32
        l1 = ellar(i)
        !write(*,*) 'l1:', l1
        TotSumGauss = 0.d0
        TotSumNGauss = 0.d0

        allocate(a3j(2*j,2*j))
        allocate(a3joC(2*j,2*j))
        allocate(bispectrum(nfields,nfields,nfields,j,j))
        allocate(bispectrum_ISWlens(nfields,nfields,nfields,j,j))
        allocate(bis(nfields**3,j))
        allocate(bis_ISWlens(nfields**3,j))
        a3joC(:,:) = 0
        do l2 = lmin, j
           bis = 0
           call get_bispectrum_sss(P,l1,l2,2,j,shape,nfields,bis)
           ! bispectrum(1,1,1,l2,:) = bis(1,:)*.5 ! .5 for account for 2 in prefactor
           call reshapeBispectrum(bis*.5,bispectrum,l2,minfields,nfields)
           bis_ISWlens = 0
           call get_bispectrum_sss(P_ISWlens,l1,l2,2,j,5,nfields,bis_ISWlens)
           call reshapeBispectrum(bis_ISWlens*.5,bispectrum_ISWlens,l2,minfields,nfields)

           min_l = max(abs(l1-l2),lmin) 
           max_l = min(j,l1+l2)
           do l3=min_l,max_l,2 !sum has to be even
              call applyInvC(bispectrum(:,:,:,l2,l3),invCll(:,:,l1),invCll(:,:,l2),invCll(:,:,l3),minfields,nfields)
              call applyInvC(bispectrum_ISWlens(:,:,:,l2,l3),invCll(:,:,l1),invCll(:,:,l2),invCll(:,:,l3),minfields,nfields)
          enddo
          ! write(*,*),bispectrum_ISWlens(1,1,1,l2,l3),bispectrum_ISWlens(1,1,1,l3,l2)

           ! bispectrum_ISWlens(1,l2,:) = bis_ISWlens(1,:)*.5 ! .5 for account for 2 in prefactor

           if (AWigner) then
              min_l = max(abs(l1-l2),lmin)  
              if (mod(l1+l2+min_l,2)/=0) then
                 min_l = min_l+1 !l3 should only lead to parity even numbers
              end if
              max_l = min(j,l1+l2)
              do l3=min_l,max_l,2 !sum has to be even

                 a3j(l2,l3) = wigner3jm0(l1,l2,l3)
              enddo
           else
              call GetThreeJs(atj(abs(l2-l1)),l1,l2,0,0)
              a3j(l2,1:2*j) = atj(1:2*j)
           endif
          
           if (nfields .gt. 1) then
             call GetThreeJs(atj(abs(l2-l1)),l1,l2,0,-2)
             a3joC(l2,1:2*j) = atj(1:2*j)
          endif
        enddo
        do l2 = lmin, j
           max_l = min(j,l1+l2)
           min_l = max(abs(l1-l2),lmin)  
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if
           do l3 = min_l,max_l,2
              call permuteBis(bispectrum_ISWlens,minfields,nfields,l2,l3)
              call permuteBis(bispectrum,minfields,nfields,l2,l3)
              ! bispectrum(1,l3,l2) = bispectrum(1,l2,l3)
              ! bispectrum_ISWlens(1,l3,l2) = bispectrum_ISWlens(1,l2,l3)
              !bispectrum(1,l3,l2) = floc(l1,l2,l3)
              !bispectrum(1,l2,l3) = floc(l1,l2,l3)
              !write(*,*) bispectrum(1,l2,l3),bispectrum(1,l3,l2),floc(l1,l2,l3),floc(l1,l3,l2)
           enddo
        enddo
        !do j = 1, intmax !l2 loop

        !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
        !$OMP PRIVATE(l2,l3,l2b,l3b,min_l,max_l,min_lb,max_lb,DB,m1,p1,q1,m2,p2,q2), &
        !$OMP PRIVATE(DSNGauss,DSNonGauss,sigsq,fnl,fnlb,fnlISW,fnlISWb,atj,atj2,tempfac,tempfacFcM),&
        !$OMP REDUCTION(+:TotSumNGauss,TotSumGauss,SumGauss, SumNGauss, TotNoise)
        do l2 = lmin, j
           !l2 = ellar(j)
           min_l = max(abs(l1-l2),lmin)
           !write(*,*) l2
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if

           max_l = min(j,l1+l2)
           !checking
           !call GetThreeJs(atj2(abs(l2-l1)),l1,l2,0,0)
           do l3=min_l,max_l, 2 !sum has to be even
              ! if (bispectrum_ISWlens(1,l2,l3) .eq. 0) then
              !    write(*,*) l2,l3,bispectrum_ISWlens(1,l2,l3)
              ! endif
              ! if (bispectrum(1,l2,l3) .eq. 0) then
              !    write(*,*) l2,l3,bispectrum(1,l2,l3)
              ! endif
              l1b=l1 !l1 = l1'

              do l2b = lmin, j
                 !l2b = ellar(k)
                 min_lb= max(abs(l1b-l2b),lmin)                      
                 !below only relevant if there would be another Wigner3J. 
                 if (mod(l1b+l2b+min_lb,2)/=0) then
                    min_lb = min_lb+1 !l3 should only lead to parity even numbers
                 end if
                 max_lb = min(j,l1b+l2b)

                 do l3b=min_lb,max_lb, 2 !min_lb,max_lb, 2 !sum has to be even
                    !4 possible permutations of the second index
                    ! DB = Cll(1,l2)*Cll(1,l2b)*4.0*a3j(l2,l3)*a3j(l2b,l3b)*FcM(l3,l1,l2)*FcM(l3b,l1,l2b)/(2*l1+1.d0)/Cll(1,l3)/Cll(1,l3b)/Cll(1,l2)/Cll(1,l2b)                  
                    ! fnlISW = bispectrum_ISWlens(1,l2,l3)
                    ! fnlISWb = bispectrum_ISWlens(1,l2b,l3b)
                    ! fnl = bispectrum(1,l2,l3)
                    ! fnlb = bispectrum(1,l2b,l3b)
                    ! tempfac = a3j(l2,l3)*prefactor(l1,l2,l3)*a3j(l2b,l3b)*prefactor(l1b,l2b,l3b)

                    ! sigsq = tempfac*( alpha**2*fnl*fnlb &
                    !      + alpha*beta*fnl*fnlISWb + alpha*beta*fnlb*fnlISW + beta**2*fnlISW*fnlISWb)
                    ! !delta (N)^2 . nine possible permutations of the first index
                    ! DSNonGauss = 9.d0*pClpp(1,l1)/Cll(1,l1)*sigsq*DB/36.!*dellar(j)*dellar(k)/36.!/tr(l1,l2,l3)/tr(l1b,l2b,l3b)
                    tempfacFcM(:,:) = 0
                    tempfac = 4.0*FcM(l3,l1,l2)*FcM(l3b,l1,l2b)/(2*l1+1.d0)*9.d0*pClpp(1,l1)/36.
                    tempfacFcM(1,1) = a3j(l2,l3)*a3j(l2b,l3b)*tempfac
                    if (nfields .gt. 1) then
                      tempfacFcM(2,1) = a3joC(l2,l3)*a3j(l2b,l3b)*tempfac
                      tempfacFcM(1,2) = a3j(l2,l3)*a3joC(l2b,l3b)*tempfac
                      tempfacFcM(2,2) = a3joC(l2,l3)*a3joC(l2b,l3b)*tempfac
                    endif
                    DSNonGauss = 0.d0
                    DSNGauss = 0.d0

                    tempfac = a3j(l2,l3)*prefactor(l1,l2,l3)*a3j(l2b,l3b)*prefactor(l1b,l2b,l3b)
                    do m2  = minfields,nfields !T,E (8 terms only)
                      do p2 = minfields,nfields !T,E
                        do q2 = minfields,nfields !T,E
                          do m1  = minfields,nfields !T,E (8 terms only)
                            do p1 = minfields,nfields !T,E
                              do q1 = minfields,nfields !T,E
                                  DB = Cllm(m1,m2,l1)*Cllm(p1,q1,l2)*Cllm(p2,q2,l2b)*tempfacFcM(q1,q2)
                                  fnlISW = bispectrum_ISWlens(m1,p1,q1,l2,l3)
                                  fnlISWb = bispectrum_ISWlens(m2,p2,q2,l2b,l3b)
                                  fnl = bispectrum(m1,p1,q1,l2,l3)
                                  fnlb = bispectrum(m2,p2,q2,l2b,l3b)
                                  sigsq = tempfac*( alpha**2*fnl*fnlb &
                                       + alpha*beta*fnl*fnlISWb + alpha*beta*fnlb*fnlISW + beta**2*fnlISW*fnlISWb)
                                  !endif
                                  if ((l1.eq.l1b) .and. (l2 .eq.l2b) .and. (l3 .eq.l3b)) then
                                    DSNGauss = sigsq*Cllm(m1,m2,l1)*Cllm(p1,q1,l2)*Cllm(p2,q2,l3)
                                  endif
                                  !delta (N)^2 . nine possible permutations of the first index
                                  DSNonGauss = DSNonGauss+ sigsq*DB!*dellar(j)*dellar(k)/36.!/tr(l1,l2,l3)/tr(l1b,l2b,l3b)
                                  !l3b = l3                                  
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo


                    !<N^2> + delta <N^2>
                    TotNoise = TotNoise + DSNonGauss
                    TotSumGauss = TotSumGauss + DSNGauss
                    TotSumNGauss = TotSumNGauss + DSNonGauss+DSNGauss
                    !delta (N)^2 Non-Gaussian covariance
                    SumGauss = SumGauss + DSNGauss
                    SumNGauss =  SumNGauss + DSNGauss + DSNonGauss

                 enddo !l3b
              enddo !l2b

           enddo !l3
        enddo !l2
        !$OMP END PARAllEl DO

        deallocate(bis)
        deallocate(bispectrum)
        deallocate(bis_ISWlens)
        deallocate(bispectrum_ISWlens)
        deallocate(a3j)
        deallocate(a3joC)
        TotSumGauss_outer = TotSumGauss_outer + TotSumGauss
        TotSumNGauss_outer = TotSumNGauss_outer + TotSumNGauss
        !write(*,*) l1,TotSumNGauss,TotSumGauss, TotSumCV/(TotSumNGauss-TotSumGauss+TotSumCV)
        !write(*,'(I4,3E16.7)') l1,TotSumNGauss_outer,TotSumGauss_outer, TotSumCV/(TotSumNGauss_outer-TotSumGauss_outer+TotSumCV)
        !write(*,*) l1,TotSumNGauss/alpha,TotSumGauss/alpha, alpha/(TotSumNGauss-TotSumGauss+alpha)
        write(*,'(I4,3E16.7)') l1,TotSumNGauss_outer/alpha,TotSumGauss_outer/alpha, alpha/(TotSumNGauss_outer-TotSumGauss_outer+alpha)

     enddo !l1
     write(12,'(I4,5E18.7)') j, sqrt(SumGauss/SumNGauss), (SumNGauss-SumGauss)/TotSumCV, TotSumNGauss_outer/alpha,TotSumGauss_outer/alpha, alpha/(TotSumNGauss_outer-TotSumGauss_outer+alpha)
     write(*,'(I4,5E18.7)') j, sqrt(SumGauss/SumNGauss), (SumNGauss-SumGauss)/TotSumCV, TotSumNGauss_outer/alpha,TotSumGauss_outer/alpha, alpha/(TotSumNGauss_outer-TotSumGauss_outer+alpha)
  enddo !lmax loop
  !write(*,'(I4,4E17.8)') lmax, (SumNGauss-SumGauss)/TotSumCV
  !write(*,'(I4,4E17.8)') ellar(intmax), SumGauss, SumNGauss, sqrt(SumGauss/SumNGauss)
  !write(*,'(A12,X,I4,X,A19,X,F11.3)') 'For lmax = ',  ellar(intmax), 'the error on fnl = ', sqrt(1./SumGauss)

  close(12)
  deallocate(Cl, Cll, pClpp,invCll,Cllm)

contains

  !subroutine 
  real(dl) function Fc(l1,l2,l3)
    integer :: l1, l2, l3
    Fc = 1.d0/2.d0*(l2*(l2+1.0)+l3*(l3+1.0)-l1*(l1+1.0))* &
         sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0))/sqrt(4.0*pi)*wigner3jm0(l1,l2,l3)
    !*Wig0(l1,l3,l2)
  end function Fc

  real(dl) function FcM(l1,l2,l3)
    integer :: l1, l2, l3
    FcM = 1.d0/2.d0*(l2*(l2+1.0)+l3*(l3+1.0)-l1*(l1+1.0))* &
         sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0))/sqrt(4.0*pi)
  end function FcM


  real(dl) function LargeArg3Js(l1,l2,l3)
    integer, intent (in) :: l1, l2, l3
    integer :: lt
    lt = l1+l2+l3
    !there are also infinities when and of the l is lt/2
    !these are not infinities, but do not know how else to deal with them
    if(l1 .eq. lt/2 .or. l2 .eq. lt/2 .or. l3 .eq. lt/2) lt = lt - 1
    LargeArg3Js = sqrt(2.d0/pi)*(-1.d0)**(lt/2.d0)/ &
         (lt*(lt-2.d0*l1)*(lt-2.d0*l2)*(lt-2.d0*l3))**(1./4.)
    !get rid of NAN's
    if(LargeArg3Js .ne. LargeArg3Js) LargeArg3Js = 0.d0

  end function LargeArg3Js

  subroutine reshapeBispectrum(bis,bispectrum,l2,minfields,nfields)
    real(dl), intent(in) :: bis(:,:)
    integer, intent(in) :: nfields,minfields
    real(dl), intent(out) :: bispectrum(:,:,:,:,:)
    integer l2
    if (nfields .eq. 1) then
      if (minfields .eq.1) then
        bispectrum(1,1,1,l2,:) = bis(1,:)
      else
        bispectrum(2,2,2,l2,:) = bis(8,:)
      endif
    else
      bispectrum(1,1,1,l2,:) = bis(1,:)
      bispectrum(1,1,2,l2,:) = bis(2,:)
      bispectrum(1,2,1,l2,:) = bis(3,:)
      bispectrum(1,2,2,l2,:) = bis(4,:)
      bispectrum(2,1,1,l2,:) = bis(5,:)
      bispectrum(2,1,2,l2,:) = bis(6,:)
      bispectrum(2,2,1,l2,:) = bis(7,:)
      bispectrum(2,2,2,l2,:) = bis(8,:)
    endif
  end subroutine reshapeBispectrum

  subroutine permuteBis(bispectrum,minfields,nfields,l2,l3)
    real(dl), intent(out) :: bispectrum(:,:,:,:,:)
    integer l2,l3,nfields,minfields
    if (nfields .eq. 1) then
        if (minfields .eq. 1) then
          bispectrum(1,1,1,l3,l2) = bispectrum(1,1,1,l2,l3)
        else
          bispectrum(2,2,2,l3,l2) = bispectrum(2,2,2,l2,l3)
        endif
    else
        bispectrum(1,1,1,l3,l2) = bispectrum(1,1,1,l2,l3)
        bispectrum(1,2,1,l3,l2) = bispectrum(1,1,2,l2,l3)
        bispectrum(1,1,2,l3,l2) = bispectrum(1,2,1,l2,l3)
        bispectrum(2,1,1,l3,l2) = bispectrum(2,1,1,l2,l3)
        bispectrum(2,1,2,l3,l2) = bispectrum(2,2,1,l2,l3)
        bispectrum(2,2,1,l3,l2) = bispectrum(2,1,2,l2,l3)
        bispectrum(1,2,2,l3,l2) = bispectrum(1,2,2,l2,l3)
        bispectrum(2,2,2,l3,l2) = bispectrum(2,2,2,l2,l3)
    endif
  end subroutine permuteBis

  subroutine applyInvC(bispectrum,invCl_l1,invCl_l2,invCl_l3,minfields,nfields)
    real(dl), intent(inout) :: bispectrum(:,:,:)
    real(dl), intent(in):: invCl_l1(:,:),invCl_l2(:,:),invCl_l3(:,:)
    integer, intent(in):: nfields,minfields
    integer m1,m2,p1,p2,q1,q2
    real(dl) tmpBis(minfields:nfields,minfields:nfields,minfields:nfields)
    tmpBis(:,:,:) = 0
    do m1 = minfields,nfields
      do m2 = minfields,nfields
        do p1 = minfields,nfields
          do p2 = minfields,nfields
            do q1 = minfields,nfields
              do q2 = minfields,nfields
                tmpBis(q1,p1,m1) = tmpBis(q1,p1,m1) +invCl_l1(q1,q2)*invCl_l2(p1,p2)*invCl_l3(m1,m2)*bispectrum(q2,p2,m2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    bispectrum(:,:,:) = tmpBis(:,:,:)
  end subroutine applyInvC

  ! real(dl) function floc(l1,l2,l3)
  !   !SW approximation 
  !   integer :: l1, l2, l3
  !   real(dl) :: amp
  !   real(dl) :: As = 2.1056d-9
  !   !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
  !   amp = (2.d0/27./pi**2)*As
  !   floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
  !        1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
  !   floc = floc*amp*2.E-7 !2.E-7  is introduced to get roughly same amplitude at l_max = 500 to full fnl_local
  ! end function floc


  ! real function prefactor(l1,l2,l3)
  !   integer, intent(in) :: l1,l2,l3

  !   prefactor = 2.0*sqrt((1./4.)*((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.))/pi)
  ! end function prefactor


  real(dl) function wigner3jm0(l1,l2,l3)
    integer, intent (in) :: l1,l2,l3
    integer :: LtX, LtY, LtZ, Lt
    LtX = -l1+l2+l3
    LtY = l1-l2+l3
    LtZ = l1+l2-l3
    Lt  = l1+l2+l3
    if(mod(Lt,2)/=0) then
       wigner3jm0  = 0.d0
    else
       wigner3jm0 = (-1)**(lt/2) *et(Lt)*rmj(Lt/2)/(rmj(LtX/2)*rmj(LtY/2)*rmj(LtZ/2)) * &
            sqrt(rmj(LtX)*rmj(LtY)*rmj(LtZ)/rmj(Lt+1))
    endif
  end function wigner3jm0

  real(dl) function rmj(x)
    integer :: x
    rmj = Sqrt(pi)*(((8.*x + 4.)*x + 1)*x + 1/30.)**(1./6.)

  end function rmj

  real(dl) function et(L)
    integer :: L
    real(dl) :: temp
    real(dl), parameter :: Euler = 2.7182818284590452353602874713526624977572470937000
    temp = L*Log(L/(L+1.))-Log(L+1.)
    et  = sqrt(Euler)*Exp(temp/2.)
  end function et

  ! real function tr(l1,l2,l3)
  !   integer, intent(in) :: l1,l2,l3
  !   if ((l1.eq.l2).and.(l2.eq.l3)) then
  !      tr  = 6.d0
  !   elseif ((l1.eq.l2).or.(l2.eq.l3).or.(l3.eq.l1)) then 
  !      tr =  2d0
  !   else
  !      tr = 1.d0
  !   endif

  ! end function tr

  !   real(dl) function fPhiISW(l1,l2,l3,CPT_l1,CTT_l3)
  !     integer :: l1, l2, l3
  !     real(dl) :: CPT_l1,CTT_l3
  !     fPhiISW = .5*((l1+1)*l1-l2*(l2+1.d0)+l3*(l3+1.d0))*CPT_l1*CTT_l3
  !   end function fPhiISW

  !   subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
  !     !Recursive evaluation of 3j symbols. Does minimal error checking on input
  !     !parameters.
  !     implicit none
  !     integer, parameter :: dl = KIND(1.d0)
  !     integer, intent(in) :: l2in,l3in, m2in,m3in
  !     real(dl), dimension(*) :: thrcof
  !     INTEGER, PARAMETER :: i8 = selected_int_kind(18)
  !     integer(i8) :: l2,l3,m2,m3
  !     integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2

  !     real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1,sumuni
  !     real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
  !     integer i,ier, index, nlim, sign2
  !     integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
  !     real(dl), parameter :: zero = 0._dl, one = 1._dl
  !     real(dl), parameter ::  tiny = 1.0d-30, srtiny=1.0d-15, huge = 1.d30,srhuge = 1.d15

  !     ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

  !     ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
  !     !                to l1max = l2+l3
  !     ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

  !     ! to achieve the numerical stability, the recursion will proceed
  !     ! simultaneously forwards and backwards, starting from l1min and l1max
  !     ! respectively.
  !     !
  !     ! lmatch is the l1-value at which forward and backward recursion are
  !     ! matched.
  !     !
  !     ! ndim is the length of the array thrcof
  !     !
  !     ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
  !     ! ier = -2 if possible 3j's exceed ndim
  !     ! ier >= 0 otherwise

  !     l2=l2in
  !     l3=l3in
  !     m2=m2in
  !     m3=m3in
  !     newfac = 0
  !     lmatch = 0
  !     m1 = -(m2+m3)

  !     ! check relative magnitude of l and m values
  !     ier = 0

  !     if (l2 < abs(m2) .or. l3 < m3) then
  !        ier = -1
  !        ! call MpiStop('error ier = -1')
  !        print*, 'error ier = -1',l2,abs(m2),l3,m3
  !        stop
  !        return
  !     end if

  !     ! limits for l1
  !     l1min = max(abs(l2-l3),abs(m1))
  !     l1max = l2+l3

  !     if (l1min >= l1max) then
  !        if (l1min/=l1max) then
  !           ier = -1

  !           !call MpiStop('error ier = -1')
  !           print*, 'error ier = -1',l1min,l1max 
  !           stop
  !           return
  !        end if

  !        ! reached if l1 can take only one value, i.e.l1min=l1max
  !        thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
  !        return

  !     end if

  !     nfin = l1max-l1min+1

  !     ! starting forward recursion from l1min taking nstep1 steps
  !     l1 = l1min
  !     thrcof(1) = srtiny
  !     sum1 = (2*l1 + 1)*tiny

  !     lstep = 1

  ! 30  lstep = lstep+1
  !     l1 = l1+1

  !     oldfac = newfac
  !     a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
  !     a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
  !     newfac = sqrt(a2*real(a1,dl))
  !     if (l1 == 1) then
  !        !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
  !        c1 = -(2*l1-1)*l1*(m3-m2)/newfac
  !     else

  !        dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
  !        denom = (l1-1)*newfac

  !        if (lstep > 2) c1old = abs(c1)
  !        c1 = -(2*l1-1)*dv/denom

  !     end if

  !     if (lstep<= 2) then

  !        ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
  !        x = srtiny*c1
  !        thrcof(2) = x
  !        sum1 = sum1+tiny*(2*l1+1)*c1*c1
  !        if(lstep==nfin) then
  !           sumuni=sum1
  !           go to 230
  !        end if
  !        goto 30

  !     end if

  !     c2 = -l1*oldfac/denom

  !     ! recursion to the next 3j-coeff x  
  !     x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
  !     thrcof(lstep) = x
  !     sumfor = sum1
  !     sum1 = sum1 + (2*l1+1)*x*x
  !     if (lstep/=nfin) then

  !        ! see if last unnormalised 3j-coeff exceeds srhuge
  !        if (abs(x) >= srhuge) then

  !           ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
  !           ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
  !           ! HAS TO BE RESCALED TO PREVENT OVERFLOW

  !           ier = ier+1
  !           do i = 1, lstep
  !              if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
  !              thrcof(i) = thrcof(i)/srhuge
  !           end do

  !           sum1 = sum1/huge
  !           sumfor = sumfor/huge
  !           x = x/srhuge

  !        end if

  !        ! as long as abs(c1) is decreasing, the recursion proceeds towards
  !        ! increasing
  !        ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
  !        ! detected, the recursion direction is reversed.

  !        if (c1old > abs(c1)) goto 30

  !     end if !lstep/=nfin

  !     ! keep three 3j-coeffs around lmatch for comparison with backward recursion

  !     lmatch = l1-1
  !     x1 = x
  !     x2 = thrcof(lstep-1)
  !     x3 = thrcof(lstep-2)
  !     nstep2 = nfin-lstep+3

  !     ! --------------------------------------------------------------------------
  !     !
  !     ! starting backward recursion from l1max taking nstep2 stpes, so that
  !     ! forward and backward recursion overlap at 3 points 
  !     ! l1 = lmatch-1, lmatch, lmatch+1

  !     nfinp1 = nfin+1
  !     nfinp2 = nfin+2
  !     nfinp3 = nfin+3
  !     l1 = l1max
  !     thrcof(nfin) = srtiny
  !     sum2 = tiny*(2*l1+1)

  !     l1 = l1+2
  !     lstep=1

  !     do
  !        lstep = lstep + 1
  !        l1= l1-1

  !        oldfac = newfac
  !        a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
  !        a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
  !        newfac = sqrt(a1*real(a2,dl))

  !        dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

  !        denom = l1*newfac
  !        c1 = -(2*l1-1)*dv/denom
  !        if (lstep <= 2) then

  !           ! if l2=l2max+1, the third term in the recursion vanishes

  !           y = srtiny*c1
  !           thrcof(nfin-1) = y
  !           sumbac = sum2
  !           sum2 = sum2 + tiny*(2*l1-3)*c1*c1

  !           cycle

  !        end if

  !        c2 = -(l1-1)*oldfac/denom

  !        ! recursion to the next 3j-coeff y
  !        y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

  !        if (lstep==nstep2) exit

  !        thrcof(nfinp1-lstep) = y
  !        sumbac = sum2
  !        sum2 = sum2+(2*l1-3)*y*y

  !        ! see if last unnormalised 3j-coeff exceeds srhuge
  !        if (abs(y) >= srhuge) then

  !           ! reached if 3j-coeff larger than srhuge so that the recursion series
  !           ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent
  !           ! overflow

  !           ier=ier+1
  !           do i = 1, lstep
  !              index=nfin-i+1
  !              if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
  !              thrcof(index) = thrcof(index)/srhuge
  !           end do

  !           sum2=sum2/huge
  !           sumbac=sumbac/huge

  !        end if

  !     end do

  !     ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
  !     ! corresponding backward recursion vals y1, y2, y3

  !     y3 = y
  !     y2 = thrcof(nfinp2-lstep)
  !     y1 = thrcof(nfinp3-lstep)

  !     ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal
  !     ! error

  !     ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
  !     nlim = nfin-nstep2+1

  !     if (abs(ratio) >= 1) then

  !        thrcof(1:nlim) = ratio*thrcof(1:nlim) 
  !        sumuni = ratio*ratio*sumfor + sumbac

  !     else

  !        nlim = nlim+1
  !        ratio = 1/ratio
  !        do n = nlim, nfin
  !           thrcof(n) = ratio*thrcof(n)
  !        end do
  !        sumuni = sumfor + ratio*ratio*sumbac

  !     end if
  !     ! normalise 3j-coeffs

  ! 230 cnorm = 1/sqrt(sumuni)

  !     ! sign convention for last 3j-coeff determines overall phase

  !     sign1 = sign(one,thrcof(nfin))
  !     sign2 = (-1)**(abs(l2+m2-l3+m3))
  !     if (sign1*sign2 <= 0) then
  !        cnorm = -cnorm
  !     end if
  !     if (abs(cnorm) >= one) then
  !        thrcof(1:nfin) = cnorm*thrcof(1:nfin)
  !        return
  !     end if

  !     thresh = tiny/abs(cnorm)

  !     do n = 1, nfin
  !        if (abs(thrcof(n)) < thresh) thrcof(n) = zero
  !        thrcof(n) = cnorm*thrcof(n)
  !     end do
  !     return 

  !   end subroutine GetThreeJs



end program bisvar
