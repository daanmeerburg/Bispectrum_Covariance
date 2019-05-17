!test program to compute components of the bispectrum covariance from lensing 

program bisvar 
  use fwigxjpf
  implicit none
  integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359

  character(80) :: Folder1, Folder2, Folder3
  character(80) :: Clfile, Cllfile

  real(dl), allocatable :: Cl(:,:), Cll(:,:)
  real(dl), pointer :: pClpp(:,:)
  !integer :: l1, l2, l3
  integer :: lmax, lmin, l1, l2, l3, l1b, l2b, l3b, el(3,6), elb(3,6)
  integer :: min_l, max_l, Lm,min_lb, max_lb
  integer :: i,j, k, l, m, n 
  !wigner 3j
  real(dl)  :: atj(0:20000),atj2(0:20000)
  real(dl), pointer :: a3j(:,:)

  real(dl) :: CMB2COBEnorm = 7428350250000.d0
  real(dl) :: DB(4,36), SumDB(4,36), SumTot, DBtot(4,36),SumTotGauss

  real(dl) :: DSNGauss, DSNonGauss, TotSumGauss, TotSumNGauss
  real(dl) :: sigsq, fnl
  !call fwig_temp_init(2*1000)

  lmax = 5000
  lmin = 2
  allocate(Cl(4,2:lmax))
  allocate(Cll(4,2:lmax))
  allocate(pClpp(3,2:lmax))
  allocate(a3j(2*lmax,2*lmax))

  Folder1 = 'SOspectra/'
  Clfile = trim(Folder1)//trim('SOspectra_lenspotentialCls.dat')
  Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')

  open(unit=17,file = Clfile, status='old')
  open(unit=18,file = Cllfile, status='old')

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
     Cll(1:4,j) = 2.*pi*Cll(1:4,j)/(real(l1,dl)*(real(l1,dl)+1.))/CMB2COBEnorm
     Cl(1:4,j) = 2.*pi*Cl(1:4,j)/(real(l1,dl)*(real(l1,dl)+1.))/CMB2COBEnorm
     pClpp(1,j) = 2.*pi*pClpp(1,j)/(real(l1,dl)*(real(l1,dl)+1.))**2
     pClpp(2:3,j) = 2.*pi*pClpp(2:3,j)/(real(l1,dl)*(real(l1,dl)+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

     !write(*,*) l1,Cll(1:4,j) 
  enddo
  close(17)
  close(18)

  !testing:
!!$  call fwig_table_init(2*lmax,9)
!!$  call fwig_temp_init(2*lmax)
!!$  l1 = 11
!!$  call deltaB1(l1,12,11,17,8,Clpp(1,:),5000,DB(1))
!!$  call deltaB2(l1,12,11,17,8,Clpp(1,:),5000,DB(2))
!!$  call deltaB3(l1,12,11,17,8,Clpp(1,:),5000,DB(3))
!!$  call deltaB4(l1,12,11,17,8,Clpp(1,:),5000,DB(4))
!!$  write(*,*) DB(1:4)
!!$  call fwig_temp_free();
!!$  call fwig_table_free();
!!$  stop

  !for fun, let us do a loop. Diagonal first.
  !lower lmax for sake of time here. On Lobster, this loop takes 52 minutes.

  !note also that I did not seperatly apply the filter that would introduce another Wigner3j
  !(is this correct?). This would lower the number of sample points. 
  lmax = 1000
  lmin = 20
  open(unit=12,file='SNratio_v1.3.txt', status = 'replace')
  call fwig_table_init(2*lmax+2,9)
  !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
  !$OMP PRIVATE(Lm,lmax,l1,l2,l3,l1b,l2b,l3b,min_l,max_l,min_lb,max_lb,DB,a3j,atj2,i,j,k,l,m,n, el, elb), &
  !$OMP PRIVATE(DSNGauss,DSNonGauss,sigsq,fnl,atj),&
  !$OMP REDUCTION(+:SumDB,Sumtot,TotSumGauss,TotSumNGauss,SumTotGauss)
  !do l1 = lmin, lmax

  do Lm = 450,500,10
     lmax = Lm
     call fwig_thread_temp_init(2*lmax)
     DB = 0.d0
     SumDB(1:4,1:36) = 0.d0
     Sumtot = 0.d0
     SumTotGauss = 0.d0
     DSNGauss = 0.d0
     DSNonGauss = 0.d0
     TotSumGauss = 0.d0
     TotSumNGauss = 0.d0
     !call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
     do l1 = lmin, lmax
        write(*,*) lmin,lmax,l1
        !call calcWigners2D(l1,lmin,lmax,a3j)
        do l2 =  max(lmin,l1), lmax
           min_l = max(abs(l1-l2),l2)
           !below only relevant if there would be another Wigner3J. 
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if
           max_l = min(lmax,l1+l2)
           call GetThreeJs(atj(abs(l2-l1)),l1,l2,0,0)
           do l3=min_l,max_l, 2 !sum has to be even
            !diagonal
            l1b=l1
            do l2b =  max(lmin,l1b), lmax
             min_lb= max(abs(l1-l2),l2)
             !below only relevant if there would be another Wigner3J. 
             if (mod(l1b+l2b+min_lb,2)/=0) then
                min_lb = min_lb+1 !l3 should only lead to parity even numbers
             end if
             max_lb = min(lmax,l1b+l2b)
             call GetThreeJs(atj2(abs(l2b-l1b)),l1b,l2b,0,0)
             do l3b=min_lb,max_lb, 2 !sum has to be even
                ! l3b=l3+100
                ! l2b=l2+100
                call assignElls(el,l1,l2,l3)
                call assignElls(elb,l1b,l2b,l3b)
                !permutations (total of 36 because only 5 ell are permutable)
                n = 1
                do i = 1,6
                   do j = 1,6
                      !call deltaB4(el(1,i),el(2,i),el(3,i),elb(2,j),elb(3,j),Cll(1,el(1,i)),Cll(1,el(2,i)),Cll(1,elb(2,j)),Clpp(1,el(1,i)),5000,DB(4,n))
                      call deltaB4(el(1,i),el(2,i),el(3,i),elb(2,j),elb(3,j),Cll(1,el(1,i)),Cll(1,el(2,i)),Cll(1,elb(2,j)),pClpp,5000,DB(4,n))
                      n = n + 1
                   enddo
                enddo
                
                !signal squared (in SW limit)
                fnl = floc(l1,l2,l3)*atj(l3)*prefactor(l1,l2,l3)

                sigsq = fnl*floc(l1b,l2b,l3b)*atj2(l3b)*prefactor(l1b,l2b,l3)
                !sigsq = fnl**2

                !delta (S/N)^2 Gaussian covariance 
                if ((l1.eq.l1b) .and. (l2 .eq.l2b) .and. (l3 .eq.l3b)) then
                  !write(*,*),l1,l2,l3,l1b,l2b,l3b,sigsq/Cll(1,l1)/Cll(1,l2)/Cll(1,l3),atj(l3),atj2(l3b)
                  DSNGauss = sigsq/Cll(1,l1)/Cll(1,l2)/Cll(1,l3)
                else
                  DSNGauss = 0
                  SumTotGauss = SumTotGauss + 1/(Cll(1,l1)*Cll(1,l2)*Cll(1,l3)) !Sum(DBtot)
                endif
                !delta (S/N)^2 Non-Gaussian covariance
                DSNonGauss = sigsq*sum(DB(4,1:36))/(Cll(1,l1)*Cll(1,l2)*Cll(1,l3)*Cll(1,l1b)*Cll(1,l2b)*Cll(1,l3b))
                !endif
                !assuming all are multiplied by Cl1Cl2Cl3 (which is true except for the last term)
                !write(*,*) l1, l2, l3, l2b, l3b, DB1, DB2, DB3

                TotSumGauss = TotSumGauss + DSNGauss
                TotSumNGauss = TotSumNGauss + DSNonGauss
                !SumDB(4,1:36) = SumDB(4,1:36) + DB(4,1:36)
                !DB(4,1:36) = DB(4,1:36)/Cll(1,l1)/Cll(1,l2)/Cll(1,l3)
                SumTot = SumTot+ sum(DB(4,1:36))/(Cll(1,l1)*Cll(1,l2)*Cll(1,l3)*Cll(1,l1b)*Cll(1,l2b)*Cll(1,l3b)) ! sum(DB(4,1:36))
                
                !write(*,'(3I4,3E17.8)') l1,l2,l3,SumDB(1:3)
              enddo !l3b
            enddo !l2b
            ! if (TotSumGauss .ne. 0) then
            !  write(*,*) TotSumGauss, TotSumNGauss
            ! endif
           enddo !l3
        enddo !l2
        !write(*,'(I4,25E17.8)') l1, SumTot, SumDB(1:4,1:6)

        
     enddo !l1
     write(12,'(I4,3E17.8)') Lm, TotSumGauss**(-1.d0/2.d0),TotSumNGauss**(-1.d0/2.d0), (TotSumGauss/TotSumNGauss)**(-1.d0/2.d0)
     write(*,'(I4,6E17.8)') Lm, TotSumGauss**(-1.d0/2.d0),TotSumNGauss**(-1.d0/2.d0), (TotSumGauss/TotSumNGauss)**(-1.d0/2.d0), SumTot, SumTotGauss
     call fwig_temp_free();
  enddo
  !$OMP END PARAllEl DO
  close(12)
  call fwig_table_free();
  deallocate(Cl, Cll, pClpp)
  deallocate(a3j)
contains
  !Eq. (29) Notes
  subroutine deltaB1(l1,l2a,l3a,l2b,l3b,Cll,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax),Cll(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0
    stepsum = 0.d0
    l_min = Max(abs(l2b-l2a),2)
    l_min = Max(abs(l3b-l3a),l_min)
    l_max = Min(abs(l2a+l2b),lmax)
    l_max = Min(abs(l3a+l3b),l_max)
    !if l2a+l2b = even/odd and  l3a+l3b = odd/even -> 0 because of conflicting wignersJs
    !hence the folliwing is always sufficient 
    if (mod(l2a+l2b+l_min,2)/=0) then
       l_min = l_min+1 
    end if
    if  (mod(l2a+l2b+l3a+l3b,2)/=0) then
       DB = 0.d0 !note that this is extra, since I am already doing this in do loop. Can be
       !removed. Might speed up code. 
    else !otherwise do the sum 
       do i = l_min, l_max, 2 !L
          stepsum = fwig6jj(2* l1 , 2* l3a , 2* l2a , 2*  i,  2*  l2b , 2*  l3b )*Clpp(i)* &
               Fc(l2b,i,l2a)*Fc(l3b,i,l3a)
          DB = DB + stepsum
          !write(*,*) DB1, stepsum, Clpp(i)
       enddo
    endif
    DB = Cll(l1)*Cll(l2a)*Cll(l3a)*(-1.d0)**(l3a+l2b)*DB
  end subroutine deltaB1

  !Eq. (31) Notes; unfortunately we have to rewrite these. I dont think we can use the same code, or can we???

  subroutine deltaB2(l1,l2a,l3a,l2b,l3b,Cll,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax),Cll(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0
    stepsum = 0.d0
    l_min = Max(abs(l2b-l3a),2)
    l_min = Max(abs(l3b-l2a),l_min)
    l_max = Min(abs(l2a+l3b),lmax)
    l_max = Min(abs(l3a+l2b),l_max)
    !if l2a+l3b = even/odd and  l3a+l2b = odd/even -> 0 because of conflicting wignersJs
    !hence the folliwing is always sufficient 
    if (mod(l2a+l3b+l_min,2)/=0) then
       l_min = l_min+1 
    end if
    if  (mod(l2a+l2b+l3a+l3b,2)/=0) then
       DB = 0.d0
    else !otherwise do the sum 
       do i = l_min, l_max, 2 !L
          stepsum = fwig6jj(2* l1 , 2* l2a , 2* l3a , 2*  i,  2*  l2b , 2*  l3b )*Clpp(i)* &
               Fc(l2b,i,l3a)*Fc(l3b,i,l2a)
          DB = DB + stepsum
          !write(*,*) DB1, stepsum, Clpp(i)
       enddo
    endif
    DB = Cll(l1)*Cll(l2a)*Cll(l3a)*(-1.d0)**(l2a+l2b)*DB
  end subroutine deltaB2

  !Eq. (33) Notes
  subroutine deltaB3(l1,l2a,l3a,l2b,l3b,Cll,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax),Cll(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0
    stepsum = 0.d0
    l_min = Max(abs(l3b-l2a),2)
    l_min = Max(abs(l3a-l2b),l_min)
    l_max = Min(abs(l2a+l3b),lmax)
    l_max = Min(abs(l3a+l2b),l_max)
    !if l2a+l3b = even/odd and  l3a+l2b = odd/even -> 0 because of conflicting wignersJs
    !hence the folliwing is always sufficient 
    if (mod(l2a+l3b+l_min,2)/=0) then
       l_min = l_min+1 
    end if
    if  (mod(l2a+l2b+l3a+l3b,2)/=0) then
       DB = 0.d0
    else !otherwise do the sum 
       do i = l_min, l_max, 2 !L
          stepsum = fwig6jj(2* l1 , 2* l2a , 2* l3a , 2*  i,  2*  l2b , 2*  l3b )*Clpp(i)* &
               Fc(l3b,i,l2a)*Fc(l3a,i,l2b)
          DB = DB + stepsum
          !write(*,*) DB1, stepsum, Clpp(i)
       enddo
    endif
    DB = Cll(l1)*Cll(l2a)*Cll(l2b)*(-1.d0)**(-l2a-l2b)*DB
  end subroutine deltaB3

  !Eq. (35) Notes
  subroutine deltaB4x(l1,l2a,l3a,l2b,l3b,Cll1,Cll2,Cll3,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp,Cll1,Cll2,Cll3
    real(dl), intent(out) :: DB
    if  (mod(l1+l2a+l3a,2)/=0 .or. mod(l1+l2b+l3b,2)/=0) then
       DB = 0.d0
    else
       DB = Clpp*Fc(l3a,l1,l2a)*Fc(l3b,l1,l2b)/(2*l1+1.d0)*Cll1*Cll2*Cll3
       !*Cll(l1)*Cll(l2a)*Cll(l2b)
    endif

  end subroutine deltaB4x

  subroutine deltaB4(l1,l2a,l3a,l2b,l3b,Cll1,Cll2,Cll3,b,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Cll1,Cll2,Cll3
    real(dl), pointer :: b(:,:)
    real(dl), intent(out) :: DB
    if  (mod(l1+l2a+l3a,2)/=0 .or. mod(l1+l2b+l3b,2)/=0) then
       DB = 0.d0
    else
       DB = b(1,l1)*Fc(l3a,l1,l2a)*Fc(l3b,l1,l2b)/(2*l1+1.d0)*Cll1*Cll2*Cll3
       !*Cll(l1)*Cll(l2a)*Cll(l2b)
    endif

  end subroutine deltaB4

  !Eq. (35) Notes
  subroutine deltaB4p(l1,l2a,l3a,l2b,l3b,a3j,Cll,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax),Cll(2:lmax), a3j(:)
    real(dl), intent(out) :: DB
    if  (mod(l1+l2a+l3a,2)/=0 .or. mod(l1+l2b+l3b,2)/=0) then
       DB = 0.d0
    else
       DB = Clpp(l1)*a3j(l2a)*a3j(l2b)*FcM(l3a,l1,l2a)*FcM(l3b,l1,l2b)/(2*l1+1.d0)*Cll(l1)*Cll(l2a)*Cll(l2b)
    endif

  end subroutine deltaB4p

  !subroutine 
  real(dl) function Fc(l1,l2,l3)
    integer :: l1, l2, l3
    Fc = 1.d0/2.d0*(l2*(l2+1.0)+l3*(l3+1.0)-l1*(l1+1.0))* &
         sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0))/sqrt(4.0*pi)*LargeArg3Js(l1,l2,l3)
    !*Wig0(l1,l3,l2)
  end function Fc

  real(dl) function FcM(l1,l2,l3)
    integer :: l1, l2, l3
    FcM = 1.d0/2.d0*(l2*(l2+1.0)+l3*(l3+1.0)-l1*(l1+1.0))* &
         sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0))/sqrt(4.0*pi)
  end function FcM

  real(dl) function Wig0(l1,l2,l3)
    integer :: l1,l2,l3
    Wig0 = fwig3jj(2* l1 , 2* l2 , 2* l3 , 2* 0, 2* 0 , 2* 0)

  end function Wig0

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
  subroutine assignElls(el,l1,l2,l3)
    integer, intent(in):: l1,l2,l3
    integer, intent(out) :: el(3,6)
    el(1,1:2) = l1
    el(1,3:4) = l2
    el(1,5:6) = l3
    el(2,1) = l1
    el(2,2) = l3
    el(2,3) = l1
    el(2,4) = l3
    el(2,5) = l1
    el(2,6) = l2
    el(3,1) = l3
    el(3,2) = l2
    el(3,3) = l3
    el(3,4) = l1
    el(3,5) = l2
    el(3,6) = l1
  end subroutine assignElls

  ! subroutine calcWigners2D(l1,lmin,lmax,a3jArr)
  !   integer, intent (in) :: l1, lmax,lmin
  !   integer ::  l2,min_l
  !   real(dl), pointer :: a3jArr(:,:)
  !   do l2 = lmin,lmax
  !      !call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
  !      call GetThreeJs(a3jArr(l2,abs(l2-l1)),l1,l2,0,0)
  !   a3jArr = a3jArr + transpose(a3jArr)
  !   end do
  ! end subroutine calcWigners2D


  real(dl) function floc(l1,l2,l3)
    !SW approximation 
    integer :: l1, l2, l3
    real(dl) :: amp
    real(dl) :: As = 2.1056d-9
    !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
    amp = (2.d0/27./pi**2)*As
    floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
         1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
    floc = floc*amp
  end function floc

  subroutine calcWigners2D(l1,lmin,lmax,a3jArr)
   integer, intent (in) :: l1, lmax,lmin
   integer ::  l2, l3,min_l
   real(dl) :: a3jArr(2*lmax,2*lmax)
   do l2 = lmin,lmax
      min_l = max(abs(l1-l2),l2)
      !below only relevant if there would be another Wigner3J.
      if (mod(l1+l2+min_l,2)/=0) then
         min_l = min_l+1 !l3 should only lead to parity even numbers
      end if
      max_l = min(lmax,l1+l2)
      !call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
      call GetThreeJs(a3jArr(abs(l2-l1),l2),l1,l2,0,0)

      a3jArr = a3jArr + transpose(a3jArr)
   end do
 end subroutine calcWigners2D

  real function prefactor(l1,l2,l3)
    integer, intent(in) :: l1,l2,l3

    prefactor = 2.0*sqrt((1./4.)*((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.))/pi)
  end function prefactor
  
  subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input
    !parameters.
    implicit none
    integer, parameter :: dl = KIND(1.d0)
    integer, intent(in) :: l2in,l3in, m2in,m3in
    real(dl), dimension(*) :: thrcof
    INTEGER, PARAMETER :: i8 = selected_int_kind(18)
    integer(i8) :: l2,l3,m2,m3
    integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2

    real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1,sumuni
    real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
    integer i,ier, index, nlim, sign2
    integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
    real(dl), parameter :: zero = 0._dl, one = 1._dl
    real(dl), parameter ::  tiny = 1.0d-30, srtiny=1.0d-15, huge = 1.d30,srhuge = 1.d15

    ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

    ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
    !                to l1max = l2+l3
    ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from l1min and l1max
    ! respectively.
    !
    ! lmatch is the l1-value at which forward and backward recursion are
    ! matched.
    !
    ! ndim is the length of the array thrcof
    !
    ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
    ! ier = -2 if possible 3j's exceed ndim
    ! ier >= 0 otherwise

    l2=l2in
    l3=l3in
    m2=m2in
    m3=m3in
    newfac = 0
    lmatch = 0
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    ier = 0

    if (l2 < abs(m2) .or. l3 < m3) then
       ier = -1
       ! call MpiStop('error ier = -1')
       print*, 'error ier = -1',l2,abs(m2),l3,m3
       stop
       return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3),abs(m1))
    l1max = l2+l3

    if (l1min >= l1max) then
       if (l1min/=l1max) then
          ier = -1

          !call MpiStop('error ier = -1')
          print*, 'error ier = -1',l1min,l1max 
          stop
          return
       end if

       ! reached if l1 can take only one value, i.e.l1min=l1max
       thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
       return

    end if

    nfin = l1max-l1min+1

    ! starting forward recursion from l1min taking nstep1 steps
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1 + 1)*tiny

    lstep = 1

30  lstep = lstep+1
    l1 = l1+1

    oldfac = newfac
    a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
    a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
    newfac = sqrt(a2*real(a1,dl))
    if (l1 == 1) then
       !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
       c1 = -(2*l1-1)*l1*(m3-m2)/newfac
    else

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
       denom = (l1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*l1-1)*dv/denom

    end if

    if (lstep<= 2) then

       ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       thrcof(2) = x
       sum1 = sum1+tiny*(2*l1+1)*c1*c1
       if(lstep==nfin) then
          sumuni=sum1
          go to 230
       end if
       goto 30

    end if

    c2 = -l1*oldfac/denom

    ! recursion to the next 3j-coeff x  
    x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
    thrcof(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (2*l1+1)*x*x
    if (lstep/=nfin) then

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(x) >= srhuge) then

          ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
          ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
          ! HAS TO BE RESCALED TO PREVENT OVERFLOW

          ier = ier+1
          do i = 1, lstep
             if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
             thrcof(i) = thrcof(i)/srhuge
          end do

          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge

       end if

       ! as long as abs(c1) is decreasing, the recursion proceeds towards
       ! increasing
       ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
       ! detected, the recursion direction is reversed.

       if (c1old > abs(c1)) goto 30

    end if !lstep/=nfin

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion

    lmatch = l1-1
    x1 = x
    x2 = thrcof(lstep-1)
    x3 = thrcof(lstep-2)
    nstep2 = nfin-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from l1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! l1 = lmatch-1, lmatch, lmatch+1

    nfinp1 = nfin+1
    nfinp2 = nfin+2
    nfinp3 = nfin+3
    l1 = l1max
    thrcof(nfin) = srtiny
    sum2 = tiny*(2*l1+1)

    l1 = l1+2
    lstep=1

    do
       lstep = lstep + 1
       l1= l1-1

       oldfac = newfac
       a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
       a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
       newfac = sqrt(a1*real(a2,dl))

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

       denom = l1*newfac
       c1 = -(2*l1-1)*dv/denom
       if (lstep <= 2) then

          ! if l2=l2max+1, the third term in the recursion vanishes

          y = srtiny*c1
          thrcof(nfin-1) = y
          sumbac = sum2
          sum2 = sum2 + tiny*(2*l1-3)*c1*c1

          cycle

       end if

       c2 = -(l1-1)*oldfac/denom

       ! recursion to the next 3j-coeff y
       y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

       if (lstep==nstep2) exit

       thrcof(nfinp1-lstep) = y
       sumbac = sum2
       sum2 = sum2+(2*l1-3)*y*y

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(y) >= srhuge) then

          ! reached if 3j-coeff larger than srhuge so that the recursion series
          ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent
          ! overflow

          ier=ier+1
          do i = 1, lstep
             index=nfin-i+1
             if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
             thrcof(index) = thrcof(index)/srhuge
          end do

          sum2=sum2/huge
          sumbac=sumbac/huge

       end if

    end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

    y3 = y
    y2 = thrcof(nfinp2-lstep)
    y1 = thrcof(nfinp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal
    ! error

    ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
    nlim = nfin-nstep2+1

    if (abs(ratio) >= 1) then

       thrcof(1:nlim) = ratio*thrcof(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

    else

       nlim = nlim+1
       ratio = 1/ratio
       do n = nlim, nfin
          thrcof(n) = ratio*thrcof(n)
       end do
       sumuni = sumfor + ratio*ratio*sumbac

    end if
    ! normalise 3j-coeffs

230 cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

    sign1 = sign(one,thrcof(nfin))
    sign2 = (-1)**(abs(l2+m2-l3+m3))
    if (sign1*sign2 <= 0) then
       cnorm = -cnorm
    end if
    if (abs(cnorm) >= one) then
       thrcof(1:nfin) = cnorm*thrcof(1:nfin)
       return
    end if

    thresh = tiny/abs(cnorm)

    do n = 1, nfin
       if (abs(thrcof(n)) < thresh) thrcof(n) = zero
       thrcof(n) = cnorm*thrcof(n)
    end do
    return 

  end subroutine GetThreeJs



endprogram bisvar
