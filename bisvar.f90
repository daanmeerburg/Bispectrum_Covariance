!test program to compute components of the bispectrum covariance from lensing 

program bisvar 
  use fwigxjpf
  implicit none
  integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359
  
  character(80) :: Folder1, Folder2, Folder3
  character(80) :: Clfile, Cllfile
  
  real(dl), allocatable :: Cl(:,:), Cll(:,:), Clpp(:,:)

  !integer :: l1, l2, l3
  integer :: lmax, lmin, l1, l2, l3, l1b, l2b, l3b
  integer :: min_l, max_l
  integer :: i,j
  
  real(dl) :: CMB2COBEnorm = 7428350250000.d0
  real(dl) :: DB(4,8), SumDB(4,8), SumTot, DBtot(4,8)
  
  !call fwig_temp_init(2*1000)
  
  lmax = 5000
  lmin = 2
  allocate(Cl(4,2:lmax))
  allocate(Cll(4,2:lmax))
  allocate(Clpp(3,2:lmax))

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
     
     read(17,*) l1, Cl(1:4,j),Clpp(1:3,j)
     read(18,*) l1, Cll(1:4,j)
     Cll(1:4,j) = 2.*pi*Cll(1:4,j)/(real(l1,dl)*(real(l1,dl)+1.))/CMB2COBEnorm
     Cl(1:4,j) = 2.*pi*Cl(1:4,j)/(real(l1,dl)*(real(l1,dl)+1.))/CMB2COBEnorm
     Clpp(1,j) = 2.*pi*Clpp(1,j)/(real(l1,dl)*(real(l1,dl)+1.))**2
     Clpp(2:3,j) = 2.*pi*Clpp(2:3,j)/(real(l1,dl)*(real(l1,dl)+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

     !write(*,*) l1,Cll(1:4,j) 
  enddo

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
  lmax = 100
  call fwig_table_init(2*lmax+2,9)
  !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
  !$OMP PRIVATE(l1,l2,l3, l2b,l3b,min_l,max_l,DB,DBtot), &
  !$OMP REDUCTION(+:SumDB,Sumtot) 
  do l1 = lmin, lmax
     call fwig_thread_temp_init(2*lmax)
     SumDB(1:4,1:8) = 0.d0
     Sumtot = 0.d0
     do l2 =  max(lmin,l1), lmax
        min_l = max(abs(l1-l2),l2)
        !below only relevant if there would be another Wigner3J. 
        !if (mod(l1+l2+min_l,2)/=0) then
        !   min_l = min_l+1 !l3 should only lead to parity even numbers
        !end if
        max_l = min(lmax,l1+l2)
        do l3=min_l,max_l !,2 sum has to be even
           !diagonal 
           l3b=l3
           l2b=l2
           l1b=l1
           !for all these, because of the Wigners inside deltaB (please check):
           !if  (mod(l2+l2b+l3+l3b,2)/=0) then
           !   DB(1:3,1) = 0.d0
           !else    
              call deltaB1(11,12,13,l2b,l3b,Clpp(1,:),5000,DB(1,1))
              call deltaB2(11,12,13,l2b,l3b,Clpp(1,:),5000,DB(2,1))
              call deltaB3(11,12,13,l2b,l3b,Clpp(1,:),5000,DB(3,1))
              call deltaB4(11,12,13,l2b,l3b,Clpp(1,:),5000,DB(4,1))
              !question: what are the permutation allowed??
              !For exmaple, are these the ones:
              call deltaB4(12,11,13,l1b,l3b,Clpp(1,:),5000,DB(4,2))
              call deltaB4(11,13,12,l3b,l2b,Clpp(1,:),5000,DB(4,3))
              call deltaB4(12,13,11,l3b,l1b,Clpp(1,:),5000,DB(4,4))
              call deltaB4(13,12,11,l2b,l1b,Clpp(1,:),5000,DB(4,5))
              call deltaB4(13,11,12,l1b,l2b,Clpp(1,:),5000,DB(4,6))
              !if that is the case, then the others are identical:
              call deltaB1(12,11,13,l1b,l3b,Clpp(1,:),5000,DB(1,2))
              call deltaB1(11,13,12,l3b,l2b,Clpp(1,:),5000,DB(1,3))
              call deltaB1(12,13,11,l3b,l1b,Clpp(1,:),5000,DB(1,4))
              call deltaB1(13,12,11,l2b,l1b,Clpp(1,:),5000,DB(1,5))
              call deltaB1(13,11,12,l1b,l2b,Clpp(1,:),5000,DB(1,6))
              call deltaB2(12,11,13,l1b,l3b,Clpp(1,:),5000,DB(2,2))
              call deltaB2(11,13,12,l3b,l2b,Clpp(1,:),5000,DB(2,3))
              call deltaB2(12,13,11,l3b,l1b,Clpp(1,:),5000,DB(2,4))
              call deltaB2(13,12,11,l2b,l1b,Clpp(1,:),5000,DB(2,5))
              call deltaB2(13,11,12,l1b,l2b,Clpp(1,:),5000,DB(2,6))
              call deltaB3(12,11,13,l1b,l3b,Clpp(1,:),5000,DB(3,2))
              call deltaB3(11,13,12,l3b,l2b,Clpp(1,:),5000,DB(3,3))
              call deltaB3(12,13,11,l3b,l1b,Clpp(1,:),5000,DB(3,4))
              call deltaB3(13,12,11,l2b,l1b,Clpp(1,:),5000,DB(3,5))
              call deltaB3(13,11,12,l1b,l2b,Clpp(1,:),5000,DB(3,6))  
           !endif
           !assuming all are multiplied by Cl1Cl2Cl3 (which is true except for the last term)
           !write(*,*) l1, l2, l3, l2b, l3b, DB1, DB2, DB3
              DBtot(1,1) = (-1.d0)**(l3+l2b)*DB(1,1)
              DBtot(1,2) = (-1.d0)**(l3-l1b)*DB(1,2)
              DBtot(1,3) = (-1.d0)**(l2-l3b)*DB(1,3)
              DBtot(1,4) = (-1.d0)**(l1-l3b)*DB(1,4)
              DBtot(1,5) = (-1.d0)**(l1-l2b)*DB(1,5)
              DBtot(1,6) = (-1.d0)**(l2-l1b)*DB(1,6)
              DBtot(2,1) = (-1.d0)**(-l2-l2b)*DB(2,1)
              DBtot(2,2) = (-1.d0)**(-l1-l1b)*DB(2,1)
              DBtot(2,3) = (-1.d0)**(-l3-l3b)*DB(2,1)
              DBtot(2,4) = (-1.d0)**(-l3-l3b)*DB(2,1)
              DBtot(2,5) = (-1.d0)**(-l2-l2b)*DB(2,1)
              DBtot(2,6) = (-1.d0)**(-l1-l1b)*DB(2,1)
              
              DBtot(3,1) = (-1.d0)**(-l2-l2b)*DB(3,1)*Cl(1,l2b)/Cl(1,l3)
              DBtot(3,2) = (-1.d0)**(-l1-l1b)*DB(3,2)*Cl(1,l1b)/Cl(1,l3)
              DBtot(3,3) = (-1.d0)**(-l3-l3b)*DB(3,3)*Cl(1,l3b)/Cl(1,l2)
              DBtot(3,4) = (-1.d0)**(-l3-l3b)*DB(3,4)*Cl(1,l3b)/Cl(1,l1)
              DBtot(3,5) = (-1.d0)**(-l2-l2b)*DB(3,5)*Cl(1,l1b)/Cl(1,l2)
              DBtot(3,6) = (-1.d0)**(-l1-l1b)*DB(3,6)*Cl(1,l2b)/Cl(1,l1)
              DBtot(4,1) = DB(4,1)*Cl(1,l2b)/Cl(1,l3)
              DBtot(4,2) = DB(4,2)*Cl(1,l1b)/Cl(1,l3)
              DBtot(4,3) = DB(4,3)*Cl(1,l3b)/Cl(1,l2)
              DBtot(4,4) = DB(4,4)*Cl(1,l3b)/Cl(1,l1)
              DBtot(4,5) = DB(4,5)*Cl(1,l1b)/Cl(1,l2)
              DBtot(4,6) = DB(4,6)*Cl(1,l2b)/Cl(1,l1)
              SumDB(1:4,1:6) = SumDB(1:4,1:6) + DBtot(1:4,1:6)
              SumTot = SumTot+ sum(DBtot) !Sum(DBtot)
           !write(*,'(3I4,3E17.8)') l1,l2,l3,SumDB(1:3)
        enddo !l3
     enddo !l2
     write(*,'(I4,25E17.8)') l1, SumTot, SumDB(1:4,1:6)
     call fwig_temp_free();
  enddo !l1
  !$OMP END PARAllEl DO
  
  call fwig_table_free();
  deallocate(Cl, Cll, Clpp)
contains
  !Eq. (29) Notes
  subroutine deltaB1(l1,l2a,l3a,l2b,l3b,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0

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
  end subroutine deltaB1

  !Eq. (31) Notes; unfortunately we have to rewrite these. I dont think we can use the same code, or can we???
  
  subroutine deltaB2(l1,l2a,l3a,l2b,l3b,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0

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
  end subroutine deltaB2

  !Eq. (33) Notes
  subroutine deltaB3(l1,l2a,l3a,l2b,l3b,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax)
    real(dl), intent(out) :: DB
    real(dl) :: stepsum
    integer :: i
    integer :: l_min, l_max
    DB  = 0.d0

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
  end subroutine deltaB3

  !Eq. (33) Notes
  subroutine deltaB4(l1,l2a,l3a,l2b,l3b,CLpp,lmax,DB)
    integer, intent(in) :: l1,l2a,l3a,l2b,l3b
    integer, intent(in) :: lmax
    real(dl), intent(in) :: Clpp(2:lmax)
    real(dl), intent(out) :: DB
    if  (mod(l1+l2a+l3a,2)/=0 .or. mod(l1+l2b+l3b,2)/=0) then
       DB = 0.d0
    else
       DB = Clpp(l1)*Fc(l3a,l1,l2a)*Fc(l3b,l1,l2b)/(2*l1+1)
    endif
    
  end subroutine deltaB4
  
  real(dl) function Fc(l1,l2,l3)
    integer :: l1, l2, l3
    Fc = 1.d0/2.d0*(l2*(l2+1.0)+l3*(l3+1.0)-l1*(l1+1.0))* &
         sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0))/sqrt(4.0*pi)*Wig0(l1,l3,l2)
  end function Fc

  real(dl) function Wig0(l1,l2,l3)
    integer :: l1,l2,l3
    Wig0 = fwig3jj(2* l1 , 2* l2 , 2* l3 , 2* 0, 2* 0 , 2* 0)
    
  end function Wig0
  

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
