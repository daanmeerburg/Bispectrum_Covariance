program FlatSky 
  implicit none
  integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359

  character(80) :: Folder1, Folder2, Folder3
  character(80) :: Clfile, Cllfile

  real(dl), pointer :: Cl(:,:), Cll(:,:)
  real(dl), pointer :: pClpp(:,:)
  integer :: l1, l2a, l3a, l2b, l3b, idphi23, idphi23b
  integer :: lmax, lmin
  integer :: i,j
  real(dl) :: phi23a, phi23b, phi31a, phi31b, phi21a, phi21b, phi2a3b, phi2b3a, phi2a2b, phi3a3b
  real(dl) :: fnl, signsq
  integer :: absl1a, absl1b, absl2l3a, absl2l3b, absl3al3b, absl2al3b, absl2al3b, absl2bl3a
  real(dl) :: l2dotl3a, l2dotl3b, l2adotl3b,l2bdotl3a,l3adotl3b,l2adotl2b
  real(dl) :: testterm


  lmax = 5000

  lmin = 2
  allocate(Cl(4,2:lmax))
  allocate(Cll(4,2:lmax))
  allocate(pClpp(3,2:lmax))

  Folder1 = 'SOspectra/'
  Clfile = trim(Folder1)//trim('SOspectra_lenspotentialCls.dat')
  Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
  !from Alex:
  Cllfile = '../SO_forecasts/CAMB/cosmo2017_10K_acc3_lensedCls.dat'
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
     Cll(1:4,j) = 2.*pi*Cll(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     Cl(1:4,j) = 2.*pi*Cl(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     pClpp(1,j) = 2.*pi*pClpp(1,j)/(l1*(l1+1.))**2
     pClpp(2:3,j) = 2.*pi*pClpp(2:3,j)/(l1*(l1+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

     !write(*,*) l1,Cll(1:4,j) 
  enddo
  close(17)
  close(18)

  !structure
  !l2/l3/l2'/l3'/dphi23/dphi23'

  do l2a, 2, lmax
     do l3a, 2, lmax
        do idphi23, 1, 100
           phi23a = (idphi23-1)*pi/99.
           absl1a = length(l2a,l3a,phi23a)
           fnl = floc(absl1a, l2a,l3a)
           do l2b,2, lmax
              do l3b, 2, lmax

                 do idphi23b, 1, 100
                    phi23b = (idphi23-1)*99
                    absl1b = length(l2b,l3b,phi23b)
                    signsq = fnl*floc(absl1b, l2b,l3b)
                    
                    !other angles
                    phi21a = angle(l2a,l3a,phi23a)
                    phi21b = angle(l2b,l3b,phi23b)

                    phi31a = pi - phi23a - phi21a
                    phi31b = pi - phi23b - phi21b
                    
                    phi2a3b = phi21a + phi31b
                    phi2b3a = phi21b + phi31a

                    phi2a2b = phi21a+phi21b
                    phi3a3b = phi31a+phi31b

                    !inner products

                    l2dotl3a = l2a*l3a*cos(phi23a)
                    l2dotl3b = l2b*l3b*cos(phi23b)

                    l2adotl3b = l2a*l3b*cos(phi2a3b)
                    l2bdotl3a = l2b*l3a*cos(phi2b3a)

                    l3adotl3b = l3a*l3b*cos(phi3a3b)
                    l2adotl2b = l2a*l2b*cos(phi2a2b)

                    !other lengths
                    absl2l3a = length(l2a,l3a,phi23a) !|l2+l3|
                    absl2l3b = length(l2b,l3b,phi23b) !|l2'+l3'|

                    absl2al2b = length(l2a,l2b,phi2a2b) !l2+l2'|
                    absl3al3b = length(l3a,l3b,phi3a3b) !|l3+l3'|
                    
                    

                    
                    testterm = Cl(1,absl1a)*Cl(1,l2a)*Cl(1,l2b)*pClpp(1,)

                    DSNonGauss = 9.d0*sigsq/36.
                    if ((l2A.eq.l2b) .and. (l3a .eq.l3b) .and. (phi23 .eq. phi23b)) then 
                       DSNGauss = sigsq/Cl(1, absl1a)/Cl(1,l2a)/Cl(1,l3a))/6.
                       !<N^2> + delta <N^2>


                    else
                       DSNGauss = 0.d0
                       !<N^2> + delta <N^2>

                    endif


                 enddo !phi23
              enddo !phi23
           enddo !l3b
        enddo !l2b
     enddo !l3a
  enddo !l2a

contains

  real(dl) function floc(l1,l2,l3)
    !SW approximation 
    integer :: l1, l2, l3
    real(dl) :: amp
    real(dl) :: As = 2.1056d-9
    !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
    amp = (2.d0/27./pi**2)*As
    floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
         1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
    floc = floc*amp*2.E-7 !2.E-7  is introduced to get roughly same amplitude at l_max = 500 to full fnl_local
  end function floc
  integer function lenghth(x,y,phi)
    real(dl), intent(in) :: x, y
    real(dl), intent(in) :: phi
    !nearest integer
    length = NINT(sqrt(x**2+y**2-2*x*y*cos(phi)))
  end function lenghth

  real function angle(x,y,phi)
    real(dl), intent(in) :: x, y
    real(dl), intent(in) :: phi

    angle = acos((x - y*cos(phi))/(1.*length(x,y,phi)))
  end function angle




end program FlatSky
