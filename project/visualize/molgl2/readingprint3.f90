program LEGGI

real(kind=10) :: x(1680), y(1680), z(1680), phi1(1680), phi2(1680), phi3(1680) 
real(kind=10) :: Lx,Ly,Lz,Volume2
integer :: k, i, st 
character(100) :: stringa, str
real(kind=10) :: xe(3)

        open(unit=30, file="readingconfig3.txt",access="sequential", status="old", action="read")

        ipos = scan(str,":",back=.true.)        

        read (str(1+ipos:),*,iostat=st) volume2

        read(30, *, iostat=st) volume2, Lx, Ly, Lz

        i=1
        
        do while(st.ge.0)

        i = i+1

!       print*, st, i
        read(30,*, iostat=st) x(i-1), y(i-1), z(i-1), phi1(i-1), phi2(i-1), phi3(i-1)
        enddo

        close(30)





print*, ".Vol:", Lx*Ly*Lz, Lx, Ly, Lz 

do k=1,1680


!if ((z(k)<38).and.(z(k)>-25)) then

xe(1) = x(k) -phi1(k)*(8+0.15)
xe(2) = y(k) -phi2(k)*(8+0.15)
xe(3) = z(k) -phi3(k)*(8+0.15)


if (floor(4*(z(k)+Lz/2)/Lz).eq.0) then
print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5	16	C[0.300000,0.500000,0.673810]"
!print*, xe(1), xe(2), xe(3), "@ 	0.25  	C[0.200000,0.200000,0.738100] " 
endif

if (floor(4*(z(k)+Lz/2)/Lz).eq.3) then
print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5  16      C[0.850000,0.500000,0.000000]"
!print*, xe(1), xe(2), xe(3), "@ 	0.25	  C[0.200000,0.200000,0.738100] "
endif

if (floor(4*(z(k)+Lz/2)/Lz).eq.2) then
print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5  16      C[0.800000,0.666666,0.273810]"
!print*, xe(1), xe(2), xe(3), "@ 	0.25	  C[0.200000,0.200000,0.738100] "
endif

if (floor(4*(z(k)+Lz/2)/Lz).eq.1) then
print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5  16      C[0.000000,0.800000,0.273810]"
!print*, xe(1), xe(2), xe(3), "@ 	0.25 	 C[0.200000,0.200000,0.738100] "
endif

!if (floor(5*(z(k)+Lz/2)/Lz).eq.20) then
!print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5  16      C[0.900000,0.100000,0.673810]"
!endif

!endif

enddo

END PROGRAM LEGGI
