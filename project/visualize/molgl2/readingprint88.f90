program LEGGI

real(kind=10) :: x(3500), y(3500), z(3500), phi1(3500), phi2(3500), phi3(3500) 
real(kind=10) :: Lx,Ly,Lz,Volume2
integer :: k, i, st, Numero
character(100) :: stringa, str
real(kind=10) :: xe(3), qu, ai, eo, qu1, ai1, eo1
integer :: summa, summa2


        open(unit=30, file="readingconfig3.txt",access="sequential", status="old", action="read")

        ipos = scan(str,":",back=.true.)        

!        read (str(1+ipos:),*,iostat=st) Volume2

!!        read(30, *, iostat=st) Volume2, Lx, Ly, Lz

        i=1
        
        do i=1,3360

!        i = i+1

!       print*, st, i
        read(30,*, iostat=st) x(i), y(i), z(i), phi1(i), phi2(i), phi3(i)
        enddo

	read(30,*,iostat=st) Lx, Ly, Lz

        close(30)

Numero = 3360

summa = 0
summa2 = 0

print*, ".Vol:", Lx*Ly*Lz, Lx, Ly, Lz 

do k=1,3360


!if ((z(k)<38).and.(z(k)>-38)) then

!xe(1) = x(k) +phi1(k)*(8+3.15)
!xe(2) = y(k) +phi2(k)*(8+3.15)
!xe(3) = z(k) +phi3(k)*(8+3.15)



!if (k>Numero/2) then
!ragg = phi1(k)*phi1(k-Numero/2) + phi2(k)*phi2(k-Numero/2) + phi3(k)*phi3(k-Numero/2)
!print*, phi3(k), phi3(k-Numero/2)
!else
!ragg = phi1(k)*phi1(k+Numero/2) + phi2(k)*phi2(k+Numero/2) + phi3(k)*phi3(k+Numero/2)
!print*,
!endif

!if (ragg>0) then

!print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@ 1.5  16      C[0.80000,0.800000,0.100000]"

!summa = summa + 1 
!print*, ragg

!else

print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.2  16      C[0.600000,0.200000,0.273810]"

!print*, ragg

!summa2 = summa2 + 1 

!endif


enddo

!print*, summa, summa2

END PROGRAM LEGGI
