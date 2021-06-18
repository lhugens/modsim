	PROGRAM changing


	implicit none

        integer, parameter :: dppp = selected_real_kind(15, 307)
        real(kind=dppp) x(3364), y(3364), z(3364)
        real(kind=dppp) phi1(3364), phi2(3364), phi3(3364)
	real(kind=dppp) volume2, Lx, Ly, Lz
        character(100) :: stringa, str
	integer k, i, st, ipos



open(unit=30, file="readingconfig3.txt",access="sequential", status="old", action="read")

ipos = scan(str,":",back=.true.)        

read (str(1+ipos:),*,iostat=st) volume2

read(30, *, iostat=st) volume2, Lx, Ly, Lz

i=1
        
do while(st.ge.0)

i = i+1

read(30,*, iostat=st) x(i-1), y(i-1), z(i-1), phi1(i-1), phi2(i-1), phi3(i-1)
        
enddo

 close(30)



print*, ".Vol:", volume2, Lx, Ly, Lz


do k=1,1680

print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@	1.5	16	C[0.000000,0.800000,0.673810]"

print*, x(k)+phi1(k)*11.15, y(k)+phi2(k)*11.15, z(k)+phi3(k)*11.15, "@ 3.15  C[0.300000,0.800000,0.373810]"

enddo


do k=1681,3360

print*, x(k), y(k), z(k), phi1(k), phi2(k), phi3(k), "@ 1.5     16      C[0.300000,0.300000,0.973810]"

print*, x(k)+phi1(k)*11.15, y(k)+phi2(k)*11.15, z(k)+phi3(k)*11.15, "@ 3.15  C[0.300000,0.800000,0.373810]"

enddo




END PROGRAM changing
