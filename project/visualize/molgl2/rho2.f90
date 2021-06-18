module globalvar 

implicit none
integer, parameter :: dp = selected_real_kind(15,307)
real(kind=dp) :: x(3,3000), phi(3,3000), lungh(3000), xp(3,3000) 

end module globalvar


program rho

use globalvar
implicit none

character(100) :: str
integer, parameter :: dpp = selected_real_kind(15,307) 
real(kind=dpp) :: Matrix(3,3), volume2, eige(3), workspace(9), NEM_ax(3), maxeige, eps, xtemp(3)
integer :: i, k, m, j, Numero, ipos, st, indice, info, ui, l, h, jii, juu 
real(kind=dpp) :: rotmatrix(3,3), phi2(3), ragg, ragg2, deltan, delta2(3), delta3(3), posx(3) 

Numero = 1680
deltan = 0

        open(unit=30, file="readingconfig3.txt",access="sequential", status="old", action="read")

        ipos = scan(str,":",back=.true.)        

        read (str(1+ipos:),*,iostat=st) volume2

        read(30, *, iostat=st) volume2

        i=1
        
        do while(st.ge.0)

        i = i+1

!       print*, st, i
        read(30,*, iostat=st) x(1,i-1), x(2,i-1), x(3,i-1), phi(1,i-1), phi(2,i-1), phi(3,i-1)
        enddo

        close(30)


do k=1,Numero/2
lungh(k) = 16
enddo

do k=Numero/2+1,Numero
lungh(k) = 6.67
enddo


do i=1,3
do k=1,3

Matrix(i,k) = 0

do m=1,Numero
do j=1,Numero
Matrix(i,k) = Matrix(i,k) + phi(i,m)*phi(k,j) 
enddo
enddo

Matrix(i,k) = Matrix(i,k)/Numero

enddo
enddo

call DSYEV("V","U",3,Matrix,3,eige,workspace,9,info)

maxeige = 0
indice = 0

do k=1,3
if (eige(k)>maxeige) then
maxeige = eige(k) 
indice = k
endif
enddo

ragg = 0

do k=1,3
NEM_ax(k) = Matrix(k,indice) 
enddo

print*, NEM_ax(1), NEM_ax(2), NEM_ax(3) 

end program rho
