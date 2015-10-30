!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Shooting method to solve boundary value problems 
!     Uses Runke-Kutta 4 as ODE solver
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 29, 2015
!-----------------------------------------------------------------------------!


program shooting
implicit none
real*8 ::h,a,b,x,y2a,y3a,y2b,A1,A2,B1,B2,guess
integer::i,j,n,ne
real*8,allocatable ::y(:)

!Solve y''' + y*y' = 0, for y(x), x in [0,10]


!Array definitions:
!y(1) = y''
!y(2) = y'
!y(3) = y

ne = 3 !number of equations
allocate(y(ne))

n=1000
a=0.0d0
b=10.0d0
h=(b-a)/dfloat(n)

!Boundary conditions
y2a = 0.0d0
y3a = 0.0d0
y2b = 1.0d0

!initial two guesses for y(1) at x=0
A1 = 1.0d0
A2 = 0.5d0

!First guess:    
!Initial condition
	y(1) = A1 !guess
	y(2) = y2a !bc
    y(3) = y3a !bc

    !Time integration with RK4 scheme to find B
	do j=1,n
		x = dfloat(j)*h
		call rk4(ne,x,h,y)  	  
	end do
    !assign estimated value for y(2) at x=10
    B1 = y(2)

!Second guess:    
!Initial condition
	y(1) = A2 !guess
	y(2) = y2a !bc
    y(3) = y3a !bc

    !Time integration with RK4 scheme
	do j=1,n
		x = dfloat(j)*h
		call rk4(ne,x,h,y)	   
	end do
    !assign estimated value for y(2) at x=10
    B2 = y(2)    


!iteration for the rest, it will converge quickly
do i=1,100

	!Initial condition
    guess = A2 + (y2b - B2)/((B2-B1)/(A2-A1)) !secant method 
    
	y(1) = guess !guess
	y(2) = y2a !bc
    y(3) = y3a !bc
 
	!Time integration with RK4 scheme
	do j=1,n
		x = dfloat(j)*h
		call rk4(ne,x,h,y)	  
	end do
    B1 = B2   	!update for the next iteration 
    B2 = y(2)  	!update for the next iteration 
    A1 = A2 	!update for the next iteration 
    A2 = guess 	!update for the next iteration 
    
    print*,i,B2
    
    !check for the final point
    if(dabs(B2-y2b).le.1.0d-6) goto 100
      
end do

100 continue

!Make final computation one more time to plot results 
open(12, file="blasius.plt")
write(12,*)'variables ="x","f1","f2","f3"'
	!Initial condition 
	y(1) = guess !guess
	y(2) = y2a !bc
    y(3) = y3a !bc
 
	write(12,*) 0.0d0,y(1),y(2),y(3)
	!Time integration with RK4 scheme
	do j=1,n
		x = dfloat(j)*h
		call rk4(ne,x,h,y)
    	write(12,*) x,y(1),y(2),y(3)   
	end do
close(12)   

    


end


!-----------------------------------------------------------------------------!
!RK4 scheme
!-----------------------------------------------------------------------------!
subroutine rk4(ne,x,h,y)
implicit none
integer::ne
integer::k
real*8 ::x,h
real*8 ::y(ne),r(ne)
real*8 ::k1(ne),k2(ne),k3(ne),k4(ne)

    call RHS(ne,y,x,r)
    
		do k=1,ne
		k1(k) = h*r(k)      
    	end do
	
    call RHS(ne,y+k1/2.0d0,x+h/2.0d0, r)
		do k=1,ne
		k2(k) = h*r(k)      
    	end do

    call RHS(ne,y+k2/2.0d0,x+h/2.0d0,r)
		do k=1,ne
		k3(k) = h*r(k)      
    	end do

    call RHS(ne,y+k3,x+h,r)
		do k=1,ne
		k4(k) = h*r(k)      
    	end do
    
 	do k=1,ne
		y(k) = y(k) + (k1(k)+2.0d0*(k2(k)+k3(k))+k4(k))/6.0d0     
    end do

end 




!-----------------------------------------------------------------------------!
!Right Hand Side (RHS) : Blasius problem
!-----------------------------------------------------------------------------!
subroutine RHS(ne,y,x,r)
implicit none
integer::ne
real*8 ::x
real*8 ::y(ne),r(ne)

r(1) =-y(1)*y(3)
r(2) = y(1)
r(3) = y(2)

end 








