
module CW3self
! This module contains the subroutines available for completing Coursework 3
! These subroutines were created by myself.
CONTAINS

SUBROUTINE Shape_fun(coord_l,fun,i,nod)
! Returns the shape functions at the ith integrating point for the local coordinates of the integrating points
 IMPLICIT NONE
 integer, INTENT (IN):: i, nod
 REAL*8,INTENT(IN):: coord_l(:,:)
 REAL*8,INTENT(OUT):: fun(:)
!
If (nod==4) then
!calculate N1-N4 with local coordinates xi and eta of gausspoints
!!!!!evt mit gauspoints(coord_l)
    fun(1)= 1.0/4.0*(1.0-coord_l(i,1))*(1.0-coord_l(i,2))
    fun(2)= 1.0/4.0*(1.0-coord_l(i,1))*(1.0+coord_l(i,2))
    fun(3)= 1.0/4.0*(1.0+coord_l(i,1))*(1.0+coord_l(i,2))
    fun(4)= 1.0/4.0*(1.0+coord_l(i,1))*(1.0-coord_l(i,2))
else
    fun(1)=1.0/4.0*(1.0-coord_l(i,1))*(1.0-coord_l(i,2))*(-coord_l(i,1)-coord_l(i,2)-1.)
    fun(2)=1.0/2.0*(1.0-coord_l(i,1))*(1.0-coord_l(i,2)**2.)
    fun(3)=1.0/4.0*(1.0-coord_l(i,1))*(1.0+coord_l(i,2))*(-coord_l(i,1)+coord_l(i,2)-1.)
    fun(4)=1.0/2.0*(1.0-coord_l(i,1)**2)*(1.0+coord_l(i,2))
    fun(5)=1.0/4.0*(1.0+coord_l(i,1))*(1.0+coord_l(i,2))*(coord_l(i,1)+coord_l(i,2)-1.)
    fun(6)=1.0/2.0*(1.0+coord_l(i,1))*(1.0-coord_l(i,2)**2.)
    fun(7)=1.0/4.0*(1.0+coord_l(i,1))*(1.0-coord_l(i,2))*(coord_l(i,1)-coord_l(i,2)-1.)
    fun(8)=1.0/2.0*(1.0-coord_l(i,1)**2)*(1.0-coord_l(i,2))
end if
RETURN
END SUBROUTINE Shape_fun

subroutine sample(coord_l, points,weights,nod)
Implicit none
integer, INTENT (IN)::  nod
REAL*8,INTENT(out):: coord_l(:,:), points(:,:),weights(:)
integer :: j

if (nod==4) then
    !creating local coordinates of gausspoints
    coord_l(1,:)=(/-(1./SQRT(3.)),-(1./SQRT(3.))/)
    coord_l(2,:)=(/-(1./SQRT(3.)),(1./SQRT(3.))/)
    coord_l(3,:)=(/(1./SQRT(3.)),(1./SQRT(3.))/)
    coord_l(4,:)=(/(1./SQRT(3.)),-(1./SQRT(3.))/)
    !Creating local coordinates of nodes
    points=coord_l/(1./SQRT(3.))
    !defining weights
    weights(:)=1.
else
    !Creating local coordinates of nodes
    points(1,:)=(/-1.,-1./)
    points(2,:)=(/-1.,0./)
    points(3,:)=(/-1.,1./)
    points(4,:)=(/0.,1./)
    points(5,:)=(/1.,1./)
    points(6,:)=(/1.,0./)
    points(7,:)=(/1.,-1./)
    points(8,:)=(/0.,-1./)
    !creating local coordinates of gausspoints
    do j=1,nod
    coord_l(j,:)=(points(j,:)*SQRT(3./5.))
    end do
    coord_l(9,:)=(/0.,0./)
    !defining weights
    weights(1)=25./81.
    weights(2)=40./81.
    weights(3)=25./81.
    weights(4)=40./81.
    weights(5)=25./81.
    weights(6)=40./81.
    weights(7)=25./81.
    weights(8)=40./81.
    weights(9)=64./81.
end if
return
end subroutine sample

subroutine shape_der(coord_l,der,i,nod) !coord_l(i,1),coord_l(i,2),der)
IMPLICIT NONE
integer, INTENT (IN):: i,nod
REAL*8,INTENT(IN):: coord_l(:,:)!coord_l(i,1),coord_l(i,2) !xi(:), eta(:)
REAL*8,INTENT(OUT)::  der(:,:)
Real*8 :: xi,eta
!calculate the derivatives der(localnr,1=xi/2=eta)
if (nod==4) then
    der(1,1)= -1.0/4.0*(1.0-coord_l(i,2))
    der(1,2)= -1./4.*(1.+coord_l(i,2))
    der(1,3)= 1./4.*(1.+coord_l(i,2))
    der(1,4)= 1./4.*(1.-coord_l(i,2))

    der(2,1)= -1./4.*(1.-coord_l(i,1))
    der(2,2)= 1./4.*(1.-coord_l(i,1))
    der(2,3)= 1./4.*(1.+coord_l(i,1))
    der(2,4)= -1./4.*(1.+coord_l(i,1))
else
xi=coord_l(i,1)
eta=coord_l(i,2)
    der(1,1)=1.0/4.0*(1.0-eta)*(2.*xi+eta)
    der(1,2)=-1.0/2.0*(1.-eta**2)
    der(1,3)=-1.0/4.0*(1.+eta)*(-2*xi+eta)
    der(1,4)=-xi*(1.+eta)
    der(1,5)=1./4.*(1.+eta)*(2*xi+eta)
    der(1,6)=1./2.*(1.-eta**2)
    der(1,7)=1./4.*(1.-eta)*(2.*xi-eta)
    der(1,8)=-xi*(1.-eta)

    der(2,1)=1.0/4.0*(1.0-xi)*(2.*eta+xi)
    der(2,2)=-eta*(1.-xi)
    der(2,3)=1.0/4.0*(1.-xi)*(-xi+2.*eta)
    der(2,4)=1.0/2.0*(1.-xi**2)
    der(2,5)=1./4.*(1.+xi)*(xi+2.*eta)
    der(2,6)=-eta*(1.+xi)
    der(2,7)=-1./4.*(1.+xi)*(xi-2.*eta)
    der(2,8)=-1./2.*(1.-xi**2)
end if
return
end subroutine

subroutine beemat(deriv, bee,nod)
Implicit none
integer, intent(in) :: nod
REAL*8,INTENT(in) :: deriv(:,:)
REAL*8,INTENT(out):: bee(:,:)
real*8 :: k=0.
integer :: i

Do i=1,nod
    bee(1,(2*i-1))=deriv(1,i)
    bee(1,(2*i))=k
    bee(2,(2*i-1))=k
    bee(2,(2*i))=deriv(2,i)
    bee(3,(2*i-1))=deriv(2,i)
    bee(3,(2*i))=deriv(1,i)
End Do
!bee(1,:)= (/deriv(1,1), k, deriv(1,2), k,deriv(1,3), k, deriv(1,4), k/)
!bee(2,:)= (/k, deriv(2,1),k,deriv(2,2),k,deriv(2,3),k,deriv(2,4)/)
!bee(3,:)= (/deriv(2,1),deriv(1,1),deriv(2,2),deriv(1,2),deriv(2,3),deriv(1,3),deriv(2,4),deriv(1,4)/)
return
end subroutine beemat

subroutine deemat(nst, prop, dee)
    Implicit none
integer, intent(in) :: nst
REAL*8,INTENT(in) :: prop(:,:)
REAL*8,INTENT(out):: dee(:,:)
real*8 :: k=0.0, l=1.0, E, v

E=prop(1,1)
v=prop(2,1)

dee(1,:)=(/l,(v/(1.-v)),k/)
dee(2,:)=(/(v/(1.-v)),l,k/)
dee(3,:)=(/k,k,((1.-2.*v)/(2.*(1.-v)))/)

dee=dee(:,:)*(E*(1.-v)/((1.+v)*(1.-2*v)))
return
end subroutine deemat

end module CW3self
