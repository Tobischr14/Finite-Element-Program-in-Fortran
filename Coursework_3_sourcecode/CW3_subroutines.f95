module CW3module
!
! This module contains the subroutines available for completing Coursework 3
! Do not modify these subroutines.
!
! To use the subroutines in a program, add this module source code to your
! project and start your program with 'USE CW3module'.
!
CONTAINS
!
FUNCTION determinant(jac) RESULT(det)
!
! This function returns the determinant of a 1x1, 2x2 or 3x3
! Jacobian matrix.
!
 IMPLICIT NONE
 REAL*8,INTENT(IN)::jac(:,:)
 REAL*8::det
 INTEGER::it
 it=UBOUND(jac,1)
 SELECT CASE(it)
 CASE(1)
   det=1.0d0
 CASE(2)
   det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
 CASE(3)
   det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
   det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
   det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
 CASE DEFAULT
   WRITE(*,*)' wrong dimension for Jacobian matrix'
 END SELECT
RETURN
END FUNCTION determinant

SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
 REAL*8,INTENT(IN OUT)::matrix(:,:)
 REAL*8::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
 INTEGER::ndim,i,k
 ndim=UBOUND(matrix,1)
 IF(ndim==2)THEN
   det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
   j11=matrix(1,1)
   matrix(1,1)=matrix(2,2)
   matrix(2,2)=j11
   matrix(1,2)=-matrix(1,2)
   matrix(2,1)=-matrix(2,1)
   matrix=matrix/det
 ELSE IF(ndim==3)THEN
   det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
   det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
   det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
   j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
   j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
   j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
   j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
   j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
   j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
   j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
   j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
   j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
   matrix(1,1)=j11
   matrix(1,2)=j12
   matrix(1,3)=j13
   matrix(2,1)=j21
   matrix(2,2)=j22
   matrix(2,3)=j23
   matrix(3,1)=j31
   matrix(3,2)=j32
   matrix(3,3)=j33
   matrix=matrix/det
 ELSE
   DO k=1,ndim
     con=matrix(k,k)
     matrix(k,k)=1.0d0
     matrix(k,:)=matrix(k,:)/con
     DO i=1,ndim
       IF(i/=k)THEN
         con=matrix(i,k)
         matrix(i,k)=0.0d0
         matrix(i,:)=matrix(i,:)-matrix(k,:)*con
       END IF
     END DO
   END DO
 END IF
RETURN
END SUBROUTINE invert


SUBROUTINE fkdiag(kdiag,g)
!
! This subroutine computes the skyline profile.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::g(:)
 INTEGER,INTENT(OUT)::kdiag(:)
 INTEGER::idof,i,iwp1,j,im,k
 idof=SIZE(g)
 DO i=1,idof
   iwp1=1
   IF(g(i)/=0)THEN
     DO j=1,idof
       IF(g(j)/=0)THEN
         im=g(i)-g(j)+1
         IF(im>iwp1)iwp1=im
       END IF
     END DO
     k=g(i)
     IF(iwp1>kdiag(k))kdiag(k)=iwp1
   END IF
 END DO
RETURN
END SUBROUTINE fkdiag

SUBROUTINE formnf(nf)
!
! This subroutine forms the nf matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::nf(:,:)
 INTEGER::i,j,m
 m=0
 DO j=1,UBOUND(nf,2)
   DO i=1,UBOUND(nf,1)
     IF(nf(i,j)/=0)THEN
       m=m+1
       nf(i,j)=m
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE formnf


SUBROUTINE fsparv(kv,km,g,kdiag)
!
! This subroutine assembles element matrices into a symmetric skyline
! global matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::g(:),kdiag(:)
 REAL*8,INTENT(IN)::km(:,:)
 REAL*8,INTENT(OUT)::kv(:)
 INTEGER::i,idof,k,j,iw,ival
 idof=UBOUND(g,1)
 DO i=1,idof
   k=g(i)
   IF(k/=0)THEN
     DO j=1,idof
       IF(g(j)/=0)THEN
         iw=k-g(j)
         IF(iw>=0)THEN
           ival=kdiag(k)-iw
           kv(ival)=kv(ival)+km(i,j)
         END IF
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE fsparv


SUBROUTINE num_to_g(num,nf,g)
!
! This subroutine finds the g vector from num and nf.
!
IMPLICIT NONE
INTEGER,INTENT(IN)::num(:),nf(:,:)
INTEGER,INTENT(OUT)::g(:)
INTEGER::i,k,nod,nodof
nod=UBOUND(num,1)
nodof=UBOUND(nf,1)
DO i=1,nod
   k=i*nodof
   g(k-nodof+1:k)=nf(:,num(i))
END DO
RETURN
END SUBROUTINE num_to_g


SUBROUTINE spabac(kv,loads,kdiag)
!
! This subroutine performs Cholesky forward and back-substitution
! on a symmetric skyline global matrix.
!
IMPLICIT NONE
REAL*8,INTENT(IN)::kv(:)
REAL*8,INTENT(IN OUT)::loads(0:)
INTEGER,INTENT(IN)::kdiag(:)
INTEGER::n,i,ki,l,m,j,it,k
REAL*8::x
n=UBOUND(kdiag,1)
loads(1)=loads(1)/kv(1)
DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1
   x=loads(i)
   IF(l/=i)THEN
      m=i-1
      DO j=l,m
         x=x-kv(ki+j)*loads(j)
      END DO
   END IF
   loads(i)=x/kv(ki+i)
END DO
DO it=2,n
   i=n+2-it
   ki=kdiag(i)-i
   x=loads(i)/kv(ki+i)
   loads(i)=x
   l=kdiag(i-1)-ki+1
   IF(l/=i)THEN
      m=i-1
      DO k=l,m
         loads(k)=loads(k)-x*kv(ki+k)
      END DO
   END IF
END DO
loads(1)=loads(1)/kv(1)
RETURN
END SUBROUTINE spabac


SUBROUTINE sparin(kv,kdiag)
!
! This subroutine performs Cholesky factorisation on a symmetric
! skyline global matrix.
!
IMPLICIT NONE
REAL*8,INTENT(IN OUT)::kv(:)
INTEGER,INTENT(IN)::kdiag(:)
INTEGER::n,i,ki,l,kj,j,ll,m,k
REAL*8::x
n=UBOUND(kdiag,1)
kv(1)=SQRT(kv(1))
DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1
   DO j=l,i
      x=kv(ki+j)
      kj=kdiag(j)-j
      IF(j/=1)THEN
         ll=kdiag(j-1)-kj+1
         ll=max(l,ll)
         IF(ll/=j)THEN
            m=j-1
            DO k=ll,m
               x=x-kv(ki+k)*kv(kj+k)
            END DO
         END IF
      END IF
      kv(ki+j)=x/kv(kj+j)
   END DO
   kv(ki+i)=SQRT(x)
END DO
RETURN
END SUBROUTINE sparin


SUBROUTINE dismsh(loads,nf,ratmax,g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*_dis.ps" displaying
! the deformed finite element mesh.
!
 IMPLICIT NONE
 REAL*8,INTENT(IN)::g_coord(:,:),loads(0:),ratmax
 INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
 CHARACTER(*),INTENT(IN)::argv
 REAL*8::width,height,scale=72,sxy,xo,yo,x,y,dismag,vmax
 REAL*8::xmin,xmax,ymin,ymax,dmax,zero=0.0d0,pt5=0.5d0,            &
   opt5=1.5d0,fpt5=5.5d0,d8=8.0d0,ept5=8.5d0,d11=11.0d0
 INTEGER::i,ii,j,jj,nn,nel,nod
 OPEN(ips,FILE=argv(1:nlen)//'_dis.ps')
!
 nn=UBOUND(nf,2)
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
 width=xmax-xmin
 height=ymax-ymin
 dmax=ratmax*width
 IF(height>width)dmax=ratmax*height
!
 vmax=zero
 DO i=1,nn
   DO j=1,2
     IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
   END DO
 END DO
 dismag=dmax/vmax
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=1,nn
   IF(g_coord(1,i)+dismag*loads(nf(1,i))<xmin)                            &
     xmin=g_coord(1,i)+dismag*loads(nf(1,i))
   IF(g_coord(1,i)+dismag*loads(nf(1,i))>xmax)                            &
     xmax=g_coord(1,i)+dismag*loads(nf(1,i))
   IF(g_coord(2,i)+dismag*loads(nf(2,i))<ymin)                            &
     ymin=g_coord(2,i)+dismag*loads(nf(2,i))
   IF(g_coord(2,i)+dismag*loads(nf(2,i))>ymax)                            &
     ymax=g_coord(2,i)+dismag*loads(nf(2,i))
!
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
!
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       portrait mode
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
 nel=UBOUND(g_num,2)
 nod=UBOUND(g_num,1)
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
!                       draw the deformed mesh
!
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
 IF(nod==5)nod=4
 IF(nod==9)nod=8
 IF(nod==10)nod=9
 IF(nod==15)nod=12
 DO i=1,nel
   ii=g_num(1,i)
   IF(ii==0)CYCLE
   x=sxy*(g_coord(1,ii)+dismag*loads(nf(1,ii))-xmin)
   y=sxy*(g_coord(2,ii)+dismag*loads(nf(2,ii))-ymin)
   WRITE(ips,'(2f9.2,a)') x, y,' m'
   DO j=2,nod
     jj=g_num(j,i)
     x=sxy*(g_coord(1,jj)+dismag*loads(nf(1,jj))-xmin)
     y=sxy*(g_coord(2,jj)+dismag*loads(nf(2,jj))-ymin)
     WRITE(ips,'(2f9.2,a)') x, y,' l'
   END DO
   WRITE(ips,'(a)')'c s'
 END DO
!
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE dismsh

SUBROUTINE mesh(g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*_msh.ps" displaying
! the undeformed finite element mesh.
!
 IMPLICIT NONE
 REAL*8,INTENT(IN)::g_coord(:,:)
 INTEGER,INTENT(IN)::g_num(:,:),ips,nlen
 CHARACTER(*),INTENT(IN)::argv
 REAL*8::xmin,xmax,ymin,ymax,width,height,scale=72,sxy,xo,yo,x,y,      &
   pt5=0.5d0,opt5=1.5d0,fpt5=5.5d0,d8=8.0d0,ept5=8.5d0,         &
   d11=11.0d0
 INTEGER::i,ii,j,jj,nn,nod,nel
 OPEN(ips,FILE=argv(1:nlen)//'_msh.ps')
!
!                       compute size of mesh
!
 nn=UBOUND(g_coord,2)
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
!
!                       draw the mesh
!
 nod=UBOUND(g_num,1)
 nel=UBOUND(g_num,2)
 IF(nod==5)nod=4
 IF(nod==9)nod=8
 IF(nod==10)nod=9
 IF(nod==15)nod=12
 DO i=1,nel
   ii=g_num(1,i)
   IF(ii==0)CYCLE
   x=sxy*(g_coord(1,ii)-xmin)
   y=sxy*(g_coord(2,ii)-ymin)
   WRITE(ips,'(2f9.2,a)')x,y,' m'
   DO j=2,nod
     jj=g_num(j,i)
     x=sxy*(g_coord(1,jj)-xmin)
     y=sxy*(g_coord(2,jj)-ymin)
     WRITE(ips,'(2f9.2,a)') x, y,' l'
   END DO
   WRITE(ips,'(a)')'c s'
 END DO
!
!                       close output file
!
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE mesh


SUBROUTINE vecmsh(loads,nf,ratmax,cutoff,g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*_vec.ps" displaying
! the nodal displacement vectors.
!
IMPLICIT NONE
 REAL*8,INTENT(IN)::g_coord(:,:),loads(0:),ratmax,cutoff
 INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
 REAL*8::width,height,scale=72,sxy,xo,yo,x1,y1,x2,y2,dismag,           &
   zero=0.0d0,pt5=0.5d0,opt5=1.5d0,fpt5=5.5d0,d8=8.d0,         &
   ept5=8.5d0,d11=11.d0,xmin,xmax,ymin,ymax,dmax,vlen,vmax
 INTEGER::i,j,k,l,nn,nels,nod,ns,i1,i2,j1,j2
 INTEGER,ALLOCATABLE::corner(:,:)
 CHARACTER(*),INTENT(IN)::argv
 LOGICAL::draw
!                       formats
 OPEN(ips,FILE=argv(1:nlen)//'_vec.ps')
!                       open output file and compute scale factors
 nn=UBOUND(nf,2)
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
 width =xmax-xmin
 height=ymax-ymin
 dmax=ratmax*width
 IF(height>width)dmax=ratmax*height
!
 vmax=zero
 DO i=1,nn
   DO j=1,2
     IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
   END DO
 END DO
 dismag=dmax/vmax
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
!
 DO i=1,nn
   IF(g_coord(1,i)+dismag*loads(nf(1,i)) < xmin)                          &
     xmin=g_coord(1,i)+dismag*loads(nf(1,i))
   IF(g_coord(1,i)+dismag*loads(nf(1,i)) > xmax)                          &
     xmax=g_coord(1,i)+dismag*loads(nf(1,i))
   IF(g_coord(2,i)+dismag*loads(nf(2,i)) < ymin)                          &
     ymin=g_coord(2,i)+dismag*loads(nf(2,i))
   IF(g_coord(2,i)+dismag*loads(nf(2,i)) > ymax)                          &
     ymax=g_coord(2,i)+dismag*loads(nf(2,i))
!
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
!
 width=xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       Portrait mode
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'/edef {exch def} bind def'
 WRITE(ips,'(a)')                                                         &
   '/arrow {/@to_y edef /@to_x edef /@from_y edef /@from_x edef'
 WRITE(ips,'(a)')'/@dx @to_x @from_x sub def /@dy @to_y @from_y sub def'
 WRITE(ips,'(a)')'/@length @dx @dx mul @dy @dy mul add sqrt def'
 WRITE(ips,'(a)')'/@angle @dy @dx atan def'
 WRITE(ips,'(a)')'gsave @from_x @from_y translate @angle rotate'
 WRITE(ips,'(a)')                                                         &
   '0 0 moveto @length 0 lineto currentpoint stroke newpath moveto'
 WRITE(ips,'(a)')'-4 -2 rlineto @length 0 moveto'
 WRITE(ips,'(a)')'-4  2 rlineto stroke grestore'
 WRITE(ips,'(a)')'} def'
 WRITE(ips,'(a)')'/*sf {'
 WRITE(ips,'(a)')'exch findfont exch'
 WRITE(ips,'(a)')                                                         &
   'dup type /arraytype eq {makefont}{scalefont} ifelse setfont'
 WRITE(ips,'(a)')'} bind def'
 WRITE(ips,'(a)')'/languagelevel where'
 WRITE(ips,'(a)')'{pop languagelevel} {1} ifelse'
 WRITE(ips,'(a)')'2 lt { % ifelse'
 WRITE(ips,'(a)')'/sf /*sf load def'
 WRITE(ips,'(a)')'} { % else'
 WRITE(ips,'(a)')'/sf /selectfont load def'
 WRITE(ips,'(a)')'} ifelse'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
 WRITE(ips,'(2f9.2,a)')xo,yo, ' translate'
 WRITE(ips,'(f9.2,a)')0.5,' setlinewidth'
!
!                       draw the displacement vectors
!
 vmax=zero
 DO i=1,nn
   vlen=loads(nf(1,i))**2+loads(nf(2,i))**2
   IF(vlen>vmax)vmax=vlen
 END DO
 vmax=SQRT(vmax)*cutoff
 DO i=1,nn
   vlen=SQRT(loads(nf(1,i))**2+loads(nf(2,i))**2)
   x1=sxy*(g_coord(1,i)-xmin)
   y1=sxy*(g_coord(2,i)-ymin)
   x2=sxy*(g_coord(1,i)+dismag*loads(nf(1,i))-xmin)
   y2=sxy*(g_coord(2,i)+dismag*loads(nf(2,i))-ymin)
   IF(vlen>vmax)THEN
     WRITE(ips,'(2f9.2,a,2f9.2,a)') x1, y1,' ', x2, y2, ' arrow'
     WRITE(ips,'(a)') 's'
   END IF
 END DO
!
!                       draw the mesh border
!
 nels=UBOUND(g_num,2)
 nod=UBOUND(g_num,1)
 IF(nod==3.OR.nod==6.OR.nod==10.OR.nod==15)ns=3
 IF(nod==4.OR.nod==5.OR.nod==8.OR.nod==9)ns=4
 ALLOCATE(corner(ns,2))
 IF(nod== 3)corner=RESHAPE((/1,2,3,2,3,1/),(/3,2/))
 IF(nod== 6)corner=RESHAPE((/1,3,5,3,5,1/),(/3,2/))
 IF(nod==10)corner=RESHAPE((/1,4,7,4,7,1/),(/3,2/))
 IF(nod==15)corner=RESHAPE((/1,5,9,5,9,1/),(/3,2/))
 IF(nod== 4)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
 IF(nod== 5)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
 IF(nod== 8)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 IF(nod== 9)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 DO i=1,nels
   DO j=1,ns
     draw=.TRUE.
     i1=g_num(corner(j,1),i)
     i2=g_num(corner(j,2),i)
     DO k=1,nels
       DO l=1,ns
         j1=g_num(corner(l,1),k)
         j2=g_num(corner(l,2),k)
         IF((i1==j2).AND.(i2==j1))THEN
           draw=.FALSE.
           EXIT
         END IF
       END DO
       IF(.NOT.draw)EXIT
     END DO
     IF(draw)THEN
       x1=sxy*(g_coord(1,i1)-xmin)
       y1=sxy*(g_coord(2,i1)-ymin)
       WRITE(ips,'(2f9.2,a)')x1, y1,' m'
       x1=sxy*(g_coord(1,i2)-xmin)
       y1=sxy*(g_coord(2,i2)-ymin)
       WRITE(ips,'(2f9.2,a)')x1, y1,' l'
       WRITE(ips,'(a)')' s'
     END IF
   END DO
 END DO
!                       close output file?
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
RETURN
END SUBROUTINE vecmsh

SUBROUTINE getname(argv,nlen)
!
! This subroutine the base name of data file.
!
 IMPLICIT NONE
 INTEGER::narg
 INTEGER,INTENT(OUT)::nlen
 INTEGER::lnblnk,iargc
 CHARACTER(*),INTENT(OUT)::argv
 LOGICAL found
 narg=iargc()
 IF(narg.lt.1)THEN
   WRITE(*,*)'Please enter the base name of data file: '
   READ(*,*) argv
  ELSE
   CALL getarg(1,argv)
 ENDIF
 nlen=lnblnk(argv)
 INQUIRE(file=argv(1:nlen)//'.dat',exist=found)
 IF(.not.found)THEN
  WRITE(*,*)'Data file not found: ',argv(1:nlen)//'.dat'
  WRITE(*,*)'Please create or check spelling.'
  STOP
 ENDIF
RETURN
END SUBROUTINE getname


end module CW3module
