! Programa funcao espectral de 2 particulas d-wave

program dwave

! ___________Definir variáveis______________________________________

implicit none
integer, dimension (2,2) :: ident
integer, dimension (2,2) :: sigmx
complex, dimension (2,2) :: sigmy
integer, dimension (2,2) :: sigmz
real*8 m2,M,U,J,W2,w,qx1i,qx1f,qy1i,qy1f
real*16 eta,mu2,deltqx,deltqy,deltaa
integer Nqx,Nqy,ii,jj,i,j3
double complex om
real*16 q,qx,qy,j2,j1,qx1,qy1,A11,A22,A12,A21,At
double complex, dimension (2,2) :: Iqd,Iq,temp1,Thet2,Iqs
double complex detI,tempa

! ____________Chamar arquivos externos e associar a variáveis______

open(1,file='inputQxQy.dat',status='old')
open(10,file='spect.out',status='unknown')
open(11,file='spect11.out',status='unknown')
open(12,file='spect12.out',status='unknown')
open(21,file='spect21.out',status='unknown')
open(22,file='spect22.out',status='unknown')


read(1,*)m2,M,U,J,W2,eta,w
read(1,*)qx1i,qx1f,qy1i,qy1f,Nqx,Nqy,deltaa

! _____________Definir constantes___________________________________

ident(1,1)=1
ident(1,2)=0
ident(2,1)=0
ident(2,2)=1

sigmx(1,1)=0
sigmx(1,2)=1
sigmx(2,1)=1
sigmx(2,2)=0

sigmy(1,1)=(0.0,0.0)
sigmy(1,2)=(0.0,-1.0)
sigmy(2,1)=(0.0,1.0)
sigmy(2,2)=(0.0,0.0)

sigmz(1,1)=1
sigmz(1,2)=0
sigmz(2,1)=0
sigmz(2,2)=-1

mu2=((1/m2)+(1/M))**(-1)
om=(w)/W2+((0.0,1.0)*eta/W2)

deltqx=(qx1f-qx1i)/Nqx
deltqy=(qy1f-qy1i)/Nqy;

do ii=1,(Nqx+1)
   do jj=1,(Nqy+1)
	qx1=qx1i+((ii-1)*deltqx)
	qy1=qy1i+((jj-1)*deltqy)
      q=(2*m2*W2*((qx1**(2))+(qy1**(2))))**(0.5)
      qx=(qx1)*((2*m2*W2)**(0.5))
      qy=(qy1)*((2*m2*W2)**(0.5))
      j2=(q**2)/(2*M*W2)
      j1=(1-(mu2/M))*(q**2)/(2*M*W2)
      tempa=(-2*M*M*W2/(m2*(q**2)))*(j2-j1+((j2-om)*(Log(om-j1)-Log(om-j2))))
      do i=1,2
         do j3=1,2
            Iqs(i,j3)=(-mu2/m2)*(Log(om-j1)-Log(om-W2))*ident(i,j3)
      	Iqd(i,j3)=tempa*(((((qx**2)-(qy**2))/(q**2))*(sigmz(i,j3)))+(((2*qx*qy)/(q**2))*(sigmx(i,j3))))
            Iq(i,j3)=Iqs(i,j3)+Iqd(i,j3)
         end do
      end do
      detI=(Iq(1,1)*Iq(2,2))-(Iq(1,2)*Iq(2,1))
      temp1(1,1)=Iq(1,1)+(detI*U)
      temp1(1,2)=Iq(1,2)-(detI*J)
      temp1(2,1)=Iq(2,1)-(detI*J)
      temp1(2,2)=Iq(2,2)+(detI*U)
      do i=1,2
         do j3=1,2
            Thet2(i,j3)=((-1)/(1+(U*((Iq(1,1))+(Iq(2,2))))+(2*J*Iq(1,2))+(((U**2)-(J**2))*detI)))*temp1(i,j3)
         end do
      end do
      A11=Log((deltaa-2*AIMAG(Thet2(1,1)))/deltaa)
      A12=Log((deltaa-2*AIMAG(Thet2(1,2)))/deltaa)
      A21=Log((deltaa-2*AIMAG(Thet2(2,1)))/deltaa)
      A22=Log((deltaa-2*AIMAG(Thet2(2,2)))/deltaa)
      At=Log((deltaa-2*AIMAG(Thet2(1,1)+Thet2(2,2)))/deltaa)
	write(10,*) qx1,qy1,At
	write(11,*) qx1,qy1,A11
	write(12,*) qx1,qy1,A12
	write(21,*) qx1,qy1,A21
	write(22,*) qx1,qy1,A22
   end do
end do

END PROGRAM dwave
