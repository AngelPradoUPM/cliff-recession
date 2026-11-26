!----------------------------------------------------------------------*
!                                                                      *
!     Equilibrium erosion of rock shores                               *
!																	   *
!                                                                      *
!     Program name - version: cliff_erosion_4.1                        *
!                                                                      *
!     Objetive: erosion of rock cliffs		                           *     
!	  																   *
!              Erosion                                                 *
!              Two directions (I<-D ; I->D)                            *
!              All forms of cliffs                                     *
!			   All angle failure                                       *
!			   Consider 3 kinds of failure:    	                       *
!                       * Rotational                                   *
!                       * Falling blocks (+ Joss Bay)                  *
!			   Talus debris in front of the cliff                      *
!              Differential rock strength                              *
!			   Sea level rise due to climate change                    *
!              Tectonic movements                                      *
!			   Uncertain angle of failure                              *
!			   Solve the tan with 2 points                             *
!			   Solve the tan with 3 points                             *
!			   Solve the integral:                                     *
!						* Left Rectangle                               *
!						* Right Rectangle                              *
!						* Mid point rule                               *
!						* Trapezoid rule                               *
!						* Simpson method                               *
!						* Gaussian quadrature                          *
!			   Stochastic time of failure:                             *
!						* Determinista                                 *
!						* Gamma                                        *
!						* Uniforme                                     *
!						* Pareto                                       *
!						* Pareto Bounded                               *
!                                                                      *
!     Input  file: datos.dat (initial data)                            *
!				   marea.dat (datos de marea)                          *
!     Output file: resultados.out                                      *
!				   Matriz_dibujos.out                                  *
!				   Matriz_roturas.out                                  *
!                                                                      *
!                                                                      *
!     Programers: R. Castedo y C. Paredes                              *
!     Email: ricardo.castedo@upm.es / carlos.paredes@upm.es            *
!                                                                      *
!     First version date: 20 / October / 2009		                   *
!     Last version date: 21 / April / 2011		                       *
!                                                                      *
!
!	  Descripcion de variables:
!
!     AD=desplazamiento medio
!	  alpha=pendiente con la que rompe el talud
!	  Am=Amplitud de marea
!	  Db= profundidad de agua en las rompientes
!	  Dbmax=profundidad maxima de ruptura
!	  Dy=lo que erosiono en eje horizontal      
!	  Dym=variable auxiliar de diferencial de y
!	  Dymt=matriz de Diferenciales de y en una marea y con un tipo de ola
!	  Dtm=Calculo del diferencial de tiempo de marea
!	  Dz=tamano del diferencial de z
!	  F=Fuerza erosiva
!	  Factor=Factor que da orientacion del acantilado
!	  H0=Altura de la ola mar adentro
!	  Hb=Altura de la ola rompiente
!	  HmaxM=Altura Maxima de Marea
!	  HminM=Altura Minima de Marea
!	  icortemin=corte donde empiezo a refinar el acantilado
!	  icortemax=corte donde termino de refinar el acantilado
!	  indiceNM=indice que determina la posicion del nivel del mar en Z
!	  indvertrot=indice de verticalidad de la rotura
!	  indvertcliff=indice de verticalidad del acantilado
!	  Iptorotura=indice que indica el punto desde el que se produce la rotura del acantilado
!	  Kh=constante hidrodinamica
!
!	  MagM=magnitud momento del terremoto
!     Mio=matriz de que almacena datos para cada incide de ola
!	  NDz=Numero de incrementos en z
!	  NDzmax=Numero maximo de incrementos de z
!	  NG=Numero de puntos que definen el perfil del acantilado original
!	  NO=numero de olas
!	  NTsim=numero de tiempos de simulacion
!	  NTmarea=numero de diferenciales para el calculo de la marea
!     Otec=sentido cambio terreno por tectonica 1:hacia arriba
!	  Pn=aspecto de la funcion de erosion
!     PR=periodo de retorno del terremoto de magnitud MagM
!	  Pzy=Creacion de la matriz de valores
!	  R=constante de Skinner
!	  SLR=incremento del nivel del mar con el tiempo
!	  sigmac=resistencia a compresión simple
!	  Tdib=tiempo en el que queremos un dibujo
!	  Trotura=tiempo en el que se produce una rotura
!	  tmarea=duracion periodo de marea (h)
!	  tsealc=tasa (m) de cambio de marea cada 100 años
!     tsim=tiempo maximo de simulacion
!	  TT=periodo de marea (inversa de la frecuencia)
!	  Utsim=unidades en que se introduce el tiempo de simulacion
!	  wt=funcion de marea (p.ej: funcion senoidal)
!	  yerosionmax=vector aux. para puntos de erosion maxima
!	  Yc=coordenada Y punto de corte rotura-acantilado
!	  Zc=coordenada Z punto de corte rotura-acantilado
!
!                                                                      *
!----------------------------------------------------------------------*
 
implicit none
parameter m=100,mm=4,mmm=1000
character*500 basura 
integer aa,aaa,controlperfil,crestatalud,DimDiby,DimPzyR,dimmatrizrotura,dimmatrizHwdraw,dimmatrizrecesion,dimrotacumulada,groundwater,i,ibusqueda,icortemin,icortemax
integer Imaxerosion[allocatable](:),l,Ialpha,Iincertidumbre
integer IDymt,IDy,IDibMatriz,Iyerosionmax,IMatrizrotura,IMatrizrecesion,Indmaxerosion,IMatrizHwdraw,Iptorotura,indiceNM,indvertrot[allocatable](:),indvertcliff,irefinomax,irefinomin
integer Ialpharot,Imrotura,Iindvertrot,IEnormalstress,IMinorppalstress,IEcohesion,IEcohesionwet,IEphi,IEphiwet,Imi,Imacizo,IRMR,Isigmacmass
integer indmaxASL1(1),io,iptoinf,IPzyRglobal,IR,Irotacumulada,IRcuna,IKh,Ipotestrato,Ipesounit,Isigmac,Isigmacwet,IPgauss,IWgauss,Iwtg,it,label(10),ilabel,iseed,Isigmacrotwet,nlabel
integer j,jdibujo,jj,k,kkentero,ll,macizo[allocatable](:),modelo,n,ndatosresis,NG,Nmaxerosion,NO,NPrefino,NPgauss,Ntmarea,NTsim,nrot
integer Ptoerosion1,Ptoerosion2,Ptoerosion3,parametrorotura,Ptocorte,PtocorteCirc2,p,PR,ptocorteaux,parametrocuna,Ptocortebase,PtocortecircMedio,Ptoultrotura
integer rotcunanotch,roturanotch,seepagelevel,selectrot,tdatos,tcambiodz,tcambiotec,tdatosrotura,tipocuna,tipotg,tipointegral,tipofractura,tiporotura,tglob,tp,Utsim,UTdib,UTrotura,Vertcirc
integer variablecontrolcuna,w
real(8), parameter :: pi= 3.141592653589793
real*8 a,AD,AD1,anchoslide,alpha[allocatable](:),alpharot[allocatable](:),Am,aux,aux1,aux2,aux3,b,CXc,CYc,curvah,CoordX,controlmin
real*8 Db(mm),Dbmax,DENO,DistPcirc2,Dtm,dummy[allocatable](:,:),dummy1[allocatable](:,:),Dy[allocatable](:),Dym,Dymt[allocatable](:,:),Dz,DibMatriz[allocatable](:,:),distY,dist,distancia(mmm),dummy2[allocatable](:,:)
real*8 Ecohesion[allocatable](:),Ecohesionwet[allocatable](:),Ecohesionrot,Enormalstress[allocatable](:),Ephi[allocatable](:),Ephiwet[allocatable](:),Ephirot,Ephirotwet,Factor,F(mm),FS,ga,gamma,grados,h,H0(mm),Hb(mm),HmaxM,HminM(mm),Hslide,Htalud,Hw,Hwdraw(mmm),hrot,incertidumbre[allocatable](:),Kh[allocatable](:),Khdef,kk,lambda,Longitud
real*8 Lcdg,malcance,mcirculo,mcorte,Mmedio,munion,mtalud,mtaludcorte,mmat,Mayorppalstress,MagM,mcliff,minestabilidad,mi[allocatable](:),Minorppalstress[allocatable](:),mrotura[allocatable](:),mrecta,Matrizrotura[allocatable](:,:),Matrizrecesion[allocatable](:,:),MatrizHwdraw[allocatable](:,:)
real*8 mcliffcuna,mroturacuna
real*8 Mio(mm,7),nslides,NUM,Otec,pesounit[allocatable](:),pesounitrot,pesounitagua,Pn,Pzy(3,m),profagua,profseepage,PzyRglobal[allocatable](:,:),potestrato[allocatable](:),potacumulada,Pgauss[allocatable](:),ptoestrato1,ptoestrato2
real*8 ptoinf,PtocortecircMedioy,promedio1,promedio2,promedio3,promedio4,promedio5,promedio6
real*8 pendientetgPcirc1,perdidacuna,retroceso,R[allocatable](:),Rcuna[allocatable](:),radio,RMR[allocatable](:),rotacumulada[allocatable](:)
real*8 s,SLR,shear,shape,shapep,shapepb,scale,sigmac[allocatable](:),sigmacwet[allocatable](:),sigmacmass[allocatable](:),sigmacrotwet[allocatable](:),sigmacwetdef
real*8 Sigmatensilerock,Sigmatensilemax,sumaincrementoy,sumaacumuladacuna,sumaacumuladarot
real*8 t,theta,tcambiosl,tcambiotec1,Tdib,tsealc,Trotura,Trotacumulada,tgrados,tmarea,tsim,TT(mm),TanF1,TanF2
real*8 vectorcuna[allocatable](:),vectorcunaaux[allocatable](:),vectorposoriginal[allocatable](:),vectorverticalrot[allocatable](:)
real*8 volacumulado,volacumuladopc,volaux,V1,V2,u1,u2,Umin,Umax,Wgauss[allocatable](:),wta,wtb,wtab2,wtg[allocatable](:),Xca,Xcb,Xmedio,Xptomedio,Xm,XDistPcirc2,YDistPcirc2,yerosionmax[allocatable](:),Yc
real*8 Yca,Ycb,Ymedio,Yptomedio,Yopc1,Yopc2,Zc
real   tic,tac,rotale

!Creacion de los ficheros de lectura y salida de datos
open(1,file='datos.dat',status='unknown')
open(2,file='timerecesion.out',status='unknown')
open(3,file='marea.dat',status='unknown')
open(4,file='Matriz_dibujo.out',status='unknown')
open(5,file='tiempos_rotura.out',status='unknown')
open(6,file='Matriz_roturas.out',status='unknown')
open(7,file='Matriz_recesion.out',status='unknown')
open(8,file='Matriz_seepage.out',status='unknown')

!Lectura de los datos del programa
read(1,*) basura
read(1,*) basura
read(1,*) modelo
read(1,*) basura
read(1,*) Utsim
read(1,*) basura
read(1,*) tsim
read(1,*) basura
read(1,*) UTdib
read(1,*) basura
read(1,*) Tdib
read(1,*) basura
read(1,*) tiporotura
read(1,*) basura
read(1,*) UTrotura
read(1,*) basura
read(1,*) Trotura
read(1,*) basura
read(1,*) iseed
read(1,*) basura
read(1,*) shape,scale
read(1,*) basura
read(1,*) shapep,Xm
read(1,*) basura
read(1,*) Umin,Umax,shapepb
read(1,*) basura
read(1,*) NO
read(1,*) basura
read(1,*) Am
read(1,*) basura
read(1,*) tmarea
read(1,*) basura
do i = 1, NO
	read(1,*) H0(i)
end do
read(1,*) basura
do i = 1, NO
	read(1,*) TT(i)
end do
read(1,*) basura
read(1,*) Ntmarea
read(1,*) basura
read(1,*) Dz
read(1,*) basura
read(1,*) grados
read(1,*) basura
read(1,*) tipotg
read(1,*) basura
read(1,*) tipointegral
read(1,*) basura
read(1,*) NPgauss
read(1,*) basura
read(1,*) Nmaxerosion
read(1,*) basura
read(1,*) tsealc
read(1,*) basura
read(1,*) MagM,PR,Otec
read(1,*) basura
read(1,*) tipofractura,rotcunanotch
read(1,*) basura
read(1,*) groundwater,profagua
read(1,*) basura !Introduccion del numero de datos diferentes de R y sigma
read(1,*) ndatosresis

allocate (alpha(ndatosresis),stat=Ialpha)
if (Ialpha.NE.0) then
	stop 'Error de dimensionamiento de angulo de rotura'
end if
allocate (incertidumbre(ndatosresis),stat=Iincertidumbre)
if (Iincertidumbre.NE.0) then
	stop 'Error de dimensionamiento de variabilidad en la rotura'
end if
read(1,*) basura !Introduccion de los valores de alpha e incertidumbre
do i = 1, ndatosresis
	read(1,*) alpha(i),incertidumbre(i)
end do

allocate (Kh(ndatosresis),stat=IKh)
if (IKh.NE.0) then
	stop 'Error de dimensionamiento de Kh'
end if
allocate (sigmac(ndatosresis),stat=Isigmac)
if (Isigmac.NE.0) then
	stop 'Error de dimensionamiento de sigmac'
end if
allocate (Enormalstress(ndatosresis),stat=IEnormalstress)
if (IEnormalstress.NE.0) then
	stop 'Error de dimensionamiento de Enormalstress'
end if
allocate (Minorppalstress(ndatosresis),stat=IMinorppalstress)
if (IMinorppalstress.NE.0) then
	stop 'Error de dimensionamiento de Minorppalstress'
end if
allocate (Ecohesion(ndatosresis),stat=IEcohesion)
if (IEcohesion.NE.0) then
	stop 'Error de dimensionamiento de Ecohesion'
end if
allocate (Ephi(ndatosresis),stat=IEphi)
if (IEphi.NE.0) then
	stop 'Error de dimensionamiento de Ephi'
end if
allocate (pesounit(ndatosresis),stat=Ipesounit)
if (Ipesounit.NE.0) then
	stop 'Error de dimensionamiento de pesounit'
end if
allocate (mi(ndatosresis),stat=Imi)
if (Imi.NE.0) then
	stop 'Error de dimensionamiento de constante m del material'
end if
allocate (macizo(ndatosresis),stat=Imacizo)
if (Imacizo.NE.0) then
	stop 'Error de dimensionamiento de macizo'
end if
allocate (RMR(ndatosresis),stat=IRMR)
if (IRMR.NE.0) then
	stop 'Error de dimensionamiento de RMR'
end if
allocate (potestrato(ndatosresis),stat=Ipotestrato)
if (Ipotestrato.NE.0) then
	stop 'Error de dimensionamiento de potestrato'
end if
read(1,*) basura !Introduccion de los valores de Kh,sigma y potencia asignada a esos valores
read(1,*) basura
do i = 1, ndatosresis
	read(1,*) Kh(i),sigmac(i),Enormalstress(i),Minorppalstress(i),Ecohesion(i),Ephi(i),pesounit(i),mi(i),macizo(i),RMR(i),potestrato(i)
end do

allocate (sigmacwet(ndatosresis),stat=Isigmacwet)
if (Isigmacwet.NE.0) then
	stop 'Error de dimensionamiento de sigmacwet'
end if
allocate (Ecohesionwet(ndatosresis),stat=IEcohesionwet)
if (IEcohesionwet.NE.0) then
	stop 'Error de dimensionamiento de Ecohesionwet'
end if
allocate (Ephiwet(ndatosresis),stat=IEphiwet)
if (IEphiwet.NE.0) then
	stop 'Error de dimensionamiento de Ephiwet'
end if

read(1,*) basura
read(1,*) tipocuna,perdidacuna

read(1,*) basura
do i = 1, ndatosresis
	read(1,*) sigmacwet(i),Ecohesionwet(i),Ephiwet(i)
end do

write(*,*)'The rock cliff erosion'
!write(*,100)

!Lectura de la matriz de valores Z e Y
read(1,*) basura
k=0
nlabel=0
10 continue
k=k+1
read(1,*,end=20) Pzy(1,k),Pzy(2,k),Pzy(3,k)
if(Pzy(3,k).eq.1.) nlabel=nlabel+1
go to 10
20 continue
NG=k-1

close(1)

select case(Utsim)
	case(0)
		!Introduccion en horas (h)
		tsim=tsim
	case(1)
		!Introduccion en dias (dias*horas(24)=horas)
		tsim=tsim*24
	case(2)
		!Introduccion en meses (mes*dias(30)*horas(24)=horas)
		tsim=tsim*30*24
	case(3)
		!Introduccion en años (año*dias(365)*horas(24)=horas)
		tsim=tsim*365*24
end select

select case(UTdib)
	case(0)
		!Introduccion en horas (h)
		Tdib=Tdib
	case(1)
		!Introduccion en dias (dias*horas(24)=horas)
		Tdib=Tdib*24
	case(2)
		!Introduccion en meses (mes*dias(30)*horas(24)=horas)
		Tdib=Tdib*30*24
	case(3)
		!Introduccion en años (año*dias(365)*horas(24)=horas)
		Tdib=Tdib*365*24
end select

select case(UTrotura)
	case(0)
		!Introduccion en horas (h)
		Trotura=Trotura
	case(1)
		!Introduccion en dias (dias*horas(24)=horas)
		Trotura=Trotura*24
	case(2)
		!Introduccion en meses (mes*dias(30)*horas(24)=horas)
		Trotura=Trotura*30*24
	case(3)
		!Introduccion en años (año*dias(365)*horas(24)=horas)
		Trotura=Trotura*365*24
end select

!Calculo de los numeros de tiempos de simulacion
NTsim = int(tsim/tmarea)

!Obtencion del factor de orientacion
if (Pzy(2,1).LT.Pzy(2,2)) then
	Factor = 1.0
else if (Pzy(2,1).GT.Pzy(2,2))then
	Factor = -1.0
end if

!Calculo de la tangente limitadora de los grados
tgrados=tan((grados*pi)/180.)

!Calculo del tiempo que tarda en cambiar el mar un dz (en años)
if (tsealc.EQ.0.) then
	tcambiodz = 0.
else
	tcambiosl = (Dz*100./tsealc) * (365*24./tmarea)
	tcambiodz = int(tcambiosl)
end if

!Calculo del tiempo que tarda en producirse un evento tectonico que produzca desplazamientos
!PR = periodo de retorno
tcambiotec1 = (PR*365*24./tmarea)
tcambiotec = int(tcambiotec1)

!Calculo del SLR (cambia un dz cada tcambiodz)
SLR = 0.0
 
!Calculo de la Altura Maxima de Marea (HmaxM)
!si uso marea sinusoidal
HmaxM = Am 
!si uso marea de datos
!HmaxM = max(w(t))
do io=1,NO
	!Calculo de Hb en funcion de io
		Hb(io) = 0.39*(9.81**0.2)*((TT(io)*H0(io)**2)**0.4)

	!Calculo de Db
		Db(io) = Hb(io)/0.78

	!Calculo de F en funcion de io
		F(io) = (Hb(io)**(13/4))*(TT(io)**(3/2))
			
	!Calculo de la Altura Minima de Marea (HminM) en funcion de io
		!si uso marea sinusoidal
		HminM(io) = -Am - Db(io)
		!si uso marea de datos
		!HminM(io) = min(w(t))-Db(io)
end do

do io= 1,NO
	Mio(io,1) = H0(io)
	Mio(io,2) = TT(io)
	Mio(io,3) = Hb(io)
	Mio(io,4) = Db(io)
	Mio(io,5) = F(io)
	Mio(io,6) = HminM(io)
	Mio(io,7) = HmaxM
end do

!Maximo valor del Db para todos los tipos de olas
Dbmax = Mio(1,4)
do io=2,NO
	Dbmax = max(Dbmax,Mio(io,4))
end do

NPrefino = int(abs(Pzy(1,1)-Pzy(1,NG))/Dz)

!Creacion de la matriz PzyR con Allocate
allocate (PzyRglobal(2,NPrefino),stat=IPzyRglobal)
if (IPzyRglobal.NE.0) then
	stop 'Error de dimensionamiento de la matriz PzyRglobal'
end if

!Coordenadas del primer punto del vector de refino
PzyRglobal(1,1) = Pzy(1,1)
PzyRglobal(2,1) = Pzy(2,1)
!Coordenadas del ultimo punto del vector de refino
PzyRglobal(1,NPrefino) = Pzy(1,NG)
PzyRglobal(2,NPrefino) = Pzy(2,NG)

!Calculo de las coordenadas del vector refinado PzyRglobal para el tiempo 0
do k = 2,NPrefino-1
	PzyRglobal(1,k) = PzyRglobal(1,k-1) - Dz
	do ibusqueda = 2,NG
		if (Pzy(1,ibusqueda).LT.PzyRglobal(1,k)) then
			if (Pzy(2,ibusqueda).EQ.Pzy(2,ibusqueda-1)) then
				PzyRglobal(2,k) = Pzy(2,ibusqueda)
			else
				PzyRglobal(2,k) = (PzyRglobal(1,k)-Pzy(1,ibusqueda-1))/(Pzy(1,ibusqueda)-Pzy(1,ibusqueda-1))&
				*(Pzy(2,ibusqueda)-Pzy(2,ibusqueda-1))+Pzy(2,ibusqueda-1)
			end if
			go to 501
		end if
	end do
	501 continue		
end do


ilabel=0
do ibusqueda=2,NG
    !Localizacion de los puntos de seguimiento
    if(Pzy(3,ibusqueda).EQ.1.)then
       do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     ilabel=ilabel+1
			 label(ilabel)=k
			 goto 100
		  end if
	   end do
    !Localizacion de la cresta del talud
	else if(Pzy(3,ibusqueda).EQ.2.)then
	   do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     crestatalud=k
			 goto 100
		  end if
	   end do
    !Localizacion del punto de surgencia o ssepage
	else if(Pzy(3,ibusqueda).EQ.3.)then
	   do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     seepagelevel=k
			 goto 100
		  end if
	   end do
	end if
	100 continue
end do

allocate (sigmacmass(ndatosresis),stat=Isigmacmass)
if (Isigmacmass.NE.0) then
	stop 'Error de dimensionamiento de sigmacmass'
end if

!Calculo de la resistencia a compresion simple del macizo rocoso
do i=1,ndatosresis
    sigmacmass(i) = 0.0
    h = 0.0
    shear = 0.0
    theta = 0.0
    Mayorppalstress = 0.0
    if ((Ephi(i).NE.0.0).and.(Ecohesion(i).NE.0.0)) then
        !salto a calcular la sigmac del macizo
        go to 551
    else if ((Enormalstress(i).EQ.0.0).and.(Minorppalstress(i).EQ.0.0)) then
        !salto a la asignacion de sigmac de probeta como la del macizo
        go to 552
    end if
    
    !Calculo de las constantes del material mmat y s
    if (macizo(i).EQ.1) then
        mmat = mi(i)*dexp((RMR(i)-100)/14)
        s = dexp((RMR(i)-100)/6)
    else if (macizo(i).EQ.2) then
        mmat = mi(i)*dexp((RMR(i)-100)/28)
        s = dexp((RMR(i)-100)/9)
    end if
    
    if (Enormalstress(i).NE.0.0) then
        !calculo con compresion normal efectiva
        h = 1.+((16.*(mmat*Enormalstress(i)+s*sigmac(i)))/(3*(mmat**2)*sigmac(i)))
        theta = (1/3.)*(90 + (180/pi)*datan(1/dsqrt(h**3-1)))
        Ephi(i) = (180/pi)*datan(1/dsqrt(4*h*dcos(theta*pi/180.)*dcos(theta*pi/180.)-1))
        shear = ((1/dtan(Ephi(i)*pi/180.))-dcos(Ephi(i)*pi/180.))*mmat*sigmac(i)/8.
        go to 550
    else if (Minorppalstress(i).NE.0.0) then
        !calculo con la fuerza a compresion principal mas pequeña, sigma 3
        Mayorppalstress = Minorppalstress(i) + dsqrt(mmat*Minorppalstress(i)*sigmac(i)+s*sigmac(i)**2)
        Enormalstress(i) = Minorppalstress (i) + (((Mayorppalstress-Minorppalstress(i))**2)/(2.*(Mayorppalstress-Minorppalstress(i))+(1/2.)*mmat*sigmac(i)))
        shear = (Enormalstress(i)-Minorppalstress(i))*dsqrt(1+((mmat*sigmac(i))/(2*(Mayorppalstress-Minorppalstress(i)))))
        Ephi(i) = 90 - (180/pi)*dasin(2*shear/(Mayorppalstress-Minorppalstress(i)))
    end if
    
    550 continue 
    Ecohesion(i) = shear - Enormalstress(i)*dtan((pi/180.)*Ephi(i))
    551 continue 
    sigmacmass(i) = 2*Ecohesion(i)*dcos(Ephi(i)*pi/180.)/(1-dsin(Ephi(i)*pi/180.))
    sigmacmass(i) = sigmacmass(i)*1000.
    !Unidades en el SI -> KPa = KN/m**2
    552 continue 
    if ((sigmac(i).NE.0.0).and.(sigmacmass(i).EQ.0.0)) then
        sigmacmass(i) = sigmac(i)*1000.
    else if (sigmac(i).EQ.0.0) then
        stop 'Error en el valor de la sigma c'
    end if
end do

allocate (sigmacrotwet(ndatosresis),stat=Isigmacrotwet)
if (Isigmacrotwet.NE.0) then
	stop 'Error de dimensionamiento de sigmacwet'
end if

!Calculo de la resistencia a compresion simple del macizo rocoso
do i=1,ndatosresis
    sigmacrotwet(i) = 0.0
    sigmacrotwet(i) = 2*Ecohesionwet(i)*dcos(Ephiwet(i)*pi/180.)/(1-dsin(Ephiwet(i)*pi/180.))
    sigmacrotwet(i) = sigmacrotwet(i)*1000.
    !Unidades en el SI -> KPa = KN/m**2
    if ((sigmacwet(i).NE.0.0).and.(sigmacrotwet(i).EQ.0.0)) then
        sigmacrotwet(i) = sigmacwet(i)*1000.
    else if (sigmacwet(i).EQ.0.0) then
        stop 'Error en el valor de la sigma c weathered'
    end if
end do
sigmacwetdef = sum(sigmacrotwet)/ndatosresis
Khdef = sum(Kh)/ndatosresis
Ephirotwet = sum(Ephiwet)/ndatosresis

!Creacion del vector R
allocate (R(NPrefino),stat=IR)
if (IR.NE.0) then
	stop 'Error de dimensionamiento del vector de resistencia de la roca'
end if

!Creacion del vector Rcuna
allocate (Rcuna(NPrefino),stat=IRcuna)
if (IRcuna.NE.0) then
	stop 'Error de dimensionamiento del vector de resistencia de la roca'
end if 

!Vector Rcuna
do i=1,NPrefino
    Rcuna(i) = Khdef*sigmacwetdef
end do

!Calculo de la constante R, para todo el acantilado, que se compone de K*Sigmac (o resist. a compresion simple)
k = 1
potacumulada = potestrato(k)
do i=1,NPrefino
	108 continue
	if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)) then !en el punto limite entre estratos coge las propiedades del de abajo
		R(i) = Kh(k)*sigmacmass(k)
		go to 107
	else
		k = k+1
		potacumulada = potacumulada + potestrato(k)
		go to 108
	end if
	107 continue
end do

!Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
do j=1,NPrefino
	if (PzyRglobal(1,j).LE.Mio(NO,7)) then
		goto 601
	end if	
end do
601 continue
irefinomax = j-1

!Calculo de la posicion más baja del vector PzyRglobal donde deja a afectar la marea
do j=1,NPrefino
	if (PzyRglobal(1,j).LT.Mio(NO,6)) then
		goto 701
	end if	
end do
701 continue
irefinomin = j-1

DimPzyR = abs(irefinomax-irefinomin)

!Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
do j=1,NPrefino
	if (PzyRglobal(1,j).LE.SLR) then
		goto 401
	end if	
end do
401 continue
indiceNM = j-1 !el ultimo valor positivo o mayor que cero
aa = NPrefino - irefinomin !coef. del numero de puntos bajo la afeccion del mar

!Contador de los dibujos
j = 1
jdibujo = 1
tdatos = int(Tdib/tmarea)
!Calculo de las veces que se dibujan los resultados
kk = tsim/Tdib
kkentero = int(kk)
DimDiby = kkentero

!Calculo del vector de los tiempos de rotura
jj = 1
!tdatosrotura = int(Trotura/tmarea)
call seed (iseed)
p = 0
trotacumulada = 0.0
select case (tiporotura)
	case(1) !caso determinista
		do while (trotacumulada.LE.tsim)
			ga = Trotura !HORAS
			trotacumulada = trotacumulada + ga
			p = p + 1
		end do
	case(2) !caso estocástico
		do while (trotacumulada.LE.tsim)
			call TroturaG(UTrotura,shape,scale,ga)
			trotacumulada = trotacumulada + ga
			p = p + 1
		end do
	case(3) !caso estocástico uniforme
		do while (trotacumulada.LE.tsim)
			call TroturaU(UTrotura,Umin,Umax,ga) !siendo Umin y Umax los valores minimo y maximo que queremos tener de intervalos de tiempo entre rotura
			trotacumulada = trotacumulada + ga
			p = p + 1
		end do
	case(4) !caso Pareto simple
		do while (trotacumulada.LE.tsim)
			call TroturaP(UTrotura,shapep,Xm,ga) !con shapep el valor de alpha y Xm el minimo valor de x
			trotacumulada = trotacumulada + ga
			p = p + 1
		end do
	case(5) !caso Pareto bounded
		do while (trotacumulada.LE.tsim)
			call TroturaPB(UTrotura,shapepb,Umin,Umax,ga) !
			trotacumulada = trotacumulada + ga
			p = p + 1
		end do
end select
if (p.EQ.1) then
	dimrotacumulada = 1
else
	dimrotacumulada = p - 1
end if

!dimensionamiento del vector de roturas acumuladas
allocate (rotacumulada(dimrotacumulada),stat=Irotacumulada)
if (Irotacumulada.NE.0) then
	stop 'Error de dimensionamiento del vector de roturas acumuladas'
end if
call seed(iseed)
select case (tiporotura)
	case(1) !caso determinista
		rotacumulada(1) = Trotura !HORAS
		write(5,'(20(F15.7,x))')Trotura,rotacumulada(1)
		do i = 2,dimrotacumulada
			rotacumulada(i) = rotacumulada(i-1) + Trotura
			write(5,'(20(F15.7,x))')Trotura,rotacumulada(i)
		end do
	case(2) !caso estocástico
		call TroturaG(UTrotura,shape,scale,ga)
		rotacumulada(1) = ga
		write(5,'(20(F15.7,x))')ga,rotacumulada(1)
		do i = 2,dimrotacumulada
			call TroturaG(UTrotura,shape,scale,ga)
			rotacumulada(i) = rotacumulada(i-1) + ga
			write(5,'(20(F15.7,x))')ga,rotacumulada(i)
		end do
	case(3) !caso estocástico uniforme
		call TroturaU(UTrotura,Umin,Umax,ga) !siendo Umin y Umax los valores minimo y maximo que queremos tener de intervalos de tiempo entre rotura
		rotacumulada(1) = ga
		write(5,'(20(F15.7,x))')ga,rotacumulada(1)
		do i = 2,dimrotacumulada
			call TroturaU(UTrotura,Umin,Umax,ga)
			rotacumulada(i) = rotacumulada(i-1) + ga
			write(5,'(20(F15.7,x))')ga,rotacumulada(i)
		end do
	case(4) !caso Pareto simple
		call TroturaP(UTrotura,shapep,Xm,ga)  !con shapep el valor de alpha y Xm el minimo valor de x
		rotacumulada(1) = ga
		write(5,'(20(F15.7,x))')ga,rotacumulada(1)
		do i = 2,dimrotacumulada
			call TroturaP(UTrotura,shapep,Xm,ga) 
			rotacumulada(i) = rotacumulada(i-1) + ga
			write(5,'(20(F15.7,x))')ga,rotacumulada(i)
		end do
	case(5) !caso Pareto bounded
		call TroturaPB(UTrotura,shapepb,Umin,Umax,ga)  !con shapep el valor de alpha y Xm el minimo valor de x
		rotacumulada(1) = ga
		write(5,'(20(F15.7,x))')ga,rotacumulada(1)
		do i = 2,dimrotacumulada
			call TroturaPB(UTrotura,shapepb,Umin,Umax,ga) 
			rotacumulada(i) = rotacumulada(i-1) + ga
			write(5,'(20(F15.7,x))')ga,rotacumulada(i)
		end do
end select

close(5) !cierre del archivo que escribe los tiempos de rotura y sus acumulados

!dimensionamiento de la matriz de resultados	
allocate (DibMatriz(NPrefino+1,DimDiby+2),stat=IDibMatriz)
if (IDibMatriz.NE.0) then
	stop 'Error de dimensionamiento de la matriz DibMatriz'
end if

!dimensionamiento de la matriz de roturas
allocate (Matrizrotura(NPrefino+1,2),stat=IMatrizrotura)
if (IMatrizrotura.NE.0) then
	stop 'Error de dimensionamiento de la matriz Matrizrotura'
end if

!dimensionamiento de la matriz de recesion
allocate (Matrizrecesion(NPrefino+2,2),stat=IMatrizrecesion)
if (IMatrizrecesion.NE.0) then
	stop 'Error de dimensionamiento de la matriz Matrizrecesion'
end if

!dimensionamiento de la matriz de difereniciales de y para cada io
allocate (Dymt(DimPzyR,io),stat=IDymt)
if (IDymt.NE.0) then
	stop 'Error de dimensionamiento de la matriz Dymt'
end if

!dimensionamiento de la matriz de difereniciales de y global
allocate (Dy(DimPzyR),stat=IDy)
if (IDy.NE.0) then
	stop 'Error de dimensionamiento de la matriz Dy'
end if

!solo creara y dimensionara los vectores de puntos/pesos de Gauss si son necesarios
if (tipointegral.EQ.6) then
	!dimensionamiento del vector de puntos de Gauss
	allocate (Pgauss(NPgauss),stat=IPgauss)
	if (IPgauss.NE.0) then
		stop 'Error de dimensionamiento del vector de puntos de Gauss'
	end if

	!dimensionamiento del vector de pesos de Gauss
	allocate (Wgauss(NPgauss),stat=IWgauss)
	if (IWgauss.NE.0) then
		stop 'Error de dimensionamiento del vector de pesos de Gauss'
	end if

	! dimensionamiento del vector de tiempos para la integral de Gauss
	allocate (wtg(NPgauss),stat=Iwtg)
	if (Iwtg.NE.0) then
		stop 'Error de dimensionamiento del vector de tiempos para la integral de Gauss'
	end if

	!llamada a los vectores de puntos / pesos de Gauss
	call gauss(NPgauss,Pgauss,Wgauss)
end if

!Creacion del vector de maximos erosivos
allocate (Imaxerosion(Nmaxerosion),stat=Indmaxerosion)
if (Indmaxerosion.NE.0) then
	stop 'Error de dimensionamiento del vector de maximos erosivos'
end if

!Creacion del vector de maximos erosivos
allocate (alpharot(ndatosresis),stat=Ialpharot)
if (Ialpharot.NE.0) then
	stop 'Error de dimensionamiento del angulo de rotura aleatorio'
end if

!Creacion del vector de maximos erosivos
allocate (indvertrot(ndatosresis),stat=Iindvertrot)
if (Iindvertrot.NE.0) then
	stop 'Error de dimensionamiento inversa de rotura'
end if

!Creacion del vector de maximos erosivos
allocate (mrotura(ndatosresis),stat=Imrotura)
if (Imrotura.NE.0) then
	stop 'Error de dimensionamiento de la pendiente de la rotura'
end if

!Creacion de la Matriz de RESULTADOS, dos primeras columnas
DibMatriz(1,1) = 0. !Tiempo inicial
DibMatriz(1,2) = 0. !
do i=1,NPrefino
	DibMatriz(i+1,1) = PzyRGlobal(1,i)
	DibMatriz(i+1,2) = PzyRGlobal(2,i)
end do

!Creacion de la Matriz de ROTURA, dos primeras columnas
Matrizrotura(1,1) = 0. !Tiempo inicial
Matrizrotura(1,2) = 0. !
do i=1,NPrefino
	Matrizrotura(i+1,1) = PzyRGlobal(1,i)
	Matrizrotura(i+1,2) = PzyRGlobal(2,i)
end do

!Creacion de la Matriz de RECESION, dos primeras columnas
Matrizrecesion(1,1) = 0. !Tiempo de erosion
Matrizrecesion(1,2) = 0. !
Matrizrecesion(2,1) = 0. !Area erosionada originalmente
Matrizrecesion(2,2) = 0. !
do i=1,NPrefino
	Matrizrecesion(i+2,1) = PzyRGlobal(1,i)
	Matrizrecesion(i+2,2) = PzyRGlobal(2,i)
end do

call cpu_time(tic)
call seed(iseed)
!para ver las pendientes de rotura
!do i=1,dimrotacumulada
!call random(rotale)
!alpharot = alpha + (rotale*2.*incertidumbre - incertidumbre)
!end do
select case (modelo)
case(1) !ESTADISTICO - el primer desarrollado
    select case (tipotg)
	    case (1) !Caso en el que usamos DOS puntos para el calculo de la tangente
	    !Calculo del bucle temporal global (Tglob->parametro de simulacion global)
		    do Tglob = 1,NTsim
    			
			    write(*,'(''Paso de tiempo: '',i8)')Tglob
			    do io = 1,NO
    							
				    !Calculo del diferencial de tiempo de marea
				    Dtm = tmarea / float(Ntmarea)
    				
				    !Inicializacion de la matriz de Diferenciales de y en una marea y con un tipo de ola
				    do k = 1,DimPzyR
					    Dymt(k,io) = 0.
				    end do
				    !Si la marea se lee de un archivo
				    !read(3,*)
				    !wt = w(t)

				    !Si la marea se genera como onda senoidal
				    select case (tipointegral)
					    case (1) !rectangulo a la izquierda
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta = SLR+(Am*dsin((t*2*pi)/tmarea))

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = Dtm*Mio(io,5)*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*Mio(io,5)*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
    						
					    case (2) !rectangulo a la derecha
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wtb = SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = Dtm*Mio(io,5)*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*Mio(io,5)*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do

					    case (3) !regla del punto medio
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = (Dtm)*Mio(io,5)*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*Mio(io,5)*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
    						
					    case (4) !regla del trapecio
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
							    wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
    						
					    case (5) !regla de Simpson
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
							    wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))
							    wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do

					    case (6) !cuadratura de Gauss
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm	
							    !obtencion de los tiempos de Gauss según los puntos en los que evaluamos
							    do i = 1,NPgauss
								    wtg(i) =SLR + (Am*dsin((((Dtm*Pgauss(i))+2*t+Dtm)*pi)/tmarea))
							    end do

							    do k = irefinomax+1,irefinomin-1,1
								    Dym = 0.
								    do i = 1,NPgauss
									    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
								    end do
								    Dym = Dym*(Dtm/2.)*Mio(io,5)
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
									    TanF2=tgrados
								    else
									    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = 0.
							    do i = 1,NPgauss
								    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
							    end do
							    Dym = Dym*(Dtm/2.)*Mio(io,5)
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
				    end select
			    end do
    			
			    !Calculo del diferencial de y para todas las mareas
			    do k = 1,DimPzyR
				    Dy(k) = 0.
				    do io = 1,NO
					    Dy(k) = Dy(k) + Dymt(k,io)
				    end do
			    end do
    		
			    !Calculo del nuevo perfil o valores de las y acumuladas
			    do k=irefinomax+1,irefinomin
				    PzyRglobal(2,k) = PzyRglobal(2,k) + Dy(k-irefinomax)
			    end do

			    !Matriz de escritura de RESULTADOS desde la tercera columna (t=1) hasta el final (t=j*tdatos)
			    if (Tglob.EQ.(j*tdatos).and.(j*Tdib).LE.tsim) then
				    !Creacion del cada vector de tiempos que introduciremos en la matriz
				    DibMatriz(1,j+2) = j*Tdib/(24.) !primera posicion del vector Y que indica el tiempo

				    !Volcamos la informacion del PzyRglobal en la matriz de dibujo
				    do i=1,NPrefino
					    DibMatriz(i+1,j+2) = PzyRglobal(2,i)
				    end do
    				
				    j = j + 1
			    end if

			    !Matriz de escritura de ROTURA desde la tercera columna (t=1) hasta el final (t=j*tdatos)	
			    if ((jj-1).EQ.dimrotacumulada) then
			        go to 999
			    else if (jj.LE.dimrotacumulada.and.Tglob.GE.(rotacumulada(jj)/tmarea)) then !.and.((jj*Trotura).LE.tsim)) then
				    !Para el calculo de los puntos maximos de erosion (puntos de rotura)
				    allocate (yerosionmax(indiceNM-irefinomax),stat=Iyerosionmax)
				    if (Iyerosionmax.NE.0) then
					    stop 'Error de dimensionamiento del vector yerosionmax'
				    end if
				    !calculo del Indice donde la erosion es maxima Above Sea water Level
				    if (Factor.EQ.1.0) then
					    do i=irefinomax+1,indiceNM
						    yerosionmax(i-irefinomax)=PzyRglobal(2,i)
					    end do
					    indmaxASL1 = minloc(yerosionmax)
					    Imaxerosion(1) = indmaxASL1(1) + irefinomax
				    else if (Factor.EQ.-1.0) then !Para I ->D
					    do i=irefinomax+1,indiceNM
						    yerosionmax(i-irefinomax)=PzyRglobal(2,i)
					    end do
					    indmaxASL1 = maxloc(yerosionmax)
					    Imaxerosion(1) = indmaxASL1(1)  + irefinomax
				    end if
				    deallocate (yerosionmax)

				    !calculo del Indice donde la erosion es maxima Under Sea water Level
				    if (Factor.EQ.1.0) then
					    !el punto mas bajo
					    ptoinf = PzyRglobal(2,irefinomin)
					    do i=irefinomin,indiceNM+2,-1
						    if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i-1)
						    else
							    iptoinf = i
							    go to 900
						    end if
					    end do
					    900 continue
					    if ((i+1).EQ.indiceNM+2) then
					        iptoinf = indiceNM+2
					    end if
					    Imaxerosion(3) = iptoinf

					    !punto en el nivel del mar
					    ptoinf = indiceNM+2
					    do i=indiceNM+2,irefinomin-1,1
						    if (PzyRglobal(2,i+1).LT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i+1)
						    else
							    iptoinf = i
							    go to 901
						    end if
					    end do
					    901 continue
					    Imaxerosion(2) = iptoinf
    					
				    else if (Factor.EQ.-1.0) then !Para I ->D
					    ptoinf = PzyRglobal(2,irefinomin)
					    do i=irefinomin,indiceNM+2,-1
						    if (PzyRglobal(2,i-1).GT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i-1)
						    else
							    iptoinf = i
							    go to 902
						    end if
					    end do
					    902 continue
					    if ((i+1).EQ.indiceNM+2) then
					        iptoinf = indiceNM+2
					    end if
					    Imaxerosion(3) = iptoinf

					    ptoinf = indiceNM+2
					    do i=indiceNM+2,irefinomin-1,1
						    if (PzyRglobal(2,i+1).GT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i+1)
						    else
							    iptoinf = i
							    go to 903
						    end if
					    end do
					    903 continue
					    Imaxerosion(2) = iptoinf

				    end if
    				
				    !comienzo de la rotura
				    dimmatrizrotura = size(Matrizrotura,dim=2)
				    allocate (dummy(NPrefino+1,dimmatrizrotura))
				    do n=1,dimmatrizrotura
					    do i=1,NPrefino+1
						    dummy(i,n) = Matrizrotura (i,n)
					    end do
				    end do

				    deallocate (Matrizrotura)
				    allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+2))

				    do n=1,dimmatrizrotura
					    do i=1,NPrefino+1
						    Matrizrotura(i,n) = dummy (i,n)
					    end do 
				    end do

				    deallocate (dummy)
    							
				    !Calculo del tiempo de rotura y posicionamiento en la matriz
				    Matrizrotura(1,dimmatrizrotura+1) = rotacumulada(jj)/(24.) !salida de datos en días
    				
				    !Matriz de dibujo de resultados cada trotura
				    do i=1,NPrefino
					    Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
				    end do
    				
				    !Calculo de la pendiente de rotura del acantilado
				    do tp=1,ndatosresis
				       call random(rotale)
				       alpharot(tp) = alpha(tp) + (rotale*2.*incertidumbre(tp) - incertidumbre(tp))
				       if (alpharot(tp).GE.90) then
					       indvertrot(tp) = 1
					       mrotura(tp) = 0.0
				       else
					       indvertrot(tp) = 0
					       mrotura(tp) = -tan((alpharot(tp)*pi)/180) * Factor
				       end if
				    end do
    				
                    select case (tipofractura)
					    case (1) !Rotura clasica
				            !ROTURA
				            tp = 1
    				        
				            !Caso 1: solo rompe desde ASL
				            if (Imaxerosion(2).EQ.Imaxerosion(3).and.Imaxerosion(2).EQ.indiceNM+2)then
					            !calculo de la parte desde ASL a max refino de oleaje
					            do i=Imaxerosion(1),irefinomax+1,-1
						            if(indvertrot(tp).EQ.1) then
							            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
						            else
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
						            end if
					            end do
					            Iptorotura = Imaxerosion(1)
				            else
					            !Caso 3: el BSL2 tiene menor o igual y que el ASL -> rotura desde BSL2
					            if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(1))).and.(PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2)))) then
						            !calculo desde BSL hasta refino de oleaje
						            do i=Imaxerosion(3),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(3))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(3))-PzyRglobal(1,Imaxerosion(3))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(3)
					            !Caso 2: rotura desde BSL2 hasta ASL y de ahi hasta refino de oleaje
					            else if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2))).and.(PzyRglobal(2,Imaxerosion(1)).LT.PzyRglobal(2,Imaxerosion(3)))) then
						            !calculo de la recta de rotura desde BSL2 a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(1))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(1),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						            do i=Imaxerosion(1),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(1)
					            !Caso 4: rompe desde BSL2 hasta BSL0, y de BSL0 a ASL
					            else 
						            !calculo de la recta de rotura desde BSL a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(2))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(2))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(2),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						            do i=Imaxerosion(1),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(1)
					            end if
				            end if

				            !Calculo del punto de interseccion entre rotura y acantilado
				            do 703 i=irefinomax,2,-1
					            !calculo de los vectores directores de cada pto del acantilado
					            u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
					            u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
					            !Calculo de la pendiente del acantilado
					            if (u1.EQ.0.0) then
						            indvertcliff = 1
						            mcliff = 0.0
					            else
						            indvertcliff = 0
						            mcliff = u2/u1
					            end if
					            !Evaluacion de los puntos de corte entre rotura y acantilado
					            if (((indvertrot(tp).EQ.indvertcliff).EQ.1).OR.(mrotura(tp).EQ.mcliff)) then
						            go to 703
					            else if (mrotura(tp).NE.mcliff) then
						            if ((mcliff.EQ.0.0).OR.(indvertcliff.EQ.1)) then
							            Yc = PzyRGlobal(2,i)
							            Zc = mrotura(tp)*Yc + PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)
						            else if ((mrotura(tp).EQ.0.0).OR.(indvertrot(tp).EQ.1)) then
							            Yc = PzyRglobal(2,Iptorotura)
							            Zc = mcliff*Yc + PzyRglobal(1,i) - mcliff*PzyRglobal(2,i)
						            else
							            Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)))/(mrotura(tp)-mcliff)
							            Zc = mrotura(tp)*Yc + PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)
						            end if
					            end if
					            !Calculo del parametro lambda que determina la posicion del 
					            !pto de corte en la recta del acantilado
					            if (indvertcliff.EQ.0.0) then
						            lambda = (Yc - PzyRglobal(2,i))/u1
					            else
						            lambda = (Zc - PzyRglobal(1,i))/u2
					            end if
					            !Calculo de la coordenada del punto de corte, en caso de producirse
					            if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
						            Ptocorte = i
						            go to 203
					            else if (lambda.EQ.1.0) then
						            Ptocorte = i-1
						            go to 203
					            else
						            go to 703
					            end if
				            703 continue
				            203 continue

				            !Definicion de la nueva geometria por encima del refino de marea hasta el punto de corte
				            !de la rotura con el acantilado
				            do i=irefinomax,ptocorte,-1
					            aux = 0.d0
					            if (indvertrot(tp).EQ.1.0) then
						            aux = PzyRglobal(2,Iptorotura)
					            else
						            aux = (PzyRglobal(1,i) + mrotura(tp)*PzyRglobal(2,Iptorotura) - PzyRglobal(1,Iptorotura))/mrotura(tp)
					            end if
					            !Caso 0: si el punto ASL es el punto donde empieza a afectar la marea y
					            !"crea" cuña de accreccion del acantilado
					            if (sign(1.0,(aux-PzyRGlobal(2,i))).EQ.(-Factor)) then
						            PzyRglobal(2,i) = aux
					            end if
				            end do
    				
				    case (2) !Rotura Joss Bay type generalizada
				            !ROTURA
    				         
				            !Caso 1: solo rompe desde ASL
				            if (Imaxerosion(2).EQ.Imaxerosion(3).and.Imaxerosion(2).EQ.indiceNM+2)then
					            !calculo de la parte desde ASL a max refino de oleaje
					            do i=Imaxerosion(1)-1,irefinomax+1,-1
					                tp=1
					                potacumulada = potestrato(tp)
					                1000 continue
					                if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)) then
						                if(indvertrot(tp).EQ.1) then
							                PzyRglobal(2,i) = PzyRglobal(2,i+1)
						                else
							                PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,i+1)-PzyRglobal(1,i+1)))/mrotura(tp)
						                end if
						                ptoestrato1 = i
						                go to 1001
						            else
						                tp = tp + 1
						                potacumulada = potacumulada + potestrato(tp)
						                go to 1000
						            end if
						            1001 continue
					            end do
					            Iptorotura = Imaxerosion(1)
					            go to 1008
				            else
					            !Caso 3: el BSL2 tiene menor o igual y que el ASL -> rotura desde BSL2
					            if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(1))).and.(PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2)))) then
						            !calculo desde BSL hasta refino de oleaje
						            do i=Imaxerosion(3)-1,irefinomax+1,-1
							            tp=1
					                    potacumulada = potestrato(tp)
					                    1002 continue
					                    if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)) then
						                    if(indvertrot(tp).EQ.1) then
							                    PzyRglobal(2,i) = PzyRglobal(2,i+1)
						                    else
							                    PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,i+1)-PzyRglobal(1,i+1)))/mrotura(tp)
						                    end if
						                    ptoestrato1 = i
						                    go to 1003
						                else
						                    tp = tp + 1
						                    potacumulada = potacumulada + potestrato(tp)
						                    go to 1002
						                end if
						                1003 continue
					                end do
						            Iptorotura = Imaxerosion(3)
						            go to 1008
					            !Caso 2: rotura desde BSL2 hasta ASL y de ahi hasta refino de oleaje
					            else if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2))).and.(PzyRglobal(2,Imaxerosion(1)).LT.PzyRglobal(2,Imaxerosion(3)))) then
						            !calculo de la recta de rotura desde BSL2 a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(1))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(1),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						            do i=Imaxerosion(1)-1,irefinomax+1,-1
							            tp=1
					                    potacumulada = potestrato(tp)
					                    1004 continue
					                    if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)) then
						                    if(indvertrot(tp).EQ.1) then
							                    PzyRglobal(2,i) = PzyRglobal(2,i+1)
						                    else
							                    PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,i+1)-PzyRglobal(1,i+1)))/mrotura(tp)
						                    end if
						                    ptoestrato1 = i
						                    go to 1005
						                else
						                    tp = tp + 1
						                    potacumulada = potacumulada + potestrato(tp)
						                    go to 1004
						                end if
						                1005 continue
					                end do
						            Iptorotura = Imaxerosion(1)
						            go to 1008
					            !Caso 4: rompe desde BSL2 hasta BSL0, y de BSL0 a ASL
					            else 
						            !calculo de la recta de rotura desde BSL a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(2))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(2))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(2),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						             do i=Imaxerosion(1)-1,irefinomax+1,-1
							            tp=1
					                    potacumulada = potestrato(tp)
					                    1006 continue
					                    if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)) then
						                    if(indvertrot(tp).EQ.1) then
							                    PzyRglobal(2,i) = PzyRglobal(2,i+1)
						                    else
							                    PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,i+1)-PzyRglobal(1,i+1)))/mrotura(tp)
						                    end if
						                    ptoestrato1 = i
						                    go to 1007
						                else
						                    tp = tp + 1
						                    potacumulada = potacumulada + potestrato(tp)
						                    go to 1006
						                end if
						                1007 continue
					                end do
						            Iptorotura = Imaxerosion(1)
						            go to 1008
					            end if
				            end if
				            1008 continue
    				        
				            tp = 1
				            potacumulada = potestrato(tp)
				            !Calculo del punto de interseccion entre rotura y acantilado
				            do 1009 i=irefinomax,2,-1
					            !calculo de los vectores directores de cada pto del acantilado
					            u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
					            u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
					            !Calculo de la pendiente del acantilado
					            if (u1.EQ.0.0) then
						            indvertcliff = 1
						            mcliff = 0.0
					            else
						            indvertcliff = 0
						            mcliff = u2/u1
					            end if
					            1010 continue
					            if (((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,i)).and.(PzyRglobal(1,i).LE.PzyRglobal(1,1)-(potacumulada-potestrato(tp)))) then
					            !Evaluacion de los puntos de corte entre rotura y acantilado
					                if (((indvertrot(tp).EQ.indvertcliff).EQ.1).OR.(mrotura(tp).EQ.mcliff)) then
						                go to 1011
					                else if (mrotura(tp).NE.mcliff) then
						                if ((mcliff.EQ.0.0).OR.(indvertcliff.EQ.1)) then
							                Yc = PzyRGlobal(2,i)
							                Zc = mrotura(tp)*Yc + PzyRglobal(1,ptoestrato1) - mrotura(tp)*PzyRglobal(2,ptoestrato1)
						                else if ((mrotura(tp).EQ.0.0).OR.(indvertrot(tp).EQ.1)) then
							                Yc = PzyRglobal(2,ptoestrato1)
							                Zc = mcliff*Yc + PzyRglobal(1,i) - mcliff*PzyRglobal(2,i)
						                else
							                Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,ptoestrato1) - mrotura(tp)*PzyRglobal(2,ptoestrato1)))/(mrotura(tp)-mcliff)
							                Zc = mrotura(tp)*Yc + PzyRglobal(1,ptoestrato1) - mrotura(tp)*PzyRglobal(2,ptoestrato1)
						                end if
					                end if
					                !Calculo del parametro lambda que determina la posicion del 
					                !pto de corte en la recta del acantilado
					                if (indvertcliff.EQ.0.0) then
						                lambda = (Yc - PzyRglobal(2,i))/u1
					                else
						                lambda = (Zc - PzyRglobal(1,i))/u2
					                end if
					                1011 continue
					                !Definicion de la nueva geometria por encima del refino de marea hasta el punto de corte
				                    !de la rotura con el acantilado
					                aux = 0.d0
					                if (indvertrot(tp).EQ.1.0) then
						                aux = PzyRglobal(2,ptoestrato1)
					                else
						                aux = (PzyRglobal(1,i) + mrotura(tp)*PzyRglobal(2,ptoestrato1) - PzyRglobal(1,ptoestrato1))/mrotura(tp)
					                end if
					                !Caso 0: si el punto ASL es el punto donde empieza a afectar la marea y
					                !"crea" cuña de accreccion del acantilado
					                if (sign(1.0,(aux-PzyRGlobal(2,i))).EQ.(-Factor)) then
						                PzyRglobal(2,i) = aux
						                ptoestrato2 = i
					                end if
					                !Calculo de la coordenada del punto de corte, en caso de producirse
					                if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
						                Ptocorte = i
						                go to 1012
					                else if (lambda.EQ.1.0) then
						                Ptocorte = i-1
						                go to 1012
					                else
						                go to 1009
					                end if
    					            
					            else if (PzyRglobal(1,1)-(potacumulada-potestrato(tp)).LE.PzyRglobal(1,i)) then
    					        
					                potacumulada = potacumulada - potestrato(tp)
					                tp = tp - 1
					                ptoestrato1 = ptoestrato2
					                go to 1010
    					            
					            else 
    					        
					                tp = tp + 1
					                potacumulada = potacumulada + potestrato(tp)
					                go to 1010
					            end if
    					            
				            1009 continue
				            1012 continue
				    end select
    				        
				    !Calculo del tiempo de rotura y posicionamiento en la matriz
				    Matrizrotura(1,dimmatrizrotura+2) = rotacumulada(jj)/(24.) !salida de datos en días
    				
				    !Matriz de dibujo de resultados cada trotura
				    do i=1,NPrefino
					    Matrizrotura(i+1,dimmatrizrotura+2) = PzyRglobal(2,i)
				    end do

				    !Creacion de la matriz de RECESION
				    dimmatrizrecesion = size(Matrizrecesion,dim=2)
				    allocate (dummy1(NPrefino+2,dimmatrizrecesion))
				    do n=1,dimmatrizrecesion
					    do i=1,NPrefino+2
						    dummy1(i,n) = Matrizrecesion (i,n)
					    end do
				    end do
				    deallocate (Matrizrecesion)
				    allocate (Matrizrecesion(NPrefino+2,dimmatrizrecesion+1))
				    do n=1,dimmatrizrecesion
					    do i=1,NPrefino+2
						    Matrizrecesion(i,n) = dummy1 (i,n)
					    end do
				    end do
				    deallocate (dummy1)

				    Matrizrecesion(1,dimmatrizrecesion+1) = rotacumulada(jj)/(24.) !salida de datos en días
				    do i= 2,NPrefino+1
					    Matrizrecesion(i+1,dimmatrizrecesion+1) = Matrizrotura(i,dimmatrizrotura+2) - Matrizrotura(i,dimmatrizrotura+1)
				    end do
				    sumaincrementoy = 0.0
				    do i=1,NPrefino
					    sumaincrementoy = sumaincrementoy + abs(Matrizrecesion(i+2,dimmatrizrecesion+1))
				    end do
				    Matrizrecesion(2,dimmatrizrecesion+1) = sumaincrementoy*Dz

                    write(2,11) rotacumulada(jj)/(24.),sumaincrementoy*Dz,(abs(Matrizrecesion(label(ilabel)+2,dimmatrizrecesion+1)),ilabel=1,nlabel)

				    jj = jj + 1
			    end if
			    999 continue

			    !Calculo de los posibles cambios en el nivel de marea por el CAMBIO CLIMATICO (CLIMATE CHANGES)
			    if (Tglob.EQ.tcambiodz.and.(tcambiodz*tmarea).LE.tsim) then
				    SLR = SLR + dz
    			    
			        do io = 1,NO
				        Mio(NO,7) = HmaxM + SLR
				        Mio(NO,6) = -Am - Mio(NO,4) + SLR
				    end do

				    !Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.Mio(NO,7)) then
						    goto 650
					    end if	
				    end do
				    650 continue
				    irefinomax = l-1

				    !Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
						    goto 651
					    end if	
				    end do
				    651 continue
				    irefinomin = l-1

				    DimPzyR = abs(irefinomax-irefinomin)

				    !Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.SLR) then
						    goto 652
					    end if	
				    end do
				    652 continue
				    indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar

				    !dimensionamiento de la matriz de diferenciales de y para cada io
				    deallocate (Dymt)
				    allocate (Dymt(DimPzyR,io),stat=IDymt)
				    if (IDymt.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dymt'
				    end if

				    !dimensionamiento de la matriz de diferenciales de y global
				    deallocate(Dy)
				    allocate (Dy(DimPzyR),stat=IDy)
				    if (IDy.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dy'
				    end if

				    tcambiodz = tcambiodz + int(tcambiosl)
			    end if
    			
			    !Calculo de los posibles cambios en el nivel de marea por la TECTONICA (TECTONICS)
			    if (Tglob.EQ.tcambiotec.and.(tcambiotec*tmarea).LE.tsim) then
    				
				    !Valido para MagM € [5.6;8.1] y AD € [0.05;8.] en metros
				    !Siendo AD = desplazamiento medio del terreno por el terremoto de magnitud MagM
				    !Wells y Coppersmith '94 -> a = -4.80 // b = 0.69
				    AD1 = -4.80 + 0.69 * MagM
				    AD = 10**(AD1)
    				
				    !Orientacion del desplazamiento; Otec = 1 -> sube el terreno; Otec = -1 -> baja el terreno
				    SLR = SLR - AD*Otec
    				
				    do io = 1,NO
				        Mio(NO,7) = HmaxM + SLR
				        Mio(NO,6) = -Am - Mio(NO,4) + SLR
				    end do

				    !Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.Mio(NO,7)) then
						    goto 653
					    end if	
				    end do
				    653 continue
				    irefinomax = l-1

				    !Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
						    goto 654
					    end if	
				    end do
				    654 continue
				    irefinomin = l-1

				    DimPzyR = abs(irefinomax-irefinomin)

				    !Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.SLR) then
						    goto 655
					    end if	
				    end do
				    655 continue
				    indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar

				    !dimensionamiento de la matriz de diferenciales de y para cada io
				    deallocate (Dymt)
				    allocate (Dymt(DimPzyR,io),stat=IDymt)
				    if (IDymt.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dymt'
				    end if

				    !dimensionamiento de la matriz de diferenciales de y global
				    deallocate(Dy)
				    allocate (Dy(DimPzyR),stat=IDy)
				    if (IDy.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dy'
				    end if

				    tcambiotec = tcambiotec + int(tcambiotec1)
			    end if
		    end do
    		
	    case (2) !Caso en el que usamos TRES puntos para el calculo de la tangente
	    !Calculo del bucle temporal global (Tglob->parametro de simulacion global)
		    do Tglob = 1,NTsim
			    write(*,'(''Paso de tiempo:'',i8)')Tglob
			    do io = 1,NO
    							
				    !Calculo del diferencial de tiempo de marea
				    Dtm = tmarea / float(Ntmarea)
    				
				    !Inicializacion de la matriz de Diferenciales de y en una marea y con un tipo de ola
				    do k = 1,DimPzyR
					    Dymt(k,io) = 0.
				    end do
    								
					    !Si la marea se lee de un archivo
					    !read(3,*)
					    !wt = w(t)

					    !Si la marea se genera como onda senoidal
				    select case (tipointegral)
					    case (1) !rectangulo a la izquierda
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta =SLR + (Am*dsin((t*2*pi)/tmarea))

							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = Dtm*(Mio(io,5))*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = Dtm*(Mio(io,5))*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*(Mio(io,5))*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1.0/dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
    						
					    case (2) !rectangulo a la derecha
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))
    					
							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = Dtm*(Mio(io,5))*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = Dtm*(Mio(io,5))*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*(Mio(io,5))*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1.0/dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do

					    case (3) !regla del punto medio
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = Dtm*(Mio(io,5))*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = Dtm*(Mio(io,5))*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*abs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = Dtm*(Mio(io,5))*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do

					    case (4) !regla del trapecio
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
							    wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))

							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*abs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
					    case (5) !regla de Simpson
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
							    wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))
							    wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*abs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do

					    case (6) !cuadratura de Gauss
						    !Calculo de la integral temporal-espacial del diferencial erosivo de y
						    t=-Dtm
						    do it = 1,Ntmarea
							    t=t+Dtm
							    !obtencion de los tiempos de Gauss según los puntos en los que evaluamos
							    do i = 1,NPgauss
								    wtg(i) =SLR + (Am*dsin((((Dtm*Pgauss(i))+2*t+Dtm)*pi)/tmarea))
							    end do

							    !Para el caso de K=1
							    k=irefinomax+1
							    Dym = 0.
							    do i = 1,NPgauss
								    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
							    end do
							    Dym = Dym*(Dtm/2.)*Mio(io,5)
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k+1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo superior
								    TanF1=(1/(2.*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k+1)-PzyRglobal(2,k+2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)

							    do k = irefinomax+2,irefinomin-1
								    Dym = 0.
								    do i = 1,NPgauss
									    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
								    end do
								    Dym = Dym*(Dtm/2.)*Mio(io,5)
								    if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.(1.0/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1)))*(2*Dz).gt.tgrados)then
									    TanF2=tgrados
								    else
									    !Formula de 3 puntos para punto medio
									    TanF1=(1.0/(2*Dz))*abs(PzyRglobal(2,k+1)-PzyRglobal(2,k-1))
									    TanF2=1.0/TanF1
								    end if
								    Dym=Dym*TanF2*Factor

								    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
							    end do
							    !Para el caso de K=NDZ(io)
							    Dym = 0.
							    do i = 1,NPgauss
								    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
							    end do
							    Dym = Dym*(Dtm/2.)*Mio(io,5)
							    if (PzyRglobal(2,k).eq.PzyRglobal(2,k-1).or.(1./dabs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2)))*(2*Dz).gt.tgrados)then
								    TanF2=tgrados
							    else
								    !Formula de 3 puntos para extremo inferior
								    TanF1=(1.0/(2*Dz))*abs(-3*PzyRglobal(2,k)+4*PzyRglobal(2,k-1)-PzyRglobal(2,k-2))
								    TanF2=1.0/TanF1
							    end if
							    Dym=Dym*TanF2*Factor
    						
							    Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
						    end do
				    end select
			    end do
    			
			    !Calculo del diferencial de y para todas las mareas
			    do k = 1,DimPzyR
				    Dy(k) = 0.
				    do io = 1,NO
					    Dy(k) = Dy(k) + Dymt(k,io)
				    end do
			    end do

			    !Calculo del nuevo perfil o valores de las y acumuladas
			    do k=irefinomax+1,irefinomin
				    PzyRglobal(2,k) = PzyRglobal(2,k) + Dy(k-irefinomax)
			    end do

			    !Matriz de escritura de ROTURA desde la tercera columna (t=1) hasta el final (t=j*tdatos)	
			    if ((jj-1).EQ.dimrotacumulada) then
			        go to 998
			    else if (jj.LE.dimrotacumulada.and.Tglob.GE.(rotacumulada(jj)/tmarea)) then !.and.((jj*Trotura).LE.tsim)) then
				    !Para el calculo de los puntos maximos de erosion (puntos de rotura)
				    allocate (yerosionmax(indiceNM-irefinomax),stat=Iyerosionmax)
				    if (Iyerosionmax.NE.0) then
					    stop 'Error de dimensionamiento del vector yerosionmax'
				    end if
				    !calculo del Indice donde la erosion es maxima Above Sea water Level
				    if (Factor.EQ.1.0) then
					    do i=irefinomax+1,indiceNM
						    yerosionmax(i-irefinomax)=PzyRglobal(2,i)
					    end do
					    indmaxASL1 = minloc(yerosionmax)
					    Imaxerosion(1) = indmaxASL1(1) + irefinomax
				    else if (Factor.EQ.-1.0) then !Para I ->D
					    do i=irefinomax+1,indiceNM
						    yerosionmax(i-irefinomax)=PzyRglobal(2,i)
					    end do
					    indmaxASL1 = maxloc(yerosionmax)
					    Imaxerosion(1) = indmaxASL1(1)  + irefinomax
				    end if
				    deallocate (yerosionmax)

				    !calculo del Indice donde la erosion es maxima Under Sea water Level
				    if (Factor.EQ.1.0) then
					    !el punto mas bajo
					    ptoinf = PzyRglobal(2,irefinomin)
					    do i=irefinomin,indiceNM+2,-1
						    if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i-1)
						    else
							    iptoinf = i
							    go to 904
						    end if
					    end do
					    904 continue
					    if ((i+1).EQ.indiceNM+2) then
					        iptoinf = indiceNM+2
					    end if
					    Imaxerosion(3) = iptoinf

					    !punto en el nivel del mar
					    ptoinf = indiceNM+2
					    do i=indiceNM+2,irefinomin-1,1
						    if (PzyRglobal(2,i+1).LT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i+1)
						    else
							    iptoinf = i
							    go to 905
						    end if
					    end do
					    905 continue
					    Imaxerosion(2) = iptoinf
    					
				    else if (Factor.EQ.-1.0) then !Para I ->D
					    ptoinf = PzyRglobal(2,irefinomin)
					    do i=irefinomin,indiceNM+2,-1
						    if (PzyRglobal(2,i-1).GT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i-1)
						    else
							    iptoinf = i
							    go to 906
						    end if
					    end do
					    906 continue
					    if ((i+1).EQ.indiceNM+2) then
					        iptoinf = indiceNM+2
					    end if
					    Imaxerosion(3) = iptoinf

					    ptoinf = indiceNM+2
					    do i=indiceNM+2,irefinomin-1,1
						    if (PzyRglobal(2,i+1).GT.PzyRglobal(2,i)) then
							    ptoinf = PzyRglobal(2,i+1)
						    else
							    iptoinf = i
							    go to 907
						    end if
					    end do
					    907 continue
					    Imaxerosion(2) = iptoinf

				    end if
    				
				    !comienzo de la rotura
				    dimmatrizrotura = size(Matrizrotura,dim=2)
				    allocate (dummy(NPrefino+1,dimmatrizrotura))
				    do n=1,dimmatrizrotura
					    do i=1,NPrefino+1
						    dummy(i,n) = Matrizrotura (i,n)
					    end do
				    end do

				    deallocate (Matrizrotura)
				    allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+2))

				    do n=1,dimmatrizrotura
					    do i=1,NPrefino+1
						    Matrizrotura(i,n) = dummy (i,n)
					    end do 
				    end do

				    deallocate (dummy)
    							
				    !Calculo del tiempo de rotura y posicionamiento en la matriz
				    Matrizrotura(1,dimmatrizrotura+1) = rotacumulada(jj)/(24.) !salida de datos en días
    				
				    !Matriz de dibujo de resultados cada trotura
				    do i=1,NPrefino
					    Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
				    end do
    				
				    !Calculo de la pendiente de rotura del acantilado
				    do tp=1,ndatosresis
				       call random(rotale)
				       alpharot(tp) = alpha(tp) + (rotale*2.*incertidumbre(tp) - incertidumbre(tp))
				       if (alpharot(tp).GE.90) then
					       indvertrot(tp) = 1
					       mrotura(tp) = 0.0
				       else
					       indvertrot(tp) = 0
					       mrotura(tp) = -tan((alpharot(tp)*pi)/180) * Factor
				       end if
				    end do
    				
				           !ROTURA
				           tp = 1
    				        
				           !Caso 1: solo rompe desde ASL
				           if (Imaxerosion(2).EQ.Imaxerosion(3).and.Imaxerosion(2).EQ.indiceNM+2)then
					            !calculo de la parte desde ASL a max refino de oleaje
					            do i=Imaxerosion(1),irefinomax+1,-1
						            if(indvertrot(tp).EQ.1) then
							            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
						            else
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
						            end if
					            end do
					            Iptorotura = Imaxerosion(1)
				           else
					            !Caso 3: el BSL2 tiene menor o igual y que el ASL -> rotura desde BSL2
					            if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(1))).and.(PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2)))) then
						            !calculo desde BSL hasta refino de oleaje
						            do i=Imaxerosion(3),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(3))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(3))-PzyRglobal(1,Imaxerosion(3))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(3)
					            !Caso 2: rotura desde BSL2 hasta ASL y de ahi hasta refino de oleaje
					            else if ((PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2))).and.(PzyRglobal(2,Imaxerosion(1)).LT.PzyRglobal(2,Imaxerosion(3)))) then
						            !calculo de la recta de rotura desde BSL2 a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(1))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(1),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						            do i=Imaxerosion(1),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(1)
					            !Caso 4: rompe desde BSL2 hasta BSL0, y de BSL0 a ASL
					            else 
						            !calculo de la recta de rotura desde BSL a ASL
						            mrecta = (PzyRglobal(1,Imaxerosion(2))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(2))-PzyRglobal(2,Imaxerosion(3)))
						            !calculo desde BSL hasta ASL
						            do i=Imaxerosion(3),Imaxerosion(2),-1
							            PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
						            end do
						            !calculo de la parte desde ASL a max refino de oleaje
						            do i=Imaxerosion(1),irefinomax+1,-1
							            if(indvertrot(tp).EQ.1) then
								            PzyRglobal(2,i) = PzyRglobal(2,Imaxerosion(1))
							            else
								            PzyRglobal(2,i) = (PzyRglobal(1,i) + (mrotura(tp)*PzyRglobal(2,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(1))))/mrotura(tp)
							            end if
						            end do
						            Iptorotura = Imaxerosion(1)
					            end if
				           end if

				            !Calculo del punto de interseccion entre rotura y acantilado
				            do 720 i=irefinomax,2,-1
					            !calculo de los vectores directores de cada pto del acantilado
					            u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
					            u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
					            !Calculo de la pendiente del acantilado
					            if (u1.EQ.0.0) then
						            indvertcliff = 1
						            mcliff = 0.0
					            else
						            indvertcliff = 0
						            mcliff = u2/u1
					            end if
					            !Evaluacion de los puntos de corte entre rotura y acantilado
					            if (((indvertrot(tp).EQ.indvertcliff).EQ.1).OR.(mrotura(tp).EQ.mcliff)) then
						            go to 720
					            else if (mrotura(tp).NE.mcliff) then
						            if ((mcliff.EQ.0.0).OR.(indvertcliff.EQ.1)) then
							            Yc = PzyRGlobal(2,i)
							            Zc = mrotura(tp)*Yc + PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)
						            else if ((mrotura(tp).EQ.0.0).OR.(indvertrot(tp).EQ.1)) then
							            Yc = PzyRglobal(2,Iptorotura)
							            Zc = mcliff*Yc + PzyRglobal(1,i) - mcliff*PzyRglobal(2,i)
						            else
							            Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)))/(mrotura(tp)-mcliff)
							            Zc = mrotura(tp)*Yc + PzyRglobal(1,Iptorotura) - mrotura(tp)*PzyRglobal(2,Iptorotura)
						            end if
					            end if
					            !Calculo del parametro lambda que determina la posicion del 
					            !pto de corte en la recta del acantilado
					            if (indvertcliff.EQ.0.0) then
						            lambda = (Yc - PzyRglobal(2,i))/u1
					            else
						            lambda = (Zc - PzyRglobal(1,i))/u2
					            end if
					            !Calculo de la coordenada del punto de corte, en caso de producirse
					            if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
						            Ptocorte = i
						            go to 721
					            else if (lambda.EQ.1.0) then
						            Ptocorte = i-1
						            go to 721
					            else
						            go to 720
					            end if
				            720 continue
				            721 continue

				            !Definicion de la nueva geometria por encima del refino de marea hasta el punto de corte
				            !de la rotura con el acantilado
				            do i=irefinomax,ptocorte,-1
					            aux = 0.d0
					            if (indvertrot(tp).EQ.1.0) then
						            aux = PzyRglobal(2,Iptorotura)
					            else
						            aux = (PzyRglobal(1,i) + mrotura(tp)*PzyRglobal(2,Iptorotura) - PzyRglobal(1,Iptorotura))/mrotura(tp)
					            end if
					            !Caso 0: si el punto ASL es el punto donde empieza a afectar la marea y
					            !"crea" cuña de accreccion del acantilado
					            if (sign(1.0,(aux-PzyRGlobal(2,i))).EQ.(-Factor)) then
						            PzyRglobal(2,i) = aux
					            end if
				            end do
    				
				    !Calculo del tiempo de rotura y posicionamiento en la matriz
				    Matrizrotura(1,dimmatrizrotura+2) = rotacumulada(jj)/(24.) !salida de datos en días
    				
				    !Matriz de dibujo de resultados cada trotura
				    do i=1,NPrefino
					    Matrizrotura(i+1,dimmatrizrotura+2) = PzyRglobal(2,i)
				    end do

				    !Creacion de la matriz de RECESION
				    dimmatrizrecesion = size(Matrizrecesion,dim=2)
				    allocate (dummy1(NPrefino+2,dimmatrizrecesion))
				    do n=1,dimmatrizrecesion
					    do i=1,NPrefino+2
						    dummy1(i,n) = Matrizrecesion (i,n)
					    end do
				    end do
				    deallocate (Matrizrecesion)
				    allocate (Matrizrecesion(NPrefino+2,dimmatrizrecesion+1))
				    do n=1,dimmatrizrecesion
					    do i=1,NPrefino+2
						    Matrizrecesion(i,n) = dummy1 (i,n)
					    end do
				    end do
				    deallocate (dummy1)

				    Matrizrecesion(1,dimmatrizrecesion+1) = rotacumulada(jj)/(24.) !salida de datos en días
				    do i= 2,NPrefino+1
					    Matrizrecesion(i+1,dimmatrizrecesion+1) = Matrizrotura(i,dimmatrizrotura+2) - Matrizrotura(i,dimmatrizrotura+1)
				    end do
				    sumaincrementoy = 0.0
				    do i=1,NPrefino
					    sumaincrementoy = sumaincrementoy + abs(Matrizrecesion(i+2,dimmatrizrecesion+1))
				    end do
				    Matrizrecesion(2,dimmatrizrecesion+1) = sumaincrementoy*Dz

                    write(2,11) rotacumulada(jj)/(24.),sumaincrementoy*Dz,(abs(Matrizrecesion(label(ilabel)+2,dimmatrizrecesion+1)),ilabel=1,nlabel)

				    jj = jj + 1
			    end if
			    998 continue

			    !Calculo de los posibles cambios en el nivel de marea por el CAMBIO CLIMATICO (CLIMATE CHANGES)
			    if (Tglob.EQ.tcambiodz.and.(tcambiodz*tmarea).LE.tsim) then
				    SLR = SLR + dz
    			    
			        do io = 1,NO
				        Mio(NO,7) = HmaxM + SLR
				        Mio(NO,6) = -Am - Mio(NO,4) + SLR
				    end do
    				
				    !Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.Mio(NO,7)) then
						    goto 453
					    end if	
				    end do
				    453 continue
				    irefinomax = l-1

				    !Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
						    goto 454
					    end if	
				    end do
				    454 continue
				    irefinomin = l-1

				    DimPzyR = abs(irefinomax-irefinomin)

				    !Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.SLR) then
						    goto 455
					    end if	
				    end do
				    455 continue
				    indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar

				    !dimensionamiento de la matriz de diferenciales de y para cada io
				    deallocate (Dymt)
				    allocate (Dymt(DimPzyR,io),stat=IDymt)
				    if (IDymt.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dymt'
				    end if

				    !dimensionamiento de la matriz de diferenciales de y global
				    deallocate(Dy)
				    allocate (Dy(DimPzyR),stat=IDy)
				    if (IDy.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dy'
				    end if

				    tcambiodz = tcambiodz + int(tcambiosl)
			    end if
    			
			    !Calculo de los posibles cambios en el nivel de marea por la TECTONICA (TECTONICS)
			    if (Tglob.EQ.tcambiotec.and.(tcambiotec*tmarea).LE.tsim) then
    				
				    !Valido para MagM € [5.6;8.1] y AD € [0.05;8.] en metros
				    !Siendo AD = desplazamiento medio del terreno por el terremoto de magnitud MagM
				    !Wells y Coppersmith '94 -> a = -4.80 // b = 0.69
				    AD1 = -4.80 + 0.69 * MagM
				    AD = 10**(AD1)
    				
				    !Orientacion del desplazamiento; Otec = 1 -> sube el terreno; Otec = -1 -> baja el terreno
				    SLR = SLR - AD*Otec
    						    
			        do io = 1,NO
				        Mio(NO,7) = HmaxM + SLR
				        Mio(NO,6) = -Am - Mio(NO,4) + SLR
				    end do

				    !Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.Mio(NO,7)) then
						    goto 456
					    end if	
				    end do
				    456 continue
				    irefinomax = l-1

				    !Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
						    goto 457
					    end if	
				    end do
				    457 continue
				    irefinomin = l-1

				    DimPzyR = abs(irefinomax-irefinomin)

				    !Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
				    do l=1,NPrefino
					    if (PzyRglobal(1,l).LE.SLR) then
						    goto 458
					    end if	
				    end do
				    458 continue
				    indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar

				    !dimensionamiento de la matriz de diferenciales de y para cada io
				    deallocate (Dymt)
				    allocate (Dymt(DimPzyR,io),stat=IDymt)
				    if (IDymt.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dymt'
				    end if

				    !dimensionamiento de la matriz de diferenciales de y global
				    deallocate(Dy)
				    allocate (Dy(DimPzyR),stat=IDy)
				    if (IDy.NE.0) then
					    stop 'Error de dimensionamiento de la matriz Dy'
				    end if

				    tcambiotec = tcambiotec + int(tcambiotec1)
			    end if
		    end do
    	
    end select
    
case(2)!PROCESO-RESPUESTA PURO - desarrollado durante la estancia en Leeds (SEP-DIC ´10)

    !Peso unitario del agua - KN/m**3:
    pesounitagua=9.81

    !Inicializacion de variables
    volacumulado=0.d0
    parametrorotura=0
    retroceso=0.5d0
    pendientetgPcirc1=dtan(1.*pi/180.)
    volacumuladopc=0.d0
    !Porque para calcular la interseccion del circulo de rotura con la superficie se lanzan verticales
    mcirculo=0.d0
    nslides=5.d0
    
    !creacion de los vectores de la cuña de derrubios
    allocate (vectorposoriginal(NPrefino))
    allocate (vectorverticalrot(NPrefino))
    do i=1,NPrefino
        vectorposoriginal(i) = PzyRglobal(2,i)
    end do
    ptocorteaux=1
    
    !Calculo de Hw para su dibujo
    Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
    Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
	dist = longitud/1000.
	do i=1,1000
	    distancia(i) = PzyRglobal(2,seepagelevel) - (i-1)*dist
	end do
	do i=1,1000
	    Hwdraw(i) = dsqrt((dabs(distancia(i)-PzyRglobal(2,seepagelevel))/longitud)*(Profseepage-Profagua)**2) + PzyRglobal(1,seepagelevel)
	end do
	!dimensionamiento de la matriz Hwdraw
    allocate (MatrizHwdraw(1000+1,2),stat=IMatrizHwdraw)
    if (IMatrizHwdraw.NE.0) then
	    stop 'Error de dimensionamiento de la matriz MatrizHwdraw'
    end if
	!Escribir Hwdraw
	MatrizHwdraw(1,1) = 0.d0
	MatrizHwdraw(1,2) = 0.d0
	do i=1,1000
	    MatrizHwdraw(i+1,1) = Hwdraw(i)
	    MatrizHwdraw(i+1,2) = distancia(i)
	end do
   
    !Calculo del bucle temporal global (Tglob->parametro de simulacion global)
    do Tglob = 1,NTsim
        write(*,'(''Paso de tiempo: '',i8)')Tglob
        !Para el calculo de los puntos maximos de erosion (puntos de rotura)
	    allocate (yerosionmax(indiceNM-irefinomax),stat=Iyerosionmax)
		if (Iyerosionmax.NE.0) then
		    stop 'Error de dimensionamiento del vector yerosionmax'
	    end if
	    !calculo del Indice donde la erosion es maxima Above Sea water Level
		if (Factor.EQ.1.0) then
		    do i=irefinomax+1,indiceNM
			    yerosionmax(i-irefinomax)=PzyRglobal(2,i)
		    end do
		    indmaxASL1 = minloc(yerosionmax)
		    Imaxerosion(1) = indmaxASL1(1) + irefinomax
		else if (Factor.EQ.-1.0) then !Para I ->D
			do i=irefinomax+1,indiceNM
		        yerosionmax(i-irefinomax)=PzyRglobal(2,i)
		    end do
			indmaxASL1 = maxloc(yerosionmax)
			Imaxerosion(1) = indmaxASL1(1)  + irefinomax
		end if
		deallocate (yerosionmax)

		!calculo del Indice donde la erosion es maxima Under Sea water Level
	    if (Factor.EQ.1.0) then
		    !el punto mas bajo
			ptoinf = PzyRglobal(2,irefinomin)
			do i=irefinomin,indiceNM+2,-1
			    if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i)) then
				    ptoinf = PzyRglobal(2,i-1)
				else
				    iptoinf = i
				    go to 2
				end if
			end do
			2 continue
		    if ((i+1).EQ.indiceNM+2) then
			    iptoinf = indiceNM+2
			end if
			Imaxerosion(3) = iptoinf

			!punto en el nivel del mar
			ptoinf = indiceNM+2
			do i=indiceNM+2,irefinomin-1,1
			    if (PzyRglobal(2,i+1).LT.PzyRglobal(2,i)) then
			        ptoinf = PzyRglobal(2,i+1)
			    else
				    iptoinf = i
				    go to 3
			    end if
			end do
		    3 continue
			Imaxerosion(2) = iptoinf
    					
	    else if (Factor.EQ.-1.0) then !Para I ->D
			ptoinf = PzyRglobal(2,irefinomin)
			do i=irefinomin,indiceNM+2,-1
			    if (PzyRglobal(2,i-1).GT.PzyRglobal(2,i)) then
		            ptoinf = PzyRglobal(2,i-1)
		        else
		            iptoinf = i
		            go to 4
			    end if
			end do
		    4 continue
	        if ((i+1).EQ.indiceNM+2) then
	            iptoinf = indiceNM+2
			end if
			Imaxerosion(3) = iptoinf
            
            !punto en el nivel del mar
            ptoinf = indiceNM+2
			do i=indiceNM+2,irefinomin-1,1
			    if (PzyRglobal(2,i+1).GT.PzyRglobal(2,i)) then
			        ptoinf = PzyRglobal(2,i+1)
			    else
			        iptoinf = i
			        go to 5
				end if
			end do
			5 continue
			Imaxerosion(2) = iptoinf
        end if
        
        !Calculo del punto mas profundo del notch
        !if (((Imaxerosion(2).EQ.indiceNM+2).and.(Imaxerosion(3).EQ.indiceNM+2))) then !.and.(dabs(PzyRglobal(2,Imaxerosion(1))-PzyRglobal(2,irefinomax)).LT.2.d0)) then
        !    roturanotch = Imaxerosion(3)
        if ((PzyRglobal(2,Imaxerosion(1)).LT.PzyRglobal(2,Imaxerosion(2))).and.(PzyRglobal(2,Imaxerosion(1)).LT.PzyRglobal(2,Imaxerosion(3))))then
            roturanotch = Imaxerosion(1)
        else if ((PzyRglobal(2,Imaxerosion(2)).LT.PzyRglobal(2,Imaxerosion(1))).and.(PzyRglobal(2,Imaxerosion(2)).LT.PzyRglobal(2,Imaxerosion(3))))then
            roturanotch = Imaxerosion(2)
        else
            roturanotch = Imaxerosion(3)
        end if
        
        IF (FACTOR.EQ.1)THEN
            IF (PzyRglobal(2,roturanotch).LE.vectorposoriginal(roturanotch))THEN
                !Continuar el bucle que calcula los circulos de rotura
                1 Continue
                !Retraemos desde el punto de rotura inicial (de metro en metro) -> segundo punto del circulo de rotura
                DistPcirc2 = PzyRglobal(2,roturanotch) - retroceso * Factor
                !Calculo interseccion circulo de rotura VS. acantilado
                do 6 i=irefinomax,2,-1
		            !calculo de los vectores directores de cada pto del acantilado
		            u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
			        u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
		            !Calculo de la pendiente del acantilado
			        if (u1.EQ.0.0) then
			            indvertcliff = 1
		                mcliff = 0.0
                    else
			            indvertcliff = 0
				        mcliff = u2/u1
			        end if
			        !Evaluacion de los puntos de corte entre circulo de rotura y acantilado
			        if (mcliff.EQ.mcirculo) then
			            go to 6
			        else if (mcirculo.EQ.0.d0) then
			            XDistPcirc2 = DistPcirc2
				        YDistPcirc2 = mcliff*XDistPcirc2 + PzyRglobal(1,i) - mcliff*PzyRglobal(2,i)
			        end if
			        !Calculo del parametro lambda que determina la posicion del 
			        !pto de corte en la recta del acantilado
			        if (indvertcliff.EQ.0.0) then
			            lambda = (XDistPcirc2 - PzyRglobal(2,i))/u1
			        else
			            lambda = (YDistPcirc2 - PzyRglobal(1,i))/u2
			        end if
			        !Calculo de la coordenada del punto de corte, en caso de producirse
			        if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
				        PtocorteCirc2 = i
				        go to 7
			        else if (lambda.EQ.1.0) then
			            PtocorteCirc2 = i-1
				        go to 7
			        else
			            go to 6
			        end if
		        6 continue
		        7 continue
        		
		        !Obtencion de la altura del talud en base al punto cresta talud y el nivel del mar
		        Htalud = abs(PzyRglobal(1,crestatalud) - PzyRglobal(1,indiceNM+1))
        		
		        !Limite de calculo de circulos de rotura -> cuando la distancia horizontal supera Htalud
		        if (DistPcirc2.LT.PzyRglobal(2,crestatalud)-Htalud) then
		            go to 69
		        end if
        		
		        !Obtenemos la pendiente de la union del punto de rotura VS. segundo punto del circulo
		        munion = (YDistPcirc2 - PzyRglobal(1,roturanotch))/(XDistPcirc2 - PzyRglobal(2,roturanotch))
        		
		        !Coordenadas del centro de rotura (SIGNO)
		        CXc = ((-PzyRglobal(2,roturanotch)/PendientetgPcirc1)-(PzyRglobal(1,roturanotch)/2.)+(YDistPcirc2/2.)-((-1/munion)*PzyRglobal(2,roturanotch)/2.)-((-1/munion)*XDistPcirc2/2.))/((-1/PendientetgPcirc1) + (1/munion))
		        CYc = (-1/PendientetgPcirc1)*(CXc - PzyRglobal(2,roturanotch)) + PzyRglobal(1,roturanotch)
        		
		        !Radio de la circunferencia
		        radio = dsqrt((CXc - PzyRglobal(2,roturanotch))**2 + (CYc - PzyRglobal(1,roturanotch))**2)
        		
		        !Tamaño del slide por la parte de arriba (b)
		        anchoslide = dabs(XDistPcirc2 - PzyRglobal(2,roturanotch))/nslides
        		
		        !Inicialización del numerador y denominador del Factor de Seguridad: FS
		        Num = 0.d0
		        Deno = 0.d0
        		
		        !Calculo del Factor de Seguridad o Safety Factor
		        do i=1,nslides
		            !Coordenada X del extremo mas proximo al inicio de la rotura
		            Xca = PzyRglobal(2,roturanotch) - (i-1)*anchoslide*Factor
		            !Coordenada X del extremo mas alejado al inicio de la rotura
		            Xcb = PzyRglobal(2,roturanotch) - i*anchoslide*Factor
		            !Coordenada X del punto medio del slide
		            Xmedio = (Xca+Xcb)/2.
        		    
		            !Coordenadas Y para el extremo mas proximo al inicio de la rotura
		            Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xca**2 + 2*CXc*Xca - CYc**2)))/2.
		            Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xca**2 + 2*CXc*Xca - CYc**2)))/2.
		            !Ver que opciones de y son las correctas -> la menor de las dos (parte inferior del circulo de rotura)
		            if (Yopc1.LT.Yopc2) then
		                Yca = Yopc1
		            else
		                Yca = Yopc2
		            end if
        		    
		            !Coordenadas Y para el extremo mas alejado al inicio de la rotura
		            Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xcb**2 + 2*CXc*Xcb - CYc**2)))/2.
		            Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xcb**2 + 2*CXc*Xcb - CYc**2)))/2.
		            !Ver que opciones de y son las correctas -> la menor de las dos (parte inferior del circulo de rotura)
		            if (Yopc1.LT.Yopc2) then
		                Ycb = Yopc1
		            else
		                Ycb = Yopc2
		            end if
        		    
		            !Coordenadas Y del punto medio del slide
		            Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xmedio**2 + 2*CXc*Xmedio - CYc**2)))/2.
		            Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xmedio**2 + 2*CXc*Xmedio - CYc**2)))/2.
		            !Ver que opciones de y son las correctas -> la menor de las dos (parte inferior del circulo de rotura)
		            if (Yopc1.LT.Yopc2) then
		                Ymedio = Yopc1
		            else
		                Ymedio = Yopc2
		            end if
        		  
		            !Pendiente punto medio del slide
		            Mmedio = (CXc-Xmedio)/(Ymedio-CYc)
        		    
		            !Calculo del angulo en el CENTRO del slide (por eso valor aboluto), sino el valor que de
		            gamma = dabs((180/pi)*datan(Mmedio))
        		    
		            !Interseccion con la topografia
		            do 8 l=roturanotch,2,-1
		                !calculo de los vectores directores de cada pto del acantilado
		                u1 = PzyRglobal(2,l-1) - PzyRglobal(2,l)
			            u2 = PzyRglobal(1,l-1) - PzyRglobal(1,l)
		                !Calculo de la pendiente del acantilado
			            if (u1.EQ.0.0) then
			                indvertcliff = 1
		                    mcliff = 0.0
                        else
			                indvertcliff = 0
				            mcliff = u2/u1
			            end if
			            !Evaluacion de los puntos de corte entre el punto medio de cada slide del circulo de rotura y acantilado
			            XPtomedio = Xmedio
				        YPtomedio = mcliff*XPtomedio + PzyRglobal(1,l) - mcliff*PzyRglobal(2,l)
			            !Calculo del parametro lambda que determina la posicion del 
			            !pto de corte en la recta del acantilado
			            if (indvertcliff.EQ.0.0) then
			                lambda = (XPtomedio - PzyRglobal(2,l))/u1
			            else
			                lambda = (YPtomedio - PzyRglobal(1,l))/u2
			            end if
			            !Calculo de la coordenada del punto de corte, en caso de producirse
			            if ((lambda.GE.0.0).AND.(lambda.LE.1.0)) then
				            go to 9
			            else
			                go to 8
			            end if
		            8 continue
		            9 continue
        		    
		            !Altura h del slide
		            Hslide = dabs(YPtomedio - Ymedio)
        		    
		            !Definicion de hw (agua subterranea)
		            if (groundwater.EQ.1) then
		                if (Xmedio.GE.PzyRglobal(2,seepagelevel)) then
		                    Hw = hslide
		                    go to 70
		                else if (Xmedio.LT.PzyRglobal(2,seepagelevel)) then
		                    Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
		                    Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
		                    distY = dabs(PzyRglobal(2,seepagelevel)-Xmedio)
		                    curvah = dsqrt((distY/longitud)*(Profseepage-Profagua)**2)
		                    Hw = dabs(curvah+PzyRglobal(1,seepagelevel))-Ymedio
		                    if (Hw.LE.0.d0) then
		                        Hw = 0.d0
		                    end if 
		                    go to 70
		                end if
		            else if (groundwater.EQ.0) then
		                Hw = 0.d0
		                go to 70
		            else if (groundwater.EQ.2) then
		                Hw = Hslide
		                go to 70
		            end if
		            70 continue 
        		    
		            !Estrato que esta al fondo del slide
		            tp = 1
		            potacumulada = potestrato (tp)
                    12 continue
		            if (((PzyRglobal(1,1)-potacumulada).LE.Ymedio).and.(Ymedio.LE.PzyRglobal(1,1)-(potacumulada-potestrato(tp)))) then
			           !Transformamos los MPa en KPa = KN/m**2
			           Ecohesionrot = Ecohesion(tp)*1000.
			           Ephirot = Ephi(tp)
			           Pesounitrot = Pesounit(tp)
			           go to 13		            
			        else
    		           tp = tp + 1
			           potacumulada = potacumulada + potestrato(tp)
			           go to 12
			        end if
			        13 continue
        			
			        !Calculo del Numerador
			        Num = Num + Ecohesionrot*(anchoslide/dcos(gamma*pi/180)) + (Hslide*Pesounitrot*anchoslide*dcos(gamma*pi/180) - Hw*pesounitagua*anchoslide/dcos(gamma*pi/180))*dtan(Ephirot*pi/180)
			        
   		            !Calculo del Denominador
		            Deno = Deno + Hslide*Pesounitrot*anchoslide*dsin((gamma*pi/180))
		           
                end do
                
                !Obtenemos el Factor de Seguridad (FS) global para todo el circulo de rotura
                !si dicho valor de FS es mayor que 1 no hay rotura: FS > 1 -> SEGURO
                FS = Num/Deno
                
                !Es estable o no, el circulo de rotura calculado
                if (FS.GE.1.d0) then
                    retroceso = retroceso + 0.5d0
                    go to 1
                else if (FS.LT.1.d0) then
                    retroceso = retroceso + 0.5d0
                    Ptoultrotura = PtocorteCirc2
                    !Pintamos en el archivo de datos caso de ser la primera rotura del instante Tglob
                    if (Parametrorotura.EQ.0) then
                        dimmatrizrotura = size(Matrizrotura,dim=2)
				        allocate (dummy(NPrefino+1,dimmatrizrotura))
				        do n=1,dimmatrizrotura
				            do i=1,NPrefino+1
					            dummy(i,n) = Matrizrotura (i,n)
					        end do
				        end do

				        deallocate (Matrizrotura)
				        allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+1))

				        do n=1,dimmatrizrotura
				            do i=1,NPrefino+1
					            Matrizrotura(i,n) = dummy (i,n)
					        end do 
			            end do

				        deallocate (dummy)
            							
				        !Calculo del tiempo de rotura y posicionamiento en la matriz
				        Matrizrotura(1,dimmatrizrotura+1) = Tglob*tmarea/(24.) !salida de datos en días
            				
				        !Matriz rotura-dibujo
				        do i=1,NPrefino
				            Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
				        end do
                    end if
                    
                    !Caso 1: rotura desde punto ASL
                    !if (roturanotch.EQ.Imaxerosion(1)) then
                        !Caso 1A: rompo desde BSL2 hasta ASL
                    !    if (PzyRglobal(2,Imaxerosion(3)).LE.PzyRglobal(2,Imaxerosion(2))) then
			        !		mrecta = (PzyRglobal(1,Imaxerosion(1))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(1))-PzyRglobal(2,Imaxerosion(3)))
			        !		do i=Imaxerosion(3),Imaxerosion(1),-1
			        !		    PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
			        !		end do
                    !    else if (PzyRglobal(2,Imaxerosion(2)).LE.PzyRglobal(2,Imaxerosion(3))) then
                            !Caso 1B: rompo desde BSL2 hasta BSL0
                    !        mrecta = (PzyRglobal(1,Imaxerosion(2))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(2))-PzyRglobal(2,Imaxerosion(3)))
			        !		do i=Imaxerosion(3),Imaxerosion(2),-1
			        !		    PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
			        !		end do
                    !    end if
                    !Caso 2: rotura desde BSL2 a BSL0 y de ahi hasta arriba
                    !else if ((roturanotch.EQ.Imaxerosion(2)).and.(roturanotch.NE.Imaxerosion(3))) then
                    !    mrecta = (PzyRglobal(1,Imaxerosion(2))-PzyRglobal(1,Imaxerosion(3)))/(PzyRglobal(2,Imaxerosion(2))-PzyRglobal(2,Imaxerosion(3)))
			        !    do i=Imaxerosion(3),Imaxerosion(2),-1
			        !	    PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Imaxerosion(3)) - PzyRglobal(1,Imaxerosion(3)))/mrecta
			        !    end do
                    !end if
                
                    !Rotura desde roturanotch
                    do i=roturanotch-1,PtocorteCirc2,-1
                        aux = 0.d0
                        aux = (2*CXc - dsqrt(4*CXc**2 + 4*(radio**2 - CXc**2 - CYc**2 + 2*CYc*PzyRglobal(1,i) - (PzyRglobal(1,i))**2)))/2.
                        if (Factor.EQ.1) then
                            if (aux.LE.PzyRglobal(2,i)) then
                                PzyRglobal(2,i) = aux
                                go to 14
                            else if (aux.GT.PzyRglobal(2,i))then
                                aux = 0.d0
                                aux = (2*CXc + dsqrt(4*CXc**2 + 4*(radio**2 - CXc**2 - CYc**2 + 2*CYc*PzyRglobal(1,i) - (PzyRglobal(1,i))**2)))/2.
                                PzyRglobal(2,i) = aux
                            end if
                        else if (Factor.EQ.-1) then
                            if (aux.GE.PzyRglobal(2,i)) then
                                PzyRglobal(2,i) = aux
                                go to 14
                            else if (aux.LT.PzyRglobal(2,i))then
                                aux = 0.d0
                                aux = (2*CXc + dsqrt(4*CXc**2 + 4*(radio**2 - CXc**2 - CYc**2 + 2*CYc*PzyRglobal(1,i) - (PzyRglobal(1,i))**2)))/2.
                                PzyRglobal(2,i) = aux
                            end if
                        end if
                        14 continue   
                    end do
                    
                    Vertcirc = 0
                    !Arreglo salto de rotura circular a perfil
                    do i=roturanotch-1,2,-1
                        if (Factor.EQ.1) then
                            if (PzyRglobal(2,i-1).LE.PzyRglobal(2,i)) then
                                if (Vertcirc.EQ.1) then
                                    go to 300
                                end if
                                !si "creo" materia me olvido de esta rotura, y vuelvo a empezar
                                if (PzyRglobal(2,i).GT.Matrizrotura(i+1,dimmatrizrotura+1)) then
                                    Parametrorotura = 0
                                    do j=1,NPrefino
				                        PzyRglobal (2,j) = Matrizrotura(j+1,dimmatrizrotura+1)
				                    end do
                                    
                                    allocate (dummy(NPrefino+1,dimmatrizrotura))
				                    do n=1,dimmatrizrotura
				                        do j=1,NPrefino+1
					                        dummy(j,n) = Matrizrotura (j,n)
					                    end do
				                    end do

				                    deallocate (Matrizrotura)
				                    allocate (Matrizrotura(NPrefino+1,dimmatrizrotura))

				                    do n=1,dimmatrizrotura
				                        do j=1,NPrefino+1
					                        Matrizrotura(j,n) = dummy (j,n)
					                    end do 
			                        end do

				                    deallocate (dummy)
                                    go to 1
                                end if
                                go to 301
                            else if (PzyRglobal(2,i-1).GT.PzyRglobal(2,i))then
                                PzyRglobal(2,i-1) = PzyRglobal(2,i)
                                Vertcirc = 1
                                !si "creo" materia me olvido de esta rotura, y vuelvo a empezar
                                if (PzyRglobal(2,i).GT.Matrizrotura(i+1,dimmatrizrotura+1)) then
                                    Parametrorotura = 0
                                    do j=1,NPrefino
				                        PzyRglobal (2,j) = Matrizrotura(j+1,dimmatrizrotura+1)
				                    end do
                                    
                                    allocate (dummy(NPrefino+1,dimmatrizrotura))
				                    do n=1,dimmatrizrotura
				                        do j=1,NPrefino+1
					                        dummy(j,n) = Matrizrotura (j,n)
					                    end do
				                    end do

				                    deallocate (Matrizrotura)
				                    allocate (Matrizrotura(NPrefino+1,dimmatrizrotura))

				                    do n=1,dimmatrizrotura
				                        do j=1,NPrefino+1
					                        Matrizrotura(j,n) = dummy (j,n)
					                    end do 
			                        end do

				                    deallocate (dummy)
                                    go to 1
                                end if
                                go to 301
                            end if
                        else if (Factor.EQ.-1) then
                            if (PzyRglobal(2,i-1).GE.PzyRglobal(2,i)) then
                                if (Vertcirc.EQ.1) then
                                    go to 300
                                end if
                                go to 301
                            else if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i))then
                                PzyRglobal(2,i-1) = PzyRglobal(2,i)
                                Vertcirc = 1
                                go to 301
                            end if
                        end if
                        301 continue
                    end do
                    300 continue
                    
                    Parametrorotura = 1
                    dimmatrizrotura = size(Matrizrotura,dim=2)
	                allocate (dummy(NPrefino+1,dimmatrizrotura))
		            do n=1,dimmatrizrotura
		                do i=1,NPrefino+1
			                dummy(i,n) = Matrizrotura (i,n)
			            end do
		            end do
		            deallocate (Matrizrotura)
		            allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+1))
		            do n=1,dimmatrizrotura
		                do i=1,NPrefino+1
			                Matrizrotura(i,n) = dummy (i,n)
			            end do 
		            end do
		            deallocate (dummy)
                    !Tiempo rotura
		            Matrizrotura(1,dimmatrizrotura+1) = Tglob*tmarea/(24.) !salida de datos en días			
		            !Matriz de dibujo de rotura
		            do i=1,NPrefino
		                Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
		            end do
            		
		            !Creacion de la matriz de RECESION
		            dimmatrizrecesion = size(Matrizrecesion,dim=2)
		            allocate (dummy1(NPrefino+2,dimmatrizrecesion))
                    do n=1,dimmatrizrecesion
		                do i=1,NPrefino+2
		                    dummy1(i,n) = Matrizrecesion (i,n)
		                end do
	                end do
		            deallocate (Matrizrecesion)
		            allocate (Matrizrecesion(NPrefino+2,dimmatrizrecesion+1))
		            do n=1,dimmatrizrecesion
		                do i=1,NPrefino+2
			                Matrizrecesion(i,n) = dummy1 (i,n)
			            end do
		            end do
		            deallocate (dummy1)
                    
                    !Tiempo de calculo
	                Matrizrecesion(1,dimmatrizrecesion+1) = Tglob*tmarea/(24.) !salida de datos en días
	                !Matriz de RECESION
	                do i= 2,NPrefino+1
    	                Matrizrecesion(i+1,dimmatrizrecesion+1) = Matrizrotura(i,dimmatrizrotura+1) - Matrizrotura(i,dimmatrizrotura)
	                end do
	                sumaincrementoy = 0.d0
		            do i=1,NPrefino
		                sumaincrementoy = sumaincrementoy + abs(Matrizrecesion(i+2,dimmatrizrecesion+1))
		            end do
		            Matrizrecesion(2,dimmatrizrecesion+1) = sumaincrementoy*Dz
		            Volacumulado = Volacumulado + sumaincrementoy*Dz
                    write(2,11) Tglob*tmarea/(24.),sumaincrementoy*Dz,Volacumulado,(abs(Matrizrecesion(label(ilabel)+2,dimmatrizrecesion+1)),ilabel=1,nlabel)
                
                    !Calculo de Hw para su dibujo
                    Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
                    Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
	                dist = longitud/1000.
	                do i=1,1000
	                    distancia(i) = PzyRglobal(2,seepagelevel) - (i-1)*dist
	                end do
	                do i=1,1000
	                    Hwdraw(i) = dsqrt((dabs(distancia(i)-PzyRglobal(2,seepagelevel))/longitud)*(Profseepage-Profagua)**2) + PzyRglobal(1,seepagelevel)
	                end do
	                !Escribir Hwdraw
	                dimmatrizHwdraw = size(MatrizHwdraw,dim=2)
	                allocate (dummy2(1000+1,dimmatrizHwdraw))
		            do n=1,dimmatrizHwdraw
		                do i=1,1000+1
			                dummy2(i,n) = MatrizHwdraw (i,n)
			            end do
		            end do
		            deallocate (MatrizHwdraw)
		            allocate (MatrizHwdraw(1000+1,dimmatrizHwdraw+2))
		            do n=1,dimmatrizHwdraw
		                do i=1,1000+1
			                MatrizHwdraw(i,n) = dummy2 (i,n)
			            end do 
		            end do
		            deallocate (dummy2)
                    !Tiempo rotura
		            MatrizHwdraw(1,dimmatrizHwdraw+1) = Tglob*tmarea/(24.) !salida de datos en días
		            MatrizHwdraw(1,dimmatrizHwdraw+2) = Tglob*tmarea/(24.) !salida de datos en días		
		            !Matriz de dibujo Hw
		            do i=1,1000
		                MatrizHwdraw(i+1,dimmatrizHwdraw+1) = Hwdraw(i)
	                    MatrizHwdraw(i+1,dimmatrizHwdraw+2) = distancia(i)
		            end do
                end if
                go to 1
                69 continue 
            END IF    
        END IF
        !Arreglo de la inestabilidad de rotura
        minestabilidad = 0.d0
        
        if (rotcunanotch.EQ.1) then
        !Rotura notch -> evaluacion de fuerzas (aproximacion geomecanica)
            Ptocorte = 0
            do 50 l=roturanotch-1,2,-1
                !calculo de los vectores directores de cada pto del acantilado
			    u1 = PzyRglobal(2,l-1) - PzyRglobal(2,l)
			    u2 = PzyRglobal(1,l-1) - PzyRglobal(1,l)
				!Calculo de la pendiente del acantilado
			    if (u1.EQ.0.0) then
			        indvertcliff = 1
			    	mcliff = 0.0
			    else
			    	indvertcliff = 0
			    	mcliff = u2/u1
			    end if
				!Evaluacion de los puntos de corte entre rotura y acantilado
			    if (mcliff.EQ.minestabilidad) then
			        go to 50
			    else
			   	    Yc = PzyRglobal(2,roturanotch)
			    	Zc = mcliff*Yc + PzyRglobal(1,l) - mcliff*PzyRglobal(2,l)
			    end if
				!Calculo del parametro lambda que determina la posicion del 
				!pto de corte en la recta del acantilado
			    if (indvertcliff.EQ.0.0) then
                    lambda = (Yc - PzyRglobal(2,l))/u1
			    else
			        lambda = (Zc - PzyRglobal(1,l))/u2
			    end if
				!Calculo de la coordenada del punto de corte, en caso de producirse
			    if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			        Ptocorte = l
			        go to 51
			    else if (lambda.EQ.1.0) then
			        Ptocorte = l-1
			    	go to 51
			    else
			        go to 50
			    end if
			    50 continue
		    51 continue
		    
		    !Sino hay punto de corte, salto porque no hay notch formado
		    if (Ptocorte.EQ.0) then
		        go to 54
		    end if
		    
		    !Vector vertical de rotura igual que cliff
		    do i=1,NPrefino
		        vectorverticalrot(i) = PzyRglobal(2,i)
		    end do
		    
		    !Plano vertical de rotura
		    nrot = 0
		    CoordX = 0.d0
		    do i = roturanotch,ptocorte,-1
		        vectorverticalrot(i) = PzyRglobal (2,roturanotch)
		        nrot = nrot + 1
		        CoordX = CoordX + dabs(PzyRglobal(2,i))
		    end do
		    
		    CoordX = CoordX/nrot
		    
		    !Localizacion del punto roturanotch en su estrato correspondiente
		    k=1
		    potacumulada = potestrato(k)
		    Sigmatensilerock = 0.d0
		    Pesounitrot = 0.d0
		    52 continue
		    if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,roturanotch)) then
			    !Transformar MPa en KPa = KN/m**2
			    Sigmatensilerock = Ecohesion(k)*0.1d0*1000.
			    Pesounitrot = Pesounit(k)
			    go to 53		            
			else
    		    k = k + 1
			    potacumulada = potacumulada + potestrato(k)
			    go to 52
			end if
			53 continue
			
			!Volumen de la posible cuña inestable
			sumaacumuladarot = 0.d0
			aux = 0.d0
			do i = roturanotch,ptocorte,-1
			    aux = abs(vectorverticalrot(i)-PzyRglobal(2,i))
			    sumaacumuladarot = sumaacumuladarot + aux
			end do
			!en KN/m**2
			sumaacumuladarot = sumaacumuladarot*dz*pesounitrot
			hrot = dabs(PzyRglobal(1,roturanotch)-PzyRglobal(1,ptocorte))
			Lcdg = dabs(PzyRglobal(2,roturanotch)) - CoordX
			Sigmatensilemax = (6.*sumaacumuladarot*Lcdg)/(hrot**2)
			
			!Se evalua la estabilidad de la cuña a vuelco
			if (Sigmatensilemax.GE.Sigmatensilerock) then
			    !rotura por vuelco
			    do i = Ptocorte,roturanotch
			        PzyRglobal(2,i) = vectorverticalrot(i)
			    end do
			end if
        
        else if (rotcunanotch.EQ.2) then
        !Rotura analitica con limite de irregularidad 0.5m
            !el punto mas bajo
		    ptoinf = PzyRglobal(2,indiceNM)
		    do i=indiceNM,irefinomax,-1
		        if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i)) then
			        ptoinf = PzyRglobal(2,i-1)
			    else
			        iptoinf = i
				    go to 22
			    end if
		    end do
	        iptoinf = indiceNM
		    22 continue
	        Ptoerosion1 = iptoinf

		    !punto en el nivel del mar
		    ptoinf = PzyRglobal(2,irefinomax)
		    do i=irefinomax,indiceNM-1,1
		        if (PzyRglobal(2,i+1).GT.PzyRglobal(2,i)) then
			        ptoinf = PzyRglobal(2,i+1)
			    else
			        iptoinf = i
				    go to 23
			    end if
		    end do
		    iptoinf = indiceNM
		    23 continue
		    Ptoerosion2 = iptoinf
    		
		    !Punto por debajo del nivel del mar
		    ptoinf = PzyRglobal(2,irefinomin)
		    do i=irefinomin,indiceNM+1,-1
		        if (PzyRglobal(2,i-1).LT.PzyRglobal(2,i)) then
			        ptoinf = PzyRglobal(2,i-1)
		        else
			        iptoinf = i
			        go to 44
			    end if
		    end do
		    iptoinf = indiceNM
		    44 continue
		    Ptoerosion3 = iptoinf
    		
		    if (PzyRglobal(2,Ptoerosion3).LE.PzyRglobal(2,Ptoerosion1)) then
		        Ptoerosion1 = Ptoerosion3
		    end if
            
            !Calculo rotura del notch en la vertical de manera discreta
            if (dabs(PzyRglobal(2,Ptoerosion1)-PzyRglobal(2,Ptoerosion2)).GT.0.5) then !.and.(PzyRglobal(2,roturanotch).LT.PzyRglobal(2,irefinomax)))then
                !Calculo del punto de interseccion entre rotura y acantilado
		        mrecta = (PzyRglobal(1,irefinomax)-PzyRglobal(1,Ptoerosion1))/(PzyRglobal(2,irefinomax)-PzyRglobal(2,Ptoerosion1))
			    do i=Ptoerosion1,roturanotch,-1
			        PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Ptoerosion1) - PzyRglobal(1,Ptoerosion1))/mrecta
			    end do
    		    
		        if (Ptoerosion1.EQ.roturanotch) then
		            do 25 l=Ptoerosion1-1,2,-1
			            !calculo de los vectores directores de cada pto del acantilado
			    	    u1 = PzyRglobal(2,l-1) - PzyRglobal(2,l)
			    	    u2 = PzyRglobal(1,l-1) - PzyRglobal(1,l)
				        !Calculo de la pendiente del acantilado
			    	    if (u1.EQ.0.0) then
			    	        indvertcliff = 1
			    		    mcliff = 0.0
			            else
			    	        indvertcliff = 0
			    		    mcliff = u2/u1
			    	    end if
				        !Evaluacion de los puntos de corte entre rotura y acantilado
			    	    if (mcliff.EQ. minestabilidad) then
			    	        go to 25
			    	    else
			   	            Yc = PzyRglobal(2,Ptoerosion1)
			    	        Zc = mcliff*Yc + PzyRglobal(1,l) - mcliff*PzyRglobal(2,l)
			    	    end if
				        !Calculo del parametro lambda que determina la posicion del 
				        !pto de corte en la recta del acantilado
			    	    if (indvertcliff.EQ.0.0) then
			    	        lambda = (Yc - PzyRglobal(2,l))/u1
			    	    else
			    	        lambda = (Zc - PzyRglobal(1,l))/u2
			    	    end if
				        !Calculo de la coordenada del punto de corte, en caso de producirse
			    	    if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			    	        Ptocorte = l
			    		    go to 24
			    	    else if (lambda.EQ.1.0) then
			    	        Ptocorte = l-1
			    		    go to 24
			    	    else
			    	        go to 25
			    	    end if
			        25 continue
			        24 continue

			        !Definicion de la nueva geometria por encima del refino de marea hasta el punto de corte
			        !de la rotura con el acantilado
			        do i=roturanotch,ptocorte,-1
			            PzyRglobal(2,i) = PzyRglobal(2,roturanotch)
			        end do
                end if
            end if
        end if
        
        54 continue !Viene de la no existencia de notch en caso geotecnico
        !Creacion de CUÑA
        
        !Pendiente de la cuña post-rotura
        mroturacuna = -tan((Ephirotwet*pi)/180.) * FACTOR
        !mroturacuna = 0.d0
        
        !cálculo del volumen de caidas después de lo que es arrastrado por el mar
        volacumuladopc = volacumulado - ((volacumulado*perdidacuna)/100.)
        if (perdidacuna.eq.100.) then
            volacumuladopc=0.d0
        end if
        selectrot=0
        controlmin=0.d0
        promedio2=0
        promedio4=0
        promedio6=0
        parametrocuna=0
        allocate (vectorcuna(NPrefino))
        allocate (vectorcunaaux(NPrefino))
        
        if ((Parametrorotura.EQ.1).and.(volacumuladopc.NE.0.d0)) then
            95 continue
            variablecontrolcuna=0
            do i=1,NPrefino
                vectorcuna(i) = PzyRglobal(2,i)
                vectorcunaaux(i) = vectorcuna(i)
            end do
            ptocorteaux=1
            !ptocorte = 0
            do j=indiceNM+1, NPrefino-1,1 !j=irefinomin+1, NPrefino-1,1
               ptocorte = 0
                !Vectores directores
                u1 = PzyRglobal(2,j+1) - PzyRglobal(2,j)
			    u2 = PzyRglobal(1,j+1) - PzyRglobal(1,j)
				!Calculo de la pendiente del acantilado
			    if (u1.EQ.0.0) then
			    	mcliffcuna = 0.0
			    else
			        mcliffcuna = u2/u1
			    end if
			    
			    do 30 i=j-1,2,-1
			        !Vectores directores
                    u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
			        u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
				    !Calculo de la pendiente del acantilado
			        if (u1.EQ.0.d0) then
			        	mcliff = 0.d0
			        else
			            mcliff = u2/u1
			        end if
			        !Evaluacion de los puntos de corte entre rotura y acantilado
			    	if (mroturacuna.EQ.mcliff) then
			    	    go to 30
			    	else if (mroturacuna.NE.mcliff) then 
			   	        if (mcliff.EQ.0.d0) then
			   	            Yc = PzyRglobal(2,i)
			    	        Zc = (mroturacuna)*Yc + PzyRglobal(1,j) - (mroturacuna)*PzyRglobal(2,j)
			    	    else
			    	        Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,j) - (mroturacuna)*PzyRglobal(2,j)))/((mroturacuna)-mcliff)
							Zc = (mroturacuna)*Yc + PzyRglobal(1,j) - (mroturacuna)*PzyRglobal(2,j)
			    	    end if
			    	end if
				    !Calculo del parametro lambda que determina la posicion del 
				    !pto de corte en la recta del acantilado
			    	if (mcliff.EQ.0.d0) then
			    	    lambda = (Zc - PzyRglobal(1,i))/u2
			    	else
			    	    lambda = (Yc - PzyRglobal(2,i))/u1
			    	end if
				    !Calculo de la coordenada del punto de corte, en caso de producirse
			    	if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			    	    Ptocorte = i
			    		go to 31
			    	else if (lambda.EQ.1.0) then
			    	    Ptocorte = i-1
			    		go to 31
			    	else
			    	    go to 30
			    	end if
			    30 continue
			    31 continue
			    
			    if (ptocorte.EQ.0) then
			        go to 32
			    end if
			    
			    do i=1,NPrefino
                    vectorcuna(i) = PzyRglobal(2,i)
                end do
			    
			    !vector cuña
                do i=j,ptocorte,-1
                    vectorcuna(i) = (PzyRglobal(1,i)+(mroturacuna)*PzyRglobal(2,j)-PzyRglobal(1,j))/(mroturacuna)
                end do
                
                !cálculo del volumen de la cuña
                sumaacumuladacuna = 0.d0
                aux = 0.d0
                do i=j,ptocorte,-1
                    aux = vectorcuna(i) - PzyRglobal(2,i)
                    sumaacumuladacuna = sumaacumuladacuna + aux
                end do
                sumaacumuladacuna = sumaacumuladacuna*dz
                
                !Verificar tamaño de cuña
                if (sumaacumuladacuna.LT.volacumuladopc) then
                    !si el volumen es < que el caido -> siguiente cuña, aunque guardo info en vector auxiliar
                    jj=j !para guardar la posición correcta
                    ptocorteaux=ptocorte
                    volaux=0.d0
                    promedio1=0.d0
                    do i=1,NPrefino,1
                        vectorcunaaux(i) = vectorcuna (i)
                    end do
                    volaux = sumaacumuladacuna
                    parametrocuna = 1
                    mtaludcorte = mroturacuna
                    promedio1 = dabs(volacumuladopc - volaux)
                    variablecontrolcuna=1
                    go to 32
                else if (sumaacumuladacuna.EQ.volacumuladopc) then
                    !si el volumen es = que el caido -> cuña correcta
                    jj=j
                    ptocorteaux=ptocorte
                    go to 33
                else if (sumaacumuladacuna.GT.volacumuladopc) then
                    !si el que menos error tiene es esta cuña es la que me quedo
                    if ((selectrot.EQ.2).and.(controlmin.EQ.2*promedio2)) then
                        jj=j !para guardar la posición correcta
                        ptocorteaux=ptocorte
                        volaux=0.d0
                        do i=1,NPrefino,1
                            vectorcunaaux(i) = vectorcuna (i)
                        end do
                        volaux = sumaacumuladacuna
                        parametrocuna = 1
                        mtaludcorte = mroturacuna
                        go to 94 !voy al final dando por buena esta cuña
                    !si el volumen es > que el caido -> cuña anterior
                    else if (parametrocuna.EQ.1) then
                        promedio2 = 0.d0
                        promedio2 = dabs(volacumuladopc-sumaacumuladacuna)
                        do i=1,NPrefino,1
                            vectorcuna(i) = vectorcunaaux(i)    
                        end do
                        sumaacumuladacuna = volaux
                        go to 33
                    else if (parametrocuna.EQ.0) then
                    !Si parametrocuna.EQ.0 tengo que irme al segundo criterio de creación de cuña
                        promedio1=1000.d0
                        promedio2=1000.d0
                        do i=1,NPrefino
                            vectorcuna(i) = PzyRglobal(2,i)
                        end do
                        go to 33
                    end if
                end if
                32 continue 
            end do
            33 continue
            
            if ((selectrot.EQ.1).and.(controlmin.EQ.promedio1)) then
                !no sigo con el calculo la pimera cuña es la de menor error
                go to 94
                !voy al final del calculo
            end if
            
            if (sumaacumuladacuna.LT.0.d0) then
                variablecontrolcuna=0
                do i=1,NPrefino
                    vectorcuna(i) = PzyRglobal(2,i)
                    vectorcunaaux(i) = vectorcuna(i)
                end do
                go to 555
            else if (sumaacumuladacuna.LT.volacumuladopc-volacumuladopc*0.05d0) then
                variablecontrolcuna=0
                do i=1,NPrefino
                    vectorcuna(i) = PzyRglobal(2,i)
                    vectorcunaaux(i) = vectorcuna(i)
                end do
                go to 555
            end if
            
            !Introducir el segundo criterio de la cuña -> SI PTOCORTE.EQ.0
            !parametrocuna.EQ.0 tengo que irme al segundo criterio de creación de cuña
            if ((ptocorte.EQ.0).or.(parametrocuna.EQ.0.)) then
                555 continue
                ptocortebase = 0
                select case (tipocuna)
                case(1)
                    !Primer tipo de cuña o tipo HUNGR / COROMINAS
                    !los parámetros del ajuste logarítmico están extraidos de COROMINAS 1996 para EARTHFLOWS
                    malcance = - 10**(-0.07*log(volacumuladopc)-0.214)*FACTOR
                    !malcance = - tan((30.*pi)/180.)*FACTOR
                    
                    !Calculo intersección angulo de ALCANCE con el terreno
                    do 34 w=Ptoultrotura+1,NPrefino-1,1
                        !Vectores directores
                        u1 = PzyRglobal(2,w+1) - PzyRglobal(2,w)
			            u2 = PzyRglobal(1,w+1) - PzyRglobal(1,w)
				        !Calculo de la pendiente del acantilado
			            if (u1.EQ.0.d0) then
			        	    mcliff = 0.d0
			            else
			                mcliff = u2/u1
			            end if
			            !Evaluacion de los puntos de corte entre pendiente HUNGR y acantilado
			    	    if (malcance.EQ.mcliff) then
			    	        go to 34
			    	    else if (mroturacuna.NE.mcliff) then 
			   	            if (mcliff.EQ.0.d0) then
			   	                Yc = PzyRglobal(2,w)
			    	            Zc = (malcance)*Yc + PzyRglobal(1,Ptoultrotura) - (malcance)*PzyRglobal(2,Ptoultrotura)
			    	        else
			    	            Yc = ((PzyRglobal(1,w)-mcliff*PzyRglobal(2,w)) - (PzyRglobal(1,Ptoultrotura) - (malcance)*PzyRglobal(2,Ptoultrotura)))/(malcance-mcliff)
							    Zc = (malcance)*Yc + PzyRglobal(1,Ptoultrotura) - (malcance)*PzyRglobal(2,Ptoultrotura)
			    	        end if
			    	    end if
				        !Calculo del parametro lambda que determina la posicion del 
				        !pto de corte en la recta del acantilado
			    	    if (mcliff.EQ.0.d0) then
			    	        lambda = (Zc - PzyRglobal(1,w))/u2
			    	    else
			    	        lambda = (Yc - PzyRglobal(2,w))/u1
			    	    end if
				        !Calculo de la coordenada del punto de corte, en caso de producirse
			    	    if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			    	        Ptocortebase = w
			    		    go to 35
			    	    else if (lambda.EQ.1.0) then
			    	        Ptocortebase = w+1
			    		    go to 35
			    	    else
			    	        go to 34
			    	    end if
                    34 continue
                    35 continue
                    
                    !por si no encuentra pto de corte con plataforma de marea salto al caso 3 o CASO EMPIRICO
                    if (Ptocortebase.EQ.0) then
                        go to 43
                    end if
                    
                    !Calculo de la Cuña <-> Volumen: con el uso de tangentes
                    aaa = 0.d0
                    jj=Ptocortebase !para guardar la posición correcta
                    36 continue
                    !Se varía la pendiente para el corte de la cuña de 1º en 1º hasta encontrar la pendiente correcta
                    aaa = aaa + 1
                    mtalud = - tan((aaa*pi)/180.)*FACTOR
                    do 92 i=Ptocortebase-2,Ptoultrotura+1,-1
!                        ptocorte=0
!                        mcorte=(abs(PzyRglobal(1,i))-abs(PzyRglobal(1,Ptocortebase)))/(abs(PzyRglobal(2,i))-abs(PzyRglobal(2,Ptocortebase)))
!                        if (abs(mcorte).LT.abs(mtalud)) then
!                            go to 37
!                        else if (abs(mcorte).GT.abs(mtalud)) then
!                            ptocorte = i+1
!                        else if (abs(mcorte).EQ.abs(mtalud)) then
!                            ptocorte = i
!                        end if
                        Ptocorte=0
                        !Vectores directores
                        u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
			            u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
				        !Calculo de la pendiente del acantilado
			            if (u1.EQ.0.d0) then
			        	    mcliff = 0.d0
			            else
			                mcliff = u2/u1
			            end if
			            !Evaluacion de los puntos de corte entre rotura y acantilado
			    	    if (mtalud.EQ.mcliff) then
			    	        go to 92
			    	    else if (mtalud.NE.mcliff) then 
			   	            if (mcliff.EQ.0.d0) then
			   	                Yc = PzyRglobal(2,i)
			    	            Zc = (mtalud)*Yc + PzyRglobal(1,Ptocortebase) - (mtalud)*PzyRglobal(2,Ptocortebase)
			    	        else
			    	            Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,Ptocortebase) - (mtalud)*PzyRglobal(2,Ptocortebase)))/((mtalud)-mcliff)
							    Zc = (mtalud)*Yc + PzyRglobal(1,Ptocortebase) - (mtalud)*PzyRglobal(2,Ptocortebase)
			    	        end if
			    	    end if
				        !Calculo del parametro lambda que determina la posicion del 
				        !pto de corte en la recta del acantilado
			    	    if (mcliff.EQ.0.d0) then
			    	        lambda = (Zc - PzyRglobal(1,i))/u2
			    	    else
			    	        lambda = (Yc - PzyRglobal(2,i))/u1
			    	    end if
				        !Calculo de la coordenada del punto de corte, en caso de producirse
			    	    if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			    	        Ptocorte = i
			    		    go to 93
			    	    else if (lambda.EQ.1.0) then
			    	        Ptocorte = i-1
			    		    go to 93
			    	    else
			    	        go to 92
			    	    end if
			            92 continue
			            93 continue
			            
			            if ((Ptocorte.EQ.0).and.(aaa.LT.90)) then
			                go to 36
                            !stop write (*,*) 'error cuña rotura'
                        else if ((Ptocorte.EQ.0).and.(aaa.GE.90)) then
                            go to 43
                        end if
                        
                        !Vector cuña
                        do aa=Ptocortebase,Ptocorte,-1
                            vectorcuna(aa) = (PzyRglobal(1,aa)+(mtalud)*PzyRglobal(2,Ptocortebase)-PzyRglobal(1,Ptocortebase))/(mtalud)
                        end do
                        !Volumen cuña
                        sumaacumuladacuna = 0.d0
                        aux = 0.d0
                        do aa=Ptocortebase,ptocorte,-1
                            aux = vectorcuna(aa) - PzyRglobal(2,aa)
                            sumaacumuladacuna = sumaacumuladacuna + aux
                        end do
                        sumaacumuladacuna = sumaacumuladacuna*dz
                        !Verificar tamaño de cuña
                        if (sumaacumuladacuna.LT.volacumuladopc) then
                            !si el volumen es < que el caido -> siguiente cuña, aunque guardo info en vector auxiliar
                            ptocorteaux=ptocorte
                            volaux=0.d0
                            do j=Ptocortebase,ptocorte,-1
                                vectorcunaaux(j) = vectorcuna (j)
                            end do
                            volaux = sumaacumuladacuna
                            parametrocuna = 1
                            mtaludcorte = mtalud !para guardarlo en archivo de salida pendiente de taludes ->transformarlo a grados
                            variablecontrolcuna=1
                            promedio3 = dabs(volacumuladopc - volaux)
                            go to 36
                        else if (sumaacumuladacuna.EQ.volacumuladopc) then
                            !si el volumen es = que el caido -> cuña correcta
                            ptocorteaux=ptocorte
                            go to 38
                        else if (sumaacumuladacuna.GT.volacumuladopc) then
                            !si el que menos error tiene es esta cuña es la que me quedo
                            if ((selectrot.EQ.4).and.(controlmin.EQ.2*promedio4)) then
                                ptocorteaux=ptocorte
                                volaux=0.d0
                                do i=1,NPrefino,1
                                    vectorcunaaux(i) = vectorcuna (i)
                                end do
                                volaux = sumaacumuladacuna
                                parametrocuna = 1
                                mtaludcorte = mtalud
                                variablecontrolcuna=1
                                go to 94 !voy al final dando por buena esta cuña
                            !si el volumen es > que el caido -> cuña anterior
                            else if (parametrocuna.EQ.0) then
                               promedio3=1000.d0
                               promedio4=1000.d0
                               !Si parametrocuna.EQ.0 tengo que irme al tercer criterio de creación de cuña
                               do w=1,NPrefino
                                    vectorcuna(w) = PzyRglobal(2,w)
                               end do
                               go to 43
                            else if (parametrocuna.EQ.1) then
                                promedio4 = 0.d0
                                promedio4 = dabs(volacumuladopc-sumaacumuladacuna)
                                !si el volumen es > que el caido -> cuña anterior
                                do j=Ptocortebase,ptocorte,-1
                                    vectorcuna(j) = vectorcunaaux(j)
                                end do
                                sumaacumuladacuna = volaux
                                go to 38
                            end if
                        end if                            
                        37 continue
                    !end do
                    38 continue
                    
                    if ((selectrot.EQ.3).and.(controlmin.EQ.promedio3)) then
                        !no sigo con el calculo la pimera cuña es la de menor error
                        go to 94
                        !voy al final del calculo
                    else if (sumaacumuladacuna.LT.volacumuladopc-volacumuladopc*0.05d0) then
                        variablecontrolcuna=0
                        do i=1,NPrefino
                            vectorcuna(i) = PzyRglobal(2,i)
                            vectorcunaaux(i) = vectorcuna(i)
                        end do
                        go to 43
                    end if
                                       
                case(2)
                    !Segundo tipo de cuña o modelo EMPIRICO (observaciones de campo)
                    !es continuación por si el segundo caso no tiene corte con terreno o el volumen de cuña es muy grande
                    43 continue
                    !Calculo del punto medio en z de la rotura en cuña
                    PtocorteCircMedioy = dabs(PzyRglobal(1,roturanotch))+(dabs(PzyRglobal(1,Ptoultrotura))-dabs(PzyRglobal(1,roturanotch)))/2.
                    do i=Ptoultrotura,roturanotch,1
                        if (PzyRglobal(1,i).LE.PtocorteCircMedioy) then
                            PtocorteCircMedio = i
                            go to 39
                        end if
                    end do
                    39 continue
                    
                    !Calculo de la Cuña <-> Volumen: con el uso de tangentes
                    aaa = 0.d0
                    ptocorteaux=PtocorteCircMedio !en este caso es al reves
                    40 continue
                    !Se varía la pendiente para el corte de la cuña de 1º en 1º hasta encontrar la pendiente correcta
                    aaa = aaa + 1
                    mtalud = - tan(((90-aaa)*pi)/180.)*FACTOR
                    do 90 i=PtocorteCircMedio+2,NPrefino,1
!                        mcorte=(abs(PzyRglobal(1,i))-abs(PzyRglobal(1,PtocorteCircMedio)))/(abs(PzyRglobal(2,i))-abs(PzyRglobal(2,PtocorteCircMedio)))
!                        if (abs(mcorte).GT.abs(mtalud)) then
!                            go to 41
!                        else if (abs(mcorte).LT.abs(mtalud)) then
!                            ptocorte = i-1
!                        else if (abs(mcorte).EQ.abs(mtalud)) then
!                            ptocorte = i
!                        end if
                        Ptocorte=0
                        !Vectores directores
                        u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
			            u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
				        !Calculo de la pendiente del acantilado
			            if (u1.EQ.0.d0) then
			        	    mcliff = 0.d0
			            else
			                mcliff = u2/u1
			            end if
			            !Evaluacion de los puntos de corte entre rotura y acantilado
			    	    if (mtalud.EQ.mcliff) then
			    	        go to 90
			    	    else if (mtalud.NE.mcliff) then 
			   	            if (mcliff.EQ.0.d0) then
			   	                Yc = PzyRglobal(2,i)
			    	            Zc = (mtalud)*Yc + PzyRglobal(1,PtocorteCircMedio) - (mtalud)*PzyRglobal(2,PtocorteCircMedio)
			    	        else
			    	            Yc = ((PzyRglobal(1,i)-mcliff*PzyRglobal(2,i)) - (PzyRglobal(1,PtocorteCircMedio) - (mtalud)*PzyRglobal(2,PtocorteCircMedio)))/((mtalud)-mcliff)
							    Zc = (mtalud)*Yc + PzyRglobal(1,PtocorteCircMedio) - (mtalud)*PzyRglobal(2,PtocorteCircMedio)
			    	        end if
			    	    end if
				        !Calculo del parametro lambda que determina la posicion del 
				        !pto de corte en la recta del acantilado
			    	    if (mcliff.EQ.0.d0) then
			    	        lambda = (Zc - PzyRglobal(1,i))/u2
			    	    else
			    	        lambda = (Yc - PzyRglobal(2,i))/u1
			    	    end if
				        !Calculo de la coordenada del punto de corte, en caso de producirse
			    	    if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			    	        Ptocorte = i
			    		    go to 91
			    	    else if (lambda.EQ.1.0) then
			    	        Ptocorte = i-1
			    		    go to 91
			    	    else
			    	        go to 90
			    	    end if
			            90 continue
			            91 continue
			            
			            if (Ptocorte.EQ.0) then
			                go to 40
                            !stop write (*,*) 'error cuña rotura'
                        end if
			            
                        !Vector cuña
                        do aa=PtocorteCircMedio,Ptocorte,1
                            vectorcuna(aa) = (PzyRglobal(1,aa)+(mtalud)*PzyRglobal(2,PtocorteCircMedio)-PzyRglobal(1,PtocorteCircMedio))/(mtalud)
                        end do
                        !Volumen cuña
                        sumaacumuladacuna = 0.d0
                        aux = 0.d0
                        do aa=PtocorteCircMedio,ptocorte,1
                            aux = vectorcuna(aa) - PzyRglobal(2,aa)
                            sumaacumuladacuna = sumaacumuladacuna + aux
                        end do
                        sumaacumuladacuna = sumaacumuladacuna*dz
                        !Verificar tamaño de cuña
                        if (sumaacumuladacuna.LT.volacumuladopc) then
                            !si el volumen es < que el caido -> siguiente cuña, aunque guardo info en vector auxiliar
                            jj=ptocorte
                            volaux=0.d0
                            do j=1,NPrefino !PtocorteCircMedio,ptocorte,1
                                vectorcunaaux(j) = vectorcuna (j)
                            end do
                            volaux = sumaacumuladacuna
                            parametrocuna = 1
                            mtaludcorte = mtalud !para guardarlo en archivo de salida pendiente de taludes ->transformarlo a grados
                            variablecontrolcuna=1
                            promedio5 = dabs(volacumuladopc - volaux)
                            go to 40
                        else if (sumaacumuladacuna.EQ.volacumuladopc) then
                            !si el volumen es = que el caido -> cuña correcta
                            jj=ptocorte
                            go to 42
                        else if (sumaacumuladacuna.GT.volacumuladopc) then
                            !si el que menos error tiene es esta cuña es la que me quedo
                            if ((selectrot.EQ.6).and.(controlmin.EQ.2*promedio6)) then
                                jj=ptocorte
                                volaux=0.d0
                                do i=1,NPrefino,1
                                    vectorcunaaux(i) = vectorcuna (i)
                                end do
                                volaux = sumaacumuladacuna
                                parametrocuna = 1
                                mtaludcorte = mtalud
                                variablecontrolcuna=1
                                go to 94 !voy al final dando por buena esta cuña
                            !si el volumen es > que el caido -> cuña anterior
                            else if (parametrocuna.EQ.0) then
                                promedio5=1000.d0
                                promedio6=1000.d0
                                !Si parametrocuna.EQ.0 error en la creación de la cuña
                                STOP 'FALLO EN LA CREACIÓN DE LA CUÑA -> CASO EMPIRICO u OBSERVACIONES DE CAMPO'
                            else if (parametrocuna.EQ.1) then
                                promedio6 = 0.d0
                                promedio6 = dabs(volacumuladopc-sumaacumuladacuna)
                                !si el volumen es > que el caido -> cuña anterior
                                do j=1,NPrefino !PtocorteCircMedio,ptocorte,1
                                    vectorcuna(j) = vectorcunaaux(j)
                                end do
                                sumaacumuladacuna = volaux
                                go to 42
                            end if
                        end if                            
                        41 continue
                    !end do
                    42 continue
                    if ((selectrot.EQ.5).and.(controlmin.EQ.promedio5)) then
                        !no sigo con el calculo la pimera cuña es la de menor error
                        go to 94
                        !voy al final del calculo
                    else if (sumaacumuladacuna.LT.volacumuladopc-volacumuladopc*0.05d0) then
                        variablecontrolcuna=0
                        do i=1,NPrefino
                            vectorcuna(i) = PzyRglobal(2,i)
                            vectorcunaaux(i) = vectorcuna(i)
                        end do
                    end if
                end select
            end if
            
            !Calculamos los valores de la menor desviacion con respecto al volacumuladopc de las distintas cuñas
            if (variablecontrolcuna.EQ.0) then
                if (promedio1.EQ.0.d0) then
                    promedio1=1000.d0
                else if (promedio2.EQ.0.d0) then
                    promedio2=1000.d0
                else if (promedio3.EQ.0.d0) then
                    promedio3=1000.d0
                else if (promedio4.EQ.0.d0) then
                    promedio4=1000.d0
                end if
                controlmin = min(promedio1,2*promedio2,promedio3,2*promedio4,promedio5,2*promedio6)
                !selecciono la cuña con menor error, y esa es con la que me quedo
                if (controlmin.EQ.promedio1) then
                    selectrot = 1
                    go to 95 !vuelvo a la primera cuña
                else if (controlmin.EQ.2*promedio2) then
                    selectrot = 2
                    go to 95 !vuelvo a la primera cuña
                else if (controlmin.EQ.promedio3) then
                    selectrot = 3
                    go to 555 !vuelvo a la primera cuña
                else if (controlmin.EQ.2*promedio4) then
                    selectrot = 4
                    go to 555 !vuelvo a la primera cuña
                else if (controlmin.EQ.promedio5) then
                    selectrot = 5
                    go to 43 !vuelvo a la primera cuña
                else if (controlmin.EQ.2*promedio6) then
                    selectrot = 6
                    go to 43 !vuelvo a la primera cuña
                end if
            end if
            
            94 continue
            !Introduzco los valores de la cuña sobre PzyRglobal
            do i=1,NPrefino,1 !jj,ptocorteaux,-1
                vectorposoriginal(i) = PzyRglobal(2,i)
                PzyRglobal(2,i) = vectorcuna(i)
            end do
        end if
        
        !Inicializo los vectores del cálculo de la cuña
        deallocate (vectorcuna)
        deallocate (vectorcunaaux)
        
		!Matriz de escritura de RESULTADOS desde la tercera columna (t=1) hasta el final (t=j*tdatos)
		if (Tglob.EQ.(jdibujo*tdatos).and.(jdibujo*Tdib).LE.tsim) then
		    !Creacion del cada vector de tiempos que introduciremos en la matriz
			DibMatriz(1,jdibujo+2) = jdibujo*Tdib/(24.) !primera posicion del vector Y que indica el tiempo
    	    !Volcamos la informacion del PzyRglobal en la matriz de dibujo
		    do i=1,NPrefino
    		    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	        end do
   		    jdibujo = jdibujo + 1
   		end if
		
		!Erosion debido a la ACCION MARINA
	    do io = 1,NO
    	    !Calculo del diferencial de tiempo de marea
		    Dtm = tmarea / float(Ntmarea)
    				
			!Inicializacion de la matriz de Diferenciales de y en una marea y con un tipo de ola
			do k = 1,DimPzyR
			    Dymt(k,io) = 0.
			end do
			!Si la marea se lee de un archivo
			!read(3,*)
			!wt = w(t)

			!Si la marea se genera como onda senoidal
			select case (tipointegral)
			case (1) !rectangulo a la izquierda
			    !Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm
					wta = SLR+(Am*dsin((t*2*pi)/tmarea))

					do k = irefinomax+1,irefinomin-1,1
					    Dym = Dtm*Mio(io,5)*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
					end do
					!Para el caso de K=NDZ(io)
					Dym = Dtm*Mio(io,5)*Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do
    						
			case (2) !rectangulo a la derecha
			    !Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm
				    wtb = SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))

					do k = irefinomax+1,irefinomin-1,1
					    Dym = Dtm*Mio(io,5)*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
					end do
					!Para el caso de K=NDZ(io)
					Dym = Dtm*Mio(io,5)*Pn(((wtb-PzyRglobal(1,k))/Mio(io,4)))
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do

			case (3) !regla del punto medio
				!Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm
					wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

					do k = irefinomax+1,irefinomin-1,1
					    Dym = (Dtm)*Mio(io,5)*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
					end do
					!Para el caso de K=NDZ(io)
					Dym = Dtm*Mio(io,5)*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do
						    
    		case (4) !regla del trapecio
				!Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm
					wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
					wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))

					do k = irefinomax+1,irefinomin-1,1
					    Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
			        end do
					!Para el caso de K=NDZ(io)
					Dym = (Dtm/2.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do
    						
			case (5) !regla de Simpson
				!Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm
					wta =SLR + (Am*dsin(((t)*2*pi)/tmarea))
					wtb =SLR + (Am*dsin(((t+Dtm)*2*pi)/tmarea))
					wtab2 =SLR + (Am*dsin(((t+(Dtm/2))*2*pi)/tmarea))

					do k = irefinomax+1,irefinomin-1,1
					    Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
					end do
					!Para el caso de K=NDZ(io)
					Dym = (Dtm/6.)*Mio(io,5)*(Pn(((wta-PzyRglobal(1,k))/Mio(io,4)))+4.*Pn(((wtab2-PzyRglobal(1,k))/Mio(io,4)))+Pn(((wtb-PzyRglobal(1,k))/Mio(io,4))))
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do

			case (6) !cuadratura de Gauss
				!Calculo de la integral temporal-espacial del diferencial erosivo de y
				t=-Dtm
				do it = 1,Ntmarea
				    t=t+Dtm	
				    !obtencion de los tiempos de Gauss según los puntos en los que evaluamos
					do i = 1,NPgauss
					    wtg(i) =SLR + (Am*dsin((((Dtm*Pgauss(i))+2*t+Dtm)*pi)/tmarea))
					end do

					do k = irefinomax+1,irefinomin-1,1
					    Dym = 0.
						do i = 1,NPgauss
						    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
						end do
						Dym = Dym*(Dtm/2.)*Mio(io,5)
						if (PzyRglobal(2,k+1).eq.PzyRglobal(2,k).or.Dz/dabs(PzyRglobal(2,k)-PzyRglobal(2,k+1)).gt.tgrados)then
						    TanF2=tgrados
						else
						    TanF2=(Dz/dabs(PzyRglobal(2,k+1)-PzyRglobal(2,k)))
						end if
						Dym=Dym*TanF2*Factor
						if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                        else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                            Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                        end if
					end do
					!Para el caso de K=NDZ(io)
					Dym = 0.
					do i = 1,NPgauss
					    Dym = Dym + Wgauss(i)*Pn(((wtg(i)-PzyRglobal(1,k))/Mio(io,4)))								
					end do
					Dym = Dym*(Dtm/2.)*Mio(io,5)
					Dym=Dym*TanF2*Factor
					if (PzyRglobal(2,k).GT.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/Rcuna(k) + Dymt(k-irefinomax,io)
                    else if (PzyRglobal(2,k).LE.vectorposoriginal(k)) then
                        Dymt(k-irefinomax,io) = Dym/R(k) + Dymt(k-irefinomax,io)
                    end if
				end do
			end select
		end do
    	
    	!Calculo del diferencial de y para todas las mareas
		do k = 1,DimPzyR
		    Dy(k) = 0.
			do io = 1,NO
			    Dy(k) = Dy(k) + Dymt(k,io)
			end do
		end do
    		
		!Calculo del nuevo perfil o valores de las y acumuladas
		do k=irefinomax+1,irefinomin
		    PzyRglobal(2,k) = PzyRglobal(2,k) + Dy(k-irefinomax)
		end do

	    !Calculo de los posibles cambios en el nivel de marea por el CAMBIO CLIMATICO (CLIMATE CHANGES)
		if (Tglob.EQ.tcambiodz.and.(tcambiodz*tmarea).LE.tsim) then
		    SLR = SLR + dz
    	    do io = 1,NO
		        Mio(NO,7) = HmaxM + SLR
				Mio(NO,6) = -Am - Mio(NO,4) + SLR
			end do

			!Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
			do l=1,NPrefino
			    if (PzyRglobal(1,l).LE.Mio(NO,7)) then
				    goto 15
				end if	
			end do
			15 continue
		    irefinomax = l-1

			!Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
			do l=1,NPrefino
			    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
				    goto 16
				end if	
			end do
			16 continue
			irefinomin = l-1

			DimPzyR = abs(irefinomax-irefinomin)
            !Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
			do l=1,NPrefino
			    if (PzyRglobal(1,l).LE.SLR) then
				    goto 17
				end if	
			end do
			17 continue
			indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar
            
            !dimensionamiento de la matriz de diferenciales de y para cada io
			deallocate (Dymt)
			allocate (Dymt(DimPzyR,io),stat=IDymt)
			if (IDymt.NE.0) then
			    stop 'Error de dimensionamiento de la matriz Dymt'
			end if

			!dimensionamiento de la matriz de diferenciales de y global
			deallocate(Dy)
			allocate (Dy(DimPzyR),stat=IDy)
			if (IDy.NE.0) then
			    stop 'Error de dimensionamiento de la matriz Dy'
			end if
            tcambiodz = tcambiodz + int(tcambiosl)
		end if
    			
		!Calculo de los posibles cambios en el nivel de marea por la TECTONICA (TECTONICS)
		if (Tglob.EQ.tcambiotec.and.(tcambiotec*tmarea).LE.tsim) then
    	    !Valido para MagM € [5.6;8.1] y AD € [0.05;8.] en metros
		    !Siendo AD = desplazamiento medio del terreno por el terremoto de magnitud MagM
		    !Wells y Coppersmith '94 -> a = -4.80 // b = 0.69
		    AD1 = -4.80 + 0.69 * MagM
		    AD = 10**(AD1)
        				
		    !Orientacion del desplazamiento; Otec = 1 -> sube el terreno; Otec = -1 -> baja el terreno
		    SLR = SLR - AD*Otec
        				
		    do io = 1,NO
		        Mio(NO,7) = HmaxM + SLR
			    Mio(NO,6) = -Am - Mio(NO,4) + SLR
		    end do

		    !Calculo de la posicion más alta del vector PzyRglobal donde empieza a afectar la marea
		    do l=1,NPrefino
		        if (PzyRglobal(1,l).LE.Mio(NO,7)) then
			        goto 18
			    end if	
		    end do
		    18 continue
		    irefinomax = l-1

			!Calculo de la posicion más baja del vector PzyRglobal donde deja de afectar la marea
			do l=1,NPrefino
			    if (PzyRglobal(1,l).LT.Mio(NO,6)) then
				    goto 19
				end if	
			end do
			19 continue
			irefinomin = l-1
            DimPzyR = abs(irefinomax-irefinomin)

			!Calculo de la posicion del vector PzyRglobal donde se hacen las Z negativas
			do l=1,NPrefino
			    if (PzyRglobal(1,l).LE.SLR) then
				    goto 21
				end if	
			end do
			21 continue
		    indiceNM = l-1 !el ultimo valor positivo o mayor que nivel del mar

			!dimensionamiento de la matriz de diferenciales de y para cada io
			deallocate (Dymt)
			allocate (Dymt(DimPzyR,io),stat=IDymt)
			if (IDymt.NE.0) then
			    stop 'Error de dimensionamiento de la matriz Dymt'
			end if

			!dimensionamiento de la matriz de diferenciales de y global
			deallocate(Dy)
			allocate (Dy(DimPzyR),stat=IDy)
			if (IDy.NE.0) then
			    stop 'Error de dimensionamiento de la matriz Dy'
			end if

			tcambiotec = tcambiotec + int(tcambiotec1)
		end if
		
		!Re-inicializacion de variables
        Volacumulado = 0.d0
        Parametrorotura = 0
        DistPcirc2 = 0.d0
        retroceso = 0.5d0
        
	end do
end select

call cpu_time(tac)

write(*,*)'tiempo CPU:  ',tac-tic,'  seg.'

11 format(F10.1,2x,F9.5,2x,10(F9.6,2x))

!Aspecto del perfil de la erosion del notch para cada periodo de marea
!Escritura de resultados para la matriz de datos cada Tdib
DibMatriz(1,2) = tac-tic !tiempo de CPU
do n = 1,NPrefino+1
	write (4,'(200(F15.7,x))') (DibMatriz(n,k),k=1,DimDiby+2)
end do

!Aspecto del perfil del acantilado antes-depues de cada caida
!Escritura de resultados para la matriz de rotura cada Trotura
if (dimmatrizrotura.EQ.size(Matrizrotura,dim=2)) then
    do n = 1,NPrefino+1
        write (6,'(200(F15.7,x))') (Matrizrotura(n,k),k=1,dimmatrizrotura+2)
    end do
else
    dimmatrizrotura = size(Matrizrotura,dim=2)
    do n = 1,NPrefino+1
	    write (6,'(200(F15.7,x))') (Matrizrotura(n,k),k=1,dimmatrizrotura)
    end do
end if
!Calculo del terreno retrocedido en cada rotura por cada punto
!Escritura de resultados para la matriz de recesion cada tiempo de rotura
if (dimmatrizrecesion.EQ.size(Matrizrecesion,dim=2)) then
    do n = 1,NPrefino+2
	    write (7,'(200(F15.7,x))') (Matrizrecesion(n,k),k=1,dimmatrizrecesion+1)
    end do
else
    dimmatrizrecesion = size(Matrizrecesion,dim=2)
    do n = 1,NPrefino+2
	    write (7,'(200(F15.7,x))') (Matrizrecesion(n,k),k=1,dimmatrizrecesion)
    end do
end if

!Calculo del agua subterranea
!Escritura de resultados para la matriz de agua subterrana cada tiempo de rotura
if (dimmatrizHwdraw.EQ.size(MatrizHwdraw,dim=2)) then
    do n = 1,1000+1
	    write (8,'(200(F15.7,x))') (MatrizHwdraw(n,k),k=1,dimmatrizHwdraw+1)
    end do
else
    dimmatrizHwdraw = size(MatrizHwdraw,dim=2)
    do n = 1,1000+1
	    write (8,'(200(F15.7,x))') (MatrizHwdraw(n,k),k=1,dimmatrizHwdraw)
    end do
end if

stop
end

! **************************************************************************
! **************************************************************************
! **************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!! FUNCIONES Y SUBRRUTINAS !!!!!!!!!!!!!!!!!!!!!!!!!!

!Subrrutina que proporciona los puntos y pesos de Gauss
subroutine gauss(NPgauss,Pgauss,Wgauss)
real*8 Pgauss(NPgauss),Wgauss(NPgauss)
!Creacion de la tabla de valores para integracion numerica con el metodo de Gauss
if (NPgauss.EQ.2) then
	Pgauss(1)=-0.5773502692
	Wgauss(1)=1.d0
	Pgauss(2)=-Pgauss(1)
	Wgauss(2)=Wgauss(1)
end if
if (NPgauss.EQ.3) then
	Pgauss(1)=-0.7745966692
	Wgauss(1)=0.5555555556
	Pgauss(2)=0.d0
	Wgauss(2)=0.8888888889
	Pgauss(3)=-Pgauss(1)
	Wgauss(3)=Wgauss(1)
end if
if (NPgauss.EQ.4) then
	Pgauss(1)=-0.8611361159
	Wgauss(1)=0.3478548451
	Pgauss(2)=-0.3399810436
	Wgauss(2)=0.6521451549
	Pgauss(3)=-Pgauss(2)
	Wgauss(3)=Wgauss(2)
	Pgauss(4)=-Pgauss(1)
	Wgauss(4)=Wgauss(1)
end if
if (NPgauss.EQ.5) then
	Pgauss(1)=-0.9061798459386639927976269
	Wgauss(1)=0.2369268850561890875142640
	Pgauss(2)=-0.5384693101056830910363144
	Wgauss(2)=0.4786286704993664680412915
	Pgauss(3)=0.d0
	Wgauss(3)=0.5688888888888888888888889
	Pgauss(4)=-Pgauss(2)
	Wgauss(4)=Wgauss(2)
	Pgauss(5)=-Pgauss(1)
	Wgauss(5)=Wgauss(1)
end if
if (NPgauss.EQ.6) then
	Pgauss(1)=-0.9324695142031520278123016
	Wgauss(1)=0.1713244923791703450402961
	Pgauss(2)=-0.6612093864662645136613996
	Wgauss(2)=0.3607615730481386075698335
	Pgauss(3)=-0.2386191860831969086305017
	Wgauss(3)=0.4679139345726910473898703
	Pgauss(4)=-Pgauss(3)
	Wgauss(4)=Wgauss(3)
	Pgauss(5)=-Pgauss(2)
	Wgauss(5)=Wgauss(2)
	Pgauss(6)=-Pgauss(1)
	Wgauss(6)=Wgauss(1)
end if
if (NPgauss.EQ.7) then
	Pgauss(1)=-0.9491079123427585245261897
	Wgauss(1)=0.1294849661688696932706114
	Pgauss(2)=-0.7415311855993944398638648
	Wgauss(2)=0.2797053914892766679014678
	Pgauss(3)=-0.4058451513773971669066064
	Wgauss(3)=0.3818300505051189449503698
	Pgauss(4)=0.d0
	Wgauss(4)=0.4179591836734693877551020
	Pgauss(5)=-Pgauss(3)
	Wgauss(5)=Wgauss(3)
	Pgauss(6)=-Pgauss(2)
	Wgauss(6)=Wgauss(2)
	Pgauss(7)=-Pgauss(1)
	Wgauss(7)=Wgauss(1)
end if
if (NPgauss.EQ.8) then
	Pgauss(1)=-0.9602898564975362
	Wgauss(1)=0.1012285362903763
	Pgauss(2)=-0.7966664774136267
	Wgauss(2)=0.2223810344533743
	Pgauss(3)=-0.5255324099163290
	Wgauss(3)=0.3137066458778873
	Pgauss(4)=-0.1834346424956498
	Wgauss(4)=0.3626837833783625
	Pgauss(5)=-Pgauss(4)
	Wgauss(5)=Wgauss(4)
	Pgauss(6)=-Pgauss(3)
	Wgauss(6)=Wgauss(3)
	Pgauss(7)=-Pgauss(2)
	Wgauss(7)=Wgauss(2)
	Pgauss(8)=-Pgauss(1)
	Wgauss(8)=Wgauss(1)
end if
if (NPgauss.EQ.9) then
	Pgauss(1)=-0.9681602395076261
	Wgauss(1)=0.08127438836157451
	Pgauss(2)=-0.8360311073266358
	Wgauss(2)=0.1806481606948574
	Pgauss(3)=-0.6133714327005904
	Wgauss(3)=0.2606106964029359
	Pgauss(4)=-0.3242534234038089
	Wgauss(4)=0.3123470770400022
	Pgauss(5)=0.d0
	Wgauss(5)=0.3302393550012596
	Pgauss(6)=-Pgauss(4)
	Wgauss(6)=Wgauss(4)
	Pgauss(7)=-Pgauss(3)
	Wgauss(7)=Wgauss(3)
	Pgauss(8)=-Pgauss(2)
	Wgauss(8)=Wgauss(2)
	Pgauss(9)=-Pgauss(1)
	Wgauss(9)=Wgauss(1)
end if 
return
end

!Subrrutina que proporciona un valor aleatorio que sigue una distribucion de gamma con shape>0
subroutine TroturaG(UTrotura,shape,scale,ga)
integer UTrotura
real*8 shape,scale,ga,non
if (shape.GT.0..and.shape.LE.0.4) then
	non = 1./shape
else if (shape.GT.0.4.and.shape.LE.4.) then
	non = 1./shape + (1./shape)*(shape - 0.4)*(1./3.6)
else 
	non = 1./dsqrt(shape)
end if
b1 = shape - 1./non
b2 = shape + 1./non
if (shape.GT.0..and.shape.LE.0.4) then
	c1 = 0.
else if (shape.GT.0.4) then
	c1 = b1*(log(b1)-1.)/2. !log en base 10 que es el log neperiano
end if
c2 = b2*(log(b2)-1.)/2.
600 continue
call random(V1)
call random(V2)
w1 = c1 + log(V1)
w2 = c2 + log(V2)
y = non*(b1*w2 - b2*w1)
if (y.LT.0.) then
	go to 600
end if
x = non*(w2-w1)
if (x.LE.log(y)) then
	ga = scale*exp(x)
	if (UTrotura.EQ.3) then
		ga = ga*365*24
	end if 
else
	go to 600
end if
return
end

!Subrrutina que proporciona un valor aleatorio que sigue una distribucion uniforme entre Umin y Umax
subroutine TroturaU(UTrotura,Umin,Umax,ga)
integer UTrotura
real*8 Umin,Umax,ga
call random(V1)
ga = (Umax-Umin)*V1+Umin
if (UTrotura.EQ.3) then
	ga = ga*365*24
end if 
return
end

!Subrrutina que proporciona un valor aleatorio que sigue una distribucion de Pareto simple
subroutine TroturaP(UTrotura,shapep,Xm,ga)
integer UTrotura
real*8 shapep,Xm,ga
call random(V1)
ga = (Xm)/(V1**(1./shapep))
if (UTrotura.EQ.3) then
	ga = ga*365*24
end if 
return
end

!Subrrutina que proporciona un valor aleatorio que sigue una distribucion de Pareto bounded
subroutine TroturaPB(UTrotura,shapepb,Umin,Umax,ga)
integer UTrotura
real*8 shapepb,Umin,Umax,ga,xx
call random(V1)
ga = (-((V1*Umax**shapepb-V1*Umin**shapepb-Umax**shapepb)/((Umax**shapepb)*(Umin**shapepb))))**(-1./shapepb)
if (UTrotura.EQ.3) then
	ga = ga*365*24
end if 
return
end

function Pn(x)
!INTRODUCIR LA FUNCION DE CHEBYCHEV
!Es una funcion definida entre [0,1] y si la diferencia entre f3 y z se sale
!de esos valores rescalados no hay afeccion del oleaje sobre el cliff o acantilado
!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------
implicit none
real*8 Pn,x,y,eqn7611

y=eqn7611(x)

Pn = (y - abs(y))/2.0

return
end

 
!-----------------------------------------------------------
REAL*8 FUNCTION eqn7611(x)
!-----------------------------------------------------------
! TableCurve C:\Documents and Settings\cparedes\Escritorio\interpolado.f90 Oct 6, 2009 5:22:20 PM 
!  
! X=  
! Y=  
! Eqn# 7611  Chebyshev Rational Order 5/6 
! r2=0.9950469458697139D0 
! r2adj=0.9865559959320805D0 
! StdErr=0.3899952850805259D0 
! Fval=146.1059150680554D0 
! a= -1.069838130676276D0 
! b= -1.441825011828062D0 
! c= -0.01221236978255141D0 
! d= 1.303412347447682D0 
! e= 0.4610744748787188D0 
! f= -0.9727210591773651D0 
! g= 0.07844716807628178D0 
! h= 0.5803220956908709D0 
! i= 0.652841133940385D0 
! j= -0.2877147080827893D0 
! k= -0.08926720346045524D0 
! l= 0.1702945908127132D0 
!-----------------------------------------------------------
  REAL*8 c 
  DIMENSION c(11+1) 
  DATA c(1)/-1.069838130676276D0/ 
  DATA c(2)/-1.441825011828062D0/ 
  DATA c(3)/-0.01221236978255141D0/ 
  DATA c(4)/1.303412347447682D0/ 
  DATA c(5)/0.4610744748787188D0/ 
  DATA c(6)/-0.9727210591773651D0/ 
  DATA c(7)/0.07844716807628178D0/ 
  DATA c(8)/0.5803220956908709D0/ 
  DATA c(9)/0.6528411339403850D0/ 
  DATA c(10)/-0.2877147080827893D0/ 
  DATA c(11)/-0.08926720346045524D0/ 
  DATA c(12)/0.1702945908127132D0/ 
  REAL*8 x,y
  REAL*8 EVALCRATL 
  y=EVALCRATL(11,x,c) 
  eqn7611=y
  RETURN
END
  
!-----------------------------------------------------------
REAL*8 FUNCTION EVALCRATL(order, x, c)
!-----------------------------------------------------------
  REAL*8 x,c(*)
  INTEGER order,i,j,k
  REAL*8 xs,n,d,tmp,t2,t1
  xs=(x-(0.5111604219506374D0))/(0.5110522089614126D0)
  i=order
  j=order
  IF (MOD(order,2).EQ.1) THEN
    i=i-1
  ELSE
    j=j-1
  ENDIF
  t2=0.0
  t1=0.0
  DO k=i+1,3,-2
    tmp=t1
    t1=2*xs*t1-t2+c(k)
    t2=tmp
  END DO
  n=xs*t1-t2+c(1)
  t2=0.0
  t1=0.0
  DO k=j+1,2,-2
    tmp=t1
    t1=2*xs*t1-t2+c(k)
    t2=tmp
  END DO
  d=xs*t1-t2+1.0
  IF (d.EQ.0.0) THEN
    EVALCRATL=0.0
  ELSE
    EVALCRATL=n/d
  ENDIF
  return
END