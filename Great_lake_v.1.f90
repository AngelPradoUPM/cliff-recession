!  Great_lake_v.1.f90 
!
!  FUNCTIONS:
!  Great_lake_v.1 - Bluff clay lakes profiles recesion (LAKE ERIE)
!
!****************************************************************************
!
!  PROGRAM: Great_lake_v.2
!
!  PURPOSE:  Bluff tills erosion over non-tidal environments.
!
!  Programers: R. Castedo                              
!  Email: ricardo.castedo@upm.es           
!                                                                      
!  First version date: 28 / September / 2011		                   
!  Last version date: 04 / November / 2011
!
!****************************************************************************
!
!  VARIABLES-PROGRAM DESCRIPTION: 
!   
!                               Failure from 0 or lakelevel in case of a new
!                               profile "creates" material
!  
!
!****************************************************************************

implicit none
parameter z=500,zz=4,zzz=1000
character*500 basura 
integer BD,Beachpoint,BaseBeach,Berm
integer changelakelevel,ciclo,controllevel,counter,crestatalud,CT
integer DB,DC,Dibujo,DimDiby,dimmatrizHwdraw,dimmatrizrotura,dimmatrizrecesion,Draw,Draw13,Draw14
integer Factor
integer groundwater
integer i,IBEACH,ibusqueda,IDibMatriz,IEcohesion,IEphi,IER,IH,IL,ilabel,IMatrizHwdraw,IMatrizrecesion,IMatrizrotura,indvertcliff,Inter,invierno
integer Ipesounit,Ipotestrato,IPzyRglobal,IRunup,iseed,Isigmac,IShearstrength,IT,ITER,ITRunup,IWF,IWH,IWn,Iyerosionmax
integer j,jdibujo
integer k
integer label(10),lakelevel,Lowerpoint
integer m
integer ndatosresis,NG,nlabel,NPrefino,nrot
integer parametrorotura,Playa,Pmonth,ptocorte,PtocorteCirc2,ptodibujo(1),ptoerosion,ptoinf,ptosup
integer rn,roturanotch
integer salto,seepagelevel
integer Tdib,Tglob,TL,tsim,tp
integer Upperpoint,UTdib,Utsim
integer Vertcirc
integer Waves,WBmax
integer yearjump
integer zerolevel
real*8, parameter :: pi= 3.141592653589793
real*8 a,anchoslide,aux,aux2,aux3[allocatable](:,:),aux4[allocatable](:,:),aux5[allocatable](:,:),aux6[allocatable](:)
real*8 b,BE,BEACH[allocatable](:),Bmat,BX
real*8 Cb,CoordX,curvah,CXc,CYc
real*8 DA,DD,DD1,Deno,DibMatriz[allocatable](:,:),dist,distancia(zzz),DistPcirc2,disty
real*8 Dprot,Drem,dummy[allocatable](:,:),dummy1[allocatable](:,:),dummy2[allocatable](:,:),DWB,Dz
real*8 e,EE,Ecohesion[allocatable](:),Ecohesionrot,Energy,Ephi[allocatable](:),Ephirot,ER[allocatable](:),Ey
real*8 Fb,FS,freshwater,Fw
real*8 gamma,GG,Gs,GR
real*8 H[allocatable](:),HB,HH,hL,hrot,Hslide,Htalud,HWB,Hw,Hwdraw(zzz)
real*8 II
real*8 KRs,Kss
real*8 L[allocatable](:),lambda,Lb,Lcdg,LimWB,LimWBmax,LL,longitud
real*8 MatrizHwdraw[allocatable](:,:),Matrizrecesion[allocatable](:,:),Matrizrotura[allocatable](:,:),mbeach,mcliff,mcirculo,minestabilidad,Mf,Mmedio,MMM,mrecta,munion
real*8 NN,n,No,nslides,Num
real*8 OO
real*8 p,pendientetgPcirc1,pesounit[allocatable](:),pesounitrot,Plat,potacumulada,potestrato[allocatable](:),profagua
real*8 Profseepage,Pzy(3,z),PzyRglobal[allocatable](:,:)
real*8 Qb
real*8 radio,retroceso,Runup[allocatable](:),Rs,Rsc
real*8 Sd,seawater,sigmac[allocatable](:),sigmatensilemax,sigmatensilerock,SFmin,Shearstrength[allocatable](:),SS,Su,sumaacumuladarot,sumaincrementoy,SZ 
real*8 T[allocatable](:),tac,TB,TER[allocatable](:),tic,TOTAL,TRunup[allocatable](:)
real*8 u1,u2,Uo
real*8 vectorverticalrot[allocatable](:),volacumulado,volacumulado1,volacumulado2,VVV
real*8 WBD,WF[allocatable](:),WFmin,WH[allocatable](:),Wn[allocatable](:)
real*8 XDistPcirc2,Xca,Xcb,Xmedio,XPtomedio,XX
real*8 Yc,Yca,Ycb,Ymedio,YDistPcirc2,yerosionmax[allocatable](:),Yopc1,Yopc2,YPtomedio,YY
real*8 Zc

!Creacion de los ficheros de lectura y salida de datos
open(1,file='data.dat',status='unknown')
open(2,file='timerecesion.out',status='unknown')
!open(3,file='marea.dat',status='unknown')
open(4,file='Matriz_dibujo.out',status='unknown')
open(6,file='Matriz_roturas.out',status='unknown')
open(7,file='Matriz_recesion.out',status='unknown')
open(8,file='Matriz_seepage.out',status='unknown')

!Input data
read(1,*) basura
read(1,*) basura !Enter run number
read(1,*) rn
read(1,*) basura !Units of simulation time (0:hours; 1:days; 2:months; 3:years)
read(1,*) Utsim
read(1,*) basura !Simulation time
read(1,*) tsim
read(1,*) basura !Units per result output (0:hours; 1:days; 2:months; 3:years)
read(1,*) UTdib
read(1,*) basura !Result output time
read(1,*) Tdib
read(1,*) basura !Vertical discretization size (m)
read(1,*) Dz
read(1,*) basura !Tidal level (m above/below[in negative] zero) 
read(1,*) lakelevel
read(1,*) basura !Sand amount (m3/m)
read(1,*) Qb
read(1,*) basura !Sand grain diameter (m)
read(1,*) Gs
read(1,*) basura !Groundwater (0:no; 1:yes; 2: fully saturated); Groundwater depth landward (m)
read(1,*) groundwater,profagua
read(1,*) basura !Number of geological layers (from top to bottom)
read(1,*) ndatosresis

!!Dimension of the used variables function of ndatosresis
!allocate (SFmin(ndatosresis),stat=ISFmin)
!if (ISFmin.NE.0) then
!	stop 'Dimension failure: SFmin'
!end if
!allocate (Kss(ndatosresis),stat=IKss)
!if (IKss.NE.0) then
!	stop 'Dimension failure: Kss'
!end if
!allocate (p(ndatosresis),stat=Ip)
!if (Ip.NE.0) then
!	stop 'Dimension failure: p'
!end if
read(1,*) basura 
!do i = 1, ndatosresis
	read(1,*) SFmin,Kss,p !Clay erosion treshold value [0.5-20 Pa]; Coef. for shear stress [Kss]; Calibration coef. [p]
!end do

!!Dimension of the used variables function of ndatosresis
!allocate (MMM(ndatosresis),stat=IMMM)
!if (IMMM.NE.0) then
!	stop 'Dimension failure: MMM'
!end if
!allocate (WFmin(ndatosresis),stat=IWFmin)
!if (IWFmin.NE.0) then
!	stop 'Dimension failure: WFmin'
!end if
read(1,*) basura 
!do i = 1, ndatosresis 
	read(1,*) WFmin,MMM !Threshold surf force [50-500 Kg/m^2]; Coef. for bluff erosion [MMM]
!end do

!Dimension of the used variables function of ndatosresis
!allocate (Rsc(ndatosresis),stat=IRsc)
!if (IRsc.NE.0) then
!	stop 'Dimension failure: Rsc'
!end if
!allocate (Cb(ndatosresis),stat=ICb)
!if (ICb.NE.0) then
!	stop 'Dimension failure: Cb'
!end if
!allocate (Bmat(ndatosresis),stat=IBmat)
!if (IBmat.NE.0) then
!	stop 'Dimension failure: Bmat'
!end if
read(1,*) basura 
!do i = 1, ndatosresis 
	read(1,*) Rsc,Cb,KRs,Bmat !Treshold run-up distance []; Coef. for run-up [Cb]; Coef. for Rs to shear [KRs]; Beach material factor [Bmat]
!end do

!Dimension of the used variables function of ndatosresis
allocate (sigmac(ndatosresis),stat=Isigmac)
if (Isigmac.NE.0) then
	stop 'Dimension failure: sigmac'
end if
allocate (Shearstrength(ndatosresis),stat=IShearstrength)
if (IShearstrength.NE.0) then
	stop 'Dimension failure: Shearstrength'
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
allocate (potestrato(ndatosresis),stat=Ipotestrato)
if (Ipotestrato.NE.0) then
	stop 'Error de dimensionamiento de potestrato'
end if
read(1,*) basura !Clay cliff erosion resistance-UCS [100-550 KPa]; Effective cohesion(KPa); Effective friction angle (º); Unit weight(KN/m**3); Strata depth (m)
do i = 1, ndatosresis
	read(1,*) sigmac(i),Shearstrength(i),Ecohesion(i),Ephi(i),pesounit(i),potestrato(i)
end do

read(1,*) basura !Debris Mobility Coeff. (m**2/day)
read(1,*) Mf

write(*,*)'Great_lake_v.1 - Bluff clay lakes profiles recesion'

!Bluff initial profile (m) [1:Study point; 2:Talus wedge; 3:Seepage point]
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
	case(1)
		!Simulation time: month
		tsim=tsim
	case(2)
		!Simulation time: years*months(12)
		tsim=tsim*12
end select

select case(UTdib)
	case(1)
		!Result output time: month
		Tdib=Tdib
	case(2)
		!Result output time: years*months(12)
		Tdib=Tdib*12
end select

!Orientation factor
if (Pzy(2,1).LT.Pzy(2,2)) then
	Factor = 1.0
else if (Pzy(2,1).GT.Pzy(2,2))then
	Factor = -1.0
end if

!Number of geometry discretization points
NPrefino = int(abs(Pzy(1,1)-Pzy(1,NG))/Dz)

!Creation of PzyRglobal -> profile points possitions matrix
allocate (PzyRglobal(2,NPrefino),stat=IPzyRglobal)
if (IPzyRglobal.NE.0) then
	stop 'Dimension failure: PzyRglobal'
end if

!First point coordenates
PzyRglobal(1,1) = Pzy(1,1)
PzyRglobal(2,1) = Pzy(2,1)
!Last point coordenates
PzyRglobal(1,NPrefino) = Pzy(1,NG)
PzyRglobal(2,NPrefino) = Pzy(2,NG)

!Coordenates matrix PzyRglobal for initial profile
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
			go to 21
		end if
	end do
	21 continue		
end do

!Localization of 1:Study point; 2:Talus wedge; 3:Seepage point; 4:Talus base point => roturanotch
ilabel=0
do ibusqueda=2,NG
    ! 1:Study point
    if(Pzy(3,ibusqueda).EQ.1.)then
       do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     ilabel=ilabel+1
			 label(ilabel)=k
			 goto 22
		  end if
	   end do
    ! 2:Talus wedge
	else if(Pzy(3,ibusqueda).EQ.2.)then
	   do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     crestatalud=k
			 goto 22
		  end if
	   end do
    ! 3:Seepage point
	else if(Pzy(3,ibusqueda).EQ.3.)then
	   do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     seepagelevel=k
			 goto 22
		  end if
	   end do
	! 4:Talus base point => roturanotch
	else if(Pzy(3,ibusqueda).EQ.4.)then
	   do k=1,NPrefino
	      if(PzyRglobal(1,k).LT.Pzy(1,ibusqueda)) then
		     roturanotch=k-1
			 goto 22
		  end if
	   end do
	end if
	22 continue
end do

!Profile location: Zero value
do j=1,NPrefino
	if (PzyRglobal(1,j).LE.0.d0) then
		goto 23
	end if	
end do
23 continue
zerolevel = j-1 !zero or first negative value

!Profile location: Tidal level (m above/below zero) 
do j=1,NPrefino
	if (PzyRglobal(1,j).LE.(PzyRglobal(1,zerolevel)+lakelevel)) then
		goto 24
	end if	
end do
24 continue
TL = j !TL or first value below it
controllevel = TL

!Drawing counter
j = 1
jdibujo = 1
!Number of drawings and matrix drawing dimension
DimDiby = int(tsim/Tdib)

!Creation of DibMatriz -> Erosive profiles matrix
allocate (DibMatriz(NPrefino+1,DimDiby+2),stat=IDibMatriz)
if (IDibMatriz.NE.0) then
	stop 'Dimension failure: DibMatriz'
end if

!Creation of Matrizrotura -> Failure profiles matrix
allocate (Matrizrotura(NPrefino+1,2),stat=IMatrizrotura)
if (IMatrizrotura.NE.0) then
	stop 'Dimension failure: Matrizrotura'
end if

!Creation of Matrizrecesion -> Diference in possitions before-after failure
allocate (Matrizrecesion(NPrefino+2,2),stat=IMatrizrecesion)
if (IMatrizrecesion.NE.0) then
	stop 'Dimension failure: Matrizrecesion'
end if

!Erosive profiles matrix, first two columns
DibMatriz(1,1) = rn !Run number
DibMatriz(1,2) = 0. !
do i=1,NPrefino
	DibMatriz(i+1,1) = PzyRGlobal(1,i)
	DibMatriz(i+1,2) = PzyRGlobal(2,i)
end do

!Failure profiles matrix, first two columns
Matrizrotura(1,1) = rn !Run number
Matrizrotura(1,2) = 0. !
do i=1,NPrefino
	Matrizrotura(i+1,1) = PzyRGlobal(1,i)
	Matrizrotura(i+1,2) = PzyRGlobal(2,i)
end do

!Recession values, first two columns
Matrizrecesion(1,1) = 1. !run number
Matrizrecesion(1,2) = 0. !
Matrizrecesion(2,1) = 0. !Area erosionada originalmente
Matrizrecesion(2,2) = 0. !
do i=1,NPrefino
	Matrizrecesion(i+2,1) = PzyRGlobal(1,i)
	Matrizrecesion(i+2,2) = PzyRGlobal(2,i)
end do

!Specific weight of fresh water for landslide stability evaluation - KN/m**3
freshwater=9.81

!Density of seawater for Surf stress (SF) evaluation - Kg/m**3
seawater=1030

!Initialization of variables
volacumulado=0.d0
volacumulado1=0.d0
volacumulado2=0.d0
changelakelevel=0 !Interannual change in lake levels -> 1: YES; 2: NO until equilibrium conditions
ciclo=0 !Number of months for a complete cicle in lake level change
Dprot=0.d0
Pmonth=0
Drem=0.d0
parametrorotura=0
retroceso=0.5d0
salto=1
counter=1
Draw13=0
Draw14=0
invierno=1 !1: winter YES; 2:winter NO
pendientetgPcirc1=dtan(1.*pi/180.) !Possitive because is going up
!Vertical lines to calculate the intersection between the circular failure and the real profile
mcirculo=0.d0
nslides=5.d0
    
!Groundwater calculus to print it
if (groundwater.EQ.0) then
    do i=1,1000
        distancia(i) = 0.d0
        Hwdraw(i) = 0.d0
    end do
    go to 64
else if (groundwater.EQ.1) then
    Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
    Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
    dist = longitud/1000.
    do i=1,1000
        distancia(i) = PzyRglobal(2,seepagelevel) - (i-1)*dist
    end do
    do i=1,1000
        Hwdraw(i) = dsqrt((dabs(distancia(i)-PzyRglobal(2,seepagelevel))/longitud)*(Profseepage-Profagua)**2) + PzyRglobal(1,seepagelevel)
    end do
end if
64 continue

!Creation of Hwdraw matrix -> Groundwater level draw
allocate (MatrizHwdraw(1000+1,2),stat=IMatrizHwdraw)
if (IMatrizHwdraw.NE.0) then
    stop 'Dimension failure: MatrizHwdraw'
end if

!Groundwater levels matrix, first two columns
MatrizHwdraw(1,1) = rn !Run number
MatrizHwdraw(1,2) = 0.d0
do i=1,1000
    MatrizHwdraw(i+1,1) = Hwdraw(i)
	MatrizHwdraw(i+1,2) = distancia(i)
end do

!WAVES CONDITIONS WINTER
Waves = 1

!Creation of H vector -> Deep water wave height
allocate (H(Waves),stat=IH)
if (IH.NE.0) then
    stop 'Dimension failure: Deep water wave height (H)'
end if
!Introduction of the Deep water wave height values
H(1) = 0.d0

!Creation of T vector -> Wave period
allocate (T(Waves),stat=IT)
if (IT.NE.0) then
    stop 'Dimension failure: Wave period (T)'
end if
!Introduction of the Wave period values
T(1) = 0.d0

!Creation of L vector -> Deep water wave length
allocate (L(Waves),stat=IL)
if (IL.NE.0) then
    stop 'Dimension failure: Deep water wave length (L)'
end if
!Calculation of the Deep water wave length values
L(1) = 1.56*(T(1)**2) !Units: [m] - CERC (2002)

!Creation of WF vector -> Percentage frecuency of each type of wave per day
allocate (WF(Waves),stat=IWF)
if (IWF.NE.0) then
    stop 'Dimension failure: Percentage frecuency of each type of wave (WF)'
end if
!Introduction of the Percentage frecuency of each type of wave 
WF(1) = 0.d0

!Creation of Wn vector -> The number of waves each hour
allocate (Wn(Waves),stat=IWn)
if (IWn.NE.0) then
    stop 'Dimension failure: The number of waves each hour (Wn)'
end if
!Calculation of the The number of waves each hour 
Wn(1) = (60*60)/T(1) !60 s/min * 60 min/hours

call cpu_time(tic)
call seed(iseed)

!INIZILIZATION OF THE TGLOB VARIABLE
Tglob = 1

!Erosion - Recession DO-WHILE LOOP that allow to modify the Tglob variable
500 IF (Tglob.LE.tsim) THEN

    !Write in commands window time steep (days)
    write(*,'(''Time steep: '',i8)')Tglob

    if (invierno.EQ.1) then
        salto = salto + counter
        if (Tglob.EQ.(3*(salto-2)+1)) then
            Tglob = Tglob+2 !From left side to right one
            Draw13 = 1 !Drawn from left side
            Draw14 = 0 !Drawn from center
        else if (((3*(salto-2)+1).LT.Tglob).and.(Tglob.LT.3*(salto-1))) then
            Tglob = Tglob+1 !From center to right side
            Draw13 = 0 !Drawn from left side
            Draw14 = 1 !Drawn from center
        else
            Draw13 = 0 !Drawn from left side
            Draw14 = 0 !Drawn from center
        end if
        
        !Counter re-inizialization
        if (counter.EQ.4) then
            counter=1
        end if
       
        !SPRING waves conditions
        Waves = 4 !Number of waves in the season SPRING
        !Creation of H vector -> Deep water wave height
        deallocate (H)
        allocate (H(Waves))
        !Creation of T vector -> Wave period
        deallocate (T)
        allocate (T(Waves))
        !Creation of L vector -> Deep water wave length
        deallocate (L)
        allocate (L(Waves))
        !Creation of WF vector -> Percentage frecuency of each type of wave per day
        deallocate (WF)
        allocate (WF(Waves))
        !Creation of Wn vector -> The number of waves each hour
        deallocate (Wn)
        allocate (Wn(Waves))
        
        !Introduction of the Deep water wave height values
        H(1) = 0.5
        H(2) = 1
        H(3) = 2
        H(4) = 3
        !Introduction of the Wave period values
        T(1) = 2.
        T(2) = 3.5
        T(3) = 4.5
        T(4) = 5.5
        !Calculation of the Deep water wave length values
        do i=1,Waves
            L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
        end do
        !Introduction of the Percentage frecuency of each type of wave 
        WF(1) = 33.03
        WF(2) = 36.33
        WF(3) = 2.23
        WF(4) = 0.1
        !Total frecuency = 71.67 %
        !Calculation of the The number of waves each hour 
        do i=1,Waves
            Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
        end do
        
        !SPRING lake level
        !TL = controllevel - int(0.4/Dz) !Increase 0.4m is the value obtained from ARMY CORPS OF ENGINEERS - DETROIT District (2011)
        
        if (Tglob.GE.tsim) then
            if (Draw13.EQ.1) then
                Draw = Tglob-2
                do j=Draw,tsim
                    !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                    if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                    !Time vectors to introduce in the matrix
		                DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	                !Transfer the position information matrix into DRAWING MATRIX
		                do i=1,NPrefino
    	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                    end do
   		                jdibujo = jdibujo + 1
                    end if
                end do
                go to 1000 !End of simulation
            else if (Draw14.EQ.1) then
                Draw = Tglob-1
                do j=Draw,tsim
                    !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                    if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                    !Time vectors to introduce in the matrix
		                DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	                !Transfer the position information matrix into DRAWING MATRIX
		                do i=1,NPrefino
    	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                    end do
   		                jdibujo = jdibujo + 1
                    end if
                end do
                go to 1000 !End of simulation
            else
                !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                if ((Tglob.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                !Time vectors to introduce in the matrix
		            DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	            !Transfer the position information matrix into DRAWING MATRIX
		            do i=1,NPrefino
    	                DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                end do
   		            jdibujo = jdibujo + 1
                end if
                go to 1000 !End of simulation
            end if
        else if (Tglob.LT.tsim) then
            if (Draw13.EQ.1) then
                Draw = Tglob-2
                do j=Draw,Tglob
                    !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                    if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                    !Time vectors to introduce in the matrix
		                DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	                !Transfer the position information matrix into DRAWING MATRIX
		                do i=1,NPrefino
    	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                    end do
   		                jdibujo = jdibujo + 1
                    end if
                end do
                go to 1001 !Loops end before -> invierno =0
            else if (Draw14.EQ.1) then
                Draw = Tglob-1
                do j=Draw,Tglob
                    !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                    if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                    !Time vectors to introduce in the matrix
		                DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	                !Transfer the position information matrix into DRAWING MATRIX
		                do i=1,NPrefino
    	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                    end do
   		                jdibujo = jdibujo + 1
                    end if
                end do
                go to 1001 !Loops end before -> invierno =0
            else
                !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                if ((Tglob.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                !Time vectors to introduce in the matrix
		            DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	            !Transfer the position information matrix into DRAWING MATRIX
		            do i=1,NPrefino
    	                DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                end do
   		            jdibujo = jdibujo + 1
                end if
                go to 1001 !Loops end before -> invierno =0
            end if
        end if
    end if
    
    !Allocate BEACH from 1 to Number of waves
    allocate (BEACH(Waves),stat=IBEACH)
    if (IBEACH.NE.0) then
        stop 'Dimension failure: Beach angle'
    end if
    do i = 1,Waves
        BEACH(i) = 0.d0
    end do
    !Allocate TER to NPrefino
    allocate (TER(NPrefino),stat=ITER)
    if (ITER.NE.0) then
        stop 'Dimension failure: TER or Total Shear erosion'
    end if
    !Allocate TRunup to NPrefino
    allocate (TRunup(NPrefino),stat=ITRunup)
    if (ITRunup.NE.0) then
        stop 'Dimension failure: TRunup or Total Run-up erosion'
    end if
    allocate (aux5(2,NPrefino))
    do j = 1,NPrefino
        TER(j) = 0.d0
        TRunup(j) = 0.d0
        aux5(2,j) = PzyRglobal(2,j)   
    end do
    
    !Variables initialization
    Ey = 0.d0
    
    !Erosion module Trenhaile (2009) - Coastal Engineering
    do i = 1,Waves !From wave 1 to Number of waves within the season
    
        !Calculate HEIGHT OF SHOALING WAVES
        DA = L(i)/2 !Depth = 0.5 * wave length - CERC (2002)
        !Integer part of the depth converted to DZ size steep 
        DB = int(DA/Dz)
        !DC is the vertical k value for the onset of shallow water (very deep)
        DC = TL + DB
     
        !Allocate WH to NPrefino
        allocate (WH(NPrefino),stat=IWH)
        if (IWH.NE.0) then
            stop 'Dimension failure: WH'
        end if
        do j = 1,NPrefino
            WH(j) = 0.d0
        end do
        
        do j = DC,TL,-1
            !Calculate n for group wave velocity
            hL = (j-TL)*Dz/L(i) !h/L0 -> h:water depth
            if (hL.GT.0.25d0) then
                LL = L(i) !wavelength in deep water
                WH(j) = H(i) !wave height in deep water
                go to 62
            else if ((hL.LE.0.25d0).AND.(hL.GT.0.05d0)) then
                DD = (dexp(2.d0*pi*hL)-dexp(-2.d0*pi*hL))/(dexp(2.d0*pi*hL)+dexp(-2.d0*pi*hL))
                LL = L(i)*dsqrt(DD) !wavelength in intermediate depths is determined using the CEM (2002) equation
                go to  61
            else if (hL.LE.0.05d0) then
                LL = T(i)*dsqrt(9.81d0*(j-TL)*Dz) !wavelength in shallow depths is determined using the CEM (2002) equation
                go to 61
            end if
            61 continue
            OO = 2.d0*2.d0*pi*(j-TL)*Dz/LL !2*k*h -> k = 2*pi/L
            NN = (dexp(OO)-dexp(-OO))/2.d0
            n = 0.5d0*(1.d0+(OO/NN)) !n for group wave velocity
            
            !Calculate wave height at and outside the breakers and breaker depth
            if ((hL.LE.0.25d0).AND.(hL.GT.0.05d0)) then
                DD1 = (dexp(2.d0*3.14*(j-TL)*Dz/LL)-dexp(-2.d0*3.14*(j-TL)*Dz/LL))/(dexp(2.d0*3.14*(j-TL)*Dz/LL)+dexp(-2.d0*3.14*(j-TL)*Dz/LL))
                WH(j) = H(i)*dsqrt(L(i)/(2.d0*T(i)*n*(9.81d0*T(i)*DD1/(2.d0*3.14)))) !H=H0*sqrt(C0/2*C*n) -> C0=L0/T ; C=(g*T/2*pi)*tanh(2*pi*h/L)
            else if (hL.LE.0.05) then
                WH(j) = H(i)*dsqrt(L(i)/(2.d0*T(i)*n*dsqrt(9.81d0*(j-TL)*Dz))) !H=H0*sqrt(C0/2*C*n) -> C0=L0/T ; C=sqrt(g*h)
            end if
            62 continue   
            if (WH(j).GE.(j-TL)*Dz/1.28d0) then !Hb >= h/1.28 -> breaking exists
                WBD = (int(WH(j)*1.28d0*10.d0))/10.d0 !WBD is actual wave breaker depth (metres) -> hb = 1.28*Hb
                HB = WH(j) !mean breaker height
                go to 63
            end if
        end do
        
        63 continue
        !Calcuate the LENGTH OF THE SURF ZONE (SZ)
        BD = int(WBD/Dz) + TL !Corresponding i-value for WBD
        SZ = PzyRglobal(2,BD) - PzyRglobal(2,TL) !length
        
        !Calculate WAVE HEIGHT IN SURF ZONE
        do j = BD-1,TL,-1
            
            GR = 0.d0
            if (SZ.EQ.0.d0) then
                go to 58        
            end if
            !Calculate SLOPE OF THE SURF ZONE
            u1 = PzyRglobal(2,BD) - PzyRglobal(2,j) !Horizontal length from BD to j
		    u2 = PzyRglobal(1,BD) - PzyRglobal(1,j) !Vertical distance from BD to j
		    GR = dabs(u2/u1)
            
            58 continue
            if (PzyRglobal(2,j).EQ.PzyRglobal(2,BD)) then
                WH(j) = HB
                go to 60
            else
                WH(j) = HB*dexp((-35.2d0*(dabs(PzyRglobal(2,BD)-PzyRglobal(2,j)))*GR)/(T(i)*dsqrt(9.81d0*WBD)))
            end if
            60 continue
        end do
        
        !Determine BEACH FACE GRADIENT (BEACH(i))(RADIANS)
        VVV =  HB/((9.81**0.5)*(Gs**0.5)*T(i)) !Sunamura´s expression (1989)
        BEACH(i) = datan(0.12/(VVV**0.5))
        
        !Determine BEACH SHAPE AND VOLUME for each wave type
        !Profile location: Tidal level + Hb -> Maximun height to find a beach berm
        LimWB = nint((PzyRglobal(1,TL)+HB)*10.d0)
        !LimWB = nint((PzyRglobal(1,TL))*10.d0)
        LimWBmax = LimWB/10.d0
        do j=1,BD
            aux2 = nint(PzyRglobal(1,j)*10.d0)
            aux2 =aux2/10.d0
	        if (aux2.LE.LimWBmax) then
	        !if (PzyRglobal(1,j).LE.(PzyRglobal(1,TL)+HB)) then
		        go to 59
	        end if	
        end do
        59 continue
        WBmax = j !WBmax or Berm Width maximun
        !LimWB = nint((PzyRglobal(1,TL)+HB)*10.d0)
        
        Beachpoint = 0
        BaseBeach = 0
        Berm = 0
        
        if (Qb.EQ.0.d0) then !No volume for beach
            Beachpoint = 0
            BaseBeach = 0
            Playa = 2 !No beach exist
            go to 30
        else if (Qb.NE.0.d0) then    
            !Determine location and volume
            do j=BD,TL,-1 !!!!!!!!!!!!!!!!!!!!ojo es TL
                if (PzyRglobal(2,j).EQ.PzyRglobal(2,j-1)) then !vertical platform
                    Plat = pi/2.d0 ! pi/2 in Rad -> 90 degrees
                    Beachpoint = 0
                    BaseBeach = 0
                    Playa = 2 !No beach exist
                    go to 25 !In the case the bottom point can varies from BD to TL
                    !go to 30 !Beach cannot create from BD
                else
                    Plat = datan(Dz/(PzyRglobal(2,j)-PzyRglobal(2,j-1))) ! Slope of the shore platform
                end if
                
                !Determine height HWB from j to WBmax
                HWB = dabs(PzyRglobal(1,WBmax)-PzyRglobal(1,j))
                !Determine horizontal distance DWB from j to WBmax
                DWB = dabs(PzyRglobal(2,WBmax)-PzyRglobal(2,j))
                
                if (BEACH(i).LE.Plat) then
                    Beachpoint = 0
                    BaseBeach = 0
                    Playa = 2 !No beach exist
                    go to 25 !In the case the bottom point can varies from BD to TL
                    !go to 30 !Beach cannot create from BD
                end if
                
                Inter = 2
                !Determine intersection point
                do 26 k=j-1,WBmax+1,-1
                    !Determine director´s vector for each point
			        u1 = PzyRglobal(2,k-1) - PzyRglobal(2,k)
			        u2 = PzyRglobal(1,k-1) - PzyRglobal(1,k)
				    !Profile slope
			        if (u1.EQ.0.d0) then
			    	    mcliff = 0.d0
			        else
			    	    mcliff = u2/u1
			        end if
				    !Cross points between BEACH and CLIFF
				    mbeach = - datan(BEACH(i))
			   	    Yc = ((PzyRglobal(1,k)-mcliff*PzyRglobal(2,k)) - (PzyRglobal(1,j) - (mbeach)*PzyRglobal(2,j)))/((mbeach)-mcliff)
					Zc = (mbeach)*Yc + PzyRglobal(1,j) - (mbeach)*PzyRglobal(2,j)
				    !Lambda paratemer that determines the presence or not of the cross point
			    	if (mcliff.EQ.0.d0) then
			    	    lambda = (Zc - PzyRglobal(1,k))/u2
			    	else
			    	    lambda = (Yc - PzyRglobal(2,k))/u1
			    	end if
				    !Coordenate of cross point if it exists
			    	if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
			   	        Inter = k !i-position of intersection point
			       	    go to 27
			    	else if (lambda.EQ.1.0) then
			    	    Inter = k-1
			    		go to 27
			    	else
			    	    go to 26
			    	end if
			    26 continue
                27 continue
                                
                if (((BEACH(i)*DWB).LT.HWB).or.(PzyRglobal(1,Inter).LT.PzyRglobal(1,WBmax))) then 
                    !Point below WBmax
                    BE = PzyRglobal(2,Inter) !X-coordenate of beachpoint
                    HWB = dabs(PzyRglobal(1,Inter)-PzyRglobal(1,j))
                    DWB = dabs(PzyRglobal(2,Inter)-PzyRglobal(2,j))
                    Beachpoint = Inter
                else if (((BEACH(i)*DWB).EQ.HWB).or.(PzyRglobal(1,Inter).EQ.PzyRglobal(1,WBmax))) then !Point equal to WBmax
                    Beachpoint = WBmax !i-position of intersection point
                    BE = PzyRglobal(2,WBmax) !X-coordenate of beachpoint
                else if (((BEACH(i)*DWB).GT.HWB).or.(PzyRglobal(1,Inter).GT.PzyRglobal(1,WBmax))) then !Point above WBmax -> beach berm and NO erosion
                    Beachpoint = WBmax !for vertical position
                    BE = PzyRglobal(2,j) - (HWB/dtan(BEACH(i))) !X-coordenate of beachpoint
                    Berm = 1 !No impact erosion
!                    Beachpoint = 0
!                    BaseBeach = 0
!                    Playa = 2
!                    Berm = 0
!                    go to 25
                end if
                
                !VOLUME determination
                a = 0.d0
                b = 0.d0
                e = 0.d0
                BX = 0.d0
                a = DWB*HWB !Rectangular maximun element
                b = 0.5d0*HWB*(PzyRglobal(2,j)-BE) !Triangle above profile
                e = 0.5d0*Dz*dabs(PzyRglobal(2,j-1)-PzyRglobal(2,j)) !First profile above slope
                TB = 0
                CT = 0
                M = 0
                do k = j-1,Beachpoint+1,-1
                    if (k.LT.(M-1))then
                        go to 28
                    else if (k.EQ.(M-1)) then
                        M = 0
                    else if (PzyRglobal(2,k).EQ.PzyRglobal(2,k-1)) then
                        !Calculate Position (CT) of top of a vertical section
                        do M = k+1,Beachpoint,-1
                            if (PzyRglobal(2,M).NE.PzyRglobal(2,k)) then
                                go to 29
                            else
                                CT = M !If it´s vertical no volumen in available
                            end if
                        end do
                        go to 28
                    end if
                    29 continue
                    BX = 0.5d0*Dz*dabs(PzyRglobal(2,k)-PzyRglobal(2,k-1)) + (j-k)*Dz*dabs(PzyRglobal(2,k)-PzyRglobal(2,k-1)) !Below the profile -> triangle + rectangle 
                    TB = TB + BX !Total volume below the profile
                    28 continue
                end do
                
                !TOTAL beach volume
                TOTAL = a - (b + e + TB)
                BaseBeach = j
                Playa = 1
                
                if (Qb.GE.TOTAL) then
                    go to 30 !Beach is ok
                else if (Qb.LT.TOTAL) then
                    Berm = 0
                    go to 25 !We need to find a beach less deeper
                end if
                25 continue
            end do
        end if
        30 continue
        
        if ((Basebeach.EQ.Beachpoint+1)) then
            !In case extremelly small beach is formed
            Playa = 2
            Berm = 0
        end if

        !NUMBER OF EACH TYPE OF WAVE at the tidal level per day - Trenhaile (2009) Fig.3
        !Non-tidal body of water, tidal duration just for one vertical value (TL) and is equal to the number of hours in a month
        No = Wn(i)*(WF(i)/100.d0)*24.d0*(30.d0-Dprot)
        
        !Allocate ER to NPrefino
        allocate (ER(NPrefino),stat=IER)
        if (IER.NE.0) then
            stop 'Dimension failure: ER or Shear erosion'
        end if
        !Allocate Runup to NPrefino
        allocate (Runup(NPrefino),stat=IRunup)
        if (IRunup.NE.0) then
            stop 'Dimension failure: Runup or Run-up erosion'
        end if
        do j = 1,NPrefino
            ER(j) = 0.d0
            Runup(j) = 0.d0
        end do
        
        !Determine SHEAR STRESS EROSION
        do j=DC,TL+1,-1
            Su = 0.d0 !Wave set-up in the surf zone Su
            Sd = 0.d0 !Wave set-down at and outside the breakers
            hL = (j-TL)*Dz/L(i) !h/L0 -> h:water depth
            
            if (hL.GT.0.25d0) then !If j is deeper than hL > 0.25 -> NO erosion
                ER(j) = 0.d0 !Erosion by shear
                go to 31
            else if ((hL.LE.0.25d0).AND.(hL.GT.0.05d0)) then !If j is between 0.25 > hL > 0.05 -> Maximum bottom orbital velocity Uo
                if (j.GT.BD) then
                    XX = 2*3.1416d0*(j-TL)*Dz/L(i)
                    YY = (dexp(XX)-dexp(-XX))/(dexp(XX)+dexp(-XX))
                    Lb = L(i)*dsqrt(YY) !Wavelength
                    HH = 2.d0*2.d0*3.1416d0*(j-TL)*Dz/Lb
                    II = (dexp(HH)-dexp(-HH))/2.d0
                    Sd = 3.1416d0*(WH(j)**2)/(4*Lb*II) !Set-down
                    GG = 2*3.1416d0*((j-TL)*Dz-Sd)/Lb
                    EE = (dexp(GG)-dexp(-GG))/2.d0
                    Uo = 3.1416d0*WH(j)/(T(i)*EE) !Maximum bottom orbital velocity
                end if
            else if (hL.LE.0.05d0) then !If j is between hL < 0.05 -> Maximum bottom orbital velocity Uo
                Lb = T(i)*dsqrt(9.81d0*(j-TL)*Dz) !Wavelength
                HH = 2.d0*2.d0*3.1416d0*(j-TL)*Dz/Lb
                II = (dexp(HH)-dexp(-HH))/2.d0
                Sd = 3.1416d0*(WH(j)**2)/(4*Lb*II) !Set-down
                Uo = (WH(j)/2.d0)*dsqrt(9.81d0/((j-TL)*Dz-Sd)) !Maximum bottom orbital velocity
            end if
            
            if ((j.LT.BD).and.(PzyRglobal(2,j).LT.PzyRglobal(2,BD))) then
                Su = 0.185767d0*(WBD/SZ)*(PzyRglobal(2,BD)-PzyRglobal(2,j)) !Set-up
                Uo = (WH(j)/2.d0)*dsqrt(9.81d0/((j-TL)*Dz+Su)) !Maximum bottom orbital velocity
            else if ((j.LT.BD).and.(PzyRglobal(2,j).EQ.PzyRglobal(2,BD))) then
                Su = 0.d0
                Uo = (WH(j)/2.d0)*dsqrt(9.81d0/((j-TL)*Dz)) !Maximum bottom orbital velocity
            end if
            
            if (Playa.EQ.1) then !Beach present
                if (TL.LT.Beachpoint) then !Lake level above beach
                    if ((j.LE.DC).and.(j.GE.BaseBeach)) then
                        Fw = 0.0034
                        SS = (1/2.d0)*Fw*seawater*(Uo**2)
                        ER(j) = No*Kss*((SS-SFmin)**p)
                    else if ((j.LT.Beachpoint).and.(j.GT.TL)) then
                        Fw = 0.0034
                        SS = (1/2.d0)*Fw*seawater*(Uo**2)
                        ER(j) = No*Kss*((SS-SFmin)**p)
                    end if
                else if (TL.GE.Beachpoint) then !Lake level within the beach or under
                    if ((j.LE.DC).and.(j.GE.BaseBeach)) then
                        Fw = 0.0034
                        SS = (1/2.d0)*Fw*seawater*(Uo**2)
                        ER(j) = No*Kss*((SS-SFmin)**p)
                    end if
                end if
            else if (Playa.EQ.2) then !No beach is present
                if (j.GT.TL) then !Security => loop just until TL-1
                    Fw = 0.0034
                    SS = (1/2.d0)*Fw*seawater*(Uo**2)
                    ER(j) = No*Kss*((SS-SFmin)**p)
                end if
            end if
            31 continue
        end do
        
        !Determine EROSION FOR WAVE IMPACT + RUN-UP
        if (Playa.EQ.1) then !Beach present
            aux2 = nint(PzyRglobal(1,Beachpoint)*10.d0)
            aux2 =aux2/10.d0
            !if (PzyRglobal(1,Beachpoint).GE.(HB+PzyRglobal(1,TL))) then !No run-up
            if (aux2.GE.(LimWBmax)) then !No run-up
                Runup(Beachpoint) = 0.d0
            !else if (PzyRglobal(1,Beachpoint).LT.(HB+PzyRglobal(1,TL))) then !Run-up exists in Beachpoint
            else if (aux2.LT.(LimWBmax)) then !Run-up exists in Beachpoint
                !Rs = ((PzyRglobal(1,TL)+HB)-PzyRglobal(1,Beachpoint))/dtan(Beach(i))
                Rs = ((LimWBmax)-PzyRglobal(1,Beachpoint))/dtan(Beach(i))
                Runup(Beachpoint) = No*Cb*Bmat*(KRs*Rs-Rsc)
            end if
            if ((TL.GT.Beachpoint).and.(TL.LT.Basebeach)) then !No impact erosion
                Ey = Ey + 0.d0
                go to 32
            else if (TL.LE.Beachpoint) then
                Runup(Beachpoint) = 0.d0
                if (berm.EQ.1) then
                    Ey = Ey + 0.d0 !No impact erosion due to berm size
                    go to 32
                end if
            end if
        else if (Playa.EQ.2) then !No beach is present
!            if ((TL.GT.roturanotch).and.(PzyRglobal(1,roturanotch).LT.(PzyRglobal(1,TL)+HB))) then !Below roturanotch and below TL+HB
!                XX = dabs((PzyRglobal(1,roturanotch)-PzyRglobal(1,TL))/(PzyRglobal(2,roturanotch)- PzyRglobal(2,TL)))
!                Rs = ((PzyRglobal(1,TL)+HB)-PzyRglobal(1,roturanotch))/XX
!                Runup(roturanotch) = No*Cb*Bmat*(KRs*Rs-Rsc)
!            else if (TL.LE.roturanotch) then !No Run-up is possible
                Runup(roturanotch) = 0.d0
!            end if
        end if
        if (SZ.EQ.0.d0) then
            Fb = 512.61*(HB/0.78)
        else
            Fb = 512.61*(HB/0.78)*dexp(-SZ*0.05)
        end if
        Energy = 0.d0
        Energy = No*MMM*(Fb-WFmin) !Erosion value for impact for each wave
        if (Energy.LT.0.d0) then
            Energy = 0.d0
        end if
        Ey = Ey + Energy !Accumulative value for impact erosion for all waves
        32 continue
        
        !Determine CUMULATIVE VALUES FOR ALL KINDS OF EROSION
        do j=1,NPrefino
            if (ER(j).LE.0.d0) then
                ER(j) = 0.d0
            end if
            if (Runup(j).LE.0.d0) then
                Runup(j) = 0.d0
            end if
            TER(j) = TER(j)+ER(j) !TER is the sum of all waves type
            TRunup(j) = TRunup(j)+Runup(j) !TRunup is the sum of all run-up possibilities
        end do
        
        !Deallocate WH for the next wave conditions
        deallocate(WH)
        !Deallocate ER for the next wave conditions
        deallocate(ER)
        !Deallocate Runup for the next wave conditions
        deallocate(Runup)  
    end do
    
    !Determine NEW VECTOR POSSITIONS after EROSIVE processes
    do j=1,NPrefino
        PzyRglobal(2,j) = PzyRglobal(2,j)-TER(j)-TRunup(j) !New possitions after TER and TRunup
    end do
    PzyRglobal(2,TL) = PzyRglobal(2,TL)-Ey !New possitions after Ey
        
    !Deallocate TER for the next erosive steep
    deallocate(TER)
    !Deallocate TRunup for the next erosive steep
    deallocate(TRunup)
    !Deallocate BEACH for the next erosive time
    deallocate(BEACH)
   
    !Determination of CRITICAL POINTS
    ptosup = crestatalud
    ptoinf = TL
    do j=crestatalud,TL-1 !Calculate the first inflection point going down
        if (PzyRglobal(2,j+1).GE.PzyRglobal(2,j)) then
            ptosup = j
        else
            ptosup = j+1 !Upper inflection point
            go to 33
        end if
    end do
    if (ptosup.EQ.(TL)) then !No inflextion point exist above lake level
        go to 34 
    end if
    33 continue

    do j=TL,crestatalud+1,-1 !Calculate the first inflection point going up
        if (PzyRglobal(2,j-1).LE.PzyRglobal(2,j)) then
            ptoinf = j-1
        else
            ptoinf = j !Lower inflection point
            go to 35
        end if
    end do
    if (ptoinf.EQ.crestatalud) then !No inflextion point exist above lake level
        go to 34 
    end if
    35 continue
    
    34 continue
    if (ptosup.EQ.ptoinf) then
        Upperpoint = ptoinf !When both points are the same
        if (ptosup.EQ.TL) then
            Upperpoint = roturanotch
            if(PzyRglobal(2,roturanotch).GT.PzyRglobal(2,ptosup)) then
                Upperpoint = ptoinf
            end if
        end if
    else if (ptosup.NE.ptoinf) then
        !Erosive points
	    allocate (yerosionmax(ptoinf-ptosup+1),stat=Iyerosionmax)
	    if (Iyerosionmax.NE.0) then
	        stop 'Dimension failure: yerosionmax'
	    end if
	    !Most erosionated point between ptoinf-ptosup
	    do i=ptosup,ptoinf
		    yerosionmax(i-ptosup+1)=PzyRglobal(2,i)
		end do
		ptodibujo = minloc(yerosionmax)
		Upperpoint = ptodibujo(1) + ptosup - 1 !i-possition of the most erosionated point between ptoinf-ptosup
	    deallocate (yerosionmax)
    end if
    
    !Determination of the Lower critical point below lake level    
    Lowerpoint = DC
    do j=DC,TL+1,-1 !Calculate the inflection point going up under water level
        if (PzyRglobal(2,j-1).LE.PzyRglobal(2,j)) then
            go to 36
        else
            Lowerpoint = j !Lower inflection point
            go to 37
        end if
        36 continue
    end do
    37 continue
    
    allocate (aux3(2,NPrefino))
    allocate (aux4(2,NPrefino))
    allocate (aux6(NPrefino))
    do j=1,NPrefino
        aux3(2,j) = PzyRglobal(2,j)
        aux4(2,j) = PzyRglobal(2,j)
        aux6(j) = PzyRglobal(2,j)
    end do
    
    !CALCULATION OF THE NEW PROFILE IF NECESSARY
    if (Lowerpoint.EQ.DC) then !No inflection point below sea level
        !No change in possition vector
        roturanotch = Upperpoint
        go to 38 
    else if (Lowerpoint.NE.DC) then
        if ((Lowerpoint.EQ.Upperpoint).and.(Lowerpoint.LE.roturanotch)) then
            !No change in possition vector
            roturanotch = Lowerpoint
            go to 38
        else if ((Lowerpoint.EQ.Upperpoint).and.(Lowerpoint.GT.roturanotch)) then !Failure from Lowerpoint to roturanotch
            if (PzyRglobal(2,roturanotch).LT.PzyRglobal(2,Lowerpoint)) then
                !Failure to the left
!                mrecta = (PzyRglobal(1,roturanotch)-PzyRglobal(1,Lowerpoint))/(PzyRglobal(2,roturanotch)-PzyRglobal(2,Lowerpoint))
!                do j=Lowerpoint,roturanotch,-1
!	                PzyRglobal(2,j) = (PzyRglobal(1,j) + mrecta*PzyRglobal(2,Lowerpoint) - PzyRglobal(1,Lowerpoint))/mrecta
!	            end do
            else if (PzyRglobal(2,roturanotch).EQ.PzyRglobal(2,Lowerpoint)) then
                !Vertical failure
!                do j=Lowerpoint,roturanotch,-1
!                    PzyRglobal(2,j) = PzyRglobal(2,Lowerpoint) 
!                end do
                roturanotch = Lowerpoint !New failure possition or cliff base possition
            else if (PzyRglobal(2,roturanotch).GT.PzyRglobal(2,Lowerpoint)) then
                !Failure to the right
!                mrecta = (PzyRglobal(1,roturanotch)-PzyRglobal(1,Lowerpoint))/(PzyRglobal(2,roturanotch)-PzyRglobal(2,Lowerpoint))
!                do j=Lowerpoint,roturanotch,-1
!	                PzyRglobal(2,j) = (PzyRglobal(1,j) + mrecta*PzyRglobal(2,Lowerpoint) - PzyRglobal(1,Lowerpoint))/mrecta
!	            end do
	            roturanotch = Lowerpoint !New failure possition or cliff base possition
            end if
        else if (Lowerpoint.NE.Upperpoint) then
            if (PzyRglobal(2,Upperpoint).LT.PzyRglobal(2,Lowerpoint)) then
                !Failure to the left
                roturanotch = Upperpoint
                mrecta = (PzyRglobal(1,Upperpoint)-PzyRglobal(1,Lowerpoint))/(PzyRglobal(2,Upperpoint)-PzyRglobal(2,Lowerpoint))
                do j=Lowerpoint,Upperpoint,-1
	                aux3(2,j) = (PzyRglobal(1,j) + mrecta*PzyRglobal(2,Lowerpoint) - PzyRglobal(1,Lowerpoint))/mrecta
	            end do
            else if (PzyRglobal(2,Upperpoint).EQ.PzyRglobal(2,Lowerpoint)) then
                !Vertical failure
                roturanotch = Lowerpoint !New failure possition or cliff base possition
                do j=Lowerpoint,Upperpoint,-1
                    aux3(2,j) = PzyRglobal(2,Lowerpoint) 
                end do
            else if (PzyRglobal(2,Upperpoint).GT.PzyRglobal(2,Lowerpoint)) then
                roturanotch = Lowerpoint !New failure possition or cliff base possition
                !Failure to the right
                mrecta = (PzyRglobal(1,Upperpoint)-PzyRglobal(1,Lowerpoint))/(PzyRglobal(2,Upperpoint)-PzyRglobal(2,Lowerpoint))
                do j=Lowerpoint,Upperpoint,-1
	                aux3(2,j) = (PzyRglobal(1,j) + mrecta*PzyRglobal(2,Lowerpoint) - PzyRglobal(1,Lowerpoint))/mrecta
	            end do
            end if
        end if
    end if
    38 continue
    
    !Correction in case the new profile "creates" material in the Platform
    if (Lowerpoint.NE.Upperpoint) then
        do j=Lowerpoint,Upperpoint,-1
            if (aux3(2,j).GT.PzyRglobal(2,j)) then
                mrecta = (PzyRglobal(1,zerolevel)-PzyRglobal(1,Lowerpoint))/(PzyRglobal(2,zerolevel)-PzyRglobal(2,Lowerpoint))
                if (Lowerpoint.EQ.zerolevel) then
                   mrecta = 1.d0
                end if
                do i= Lowerpoint,zerolevel,-1
                    PzyRglobal(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,Lowerpoint) - PzyRglobal(1,Lowerpoint))/mrecta
                    if ((aux6(i).LT.PzyRglobal(2,i)).and.(i.NE.zerolevel).and.(i.NE.Lowerpoint)) then
                        do k=Lowerpoint,zerolevel,-1
                            PzyRglobal(2,k) = aux6(k)
                            go to 900
                        end do
                    end if
                end do
                go to 900
            else
                PzyRglobal(2,j) = aux3(2,j)
            end if
        end do
    end if
    900 continue
    deallocate (aux3)
    deallocate (aux6)
    
    !In case ptosup and ptoinf are not the same failure between them
    if ((Upperpoint.NE.ptoinf).and.(PzyRglobal(2,ptoinf).GT.PzyRglobal(2,Upperpoint))) then
        mrecta = (PzyRglobal(1,Upperpoint)-PzyRglobal(1,ptoinf))/(PzyRglobal(2,Upperpoint)-PzyRglobal(2,ptoinf))
        do i= ptoinf,Upperpoint,-1
            aux4(2,i) = (PzyRglobal(1,i) + mrecta*PzyRglobal(2,ptoinf) - PzyRglobal(1,ptoinf))/mrecta    
        end do
        !if the new profile "creates" material -> NO new profile appears and roturanotch is ptoinf
        Dibujo = 1
        if (aux5(2,Upperpoint).EQ.aux4(2,Upperpoint)) then
            Dibujo = 0
            roturanotch = ptoinf
            go to 901
        end if
        do j=ptoinf,Upperpoint,-1
            if (aux4(2,j).GT.PzyRglobal(2,j)) then
                Dibujo = 0
                roturanotch = ptoinf
                go to 901
            end if
        end do
        if (Dibujo.EQ.1) then
            do i= ptoinf,Upperpoint,-1
                PzyRglobal(2,i) = aux4(2,i)    
            end do
        end if    
    end if
    901 continue
    deallocate (aux4)
    deallocate (aux5)
    
    !LANDSLIDE FAILURE EVALUATION
    1 Continue
    !Second failure point (meter to meter)
    DistPcirc2 = PzyRglobal(2,roturanotch) - retroceso
    !Determine intersection circular failure VS. cliff
    do 39 i=roturanotch-1,2,-1
        !Determine director´s vector for each point
        u1 = PzyRglobal(2,i-1) - PzyRglobal(2,i)
	    u2 = PzyRglobal(1,i-1) - PzyRglobal(1,i)
	    !Profile slope
	    if (u1.EQ.0.0) then
	        indvertcliff = 1
	        mcliff = 0.0
        else
		    indvertcliff = 0
		    mcliff = u2/u1
        end if
	    !Intersection points between circular failure and cliff
	    if (mcliff.EQ.mcirculo) then
		    go to 39
	    else if (mcirculo.EQ.0.d0) then
		    XDistPcirc2 = DistPcirc2
			YDistPcirc2 = mcliff*XDistPcirc2 + PzyRglobal(1,i) - mcliff*PzyRglobal(2,i)
		end if
	    !Lambda parameter that determine intersection point within the cliff
		if (indvertcliff.EQ.0.0) then
		    lambda = (XDistPcirc2 - PzyRglobal(2,i))/u1
	    else
		    lambda = (YDistPcirc2 - PzyRglobal(1,i))/u2
		end if
	    !Coordenate of cross point if it exists
		if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
		    PtocorteCirc2 = i
			go to 40
	    else if (lambda.EQ.1.0) then
		    PtocorteCirc2 = i-1
			go to 40
		else
		    go to 39
		end if
	39 continue
	40 continue
        		
    !Cliff height from talus wedge and talus base point - roturanotch
	Htalud = abs(PzyRglobal(1,crestatalud) - PzyRglobal(1,roturanotch))
        		
	!Circular failure calculus limit -> when horizontal distance exceeds Htalud
	if (DistPcirc2.LT.PzyRglobal(2,crestatalud)-Htalud) then
	    go to 41 !No more circular failures can exist
	end if
        		
    !Determine slope between bottom failure point VS. upper failure point
	munion = (YDistPcirc2 - PzyRglobal(1,roturanotch))/(XDistPcirc2 - PzyRglobal(2,roturanotch))
        		
    !Failure center coordenates (SIGN)
	CXc = ((-PzyRglobal(2,roturanotch)/PendientetgPcirc1)-(PzyRglobal(1,roturanotch)/2.)+(YDistPcirc2/2.)-((-1/munion)*PzyRglobal(2,roturanotch)/2.)-((-1/munion)*XDistPcirc2/2.))/((-1/PendientetgPcirc1) + (1/munion))
	CYc = (-1/PendientetgPcirc1)*(CXc - PzyRglobal(2,roturanotch)) + PzyRglobal(1,roturanotch)
        		
    !Radius of a circunference
	radio = dsqrt((CXc - PzyRglobal(2,roturanotch))**2 + (CYc - PzyRglobal(1,roturanotch))**2)
        		
    !Slide size upper part (b)
	anchoslide = dabs(XDistPcirc2 - PzyRglobal(2,roturanotch))/nslides
        		
	!Inizialization of the numerator and denominator for Factor of Safety: FS
	Num = 0.d0
	Deno = 0.d0
        		
    !Safety Factor
	do i=1,nslides
        !X-Coordenate of Right extreme
		Xca = PzyRglobal(2,roturanotch) - (i-1)*anchoslide*Factor
		!X-Coordenate of Left extreme
		Xcb = PzyRglobal(2,roturanotch) - i*anchoslide*Factor
		!X-Coordenate of Medium slide point
		Xmedio = (Xca+Xcb)/2.
        		    
		!Y-Coordenate of Right extreme
		Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xca**2 + 2*CXc*Xca - CYc**2)))/2.
		Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xca**2 + 2*CXc*Xca - CYc**2)))/2.
		!Choose the correct Y value -> the lowest one (bottom part of the circular failure)
		if (Yopc1.LT.Yopc2) then
		    Yca = Yopc1
		else
		    Yca = Yopc2
		end if
        		    
		!Y-Coordenate of Left extreme
		Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xcb**2 + 2*CXc*Xcb - CYc**2)))/2.
		Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xcb**2 + 2*CXc*Xcb - CYc**2)))/2.
		!Choose the correct Y value -> the lowest one (bottom part of the circular failure)
		if (Yopc1.LT.Yopc2) then
		    Ycb = Yopc1
		else
		    Ycb = Yopc2
		end if
        		    
		!Y-Coordenate of Medium slide point
		Yopc1 = (2*CYc + dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xmedio**2 + 2*CXc*Xmedio - CYc**2)))/2.
		Yopc2 = (2*CYc - dsqrt(4*CYc**2 + 4*(radio**2 - CXc**2 - Xmedio**2 + 2*CXc*Xmedio - CYc**2)))/2.
		!Choose the correct Y value -> the lowest one (bottom part of the circular failure)
		if (Yopc1.LT.Yopc2) then
		    Ymedio = Yopc1
		else
		    Ymedio = Yopc2
		end if
        		  
		!Slope of medium slide point
		Mmedio = (CXc-Xmedio)/(Ymedio-CYc)
        		    
		!Center slide angle
		gamma = dabs((180/pi)*datan(Mmedio))
        		    
		!Determine intersection between each slide and the cliff profile
		do 42 j=roturanotch,2,-1
		    !Determine director´s vector for each point
		    u1 = PzyRglobal(2,j-1) - PzyRglobal(2,j)
			u2 = PzyRglobal(1,j-1) - PzyRglobal(1,j)
		    !Determine cliff slope
			if (u1.EQ.0.0) then
			    indvertcliff = 1
		        mcliff = 0.0
            else
			    indvertcliff = 0
				mcliff = u2/u1
			end if
			!Determine cross point between medium slide point and cliff
			XPtomedio = Xmedio
			YPtomedio = mcliff*XPtomedio + PzyRglobal(1,j) - mcliff*PzyRglobal(2,j)
			!Lambda parameter that determine intersection point within the cliff
			if (indvertcliff.EQ.0.0) then
			    lambda = (XPtomedio - PzyRglobal(2,j))/u1
			else
			    lambda = (YPtomedio - PzyRglobal(1,j))/u2
			end if
			!Coordenate of cross point if it exists
			if ((lambda.GE.0.0).AND.(lambda.LE.1.0)) then
			    go to 43
			else
			    go to 42
			end if
		42 continue
		43 continue
        		    
		!Slide height
		Hslide = dabs(YPtomedio - Ymedio)
        		    
		!Groundwater definition -> Hw
		if (groundwater.EQ.1) then
		    if (Xmedio.GE.PzyRglobal(2,seepagelevel)) then
		        Hw = hslide
		        go to 44
		    else if (Xmedio.LT.PzyRglobal(2,seepagelevel)) then
		        Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
		        Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
		        distY = dabs(PzyRglobal(2,seepagelevel)-Xmedio)
		        curvah = dsqrt((distY/longitud)*(Profseepage-Profagua)**2)
		        Hw = dabs(curvah+PzyRglobal(1,seepagelevel))-Ymedio
		        if (Hw.LE.0.d0) then
		            Hw = 0.d0
		        end if 
		        go to 44
		    end if
		else if (groundwater.EQ.0) then
		    Hw = 0.d0
		    go to 44
		else if (groundwater.EQ.2) then
		    Hw = Hslide
		    go to 44
		end if
        44 continue 
        		    
		!Geological layer at the bottom
		tp = 1
		potacumulada = potestrato (tp)
        45 continue
		if (((PzyRglobal(1,1)-potacumulada).LE.Ymedio).and.(Ymedio.LE.PzyRglobal(1,1)-(potacumulada-potestrato(tp)))) then
		    !Units in KPa = KN/m**2
			Ecohesionrot = Ecohesion(tp)
			Ephirot = Ephi(tp)
			Pesounitrot = Pesounit(tp)
			go to 46		            
		else
    		tp = tp + 1
			potacumulada = potacumulada + potestrato(tp)
			go to 45
		end if
		46 continue
        			
	    !Determine the Numerator
		Num = Num + Ecohesionrot*(anchoslide/dcos(gamma*pi/180)) + (Hslide*Pesounitrot*anchoslide*dcos(gamma*pi/180) - Hw*freshwater*anchoslide/dcos(gamma*pi/180))*dtan(Ephirot*pi/180)
			        
   		!Determine the Denominator
		Deno = Deno + Hslide*Pesounitrot*anchoslide*dsin((gamma*pi/180))
    end do
                
    !Determine the global SAFETY FACTOR (FS) for each circular failure
    !If the value is < 1 -> FAILURE OCCURS
    FS = Num/Deno
                
    !Calculated circular failure is stable or instable?
    if (FS.GE.1.d0) then
        retroceso = retroceso + 0.5d0
        go to 1
    else if (FS.LT.1.d0) then
        retroceso = retroceso + 0.5d0
        !Save the possitions vector previous failure in Failure Matrix in case its the first failure of the simulation steep
        if (Parametrorotura.EQ.0) then
            dimmatrizrotura = size(Matrizrotura,dim=2)
			allocate (dummy(NPrefino+1,dimmatrizrotura))
			do j=1,dimmatrizrotura
			    do i=1,NPrefino+1
				    dummy(i,j) = Matrizrotura (i,j)
				end do
			end do
		    
		    deallocate (Matrizrotura)
			allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+1))
		    do j=1,dimmatrizrotura
			    do i=1,NPrefino+1
				    Matrizrotura(i,j) = dummy (i,j)
				end do 
			end do

		    deallocate (dummy)
            !Determine failue time and matrix possition
			Matrizrotura(1,dimmatrizrotura+1) = Tglob !Output time in months
		    !Failure matrix
			do i=1,NPrefino
			    Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
			end do
        end if
                
        !Circular failure from roturanotch to PtocorteCirc2
        do i=roturanotch-1,PtocorteCirc2,-1
            aux = 0.d0
            aux = (2*CXc - dsqrt(4*CXc**2 + 4*(radio**2 - CXc**2 - CYc**2 + 2*CYc*PzyRglobal(1,i) - (PzyRglobal(1,i))**2)))/2.
            if (aux.LE.PzyRglobal(2,i)) then
                PzyRglobal(2,i) = aux
                go to 47
            else if (aux.GT.PzyRglobal(2,i))then
                aux = 0.d0
                aux = (2*CXc + dsqrt(4*CXc**2 + 4*(radio**2 - CXc**2 - CYc**2 + 2*CYc*PzyRglobal(1,i) - (PzyRglobal(1,i))**2)))/2.
                PzyRglobal(2,i) = aux
            end if
            47 continue   
        end do
                    
        Vertcirc = 0
        !Fix geometry in case circular failure and cliff have strange geometry
        do i=roturanotch-1,2,-1
            if (PzyRglobal(2,i-1).LE.PzyRglobal(2,i)) then
                if (Vertcirc.EQ.1) then
                    !go to 48
                end if
                !If material is created we avoid this circular failure
                if (PzyRglobal(2,i).GT.Matrizrotura(i+1,dimmatrizrotura+1)) then
                    Parametrorotura = 0
                    do j=1,NPrefino
				        PzyRglobal (2,j) = Matrizrotura(j+1,dimmatrizrotura+1)
				    end do
                                    
                    allocate (dummy(NPrefino+1,dimmatrizrotura))
				    do k=1,dimmatrizrotura
				        do j=1,NPrefino+1
					        dummy(j,k) = Matrizrotura (j,k)
					    end do
				    end do
				    deallocate (Matrizrotura)
		            allocate (Matrizrotura(NPrefino+1,dimmatrizrotura))

				    do k=1,dimmatrizrotura
				        do j=1,NPrefino+1
					        Matrizrotura(j,k) = dummy (j,k)
					    end do 
			        end do
				    deallocate (dummy)
                    go to 1
                end if
                go to 49
            else if (PzyRglobal(2,i-1).GT.PzyRglobal(2,i))then
                PzyRglobal(2,i-1) = PzyRglobal(2,i)
                Vertcirc = 1
                !If material is created we avoid this circular failure
                if (PzyRglobal(2,i).GT.Matrizrotura(i+1,dimmatrizrotura+1)) then
                    Parametrorotura = 0
                    do j=1,NPrefino
				        PzyRglobal (2,j) = Matrizrotura(j+1,dimmatrizrotura+1)
				    end do                
                    allocate (dummy(NPrefino+1,dimmatrizrotura))
				    do k=1,dimmatrizrotura
				        do j=1,NPrefino+1
					        dummy(j,k) = Matrizrotura (j,k)
					    end do
		            end do

				    deallocate (Matrizrotura)
				    allocate (Matrizrotura(NPrefino+1,dimmatrizrotura))
                    do k=1,dimmatrizrotura
                        do j=1,NPrefino+1
                            Matrizrotura(j,k) = dummy (j,k)
                        end do 
		            end do
                    deallocate (dummy)
                    go to 1
                end if
                go to 49
            end if
            49 continue
        end do
        48 continue
        
        !Save the possitions vector after failure in FAILURE MATRIX          
        Parametrorotura = 1
        dimmatrizrotura = size(Matrizrotura,dim=2)
	    allocate (dummy(NPrefino+1,dimmatrizrotura))
		do k=1,dimmatrizrotura
		    do i=1,NPrefino+1
			    dummy(i,k) = Matrizrotura (i,k)
			end do
		end do
		deallocate (Matrizrotura)
		allocate (Matrizrotura(NPrefino+1,dimmatrizrotura+1))
		do k=1,dimmatrizrotura
		    do i=1,NPrefino+1
			    Matrizrotura(i,k) = dummy (i,k)
			end do 
		end do
		deallocate (dummy)
        !Failure time
		Matrizrotura(1,dimmatrizrotura+1) = Tglob !Output time in months
		!Failure matrix
		do i=1,NPrefino
		    Matrizrotura(i+1,dimmatrizrotura+1) = PzyRglobal(2,i)
		end do
            		
		!RECESSION MATRIX
		dimmatrizrecesion = size(Matrizrecesion,dim=2)
		allocate (dummy1(NPrefino+2,dimmatrizrecesion))
        do k=1,dimmatrizrecesion
		    do i=1,NPrefino+2
		        dummy1(i,k) = Matrizrecesion (i,k)
		    end do
	    end do
		deallocate (Matrizrecesion)
		allocate (Matrizrecesion(NPrefino+2,dimmatrizrecesion+1))
		do k=1,dimmatrizrecesion
		    do i=1,NPrefino+2
			    Matrizrecesion(i,k) = dummy1 (i,k)
	        end do
	    end do
		deallocate (dummy1)
        !Failure time
	    Matrizrecesion(1,dimmatrizrecesion+1) = Tglob !Output time in months
	    !Recession matrix
	    do i= 2,NPrefino+1
	        Matrizrecesion(i+1,dimmatrizrecesion+1) = Matrizrotura(i,dimmatrizrotura+1) - Matrizrotura(i,dimmatrizrotura)
        end do
        sumaincrementoy = 0.d0
		do i=1,NPrefino
		    sumaincrementoy = sumaincrementoy + abs(Matrizrecesion(i+2,dimmatrizrecesion+1))
		end do
		Matrizrecesion(2,dimmatrizrecesion+1) = sumaincrementoy*Dz
		Volacumulado1 = Volacumulado1 + sumaincrementoy*Dz
        write(2,11) Tglob,sumaincrementoy*Dz,Volacumulado1,(abs(Matrizrecesion(label(ilabel)+2,dimmatrizrecesion+1)),ilabel=1,nlabel)
                
        !GROUNDWATER LEVELS MATRIX
        if (groundwater.EQ.0) then
            go to 65
        else if (groundwater.EQ.1) then
            Profseepage = dabs(PzyRglobal(1,seepagelevel)-PzyRglobal(1,1))
            Longitud = dabs(PzyRglobal(2,seepagelevel)-PzyRglobal(2,1))
	        dist = longitud/1000.
	        do i=1,1000
	            distancia(i) = PzyRglobal(2,seepagelevel) - (i-1)*dist
	        end do
	        do i=1,1000
                Hwdraw(i) = dsqrt((dabs(distancia(i)-PzyRglobal(2,seepagelevel))/longitud)*(Profseepage-Profagua)**2) + PzyRglobal(1,seepagelevel)
	        end do
	        !Groundwater levels to print it
	        dimmatrizHwdraw = size(MatrizHwdraw,dim=2)
	        allocate (dummy2(1000+1,dimmatrizHwdraw))
		    do j=1,dimmatrizHwdraw
                do i=1,1000+1
	            dummy2(i,j) = MatrizHwdraw (i,j)
    		    end do
		    end do
		    deallocate (MatrizHwdraw)
		    allocate (MatrizHwdraw(1000+1,dimmatrizHwdraw+2))
		    do j=1,dimmatrizHwdraw
		        do i=1,1000+1
			        MatrizHwdraw(i,j) = dummy2 (i,j)
	            end do 
	        end do
		    deallocate (dummy2)
            !Failure time
	        MatrizHwdraw(1,dimmatrizHwdraw+1) = Tglob !Output time in months
		    MatrizHwdraw(1,dimmatrizHwdraw+2) = Tglob !Output time in months
		    !Groundwater print matrix -> Hw
		    do i=1,1000
		        MatrizHwdraw(i+1,dimmatrizHwdraw+1) = Hwdraw(i)
	            MatrizHwdraw(i+1,dimmatrizHwdraw+2) = distancia(i)
            end do
        end if
        65 continue
    end if
    go to 1
    41 continue
    
    if (Parametrorotura.EQ.1) then !Landslide occurs -> NO toppling
        Volacumulado2 = 0.d0
        !Total DEBRIS volume
        Volacumulado = Volacumulado1
        go to 54
    else if (Parametrorotura.EQ.0) then
        Volacumulado1 = 0.d0
        !TOPPLING FAILURE EVALUATION
        minestabilidad = 0.d0
        Ptocorte = 0
        do 50 j=roturanotch-1,2,-1
            !Determine director´s vector for each cliff point
		    u1 = PzyRglobal(2,j-1) - PzyRglobal(2,j)
		    u2 = PzyRglobal(1,j-1) - PzyRglobal(1,j)
		    !Determine cliff slope
		    if (u1.EQ.0.0) then
	            indvertcliff = 1
			    mcliff = 0.0
	        else
	            indvertcliff = 0
			    mcliff = u2/u1
	        end if
	        !Determine cross points between vertical failure and cliff
		    if (mcliff.EQ.minestabilidad) then
		        go to 50
	        else
	   	        Yc = PzyRglobal(2,roturanotch)
	    	    Zc = mcliff*Yc + PzyRglobal(1,j) - mcliff*PzyRglobal(2,j)
	        end if
	        !Lambda parameter that determine intersection point within the cliff
		    if (indvertcliff.EQ.0.0) then
                lambda = (Yc - PzyRglobal(2,j))/u1
	        else
	            lambda = (Zc - PzyRglobal(1,j))/u2
	        end if
		    !Coordenate of cross point if it exists
	        if ((lambda.GE.0.0).AND.(lambda.LT.1.0)) then
	            Ptocorte = j
	            go to 51
	        else if (lambda.EQ.1.0) then
			    Ptocorte = j-1
			    go to 51
	        else
			    go to 50
	        end if
        50 continue
	    51 continue	    
        !If cutoff doesn´t exist -> no notch => next erosive event
        if (Ptocorte.EQ.0) then
	        go to 54
        end if
        !Vertical vector equal to cliff
        allocate (vectorverticalrot(NPrefino))
	    do i=1,NPrefino
	        vectorverticalrot(i) = PzyRglobal(2,i)
        end do
        !Vetical plane of failure
        nrot = 0
        CoordX = 0.d0
        do i = roturanotch,ptocorte,-1
	        vectorverticalrot(i) = PzyRglobal (2,roturanotch)
		    nrot = nrot + 1
		    CoordX = CoordX + dabs(PzyRglobal(2,i))
        end do		    
        CoordX = CoordX/nrot
	    !Localization of the point roturanotch in its geological layer
        k=1
	    potacumulada = potestrato(k)
	    Sigmatensilerock = 0.d0
	    Pesounitrot = 0.d0
        52 continue
        if ((PzyRglobal(1,1)-potacumulada).LE.PzyRglobal(1,roturanotch)) then
	        !Units in KPa = KN/m**2
            !!!!Sigmatensilerock = Ecohesion(k)*0.1d0
		    Sigmatensilerock = Shearstrength(k)
	        Pesounitrot = Pesounit(k)
		    go to 53		            
	    else
            k = k + 1
		    potacumulada = potacumulada + potestrato(k)
		    go to 52
	    end if
	    53 continue
    			
	    !Volume of the possible unstable wedge
        sumaacumuladarot = 0.d0
	    aux = 0.d0
	    do i = roturanotch,ptocorte,-1
	        aux = abs(vectorverticalrot(i)-PzyRglobal(2,i))
	        sumaacumuladarot = sumaacumuladarot + aux
	    end do
	    !Units in KN/m**2
	    Volacumulado2 = sumaacumuladarot*dz
	    sumaacumuladarot = Volacumulado2*pesounitrot
	    hrot = dabs(PzyRglobal(1,roturanotch)-PzyRglobal(1,ptocorte))
	    Lcdg = dabs(PzyRglobal(2,roturanotch)) - CoordX
	    Sigmatensilemax = (6.*sumaacumuladarot*Lcdg)/(hrot**2)
    			
	    !SAFETY FACTOR VS. TOPPLING FAILURE
	    if (Sigmatensilemax.GE.Sigmatensilerock) then
		    do i = Ptocorte,roturanotch
		        PzyRglobal(2,i) = vectorverticalrot(i)
	        end do
	        !Total DEBRIS volume
	        Volacumulado = Volacumulado2
	        deallocate (vectorverticalrot)
	        go to 54
	    else if (Sigmatensilemax.LT.Sigmatensilerock) then
	        !No DEBRIS
	        Volacumulado = 0.d0
	        deallocate (vectorverticalrot)
	        go to 54
	    end if
	end if
    54 continue !It comes from the absence of notch -> ptocorte = 0


    !RESULTS MATRIX from third collum (t=1) until the end (t=j*tdatos)
    if ((Tglob.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	    !Time vectors to introduce in the matrix
		DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
    	!Transfer the position information matrix into DRAWING MATRIX
		do i=1,NPrefino
    	    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	    end do
   		jdibujo = jdibujo + 1
    end if
    
    !DETERMINE TIME WITHOUT EROSION DUE TO DEBRIS MATERIALS
    if (Volacumulado.EQ.0.d0) then
        Pmonth = 0
        Dprot = 0.d0
        Drem = 0.d0
        go to 55
    end if
    Drem = Volacumulado/Mf !Result in days
    Pmonth = int(Drem/30.) !If the value is greater than 1, Pmonth is the time simulation jump (months)
    Dprot = dabs(Pmonth*30-Drem) !The value that modifies in No the number of days
    55 continue
    
    !DETERMINE YEAR SEASON TO INTRODUCE THE WAVE VALUES
    if (Pmonth.EQ.0) then
        !YEARS WITHOUT CHANGE IN SEA LEVEL until equilibrium conditions are reached
!        if (Tglob.GE.(12*5)) then !5 years without lakelevel change
!            changelakelevel = 1
!        end if
        
        if (changelakelevel.EQ.1) then
            !15 years cicle to change the lake level
            if (Tglob.EQ.((5*12)+9+ciclo)) then
                !1st change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+18+ciclo)) then
                !2nd change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+27+ciclo)) then
                !3rd change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+36+ciclo)) then
                !4th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+45+ciclo)) then
                !5th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+54+ciclo)) then
                !6th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+63+ciclo)) then
                !7th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+72+ciclo)) then
                !8th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+81+ciclo)) then
                !9th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+90+ciclo)) then !Half of the cicle 7.5 years
                !10th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+99+ciclo)) then
                !11th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+108+ciclo)) then
                !12th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+117+ciclo)) then
                !13th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+126+ciclo)) then
                !14 change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+135+ciclo)) then
                !15th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+144+ciclo)) then
                !16th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+153+ciclo)) then
                !17th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+162+ciclo)) then
                !18th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+171+ciclo)) then
                !19th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 56
            else if (Tglob.EQ.((5*12)+180+ciclo)) then !Cicle´s end 15 years
                !20th change
                controllevel = controllevel-1 !TL rise in 0.1
                ciclo = ciclo + 180 !180 = 15(years)*12(months)
                go to 56
            end if
        end if
        56 continue
        
        if (Tglob.EQ.3*salto) then
            !SUMMER waves conditions
            Waves = 4 !Number of waves in the season SUMMER
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.5
            H(2) = 1
            H(3) = 2
            H(4) = 3
            !Introduction of the Wave period values
            T(1) = 2.
            T(2) = 3.5
            T(3) = 4.5
            T(4) = 6.5
            !Calculation of the Deep water wave length values
            do i=1,Waves
                L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
            end do
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 31.27
            WF(2) = 41.73
            WF(3) = 1.9
            WF(4) = 0.33
            !Total frecuency = 75.23 %
            !Calculation of the The number of waves each hour 
            do i=1,Waves
                Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
            end do
            
            !SUMMER lake level
            !TL = controllevel - int(0.3/Dz) !Increse 0.3m from WINTER -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)

            go to 1001 !Loops end before -> invierno =0
        else if (Tglob.EQ.3*(salto+1)) then
            !FALL waves conditions
            Waves = 7 !Number of waves in the season FALL
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.5
            H(2) = 1
            H(3) = 2
            H(4) = 3
            H(5) = 4
            H(6) = 5
            H(7) = 6
            !Introduction of the Wave period values
            T(1) = 2.
            T(2) = 3.8
            T(3) = 4.5
            T(4) = 5.5
            T(5) = 6.5
            T(6) = 7.
            T(7) = 7.5
            !Calculation of the Deep water wave length values
            do i=1,Waves
                L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
            end do
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 8.87
            WF(2) = 36.27
            WF(3) = 16.57
            WF(4) = 9.40
            WF(5) = 1.5
            WF(6) = 0.2
            WF(7) = 0.1
            !Total frecuency = 72.9 %
            !Calculation of the The number of waves each hour 
            do i=1,Waves
                Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
            end do

            !FALL lake level
            !TL = controllevel !Equal level from WINTER -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)
            
            go to 1001 !Loops end before -> invierno =0
        else if (Tglob.EQ.3*(salto+2)) then
            invierno = 1 !Winter conditions
            counter = counter + 3 !Winter conditions
            !WINTER waves conditions
            Waves = 1
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.d0
            !Introduction of the Wave period values
            T(1) = 0.d0
            !Calculation of the Deep water wave length values
            L(1) = 1.56*(T(1)**2) !Units: [m] - CERC (2002)
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 0.d0
            !Calculation of the The number of waves each hour 
            Wn(1) = (60*60)/T(1) !60 s/min * 60 min/hours
            
            !WINTER lake level
            !TL = controllevel !No variation from FALL -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)
            
            go to 1002 !Loop end after re-inizialization of variables
        end if
        
    else if (Pmonth.NE.0) then
    
        Draw = Tglob
        Tglob = Tglob + Pmonth !Simulation jump due to debris toe protection
        !YEARS WITHOUT CHANGE IN SEA LEVEL until equilibrium conditions are reached
!        if (Tglob.GE.(12*5)) then !5 years without lakelevel change
!            changelakelevel = 1
!        end if
        if (Tglob.GT.Draw+12) then
            yearjump = int((Tglob-Draw)/12)
            salto = salto + 4*yearjump !Winter conditions
        end if
        
        if (changelakelevel.EQ.1) then
            !15 years cicle to change the lake level 
            if ((Draw.LE.((5*12)+9+ciclo)).and.(((5*12)+9+ciclo).LE.Tglob)) then
                !1st change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+18+ciclo)).and.(((5*12)+18+ciclo).LE.Tglob)) then
                !2nd change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+27+ciclo)).and.(((5*12)+27+ciclo).LE.Tglob)) then
                !3rd change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+36+ciclo)).and.(((5*12)+36+ciclo).LE.Tglob)) then
                !4th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+45+ciclo)).and.(((5*12)+45+ciclo).LE.Tglob)) then
                !5th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+54+ciclo)).and.(((5*12)+54+ciclo).LE.Tglob)) then
                !6th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+63+ciclo)).and.(((5*12)+63+ciclo).LE.Tglob)) then
                !7th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+72+ciclo)).and.(((5*12)+72+ciclo).LE.Tglob)) then
                !8th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+81+ciclo)).and.(((5*12)+81+ciclo).LE.Tglob)) then
                !9th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+90+ciclo)).and.(((5*12)+90+ciclo).LE.Tglob)) then !Half of the cicle 7.5 years
                !10th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+99+ciclo)).and.(((5*12)+99+ciclo).LE.Tglob)) then
                !11th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+108+ciclo)).and.(((5*12)+108+ciclo).LE.Tglob)) then
                !12th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+117+ciclo)).and.(((5*12)+117+ciclo).LE.Tglob)) then
                !13th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+126+ciclo)).and.(((5*12)+126+ciclo).LE.Tglob)) then
                !14 change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+135+ciclo)).and.(((5*12)+135+ciclo).LE.Tglob)) then
                !15th change
                controllevel = controllevel+1 !TL decrease in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+144+ciclo)).and.(((5*12)+144+ciclo).LE.Tglob)) then
                !16th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+153+ciclo)).and.(((5*12)+153+ciclo).LE.Tglob)) then
                !17th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+162+ciclo)).and.(((5*12)+162+ciclo).LE.Tglob)) then
                !18th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+171+ciclo)).and.(((5*12)+171+ciclo).LE.Tglob)) then
                !19th change
                controllevel = controllevel-1 !TL rise in 0.1
                go to 57
            else if ((Draw.LE.((5*12)+180+ciclo)).and.(((5*12)+180+ciclo).LE.Tglob)) then !Cicle´s end 15 years
                !20th change
                controllevel = controllevel-1 !TL rise in 0.1
                ciclo = ciclo + 180 !180 = 15(years)*12(months)
                go to 57
            end if
        end if
        57 continue
        
        if (Tglob.GE.tsim) then
            do j=Draw,tsim
                !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                !Time vectors to introduce in the matrix
		            DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
                    !Transfer the position information matrix into DRAWING MATRIX
	                do i=1,NPrefino
 	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                end do
   		            jdibujo = jdibujo + 1
                end if
            end do
            go to 1000 !End of simulation
        else if (Tglob.LT.tsim) then
            do j=Draw,Tglob
                !Results matrix from third collum (t=1) until the end (t=j*tdatos)
                if ((j.EQ.(jdibujo*Tdib)).and.((jdibujo*Tdib).LE.tsim)) then
	                !Time vectors to introduce in the matrix
		            DibMatriz(1,jdibujo+2) = jdibujo*Tdib !Frist matrix position -> Time (months)
                    !Transfer the position information matrix into DRAWING MATRIX
	                do i=1,NPrefino
 	                    DibMatriz(i+1,jdibujo+2) = PzyRglobal(2,i)
	                end do
   		            jdibujo = jdibujo + 1
                end if
            end do
        end if
        
        if ((3*salto.LE.Tglob).and.(Tglob.LT.3*(salto+1))) then
            !SUMMER waves conditions
            Waves = 4 !Number of waves in the season SUMMER
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.5
            H(2) = 1
            H(3) = 2
            H(4) = 3
            !Introduction of the Wave period values
            T(1) = 2.
            T(2) = 3.5
            T(3) = 4.5
            T(4) = 6.5
            !Calculation of the Deep water wave length values
            do i=1,Waves
                L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
            end do
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 31.27
            WF(2) = 41.73
            WF(3) = 1.9
            WF(4) = 0.33
            !Total frecuency = 75.23 %
            !Calculation of the The number of waves each hour 
            do i=1,Waves
                Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
            end do
            
            !SUMMER lake level
            !TL = controllevel - int(0.3/Dz) !Increase 0.3m from WINTER -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)
       
            go to 1001 !Loops end before -> invierno =0
        else if ((3*(salto+1).LE.Tglob).and.(Tglob.LT.3*(salto+2))) then
            !FALL waves conditions
            Waves = 7 !Number of waves in the season FALL
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.5
            H(2) = 1
            H(3) = 2
            H(4) = 3
            H(5) = 4
            H(6) = 5
            H(7) = 6
            !Introduction of the Wave period values
            T(1) = 2.
            T(2) = 3.8
            T(3) = 4.5
            T(4) = 5.5
            T(5) = 6.5
            T(6) = 7.
            T(7) = 7.5
            !Calculation of the Deep water wave length values
            do i=1,Waves
                L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
            end do
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 8.87
            WF(2) = 36.27
            WF(3) = 16.57
            WF(4) = 9.40
            WF(5) = 1.5
            WF(6) = 0.2
            WF(7) = 0.1
            !Total frecuency = 72.9 %
            !Calculation of the The number of waves each hour 
            do i=1,Waves
                Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
            end do

            !FALL lake level
            !TL = controllevel !Equal from WINTER -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)
           
            go to 1001 !Loops end before -> invierno =0
        else if ((3*(salto+2).LE.Tglob).and.(Tglob.LT.3*(salto+3))) then
            invierno = 1 !Winter conditions
            counter = counter + 3 !Winter conditions
            !WINTER waves conditions
            Waves = 1
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.d0
            !Introduction of the Wave period values
            T(1) = 0.d0
            !Calculation of the Deep water wave length values
            L(1) = 1.56*(T(1)**2) !Units: [m] - CERC (2002)
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 0.d0
            !Calculation of the The number of waves each hour 
            Wn(1) = (60*60)/T(1) !60 s/min * 60 min/hours
            
            !WINTER lake level
            !TL = controllevel !No variation from FALL -> ARMY CORPS OF ENGINEERS - DETROIT District (2011)
            
            go to 1002 !Loop end after re-inizialization of variables
        else if (Tglob.EQ.3*(salto+3)) then
            counter = 1 !Winter conditions
            salto = salto + 4 !Winter conditions
            !SPRING waves conditions
            Waves = 4 !Number of waves in the season SPRING
            !Creation of H vector -> Deep water wave height
            deallocate (H)
            allocate (H(Waves))
            !Creation of T vector -> Wave period
            deallocate (T)
            allocate (T(Waves))
            !Creation of L vector -> Deep water wave length
            deallocate (L)
            allocate (L(Waves))
            !Creation of WF vector -> Percentage frecuency of each type of wave per day
            deallocate (WF)
            allocate (WF(Waves))
            !Creation of Wn vector -> The number of waves each hour
            deallocate (Wn)
            allocate (Wn(Waves))
            
            !Introduction of the Deep water wave height values
            H(1) = 0.5
            H(2) = 1
            H(3) = 2
            H(4) = 3
            !Introduction of the Wave period values
            T(1) = 2.
            T(2) = 3.5
            T(3) = 4.5
            T(4) = 5.5
            !Calculation of the Deep water wave length values
            do i=1,Waves
                L(i) = 1.56*(T(i)**2) !Units: [m] - CERC (2002)
            end do
            !Introduction of the Percentage frecuency of each type of wave 
            WF(1) = 33.03
            WF(2) = 36.33
            WF(3) = 2.23
            WF(4) = 0.1
            !Total frecuency = 71.67 %
            !Calculation of the The number of waves each hour 
            do i=1,Waves
                Wn(i) = (60*60)/T(i) !60 s/min * 60 min/hours
            end do
            
            !SPRING lake level
            !TL = controllevel - int(0.4/Dz) !Increase 0.4m is the value obtained from ARMY CORPS OF ENGINEERS - DETROIT District (2011)
            
            go to 1001 !Loops end before -> invierno =0
        end if
    end if
    
    1001 continue
    invierno = 0
    1002 continue
    Pmonth = 0
    Drem = 0
	!RE-INIZIALIZATION OF VARIABLES
    Volacumulado1 = 0.d0
    Volacumulado2 = 0.d0
    Volacumulado = 0.d0    
    Parametrorotura = 0
    DistPcirc2 = 0.d0
    retroceso = 0.5d0
    
    Tglob = Tglob + 1
    GO TO 500
END IF
1000 continue

call cpu_time(tac)
write(*,*)'tiempo CPU:  ',tac-tic,'  seg.'
11 format(F10.1,2x,F9.5,2x,10(F9.6,2x))

!CLIFF PROFILE FOR EACH Tdib TIME
DibMatriz(1,2) = tac-tic !CPU time
do i = 1,NPrefino+1
	write (4,'(200(F15.7,x))') (DibMatriz(i,k),k=1,DimDiby+2)
end do

!CLIFF PROFILE BEFORE-AFTER FAILURE EVENT
if (dimmatrizrotura.EQ.size(Matrizrotura,dim=2)) then
    do i = 1,NPrefino+1
        write (6,'(200(F15.7,x))') (Matrizrotura(i,k),k=1,dimmatrizrotura)
    end do
else
    dimmatrizrotura = size(Matrizrotura,dim=2)
    do i = 1,NPrefino+1
	    write (6,'(200(F15.7,x))') (Matrizrotura(i,k),k=1,dimmatrizrotura)
    end do
end if

!RECESION VALUES FOR EACH POINT AT EACH FAILURE EVENT
if (dimmatrizrecesion.EQ.size(Matrizrecesion,dim=2)) then
    do i = 1,NPrefino+2
	    write (7,'(200(F15.7,x))') (Matrizrecesion(i,k),k=1,dimmatrizrecesion+1)
    end do
else
    dimmatrizrecesion = size(Matrizrecesion,dim=2)
    do i = 1,NPrefino+2
	    write (7,'(200(F15.7,x))') (Matrizrecesion(i,k),k=1,dimmatrizrecesion)
    end do
end if

!GROUNDWATER VALUES FOR EACH FAILURE EVENT
if (dimmatrizHwdraw.EQ.size(MatrizHwdraw,dim=2)) then
    do i = 1,1000+1
	    write (8,'(200(F15.7,x))') (MatrizHwdraw(i,k),k=1,dimmatrizHwdraw+1)
    end do
else
    dimmatrizHwdraw = size(MatrizHwdraw,dim=2)
    do i = 1,1000+1
	    write (8,'(200(F15.7,x))') (MatrizHwdraw(i,k),k=1,dimmatrizHwdraw)
    end do
end if

stop
end