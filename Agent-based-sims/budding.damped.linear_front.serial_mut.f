      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Code compiled with intel compiler (ifort)
      !!
      !!
      !!  !------------------------------------------------------------!
      !!  !---   Completely overdamped MD for dimers with growth    ---!
      !!  !------------------------------------------------------------!
      !!
      !!  BC's = PBC in x
      !!
      !!  Cells (1) grow in growth layer, (2) are pushed in 
      !!    propagation layer, (3) are held fixed in boundary 
      !!    layer, & (4) are removed beyond boundary layer
      !!
      !!  Depths defined as distance to closest cell in front
      !!  
      !!  Options - Restart: T/F = do don't use restart file
      !!              Movie: T/F = do/don't output movie
      !!             Bottom: T/F = keep/discard bottom
      !!
      !!  F = b*m*dr/dt (m=1 implicit)   
      !!  T = b*I*dth/dt (I=inertia, th=orientation angle)   
      !!
      !!
      !!  !------------------------------------------------------------!
      !!  !-------               Mutation algorithm             -------!
      !!  !------------------------------------------------------------!
      !!
      !!  Mutations introducted serially to propagating front
      !! 
      !!  Mutations stopped when (1) mutant dies away (nmuts=0)
      !!    or takes over from (nmuts=nlayer)
      !!  
      !!  Two retart files for `burned in' config (file2) & 
      !!    for mutation config (file3)
      !!
      !!
      !!  !------------------------------------------------------------!
      !!  !-------               Inputs and outputs             -------!
      !!  !------------------------------------------------------------!
      !!
      !!  Input parameters:
      !!      Ratio of cell length to width at birth (alpha0)
      !!      Ratio of cell length to width at division (alphamax)
      !!      Length of box (Lx)
      !!      Range of attractive force (att)
      !!      Growth rate difference of mutant to WT cells (s)
      !!      Growth rate of WT cells (rateWT)
      !!      Friction coefficient for overdamped dynamics (b)
      !!      # steps to initialize before starting mutations (steps_burn)
      !!      # steps in between mutations (steps_decorr)
      !!      # mutations to introduce serially (trials)
      !!      # steps in between calulating distance to front (layerskip)
      !!      # steps in between outputting data (dataskip)
      !!      # steps in between outputting production file (prodskip)
      !!      # steps in between outputting restart file (restskip)
      !!      Time step (dt)
      !!      Numerical parameter for finding distance to front (layerwidth)
      !!      Depth of region where cells actively grow (layerdepth)
      !!      Depth of region where cells are pushed by forces (propdepth)
      !!      Depth of region where cell positions are fixed (bounddepth)
      !!      Noise on growth rate to desynronize cell cycles (desync)
      !!      Random seed (seed0)
      !!      Movie file (file1)
      !!      Logical parameter to output movie  (movie)
      !!      Logical parameter to initialize with restart file (restart)
      !!      Logical parameter to include bottom of colony (bottom)
      !!
      !!  Output files: file1: movie file
      !!                file2: restart file for burned config
      !!                file3: restart file for mutants 
      !!                file4: mutant statistics
      !!                file5: mutant config
      !!                file7: # cells in various layers
      !!                file8: width of mutant clone in growth layer
      !!                file9: binned front-line of colony
      !!                      
      !!
      !!     Carl Schreck
      !!     11/28/2022
      !!
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program main

      implicit none
      integer Ntot,Ntot2,LMAX
      parameter(Ntot=2**16)
      parameter(Ntot2=2**24)
      parameter(LMAX=5000)
      double precision x(Ntot),vx(Ntot),ax(Ntot),bx(Ntot),fx(Ntot),s
      double precision y(Ntot),vy(Ntot),ay(Ntot),by(Ntot),fy(Ntot)
      double precision th(Ntot),vth(Ntot),ath(Ntot),bth(Ntot),fth(Ntot)
      double precision xp(Ntot),yp(Ntot),D(Ntot),alpha(Ntot),rate(Ntot)
      double precision inert(Ntot),depth(Ntot),tdiv,alpha0,rate0(Ntot)
      double precision att,Lx,exp,desync,kinetic,KE,V,ran2,cc,ss,corr,b
      double precision layerwidth,layerdepth,tmp,dr(2),dk(2),xa(2),ya(2)
      double precision ddsq,dd,maxdis,width,propdepth,bounddepth,rateWT
      double precision rateMUT,dt,propdist,bounddist,alphamax,dB(Ntot)
      double precision xB(Ntot),yB(Ntot),thB(Ntot),vxB(Ntot),vyB(Ntot)
      double precision vthB(Ntot),axB(Ntot),ayB(Ntot),bxB(Ntot),offset
      double precision byB(Ntot),bthB(Ntot),depthB(Ntot),rateB(Ntot)
      double precision rate0B(Ntot),alphaB(Ntot),athB(Ntot),xrem(Ntot2)
      double precision yrem(Ntot2),threm(Ntot2),alpharem(Ntot2),term
      double precision drem(Ntot2),amutstot,mutdepth,rf(2),rb(2),edge
      double precision frontline(LMAX),frontcelltype(LMAX),frontwidth
      integer N,seed,steps,i,j,k,countn,nl(12*Ntot,2),restartexist,div
      integer dataskip,prodskip,layerskip,restskip,trials,steps_burn
      integer forcelist(Ntot),proplist(Ntot),nrem,nprop,nsum,seedstart
      integer steps_decorr,kmut,nmuts,initmut,ksum,imut,divcells(Ntot)
      integer nlayer,proplistB(Ntot),NB,seedB,seed0,ksumB,celltype(Ntot)
      integer celltypeB(Ntot),nmutstot,nsumB,kmut_start,k_start,nlayerB
      integer maxbins,bin
      logical restart,movie,bottom,restartmuts
      character file1*80,file2*80,file3*80,file4*80
      character file5*80,file7*80,file8*80,file9*80
      common /f1com/ exp,alpha
      common /f2com/ nl,countn 
      common /f3com/ proplist
      common /f4com/ bottom
      common /f5com/ alphamax,alpha0

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) alphamax
      read(*,*) Lx
      read(*,*) att
      read(*,*) s

      ! read rates
      read(*,*) rateWT
      read(*,*) b

      ! read steps
      read(*,*) steps_burn
      read(*,*) steps_decorr
      read(*,*) trials

      ! read skip steps
      read(*,*) layerskip
      read(*,*) dataskip
      read(*,*) prodskip
      read(*,*) restskip
      read(*,*) dt
 
      ! read growth layer parameters
      read(*,*) layerwidth      
      read(*,*) layerdepth

      ! read layer parameters for force calc
      read(*,*) propdepth
      read(*,*) bounddepth

      ! read run parameters
      read(*,*) desync
      read(*,*) seed0
      
      ! read output files
      read(*,*) file1
      read(*,*) file2
      read(*,*) file3
      read(*,*) file4
      read(*,*) file5
      read(*,*) file7
      read(*,*) file8
      read(*,*) file9

      ! read options
      read(*,*) movie
      read(*,*) restart
      read(*,*) bottom

      ! parameters
      exp=2d0     ! 2=LS, 2.5=Hertzian, >2.9=RLJ
      width=0.2d0 ! width of neighborlist 
      term=0.75d0 ! terminate when fraction of growth layer is mutants
      edge=1d0    ! distance from edge (front) to qualify as edge

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+dt*rateWT)*dt

      ! growth rate of mutatnts
      rateMUT=(1d0+s)*rateWT

      ! distances of propagation/boundary layer from front
      propdist=layerdepth+propdepth
      bounddist=layerdepth+propdepth+bounddepth

      ! initialize system from scratch or restart file
      call initialize(file1,file2,file3,file4,file5,file7,file8,file9,
     +    restart,movie,Lx,b,att,desync,rateWT,seed0,kmut_start,k_start,
     +    NB,seedB,dB,xB,yB,thB,vxB,vyB,vthB,axB,ayB,athB,bxB,byB,bthB,
     +    rateB,rate0B,alphaB,depthB,proplistB,celltypeB,nsumB,ksumB,N,
     +    seed,d,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,rate,rate0,alpha,
     +    depth,proplist,celltype,nsum,ksum,nmuts,nlayer,nmutstot,
     +    amutstot,initmut,drem,xrem,yrem,threm,alpharem,mutdepth)

      ! calculate inertia of each cell
      call calc_inert(N,inert,D)
          
      ! set restart values
      if(kmut_start.eq.0) then
         restartmuts=.FALSE.
         kmut_start=1
      else
         restartmuts=.TRUE.
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!        burn in front        !!!!!!!!!!!!!!!!!
      if(.not.restartmuts) then
         do k=1,steps_burn
            ksum=ksum+1

            ! grow/divide cells
            call grow(dt,N,nsum,depth,layerdepth,rate,rate0,Lx,
     +           width,att,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           xp,yp,seed,desync,div,divcells,celltype)

            ! calc propagation list
            call calc_proplist(N,nprop,depth,proplist,
     +           vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
            
            ! remove cells & make neighbor list
            call checklist(N,x,y,xp,yp,maxdis)
            if(maxdis.gt.width*d(1)) then
               call remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +              alpha,depth,rate,rate0,proplist,bounddist,celltype) 
               call makelist(N,x,y,d,Lx,xp,yp,width,att)
            endif

            ! calculate inertia of each cell
            call calc_inert(N,inert,D)
            
            ! Gear precictor-corrector
            call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
            call force(N,x,y,th,d,V,fx,fy,fth,Lx,att)           
            call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           fx,fy,fth,inert,b)
            KE=kinetic(N,vx,vy,vth,inert)

            ! calc distance to front     
            if(mod(k,layerskip).eq.0) then
               call calc_depth(N,x,y,d,Lx,layerwidth,depth)
            endif
            
            ! calc propagation list
            call calc_nlayer(N,nlayer,nmuts,layerdepth,depth,celltype) 
 
            ! output data to screen & widths
            if(mod(k,dataskip).eq.0) then
               write(7,'(8I,2ES16.7)')k,0,k,nsum,N,nprop,
     +              nlayer,0,V/dble(nprop),KE/dble(nprop)
            endif         

            ! save config
            if(movie.and.mod(k,prodskip).eq.0) then
               write(1,*) 2*N
               do i=1,N
                 cc=dcos(th(i))
                 ss=dsin(th(i))
                 dd=alpha(i)-1d0
                 dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
                 dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
                 do j=1,2
                    xa(j)=x(i)+dr(j)*cc
                    ya(j)=y(i)+dr(j)*ss
                 enddo
                 dk(1)=D(i)
                 dk(2)=dd*D(i)
                 write(1,'(4F,I)')xa(1),ya(1),dk(1),depth(i),celltype(i)
                 write(1,'(4F,I)')xa(2),ya(2),dk(2),depth(i),celltype(i)
               enddo
               flush(1)
            endif         
         enddo

         ! burn in config
         NB=N
         seedB=seed
         ksumB=ksum
         nsumB=nsum
         nlayerB=nlayer
         do i=1,N
            dB(i)=d(i)
            xB(i)=x(i)
            yB(i)=y(i)
            thB(i)=th(i)
            vxB(i)=vx(i)
            vyB(i)=vy(i)
            vthB(i)=vth(i)
            axB(i)=ax(i)
            ayB(i)=ay(i)
            athB(i)=ath(i)
            bxB(i)=bx(i)
            byB(i)=by(i)
            bthB(i)=bth(i)                     
            rateB(i)=rate(i)                   
            rate0B(i)=rate0(i)
            alphaB(i)=alpha(i)
            depthB(i)=depth(i)
            proplistB(i)=proplist(i)
            celltypeB(i)=celltype(i)
         enddo
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!    loop over mutations      !!!!!!!!!!!!!!!!!
      do kmut=kmut_start,trials
         if(restartmuts) then
            k=k_start
            restartmuts=.FALSE.
         else
            k=0
            nmuts=0
            nmutstot=0
            amutstot=0d0
            initmut=0

            ! start with burned-in configuration
            N=NB
            ksum=ksumB
            nsum=nsumB
            nlayer=nlayerB
            do i=1,N
               d(i)=dB(i)
               x(i)=xB(i)
               y(i)=yB(i)
               th(i)=thB(i)
               vx(i)=vxB(i)
               vy(i)=vyB(i)
               vth(i)=vthB(i)
               ax(i)=axB(i)
               ay(i)=ayB(i)
               ath(i)=athB(i)
               bx(i)=bxB(i)
               by(i)=byB(i)
               bth(i)=bthB(i)                     
               rate(i)=rateB(i)                   
               rate0(i)=rate0B(i)
               alpha(i)=alphaB(i)
               depth(i)=depthB(i)
               proplist(i)=proplistB(i)
               celltype(i)=celltypeB(i)
            enddo

            ! burn seeds
            seed=seed0
            do while (seed.ne.seedB)
               tmp=ran2(seed)
            enddo
         endif

         ! do k=kstart_mut+1,steps_mut
         do while((nmuts.gt.0.and.nmuts.lt.term*nlayer).or.initmut.eq.0)
            k=k+1
            ksum=ksum+1

            ! grow/divide cells
            call grow(dt,N,nsum,depth,layerdepth,rate,rate0,Lx,
     +           width,att,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           xp,yp,seed,desync,div,divcells,celltype)

            ! mutate cells
            if(initmut.eq.0.and.k.gt.steps_decorr) then
               if(div.gt.0) then
                  initmut=1        
       
                  ! choose random cell to mutate
                  imut=divcells(ceiling(div*ran2(seed)))      
                  rate0(imut)=ratemut
                  rate(imut)=rate0(imut)*
     +                 (1d0+(ran2(seed)-0.5d0)*desync)
                  celltype(imut)=kmut

                  ! update depth of mutated cell
                  call calc_depth(N,x,y,d,Lx,layerwidth,depth)
                  mutdepth=depth(imut)                  
               endif
            endif

            ! calc propagation list
            call calc_proplist(N,nprop,depth,proplist,
     +           vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)

            ! remove cells & make neighbor list
            call checklist(N,x,y,xp,yp,maxdis)
            if(maxdis.gt.width*d(1)) then
               call removemut(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +              alpha,depth,rate,rate0,proplist,bounddist,celltype,
     +              nmutstot,amutstot,drem,xrem,yrem,threm,alpharem) 
               call makelist(N,x,y,d,Lx,xp,yp,width,att)
            endif

            ! calculate inertia of each cell
            call calc_inert(N,inert,D)
 
            ! Gear precictor-corrector
            call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
            call force(N,x,y,th,d,V,fx,fy,fth,Lx,att)           
            call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           fx,fy,fth,inert,b)
            KE=kinetic(N,vx,vy,vth,inert)

            ! calc distance to front     
            if(mod(k,layerskip).eq.0) then
               call calc_depth(N,x,y,d,Lx,layerwidth,depth)
            endif

            ! calc propagation list
            call calc_nlayer(N,nlayer,nmuts,layerdepth,depth,celltype)

            ! burn in configuration
            if(initmut.eq.0.and.k.eq.steps_decorr) then
               NB=N
               seedB=seed
               ksumB=ksum
               nsumB=nsum
               nlayerB=nlayer
               do i=1,N
                  dB(i)=d(i)
                  xB(i)=x(i)
                  yB(i)=y(i)
                  thB(i)=th(i)
                  vxB(i)=vx(i)
                  vyB(i)=vy(i)
                  vthB(i)=vth(i)
                  axB(i)=ax(i)
                  ayB(i)=ay(i)
                  athB(i)=ath(i)
                  bxB(i)=bx(i)
                  byB(i)=by(i)
                  bthB(i)=bth(i)                     
                  rateB(i)=rate(i)                   
                  rate0B(i)=rate0(i)
                  alphaB(i)=alpha(i)
                  depthB(i)=depth(i)
                  proplistB(i)=proplist(i)
                  celltypeB(i)=celltype(i)
               enddo

               ! output burned config to restart file
               open(unit=2,file=TRIM(file2))
               write(2,*) N, seed
               write(2,*) ksum, nsum
               do i=1,N
                  write(2,'(17E26.18,2I)') xB(i),yB(i),thB(i),
     +                 vxB(i),vyB(i),vthB(i),axB(i),ayB(i),athB(i),
     +                 bxB(i),byB(i),bthB(i),dB(i),alphaB(i),depthB(i),
     +                 rateB(i),rate0B(i),proplistB(i),celltypeB(i)
               enddo
               flush(2)
               close(2)
            endif

            ! output data to screen & widths
            if(mod(k,dataskip).eq.0) then
              write(7,'(8I,2ES16.7)')ksum,kmut,k,nsum,N,nprop,
     +              nlayer,nmuts,V/dble(nprop),KE/dble(nprop)
              flush(7)

               call calc_wid(N,x,y,depth,layerdepth,celltype,edge,rf,rb)
               write(8,'(I,4ES)') kmut, rf(1), rf(2), rb(1), rb(2)
               flush(8)

               frontwidth=Lx/dnint(Lx/3.2d0)
               call calc_front_line(N,Lx,x,y,d,depth,frontwidth,
     +              celltype,frontline,frontcelltype)
               maxbins=nint(Lx/frontwidth)
               offset=(1d0+dble(maxbins))*frontwidth/2d0
               do bin=1,maxbins
                  write(9,'(I,3ES)') kmut, frontwidth*dble(bin)-offset, 
     +                 frontline(bin), frontcelltype(bin)
               enddo
               flush(9)
            endif         

            ! save config
            if(movie.and.mod(k,prodskip).eq.0) then
               write(1,*) 2*N
               do i=1,N
                  cc=dcos(th(i))
                  ss=dsin(th(i))
                  dd=alpha(i)-1d0
                  dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
                  dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
                  do j=1,2
                     xa(j)=x(i)+dr(j)*cc
                     ya(j)=y(i)+dr(j)*ss
                  enddo
                  dk(1)=D(i)
                  dk(2)=dd*D(i)
                 write(1,'(4F,I)')xa(1),ya(1),dk(1),depth(i),celltype(i)
                 write(1,'(4F,I)')xa(2),ya(2),dk(2),depth(i),celltype(i)
               enddo
               flush(1)
            endif         

            ! save mutant restart file
            if(restart.and.mod(k,restskip).eq.0) then
               open(unit=3,file=TRIM(file3))
               write(3,*) N, seed               
               write(3,*) kmut, k
               write(3,*) ksum, nsum
               write(3,*) nlayer, initmut, nmuts
               write(3,*) mutdepth, nmutstot, amutstot
               do i=1,N
                  write(3,'(17E26.18,2I)') x(i),y(i),th(i),
     +                 vx(i),vy(i),vth(i),ax(i),ay(i),ath(i),
     +                 bx(i),by(i),bth(i),d(i),alpha(i),depth(i),
     +                 rate(i),rate0(i),proplist(i),celltype(i)
               enddo
               do i=1,nmutstot
                  write(3,'(5E26.18)') xrem(i),yrem(i),
     +                 threm(i),drem(i),alpharem(i)
               enddo
               flush(3)
               close(3)
            endif
         enddo

         ! collect the remainder of the mutated cells
         do i=1,N
            if(celltype(i).ne.0) then
               nmutstot=nmutstot+1
               amutstot=amutstot+1d0+(alpha(i)-1d0)**2
               drem(nmutstot)=d(i)
               xrem(nmutstot)=x(i)
               yrem(nmutstot)=y(i)
               threm(nmutstot)=th(i)
               alpharem(nmutstot)=alpha(i)
            endif
         enddo
 
         ! output file to restart new mutant run
         open(unit=3,file=TRIM(file3))
         write(3,*) NB, seedB
         write(3,*) kmut+1, 0
         write(3,*) ksumB, nsumB
         write(3,*) nlayerB, 0, 0
         write(3,*) 0d0, 0, 0d0
         do i=1,NB
            write(3,'(17E26.18,2I)') xB(i),yB(i),thB(i),
     +           vxB(i),vyB(i),vthB(i),axB(i),ayB(i),athB(i),
     +           bxB(i),byB(i),bthB(i),dB(i),alphaB(i),depthB(i),
     +           rateB(i),rate0B(i),proplistB(i),celltypeB(i)
         enddo
         flush(3)
         close(3)

         ! output mutant config
         write(4,'(I,F,3I,F)') kmut, mutdepth, 
     +        nmuts, nlayer, nmutstot, amutstot
         flush(4)
 
         ! output config file containing mutant positions
         write(5,*) 2*nmutstot
         do i=1,nmutstot
            cc=dcos(threm(i))
            ss=dsin(threm(i))
            dd=alpharem(i)-1d0
            dk(1)=drem(i)
            dk(2)=dd*drem(i)
            dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*drem(i)/2d0
            dr(2)=-(1d0+dd)/(1d0+dd**2)*drem(i)/2d0
            do j=1,2
               xa(j)=xrem(i)+dr(j)*cc
               ya(j)=yrem(i)+dr(j)*ss
               write(5,'(3F)') xa(j), ya(j), dk(j)
            enddo
         enddo
         flush(5)
      enddo

      end ! end main


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!  initialize cell positions & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine initialize(file1,file2,file3,file4,file5,file7,file8,
     +     file9,restart,movie,Lx,b,att,desync,rateWT,seed0,kmut_start,
     +     k_start,NB,seedB,dB,xB,yB,thB,vxB,vyB,vthB,axB,ayB,athB,
     +     bxB,byB,bthB,rateB,rate0B,alphaB,depthB,proplistB,celltypeB,
     +     nsumB,ksumB,N,seed,d,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     rate,rate0,alpha,depth,proplist,celltype,nsum,ksum,nmuts,
     +     nlayer,nmutstot,amutstot,initmut,drem,xrem,yrem,threm,
     +     alpharem,mutdepth)

      parameter(Ntot=2**16,Ntot2=2**24)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision b,Lx,x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),depth(Ntot),d(Ntot),alpha(Ntot),V
      double precision rate(Ntot),exp,att,dd,ddsq,rateWT,ran2,tmp
      double precision rate0(Ntot),alpha0,alphamax,desync,xB(Ntot)
      double precision yB(Ntot),thB(Ntot),vxB(Ntot),vyB(Ntot),vthB(Ntot)
      double precision axB(Ntot),ayB(Ntot),athB(Ntot),bxB(Ntot),dB(Ntot)
      double precision byB(Ntot),bthB(Ntot),depthB(Ntot),alphaB(Ntot)
      double precision rateB(Ntot),rate0B(Ntot),amutstot,drem(Ntot2)
      double precision xrem(Ntot2),yrem(Ntot2),threm(Ntot2)
      double precision alpharem(Ntot2),mutdepth
      integer N,NB,seed,seedB,seedstart,i,proplist(Ntot),celltype(Ntot)
      integer proplistB(Ntot),celltypeB(Ntot),kmut_start,k_start,ksum
      integer nmuts,nmutstot,initmut,nlayer,ksumB,nsumB,nsum,seed0
      character file1*80,file2*80,file3*80,file4*80
      character file5*80,file7*80,file8*80,file9*80
      logical restart,movie,restartex2,restartex3
      common /f2com/ nl,countn 
      common /f4com/ bottom
      common /f5com/ alphamax,alpha0

      ! check if restart file exists
      inquire(file=file2,exist=restartex2)
      inquire(file=file3,exist=restartex3)

      if(restart.and.restartex2.and.restartex3) then  
         ! open files
         open(unit=2,file=TRIM(file2))
         open(unit=3,file=TRIM(file3))
         open(unit=4,file=TRIM(file4),ACCESS="APPEND")
         open(unit=5,file=TRIM(file5),ACCESS="APPEND")
         open(unit=7,file=TRIM(file7),ACCESS="APPEND")
         open(unit=8,file=TRIM(file8),ACCESS="APPEND")
         open(unit=9,file=TRIM(file9),ACCESS="APPEND")
         if(movie) open(unit=1,file=TRIM(file1),ACCESS="APPEND")

         ! read restart file
         read(2,*) NB, seedB
         read(2,*) ksumB, nsumB
         do i=1,NB
            read(2,*) xB(i),yB(i),thB(i),vxB(i),vyB(i),vthB(i),axB(i),
     +           ayB(i),athB(i),bxB(i),byB(i),bthB(i),dB(i),alphaB(i),
     +           depthB(i),rateB(i),rate0B(i),proplistB(i),celltypeB(i)
         enddo
         close(2)

         ! read restart file
         read(3,*) N, seedstart
         read(3,*) kmut_start, k_start
         read(3,*) ksum, nsum
         read(3,*) nlayer, initmut, nmuts
         read(3,*) mutdepth, nmutstot, amutstot
         do i=1,N
            read(3,*) x(i),y(i),th(i),vx(i),vy(i),vth(i),ax(i),
     +           ay(i),ath(i),bx(i),by(i),bth(i),d(i),alpha(i),
     +           depth(i),rate(i),rate0(i),proplist(i),celltype(i)
         enddo
         do i=1,nmutstot
            read(3,*) xrem(i),yrem(i),threm(i),drem(i),alpharem(i)
         enddo
         close(3)

         ! burn seeds
         seed=seed0
         do while (seed.ne.seedstart)
            tmp=ran2(seed)
         enddo
      else ! no restart file exists
         N=floor(Lx)
         seed=seed0

         ! set initial # of cells & mutations 
         nlayer=N
         nsum=N
         initmut=0
         nmuts=0
         nmutstot=0
         amutstot=0d0
         mutdepth=0d0

         ! set initial indices
         kmut_start=0
         k_start=0
         ksum=0

         ! open files
         open(unit=4,file=TRIM(file4))
         open(unit=5,file=TRIM(file5))
         open(unit=7,file=TRIM(file7))
         open(unit=8,file=TRIM(file8))
         open(unit=9,file=TRIM(file9))
         if(movie) open(unit=1,file=TRIM(file1))

         ! random initial config
         do i=1,N     
            d(i)=1d0
            x(i)=dble(i)-dble(N+1)/2d0
            y(i)=0d0
            th(i)=(ran2(seed)-0.5d0)*2d0*pi
            depth(i)=0d0
            rate0(i)=rateWT
            proplist(i)=1
            celltype(i)=0
         enddo
     
         ! assign initial aspect ratios & rates
         do i=1,N     
c            alpha(i)=alpha0*(1d0+dble(i-1)/2d0)
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
         enddo
         
         ! assign initial aspect ratios & rates
         do i=1,N     
            alpha(i)=alpha0*(1d0+ran2(seed))
         enddo         
         
         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         call makelist(N,x,y,d,Lx,xp,yp,width,att)
         call force(N,x,y,th,d,V,fx,fy,fth,Lx,att)  
         do i=1,N
            vx(i)=b*fx(i)
            vy(i)=b*fy(i)
            vth(i)=b*fth(i)/inert(i)
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0         
         enddo
      endif

      return
      end ! end initialize



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!          grow & divide cells         !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine grow(dt,N,nsum,depth,layerdepth,rate,rate0,Lx,
     +     width,att,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     xp,yp,seed,desync,div,divcells,celltype)
 
      parameter(Ntot=2**16)
      double precision exp,Lx,x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot),dt
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),att,ran2,rate0(Ntot),alpha0,alphamax
      double precision desync,layerdepth,growfrac
      integer N,nsum,seed,i,div,divcells(Ntot),celltype(Ntot)
      common /f1com/ exp,alpha
      common /f5com/ alphamax,alpha0

      div=0
      do i=1,N
         ! grow cell i
         if(depth(i).lt.layerdepth) then
            growfrac=1d0+dt*rate(i)
            alpha(i)=1d0+dsqrt((1d0+(alpha(i)-1d0)**2)*growfrac-1d0)
         endif

         ! divide cell i
         if(alpha(i).gt.alphamax) then
            divcells(div+1)=i
            divcells(div+2)=N
            div=div+2

            ! divide into 2 - N=current cels, nsum=total 
            N=N+1
            nsum=nsum+1

            ! divide into 2 - 1st assigned index N+1
            D(N)=D(i)
            x(N)=x(i)+alpha0/2d0*dcos(th(i))
            y(N)=y(i)+alpha0/2d0*dsin(th(i))
            th(N)=th(i)
            rate0(N)=rate0(i)
            rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(N) 
            celltype(N)=celltype(i)
            alpha(N)=alpha0
            vx(N)=vx(i)
            vy(N)=vy(i)
            vth(N)=vth(i)               
            ax(N)=ax(i)
            ay(N)=ay(i)
            ath(N)=ath(i)               
            bx(N)=bx(i)
            by(N)=by(i)
            bth(N)=bth(i)

            ! divide into 2 - 1st assigned index i
            x(i)=x(i)-alpha0/2d0*dcos(th(i))
            y(i)=y(i)-alpha0/2d0*dsin(th(i))
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
            alpha(i)=alpha0

            ! update neighbor list
            call makelistind(N,N,x,y,d,Lx,xp,yp,width,att)

            ! update depth
            depth(N)=max(0d0,depth(i)-alpha0/2d0*dsin(th(i)))
            depth(i)=max(0d0,depth(i)+alpha0/2d0*dsin(th(i)))
         endif
      enddo

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!  remove cells that fall behind propagation layer   !!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +     alpha,depth,rate,rate0,proplist,bounddist,celltype)

      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision rate0(Ntot),alpha(Ntot),bounddist
      integer N,nrem,i,j,proplist(Ntot),celltype(Ntot)

      nrem=0
      do i=N,1,-1
         if(depth(i).gt.bounddist) then
            nrem=nrem+1
            do j=i+1,N
               x(j-1)=x(j)
               y(j-1)=y(j)
               th(j-1)=th(j)
               vx(j-1)=vx(j)
               vy(j-1)=vy(j)
               vth(j-1)=vth(j)
               ax(j-1)=ax(j)
               ay(j-1)=ay(j)
               ath(j-1)=ath(j)
               bx(j-1)=bx(j)
               by(j-1)=by(j)
               bth(j-1)=bth(j)
               d(j-1)=d(j)
               alpha(j-1)=alpha(j)
               depth(j-1)=depth(j)
               rate(j-1)=rate(j)
               rate0(j-1)=rate0(j)
               proplist(j-1)=proplist(j)                     
               celltype(j-1)=celltype(j)
            enddo
         endif
      enddo
      N=N-nrem
         
      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!  remove cells, save mutants   !!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine removemut(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +     alpha,depth,rate,rate0,proplist,bounddist,celltype,
     +     nmutstot,amutstot,drem,xrem,yrem,threm,alpharem)

      parameter(Ntot=2**16,Ntot2=2**24)
      double precision rate(Ntot),rate0(Ntot),x(Ntot),y(Ntot),th(Ntot)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),depth(Ntot)
      double precision d(Ntot),alpha(Ntot),alpharem(Ntot2),drem(Ntot2)
      double precision threm(Ntot2),xrem(Ntot2),yrem(Ntot2)
      double precision amutstot,bounddist
      integer N,nrem,i,j,proplist(Ntot),celltype(Ntot),nmutstot

      nrem=0
      do i=N,1,-1
         if(depth(i).gt.bounddist) then
            nrem=nrem+1
            if(celltype(i).ne.0) then
               nmutstot=nmutstot+1
               amutstot=amutstot+1d0+(alpha(i)-1d0)**2
               drem(nmutstot)=d(i)
               xrem(nmutstot)=x(i)
               yrem(nmutstot)=y(i)
               threm(nmutstot)=th(i)
               alpharem(nmutstot)=alpha(i)
            endif

            do j=i+1,N
               x(j-1)=x(j)
               y(j-1)=y(j)
               th(j-1)=th(j)
               vx(j-1)=vx(j)
               vy(j-1)=vy(j)
               vth(j-1)=vth(j)
               ax(j-1)=ax(j)
               ay(j-1)=ay(j)
               ath(j-1)=ath(j)
               bx(j-1)=bx(j)
               by(j-1)=by(j)
               bth(j-1)=bth(j)
               d(j-1)=d(j)
               alpha(j-1)=alpha(j)
               depth(j-1)=depth(j)
               rate(j-1)=rate(j)
               rate0(j-1)=rate0(j)
               proplist(j-1)=proplist(j)                     
               celltype(j-1)=celltype(j)
            enddo
         endif
      enddo
      N=N-nrem
         
      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!  calc propagation list  !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
 
      parameter(Ntot=2**16)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),depth(Ntot)
      double precision propdist
      integer N,nprop,i,proplist(Ntot)

      nprop=0
      do i=1,N
         if(depth(i).lt.propdist) then
            proplist(i)=1
            nprop=nprop+1
         else
            proplist(i)=0
            vx(i)=0d0
            vy(i)=0d0
            vth(i)=0d0
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0
         endif
      enddo

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  calc number of all cells & mutants in growth layer  !!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_nlayer(N,nlayer,nmuts,layerdepth,depth,celltype)
 
      parameter(Ntot=2**16)
      double precision depth(Ntot),layerdepth
      integer N,nlayer,nmuts,i,celltype(Ntot)

      nlayer=0
      nmuts=0
      do i=1,N
         if(depth(i).lt.layerdepth) then
            nlayer=nlayer+1
            if(celltype(i).gt.0) then
               nmuts=nmuts+1
            endif
         endif
      enddo

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_inert(N,inert,D)
      parameter(Ntot=2**16)
      double precision inert(Ntot),D(Ntot),alpha(Ntot),dd,ddsq,exp
      integer i,N
      common /f1com/ exp,alpha
      
      do i=1,N
         dd=alpha(i)-1d0
         ddsq=dd*dd
         inert(i)=((1d0+ddsq**2)/(1d0+ddsq)+2d0*ddsq*
     +        (1d0+dd)**2/(1d0+ddsq)**2)*d(i)**2/8d0
      enddo

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine checklist(N,x,y,xp,yp,maxdis)
      parameter(Ntot=2**16)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),df
      integer N

      df=2d0

      maxdis=0d0
      do i=1,N
	maxdis=max(dabs(x(i)-xp(i)),maxdis)
	maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(df*maxdis*maxdis)

      return
      end ! end checklist


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   make neighbor list   !!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelist(N,x,y,d,Lx,xp,yp,width,att)
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),alphamax
      double precision alpha0,att,Lx,Ly,xij,yij,dijlist
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      countn=0      
      do i=2,N
         do j=1,i-1
            xij=x(i)-x(j)
            xij=xij-dnint(xij/Lx)*Lx           
            dij=alphamax*d(i) ! max distance - aspect ratio = 2
            dijlist=dij+(width+att)*d(1)
            if(dabs(xij).lt.dijlist) then
               yij=y(i)-y(j)
               rijsq=xij*xij+yij*yij
               if(rijsq.lt.dijlist**2) then
                  countn=countn+1
                  nl(countn,1)=i
                  nl(countn,2)=j
               endif
            endif
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo      
      
      return
      end ! end makelist
      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!   make neighbor list only for cell i   !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelistind(i,N,x,y,d,Lx,xp,yp,width,att)
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),att,Lx,Ly
      double precision xij,yij,dijlist,alphamax,alpha0
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      do j=1,i-1         
         xij=x(i)-x(j)
         xij=xij-dnint(xij/Lx)*Lx            
         dij=alphamax*d(i) ! max distance = aspect ratio
         dijlist=dij+(width+att)*d(1)
         if(dabs(xij).lt.dijlist) then
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            if(rijsq.lt.dijlist**2) then
               countn=countn+1
               nl(countn,1)=i
               nl(countn,2)=j
            end if
         endif
      enddo
      
      xp(i)=x(i)
      yp(i)=y(i)
      
      return
      end ! end makelistind
      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!    predicts new positions and velocities    !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)     
      parameter(Ntot=2**16)
      integer N,i,proplist(Ntot)
      double precision x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),dt,c1,c2,c3
      common /f3com/ proplist

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/3d0

      do i=1,N
         if(proplist(i).eq.1) then 
            x(i) = x(i) + c1*vx(i) + c2*ax(i) + c3*bx(i)
            y(i) = y(i) + c1*vy(i) + c2*ay(i) + c3*by(i)
            th(i) = th(i) + c1*vth(i) + c2*ath(i) + c3*bth(i)         
            vx(i) = vx(i) + c1*ax(i) + c2*bx(i)
            vy(i) = vy(i) + c1*ay(i) + c2*by(i)     
            vth(i) = vth(i) + c1*ath(i) + c2*bth(i)     
            ax(i) = ax(i) + c1*bx(i)
            ay(i) = ay(i) + c1*by(i)
            ath(i) = ath(i) + c1*bth(i)
         endif
      enddo

      end ! end prediction step



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   corrects prediction   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     fx,fy,fth,inert,b)
      parameter(Ntot=2**16)
      integer i,N,proplist(Ntot)
      double precision b,dt,x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),c1,c2,c3,cg0,cg2,cg3
      double precision gear0,gear2,gear3,corrx,corry,corrth
      common /f3com/ proplist

      gear0 = 3d0/8d0
      gear2 = 3d0/4d0
      gear3 = 1d0/6d0

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/2d0

      cg0 = gear0*c1
      cg2 = gear2*c1/c2
      cg3 = gear3*c1/c3

      do i=1,N
         if(proplist(i).eq.1) then 
            vxi = b*fx(i)
            vyi = b*fy(i)
            vthi = b*fth(i)/inert(i)
            corrx = vxi - vx(i)
            corry = vyi - vy(i)
            corrth = vthi - vth(i)
            x(i) = x(i) + cg0*corrx
            y(i) = y(i) + cg0*corry
            th(i) = th(i) + cg0*corrth        
            vx(i) = vxi
            vy(i) = vyi
            vth(i) = vthi
            ax(i) = ax(i) + cg2*corrx
            ay(i) = ay(i) + cg2*corry
            ath(i) = ath(i) + cg2*corrth
            bx(i) = bx(i) + cg3*corrx
            by(i) = by(i) + cg3*corry
            bth(i) = bth(i) + cg3*corrth
         endif
      enddo

      end ! end correction step

            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!           force            !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine force(N,x,y,th,d,V,fx,fy,fth,Lx,att) ! dimer force
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),D(Ntot),Lx
      double precision radi_up(Ntot),fx(Ntot),fy(Ntot),fth(Ntot),V,Vij
      double precision f_x,f_y,fc,fr,LJ,LJ0,exp,dij,rij,xij,yij,dij_up
      double precision fact,att,fthi,fthj,rijsq,c(Ntot),s(Ntot),dd,dd2
      double precision xa(Ntot,2),ya(Ntot,2),dk(Ntot,2),dr(Ntot,2),di1j1
      integer countn,nl(12*Ntot,2),N,ki,kj,jj,up,down,proplist(Ntot)
      common /f1com/ exp,alpha
      common /f2com/ nl,countn
      common /f3com/ proplist
      
      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
      enddo
      V=0d0
 
      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i))
         s(i)=dsin(th(i))
         dd=alpha(i)-1d0
         dd2=(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dr(i,1)=dd2*dd**2
         dr(i,2)=-dd2
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)
         dk(i,2)=dd*D(i)
         radi_up(i)=(dk(i,2)-2d0*dr(i,2))/2d0
      enddo

      ! inter-particle interactions      
      do k=1,countn
         i=nl(k,1)
         j=nl(k,2)         
         if(proplist(i).eq.1.or.proplist(j).eq.1) then 
            dij_up=radi_up(i)+radi_up(j)
            xij=x(i)-x(j)
            xij=xij-dnint(xij/Lx)*Lx
            if(dabs(xij).lt.dij_up) then 
               yij=y(i)-y(j)        
               rijsq=xij**2+yij**2
               if(rijsq.lt.dij_up*dij_up) then
                  di1j1=(dk(i,1)+dk(j,1))/2d0
                  do ki=1,2
                     do kj=1,2
                        dij=(dk(i,ki)+dk(j,kj))/2d0
                        xij=xa(i,ki)-xa(j,kj)
                        xij=xij-dnint(xij/Lx)*Lx
                        yij=ya(i,ki)-ya(j,kj)
                        rijsq=xij**2+yij**2
                        if(rijsq.lt.(dij+att)**2) then
                           rij=dsqrt(rijsq)
                           if(exp.eq.2d0) then
                              fc=(1d0-rij/dij)/dij     
                              Vij=(1d0-rij/dij)**2/exp
     +                             -(att/dij)**2/exp
                              fact=(dij/di1j1)**2
                           elseif(exp.lt.2.9) then
                              fc=(1d0-rij/dij)/dij     
                              Vij=(1d0-rij/dij)**2/exp
     +                             -(att/dij)**2/exp
                              fact=(dij/di1j1)**exp
                           else
                              LJ=(dij/rij)*(dij/rij)
                              LJ=LJ*LJ*LJ
                              LJ0=(dij/(dij+att))**6
                              fc=1d0/rij*LJ*(LJ-1d0)
                              Vij=(LJ-1d0)**2-(LJ0-1d0)**2
                              fact=(dij/di1j1)**2
                           endif                    
                           fr=fc/rij*fact
                           f_x=fr*xij
                           f_y=fr*yij
                           if(proplist(i).eq.1) then
                              fx(i)=fx(i)+f_x
                              fy(i)=fy(i)+f_y
                              fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)
                           endif
                           if(proplist(j).eq.1) then
                              fx(j)=fx(j)-f_x
                              fy(j)=fy(j)-f_y
                              fth(j)=fth(j)-dr(j,kj)*(c(j)*f_y-s(j)*f_x)
                           endif
                           V=V+Vij*fact                     
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
      enddo
 
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo         
         V=V/72d0
      endif
      
      return							
      end ! end force calc
     

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!    calc kinetic energy    !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      function kinetic(N,vx,vy,vth,inert)
      parameter(Ntot=2**16)
      integer i,N,proplist(Ntot)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),inert(Ntot),kinetic
      common /f3com/ proplist

      kinetic=0d0
      do i=1,N
         if(proplist(i).eq.1) then         
            kinetic=kinetic+vx(i)**2+vy(i)**2+inert(i)*vth(i)**2   
         endif
      enddo   
      kinetic=kinetic/2d0

      end ! end kinetic energy calc


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_depth(N,x,y,d,Lx,layerwidth,depth)
      parameter(Ntot=2**16,LMAX=5000,DEPTHMAX=1000)
      double precision x(Ntot),y(Ntot),d(Ntot),Lx,layerwidth,xi,xij,yij
      double precision ymin,ymax,depth(Ntot),depthD(Ntot),depthU(Ntot)
      double precision offset,dx,drsq
      integer N,numbins,i,j,bin,dbin,jj
      integer bini(Ntot),npart(LMAX),ibin(LMAX,DEPTHMAX)
      integer frontUind(Ntot),frontDind(Ntot)
      logical bottom
      common /f4com/ bottom

      numbins=nint(Lx/layerwidth)
      
      ! calc vertical depths
      do bin=1,numbins
         npart(bin)=0
      enddo         
      offset=floor(Lx/2d0/layerwidth)+1
      do i=1,N
         xi=x(i)-dnint(x(i)/Lx)*Lx
         bin=floor(xi/layerwidth)+offset
         bini(i)=bin              ! x bin of cell i
         npart(bin)=npart(bin)+1  ! # cells in bin
         ibin(bin,npart(bin))=i   ! ibin = cell index in bin      
      enddo
      do i=1,N
         ymin=0d0
         ymax=0d0
         do dbin=-1,1
            bin=mod(numbins+bini(i)+dbin-1,numbins)+1
            do jj=1,npart(bin)                  
               j=ibin(bin,jj)
               xij=x(i)-x(j)
               xij=xij-dnint(xij/Lx)*Lx
               if(dabs(xij)<layerwidth) then
                  if(y(j).gt.ymax) then
                     ymax=y(j)
                  elseif(y(j).lt.ymin) then
                     ymin=y(j)
                  endif
               endif
            enddo
         enddo
         depthU(i)=ymax-y(i)
         depthD(i)=y(i)-ymin 
      enddo

      ! assign cells near front to be at front
      do i=1,N
         if(depthU(i).lt.D(i)) then 
            frontUind(i)=1
            depthU(i)=0d0
         else
            frontUind(i)=0
         endif
      enddo
      if(bottom) then
         do i=1,N
            if(depthD(i).lt.D(i)) then
               frontDind(i)=1
               depthD(i)=0d0
            else
               frontDind(i)=0
            endif
         enddo
      endif

      ! calc new depths from nearest cell at front
      do i=1,N
         dbin=ceiling(depthU(i)/layerwidth)
         do bin=bini(i)-dbin,bini(i)+dbin
            do jj=1,npart(bin)
               j=ibin(bin,jj)
               if(frontUind(j).eq.1) then
                  dx=x(i)-x(j)
                  dx=dx-dnint(dx/Lx)*Lx
                  if(dabs(dx).lt.depthU(i)) then
                     dy=y(i)-y(j)
                     drsq=dx*dx+dy*dy
                     if(drsq.lt.depthU(i)**2) then
                        depthU(i)=dsqrt(drsq)
                     endif
                  endif
               endif
            enddo
         enddo
      enddo         
      if(bottom) then
         do i=1,N
            dbin=ceiling(depthU(i)/layerwidth)
            do bin=bini(i)-dbin,bini(i)+dbin
               do jj=1,npart(bin)
                  j=ibin(bin,jj)
                  if(frontUind(j).eq.1) then
                     dx=x(i)-x(j)
                     dx=dx-dnint(dx/Lx)*Lx
                     if(dabs(dx).lt.depthD(i)) then
                        dy=y(i)-y(j)
                        drsq=dx*dx+dy*dy
                        if(drsq.lt.depthd(i)**2) then
                           depthD(i)=dsqrt(drsq)
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      ! assign depth to closest boundary
      if(bottom) then
         do i=1,N
            depth(i)=min(depthU(i),depthD(i))  
         enddo         
      else
         do i=1,N
            depth(i)=depthU(i)
         enddo
      endif

      return
      end ! end depth calc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!    calc width of mutant sector in boundary layer    !!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_wid(N,x,y,depth,layerdepth,celltype,edge,rf,rb)

      parameter(Ntot=2**16,LMAX=5000,DEPTHMAX=1000)
      double precision x(Ntot),y(Ntot),layerdepth,rf(2),rb(2)
      double precision xminF,xmaxF,xminB,xmaxB,depth(Ntot),edge
      integer N,i,cellsF,cellsB,celltype(Ntot)

      cellsF=0
      cellsB=0
      xminF=0d0
      xmaxF=0d0
      xminB=0d0
      xmaxB=0d0
      xminFy=0d0
      xmaxFy=0d0
      xminBy=0d0
      xmaxBy=0d0
      do i=1,N
         if(celltype(i).gt.0) then
            if(depth(i).lt.edge) then
               if(cellsF.eq.0) then
                  cellsF=1
                  xminF=x(i)
                  xmaxF=x(i)
                  xminFy=y(i)
                  xmaxFy=y(i)
               else
                  if(x(i).lt.xminF) then
                     xminF=x(i)
                     xminFy=y(i)
                  else if(x(i).gt.xmaxF) then
                     xmaxF=x(i)
                     xmaxFy=y(i)
                  endif
               endif
            else if(layerdepth-depth(i).lt.edge) then
               if(cellsB.eq.0) then
                  cellsB=1
                  xminB=x(i)
                  xmaxB=x(i)
                  xminBy=y(i)
                  xmaxBy=y(i)
               else
                  if(x(i).lt.xminB) then
                     xminB=x(i)
                     xminBy=y(i)
                  else if(x(i).gt.xmaxB) then
                     xmaxB=x(i)
                     xmaxBy=y(i)
                  endif
               endif
            endif
         endif
      enddo

      if(cellsF.eq.0) then
         rF(1)=0d0
         rF(2)=0d0
      else
         rF(1)=(xmaxFy+xminFy)/2d0
         rF(2)=xmaxF-xminF+1d0
      endif

      if(cellsB.eq.0) then
         rB(1)=0d0
         rB(2)=0d0
      else
         rB(1)=(xmaxBy+xminBy)/2d0
         rB(2)=xmaxB-xminB+1d0
      endif

      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!            calc binned front-line            !!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_front_line(N,Lx,x,y,d,depth,frontwidth,
     +     celltype,frontline,frontcelltype)

      parameter(Ntot=2**16,LMAX=5000)
      double precision Lx,xbox,frontwidth,x(Ntot),y(Ntot),d(Ntot)
      double precision depth(Ntot),frontline(LMAX),frontcelltype(LMAX)
      integer N,i,celltype(Ntot),bin,maxbins,count(LMAX)

      maxbins=nint(Lx/frontwidth)

      do bin=1,maxbins
         frontline(bin)=0d0         
         frontcelltype(bin)=0d0
         count(bin)=0
      enddo

      do i=1,N
         if(depth(i).lt.D(i)) then
            xbox=x(i)-(dnint(x(i)/Lx)-0.5d0)*Lx           
            bin=ceiling(xbox/frontwidth)
            frontline(bin)=frontline(bin)+y(i)     
            if(celltype(i).gt.0) then
               frontcelltype(bin)=frontcelltype(bin)+1d0
            endif
            count(bin)=count(bin)+1
         endif
      enddo

      do bin=1,maxbins
         frontline(bin)=frontline(bin)/dble(count(bin))
         frontcelltype(bin)=frontcelltype(bin)/dble(count(bin))
      enddo

      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!    random number generator    !!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END ! end ran2
