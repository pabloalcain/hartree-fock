      Program nrhf
************************************************************************
*  Program to calculate nonrelativistic hartree-fock wave functions
*  for closed-shell atoms
*
*  prepared for Physics 333
*  version 1:    March 1990
*  
*            W. R. Johnson
*            Department of Physics, Notre Dame University
*            Notre Dame, IN 46556 
*
************************************************************************
      implicit real*8(a-h,o-z)
        call datain(h1,gi,*901)
        call hartnr(h1,*901)
        call pickup(*901)
        call output(gi,*901)
        stop
 901    write(6,*) ' Error exit in HF Program'
        stop
      end

      subroutine datain(h1,gi,*) 
**********************************************************************
*
*    Read in atomic configuration data 
*    and set up the radial grid
*
*  Input data
*  Card 1:  idn,jx,jz,jc  (a4,2i4)
*                idn  : atomic symbol
*                jx   : number of shell data cards to follow
*                jz   : nuclear charge
*                jc   : index of last core shell
*                       jc = 0 ==> closed core only
*                       jc = k ==> k+1, ..., jx are valence orbitals
*                     
*
* Cards 2 .. jx + 1:   ns(j),ls(j),js(j) (3i4)
*                ns(j) : principal quantum number of shell
*                ls(j) : angular momentum quantum number of shell
*                js(j) : occupation number of shell
*                      : if js(j) = 0 ==> shell filled
*                                     ==> js(j) = 4*ls(j)+2
* Card jx+2 : r0, h  (*)
*
*                 r(i) = r0 [exp((i-1)*h) -1] , i = 1, ..., NGP
*                  if r0 = 0.0 ==>  r0 = 0.0005
*                  if h  = 0.0 ==>  h = 0.03
* 
* Card jx+3 : h1  (*)
*         
*                 V(r) = h1*Vold(r) + (1.0-h1)*Vnew(r)
*                        h1 = 0.0 ==> no damping
*                        h1 = 1.0 ==> potential unchanged
*
* Card jx+4 : gI  (*)
*                 gI = Nuclear gyromagnetic ratio
*                 gI = 0 ==> gI = 1
*
*********************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500,NOR=30)
      character*4 idn,lab,label
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/charge/z(NGP)
      common/orblab/label
      data r0def /0.0005d0/, hdef/0.03d0/

***  read first data card

      read(5,1010) idn,jx,jz,jc
      write(6,1000) idn
 1000 format('   Input data for ',a/)
      write(6,1010) idn,jx,jz,jc
 1010 format(a,3i4)
      if(jx.gt.NOR) then
        write(6,1011) ' Too many input orbitals !!  shell dim = ',NOR
 1011   format(a,i4)
        go to 901
      endif
      if(jc.eq.0) jc = jx

***  read orbital data cards and set up initial energies

      ncore = 0
      do 100 ja = 1,jx
         read(5,1020) ns(ja),ls(ja),js(ja)
 1020    format(3i4)
         call enclab(ns(ja),ls(ja),*901)
         lab(ja)=label
         if(ja.le.jc) then
            js(ja) = 4*ls(ja) + 2
            ncore = ncore + js(ja)
            wi(ja)=-jz**2/(2.0d0*ns(ja)**2)
         else
            js(ja) = 1
            wi(ja) = -(jz-ncore)**2/(2.0d0*ns(ja)**2)
         endif
         write(6,1030) lab(ja),js(ja),wi(ja)
 1030    format(4x,a,4x,i4,4x,f16.5)
 100  continue

***  read grid parameters -- set up radial grid and initial z(i)

      read(5,*) r0,h
      if(r0.eq.0.0)  r0 = r0def
      if(h.eq.0.0)   h = hdef
      write(6,1040) r0,h
 1040 format(/' r0 =',f10.5,4x,' h=',f10.5/)
      r(1)=0d0
      rp(1)=r0
      rpor(1)=0d0
      z(1) = jz
      do 200 i=2,NGP
         rp(i) = r0*dexp((i-1)*h)
         r(i) = rp(i) - r0
         rpor(i) = rp(i)/r(i)
         z(i) = jz
 200  continue

*** read hartree damping factor
      read(5,*) h1
      write(6,1050) h1
 1050 format(' damping factor =',f8.3/)

*** read gyromagnetic ratio
      read(5,*) gi
      write(6,1060) gi
 1060 format(' gI = ',f10.7)
      return
 901  return 1
      end

      subroutine hartnr(h1,*)
**********************************************************************
*
*   Routine to obtain nonrelativistic Hartree wavefunctions
*
*   V(r) = U(r) - Z/r
*   U(r) = Udir(r) - Uexc(r)
*   Udir(r) = Sum[ js(j)* y0(j,r) ] 
*   Uexc(r) = Udir(r)/N , 
*   N = # of occupied orbitals = Sum [js(j)]
*   
*   Iteration damping factor h1 : 
*                V(r) = h1*Vold(r) + (1.0-h1)*Vnew(r)
*
**********************************************************************
      implicit real*8(a-h,o-z)
      parameter(NOR=30,NGP=500)
      character*4 idn,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/charge/z(NGP)
      common/shells/p(NGP,NOR),q(NGP,NOR)
      dimension rho(NGP),y(NGP),u(NGP),zold(NGP)
      data eps /1.0e-3/, itmax / 50/

      occ = 0.0
      do 100 ja = 1,jc
        occ = occ + js(ja)
 100  continue      
      exc = 1.0d0/occ

      write(6,1000) idn
 1000 format(/' Start Hartree Iteration for ',a/)

      do 300 it = 1,itmax

***  rho(i) = charge density

         do 120 i = 1,NGP
            rho(i) = 0.0
 120     continue

***  mx is max of ms(ja)
         mx = 0
         do 200  ja = 1,jc
            wf(ja) = wi(ja)
            call master(p(1,ja),q(1,ja),z,wf(ja),ns(ja),ls(ja),
     &                    ms(ja),ks(ja),*901)
            ma = ms(ja)
            if(ma.gt.mx) mx = ma
            do 140 i = 1,ma
               rho(i) = rho(i) + js(ja)*p(i,ja)*p(i,ja)
 140        continue
 200     continue

***  y(r) = Udir(r) = Sum[ js(j) y0(j,r) ]
         call yfun(rho,y,0,mx,*901)
      
***  mixture of old and new potentials to damp oscillations
*                V(r) = h1*Vold(r) + (1.0-h1)*Vnew(r)

         do 210 i = 1,NGP
            zold(i) = z(i)
            zn = jz - (1.0d0-exc)*y(i)*r(i)
            z(i) = h1*zold(i) + (1.0d0 - h1)*zn
 210     continue

***  use first order perturbation theory to guess new energies
*    dw = 1st order change in E(n,l) due to change of potential
*    del = max of relative change in energy
*
         del = 0.0d0
         do 250 ja = 1,jc
            ma = ms(ja)
            do 220 i = 1,ma
                u(i) = p(i,ja)**2*(z(i)-zold(i))*rpor(i)
 220        continue
            dw = -rint(u,1,ma,7,h)
            del = max(del,abs(dw/wf(ja)))
            write(6,1040) lab(ja),wi(ja),wf(ja),ks(ja)
 1040       format('  E(',a,') =',2f16.5,6x,'k =',i3)
            wi(ja) = wf(ja) + dw
 250     continue

         write(6,1045) it,del
 1045    format('  loop ',i3,'   rel-err =',1pe9.1/)
***  skip out of loop if del < eps  (1.e-9)
         if(del.lt.eps) go to 350

 300  continue
      write(6,1050)
 1050 format(' Warning hartree: Iteration has not converged' )

 350  continue

***  iteration ended, replace z(i) by value consistent with p(i),q(i)

      do 360 i = 1,NGP
         z(i) = zold(i)
 360  continue

***  replace wi(j) by wf(j)

      do 380 ja = 1,jc
         wi(ja) = wf(ja)
 380  continue

      return
 901  return 1
      end

      subroutine pickup(*)
************************************************************************
*
*   Control the flow of the h-f program by selecting
*   orbital with maximum error and solving the Schrodinger
*   equation for the corresponding orbital
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NOR=30,NFC=100)
      character*4 idn,lab
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      dimension dw(NOR),ddw(NOR),erorb(NOR)
      data emax0/1d30/, del0/1d-3/, eps/1.0d-8/
******   initialization   ***

      itmax = NFC*jc
      k = 0
      ifn = 0
      emax = emax0
      del = del0
      write(6,1000) idn
 1000 format(/'    Hartree-Fock Iteration for ',a/)
      do 100 ja = 1,jc
         k = k + 1
         temp = wf(ja)
         call exchan(ja,del,ifn,*901)
         ddw(ja) = wf(ja) - temp
         write(6,1010) lab(ja),wf(ja),k
 1010    format(2x,a,f16.8,4x,'**********',4x,'initial ',i4)
         k = k + 1
         temp = wf(ja)
         call exchan(ja,del,ifn,*901)
         dw(ja) = wf(ja) - temp
         write(6,1010) lab(ja),wf(ja),k
 100  continue

 200  continue
      ifn = 0
      do 500 kk = k+1,itmax

***   choose orbital with largest error   ***
         do 300 ja = 1,jc
            erorb(ja) = abs((ddw(ja)-dw(ja))/wf(ja))
 300     continue
         ja = 1
         emax = erorb(1)
         do 400 jb = 2,jc
            if(erorb(jb).gt.emax) then
               ja = jb
               emax = erorb(jb)
            endif
 400     continue
         if(emax.le.eps) then
            del = 1d-12
            k = kk
            go to 550
         elseif(emax.lt.1d-7) then
            del = 1d-9
         elseif(emax.lt.1d-5) then
            del = 1d-7
         elseif(emax.lt.1d-3) then
            del = 1d-5
         endif
         ddw(ja)=dw(ja)
         temp = wf(ja)
         call exchan(ja,del,ifn,*901)
         dw(ja) = wf(ja) - temp
         write(6,1020) lab(ja),wf(ja),emax,kk
 1020    format(2x,a,f16.8,4x,1pe10.1,4x,'iterate ',i4)
 500  continue 

***  too many iterations

      write(6,*) ' Error in pickup: too many iterations'
      return 1

***   final convergence test   ***

 550  k = k - 1
      elim = emax
      ifn = 1
      do 600 ja = 1,jc
         k = k + 1
         ddw(ja) = dw(ja)
         temp = wf(ja) 
         call exchan(ja,del,ifn,*901)
         dw(ja) = wf(ja) - temp
         emax = abs((ddw(ja)-dw(ja))/wf(ja))
         if(emax.gt.elim) elim = emax
         write(6,1030) lab(ja),wf(ja),emax,k
 1030    format(2x,a,f16.8,4x,1pe10.1,4x,'final   ',i4)
 600  continue
      if(elim.gt.eps) go to 200
      return
 901  return 1
      end

      subroutine exchan(ja,del,ifn,*)
      implicit real*8(a-h,o-z)
      parameter(NOR=30,NGP=500)
      character*4 idn,lab,lba,lbb
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/charge/z(NGP)
      common/shells/ph(NGP,NOR),qh(NGP,NOR)
      common/wavefn/pf(NGP,NOR),qf(NGP,NOR)
      dimension ind(NOR)
      dimension vdir(NGP),p(NGP),q(NGP)
      dimension u(NGP),v(NGP),w(NGP)
      dimension rv(NGP,NOR)
      data itmax /30/, ind/NOR*0/, init /0/, eps/1d-12/

***   set up direct potential first time thru
*     and set pf(i,ja) = ph(i,ja)
*             qf(i,ja) = qh(i,ja)
*
      save init, vdir
      if(init.eq.0) then
         init = 1
         mx = 0
         do 100 i = 1,NGP
            u(i) = 0.0d0
 100     continue
         do 150 jb = 1,jc
            mb = ms(jb)
            if(mb.gt.mx) mx = mb
            do 120 i = 1,mb
               u(i) = u(i) + js(jb)*ph(i,jb)*ph(i,jb)
 120        continue
            do 130 i = 1,NGP
               pf(i,jb) = ph(i,jb)
               qf(i,jb) = qh(i,jb)
 130        continue
 150     continue
         call yfun(u,vdir,0,mx,*901)
      endif

***   solve for orbital ja
      la = ls(ja)
      ma = ms(ja)
      lba = lab(ja)

***   v(r) = 2 [ Vdir(r)-U(r) ]*pa(r)

      do 200 i = 1,NGP
            v(i) = 0.0d0
 200  continue

      do 250 i = 2,ma
         v(i) = 2.0d0*(vdir(i) + (z(i) - jz)/r(i) )*pf(i,ja)
 250  continue

***    Exchange contributions from subshell jb

      do 500 jb = 1,jc
         lb = ls(jb)
         mb = ms(jb) 
         lbb = lab(jb)
         lmin = iabs(la-lb)
         lmax = la + lb
         mm = min0(ma,mb)
         do 300 i = 1,mm
            u(i) = pf(i,jb)*pf(i,ja)
 300     continue

***    Check scalar product
*
*        if((lb.eq.la).and.(jb.ne.ja)) then
*           do 320 i=1,mm
*              w(i) = u(i)*rp(i)
*320        continue
*           scal = rint(w,1,mm,7,h)
*           write(6,1000) lba,lbb,scal
*1000       format(' <',a,'|',a,'> =  ',1pd8.1)
*        endif

***   Loop thru possible l values

         do 400 l = lmin,lmax
            if(mod(la+lb+l,2).ne.0) go to 400
            call yfun(u,w,l,mm,*901)
            cof = js(jb)*clam(la,lb,l)**2
            do 350 i = 1,mb
               v(i) = v(i) - cof*w(i)*pf(i,jb)
 350        continue
 400     continue
 500  continue

***    Set up damping of inhomogeneous term

      if(ind(ja).eq.0) then
         ind(ja) = 1
         do 510 i=1,ma 
            rv(i,ja) = v(i)
 510     continue
      elseif(ifn.eq.0) then
         do 520 i=1,ma
            rv(i,ja) = 0.5d0*(v(i)+rv(i,ja))
            v(i) = rv(i,ja)
 520     continue
      endif 

***   Estimate correction to energy

      do 550 i = 1,ma
         u(i) = ph(i,ja)*v(i)*rp(i)
         w(i) = ph(i,ja)*pf(i,ja)*rp(i)
 550  continue
 
***   Correct orbital energy in perturbation theory
      de =  0.5d0*rint(u,1,ma,7,h)/rint(w,1,ma,7,h)
      wa = wi(ja) + de 
         
***   Wave function for corrected orbital  

      do 800 iter = 1,itmax
         temp = wa
         call solve(p,q,z,v,wa,ja,del,*901)
         da = abs((wa-temp)/wa)
         if(da.lt.eps) go to 820
 800  continue
      write(6,*) ' Error in fock:  norm loop fails to converge'
      return 1

 820  continue

***  modify direct potential except on final loop

      if(ifn.eq.0) then
         do 850 i = 1,ma
            u(i) = p(i)*p(i) - pf(i,ja)*pf(i,ja)
 850     continue
         call yfun(u,w,0,ma,*901)
         do 860 i=1,NGP
            vdir(i) = vdir(i) + js(ja)*w(i)
 860     continue
      endif

***  store energy and wave function

      wf(ja) = wa
      do 870 i = 1,ma
         pf(i,ja) = p(i)
         qf(i,ja) = q(i)
 870  continue              

      return
 901  return 1
      end

      subroutine output(gi,*)
************************************************************************
*
*   routine to print out information about atom in Hartree approximation
*   (1)  orbital energies
*   (2)  scalar products
*   (3)  averages of r**k  and r**(k+1)
*   (4)  hyperfine constants of valence electrons
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NOR=30,NGP=500,NSCAL=10)
      character*4 idn,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/wavefn/p(NGP,NOR),q(NGP,NOR)
      dimension u(NGP),v(NGP),w(NGP),rhs(NGP)
      dimension over(NSCAL),icnt(NSCAL),jov(NSCAL,NSCAL)

      alpha = 137.0359895d0
      ryd = 219474.624d0/2.0d0
      cl = 29979245800.0d0
      ryc = ryd*cl/1.0d6
      ratmas = 1.0d0/1836.152701d0
      hfc = ratmas*ryc/alpha**2

***  calculate total energy of core

      do 10 i = 1,NGP
         u(i) = 0.0
 10   continue
      mx = 0
      do 30 ja = 1,jc
         ma = ms(ja)
         if(ma.gt.mx) mx = ma
         do 20 i = 1,ma
            u(i) = u(i) + js(ja)*p(i,ja)**2
 20      continue
 30   continue
      call yfun(u,v,0,mx,*901)
      etot = 0.0
      do 90 ja = 1,jc
          mk = ms(ja)
          la = ls(ja)
          do 40 i = 1,NGP
              rhs(i) = 0.0
 40       continue
          do 80 jb = 1,jc
             mb = ms(jb)
             lb = ls(jb)
             mm = min0(ma,mb)
             do 50 i = 1,mm
                 u(i) = p(i,jb)*p(i,ja)
 50          continue
             low = iabs(la - lb)
             lhi = la + lb
             do 70 ll=low,lhi
                lt = la + lb + ll
                if(mod(lt,2).ne.0) go to 70
                call yfun(u,w,ll,mm,*901)
                adg = 0.5d0*js(jb)*clam(la,lb,ll)**2
                do 60 i = 1,mb
                   rhs(i) = rhs(i) + adg*w(i)*p(i,jb)
 60              continue
 70           continue
 80       continue
          do 85 i = 1,ma
             u(i) = ( p(i,ja)**2*v(i) - p(i,ja)*rhs(i) )*rp(i)
 85       continue
          term = rint(u,1,ma,7,h)
          etot = etot + js(ja)*( wf(ja) - 0.5d0*term )
 90   continue

***   calculate valence wave functions and energies

      if(jx.gt.jc) then
         call valenc(*901)
      endif

***   print out orbital energies

      write(6,1000) idn,jz
 1000 format(/'  Hartree-Fock Energy Levels for ',a,'  Z =',i3/)
      write(6,1010)
 1010 format('  Shell','  #el',6x,'  Energy')

      do 150 ja = 1,jx
         write(6,1020) lab(ja),js(ja),wf(ja)
 1020    format(3x,a,i4,f16.6)
 150  continue

***  print out total energy

      write(6,1025) etot
 1025 format(/3x,'Ehf core',f18.8)

***   calculate and print scalar products

      write(6,1100)
 1100 format(/'  *****  scalar products   *****')

      ic = 0
      do 220 lin = 0,4
              ic = ic + 1
              ins = 0
              do 200 j = 1,jx
                  if(ls(j).ne.lin) go to 200
                  ins = ins + 1
                  jov(ins,ic) = j
 200          continue
              icnt(ic) = ins
 220  continue

      ic = 0
      do 260 lin = 0,4
         ic = ic + 1
         jns = 0
         do 240 j = 1,jx
             if(ls(j).ne.lin) go to 240
             jns = jns + 1
             ins = 0
             do 230 jl = 1,j
                if(ls(jl).ne.lin) go to 230
                ins = ins + 1
                mm=min0(ms(j),ms(jl))
                do 225 i=1,mm
                   u(i)=p(i,j)*p(i,jl)*rp(i)
 225            continue
                over(ins) = rint(u,1,mm,7,h)
                if(jl.eq.j) over(ins) = over(ins) - 1.0
                if(jns.eq.1) then
                   write(6,1200) (lab(jov(ii,ic)),ii = 1,icnt(ic))
 1200              format(/11x,8('    ',a,'> '))
                endif
 230        continue
            write(6,1210) lab(j),(over(ii),ii=1,ins)
 1210       format('  <',a,'|   ',1p,8e10.1)
 240     continue
 260  continue

  500 continue
      write(6,1500)
 1500 format(/)
      open(unit=1,file='fort.1',form='unformatted',status='unknown')
      rewind 1
      write(1) idn,jz,jx,jc
      write(1) r,rp,rpor,h
      do 600 ia = 1,jx
          write(1) ns(ia),ls(ia),ms(ia),lab(ia),wf(ia)
          write(1) (p(i,ia),i=1,NGP)
 600  continue
      rewind 1
      return
 901  return 1
      end

      subroutine valenc(*)
****************************************************************
*
*     determine valence orbitals
*
****************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500,NOR=30)
      character*4 idn,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/wavefn/p(NGP,NOR),q(NGP,NOR)
      dimension u(NGP),v(NGP),y(NGP),rhs(NGP)
      dimension z(NGP),zex(NGP),phs(NGP),qhs(NGP)
      data itmax /40/, eps /1.0d-9/, enorm/1.0d-12/, ept /1.0e30/

***   set up direct potential and starting df functions
      
      pi = dacos(-1.0d0)
      th = 1.0d0/3.0d0
      cex = (3.0d0/(32.0d0*pi**2))**th
      do 100  i = 1,NGP
         u(i) = 0.0
         zex(i) = 0.0
 100  continue
      mx=0
      do 120  ja = 1,jc
         ma = ms(ja)
         if(ma.gt.mx) mx = ma
         do 110  i = 1,ma
            u(i) = u(i) + js(ja)*p(i,ja)**2
 110     continue
 120  continue
      do 130 i = 2,mx
         zex(i) = cex*(r(i)*u(i))**th
 130  continue
      zex(1) = 0.0d0

      call yfun(u,y,0,mx,*901)

      do 140 i = 1,NGP
         z(i) = jz - y(i)*r(i) + zex(i)
 140  continue

***  loop over all valence states

      do 500 jv = jc+1,jx
         wf(jv) = wi(jv)
         call master(phs,qhs,z,wf(jv),ns(jv),ls(jv),
     &                    ms(jv),ks(jv),*901)
         write(6,1000) lab(jv),lab(jv),wf(jv),ept,ks(jv)
 1000    format(/5x,'Valence Iteration Loop for ',a,' shell'//
     1          2x,a,f16.8,4x,f10.1,4x,'initial ',i4)
*        write(6,1000) wi(jv),wf(jv),ks(jv)
*1000    format(/5x,'Valence electron iteration'/
*    1      5x,'wi=',f14.7,4x,'wf=',f14.7,4x,'# iter=',i3)

         wi(jv) = wf(jv)
         mv = ms(jv)
         lv = ls(jv)
         do 150 i = 1,mv
            p(i,jv) = phs(i)
            q(i,jv) = qhs(i)
 150     continue
         if(mv.lt.NGP) then
            mv1 = mv + 1
            do 155 i = mv1,NGP
               p(i,jv) = 0.0
               q(i,jv) = 0.0
 155        continue
         endif

***  exchange iteration loop
         err = 1.0
         del  = 1.0d-3
         do 400 it = 1,itmax

***   set accuracy for norm calculation in routine solve
            if(err.lt.1.0d-8)  then
               del = 1.0d-11
            elseif(err.lt.1.0d-6) then
               del = 1.0d-9
            elseif(err.lt.1.0d-4) then
               del = 1.0d-7
            endif

***   calculate inner shell overlap integrals
*
*           do 170 ja = 1,jc
*              if(ls(ja).ne.lv) go to 170
*              mm = min0(ms(ja),mv)
*              do 160 i = 1,mm
*                 u(i) = p(i,ja)*p(i,jv)*rp(i)
*160           continue
*              over = rint(u,1,mm,7,h)
*              write(6,1010) lab(ja),lab(jv),over
*1010          format('     <',a,'|',a,'> =', 1pe10.2)
*170        continue
*
***   calculate valence norm
*
*           do 180 i = 1,mv
*              u(i) = p(i,jv)**2*rp(i)
*180        continue 
*           over = 1.0d0 - rint(u,1,mv,7,h)
*           write(6,1020) lab(jv),lab(jv),over
*1020       format(' 1.0-<',a,'|',a,'> =', 1pe10.2/)
*
***   setup rhs

            v(1) = 0.d0
            do 200 i = 2,mv
               v(i) = 2.0d0*zex(i)*p(i,jv)/r(i)
 200        continue
            if(mv.lt.NGP) then
               mv1 = mv + 1
               do 205 i = mv1,NGP
                  v(i) = 0.0
 205           continue
            endif
            do 300 jb = 1,jc
               lb = ls(jb)
               mb = ms(jb)
               mm = min0(mb,mv)
               do 240 i = 1,mm
                  u(i) = p(i,jb)*p(i,jv)
 240           continue
               low = iabs(lb-lv)
               lhi = lb + lv
               do 260 ll = low,lhi
                  lt = lb + lv + ll
                  if(mod(lt,2).ne.0) go to 260
                  call yfun(u,y,ll,mm,*901)
                  cl = js(jb)*clam(lb,lv,ll)**2
                  do 250 i = 1,mb
                     v(i) = v(i) - cl*y(i)*p(i,jb)
 250              continue
 260           continue
 300        continue

***   Apply 50%-50% damping

            if(it.eq.1) then
               do 320 i = 1,mv
                  rhs(i) = v(i)
 320           continue
            else
               do 330 i = 1,mv
                  rhs(i) = 0.5*( rhs(i) + v(i) )
                  v(i) = rhs(i)
 330           continue
            endif

***   Estimate correction to energy

            do 350 i = 1,mv
               u(i) = phs(i)*v(i)*rp(i)
               y(i) = phs(i)*p(i,jv)*rp(i)
 350        continue

***   Correct orbital energy in perturbation theory

            de = 0.5d0*rint(u,1,mv,7,h)/rint(y,1,mv,7,h)
            wv = wi(jv) + de
            temp = wf(jv)
***   Determine corrected wave function for orbital

***   Loop for proper norm
            do 380 inor = 1,itmax
               tempn = wv
               call solve(p(1,jv),q(1,jv),z,v,wv,jv,del,*901)
               wf(jv) = wv
               ernorm = abs(wf(jv)-tempn)
               if(ernorm.lt.enorm) go to 390
 380        continue
            write(6,*) ' Error in valence: too many norm iterations'
            return 1
 390        continue

            err = abs(( wf(jv) - temp )/wf(jv))

            write(6,1030) lab(jv),wf(jv),err,it
 1030       format(2x,a,f16.8,4x,1pe10.1,4x,'iterate ',i4)

*           write(6,1030) it, wf(jv), err
*1030       format(4x,'iteration',i3,4x,' wf(au) =',1pe14.6,4x,
*    &                'err =',e10.2)

            if(err.le.eps) go to 450
 400     continue
         write(6,*) ' Error in valence: too many exchange iterations'
         return 1
 450     continue

 500  continue
      return
 901  return 1
      end

      subroutine solve(p,q,z,rv,wa,ja,del,*)
************************************************************************
*
*  Routine to solve the inhomogeneous Schrodinger equation 
*          using Green's functions
*
*  p'' + 2( wa - z(r)/r - la(la+1)/2r**2)p = rv(r)
*  q = p'
*
*  p1'' + 2( wa - z(r)/r - la(la+1)/2r**2)p1 = 0 : reg at origin
*  q1 = p1'
*
*  p2'' + 2( wa - z(r)/r - la(la+1)/2r**2)p2 = 0 : reg at infinite
*  q2 = p2'
*  
*  Solution: p(r) = [   p1(r) Int{r,inf} p2(r)*rv(r)
*                     + p2(r) Int{0,r} p1(r)*rv(r)     ]/wron
*
*            q(r) = [   q1(r) Int{r,inf} p2(r)*rv(r)
*                     + q2(r) Int{0,r} p1(r)*rv(r)     ]/wron
*
*            wron = [ q2(r)*p1(r) - p2(r)*q1(r) ]  = const
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500,NOR=30)
      character*4 idn,lab
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/wavefn/pf(NGP,NOR),qf(NGP,NOR)
      dimension p(NGP),q(NGP),z(NGP),rv(NGP)
      dimension p1(NGP),q1(NGP),p2(NGP),q2(NGP)
      dimension u(NGP),w(NGP),x(NGP),y(NGP)
     
      la = ls(ja)
      ma = ms(ja)

***  stage 1 : solve equation with input energy
      call outsch (p1,q1,z,wa,la,ma)
      call insch(p2,q2,z,wa,la,ma,2,*901)
      call wron(wr,p1,q1,p2,q2,ma)
      do 100 i = 1,ma
         u(i) = p1(i)*rv(i)*rp(i)
         w(i) = p2(i)*rv(i)*rp(i)
 100  continue
      call yint(u,w,x,y,ma,h)
      do 200 i = 1,ma
         p(i) = (p1(i)*y(i)+p2(i)*x(i))/wr
         q(i) = (q1(i)*y(i)+q2(i)*x(i))/wr 
         u(i) = p(i)*p(i)*rp(i)
 200  continue

***  stage 2 : use perturbation theory to modify the energy 
*              to normalize wave function

      anorm = rint(u,1,ma,7,h)
*     write(6,1000) lab(ja),lab(ja),anorm,wa
*1000 format(' <',a,'|',a,'> = ',1pd18.10,'   w =',1pd18.10)
      if(abs(anorm - 1.0d0).lt.del) return
      do 250 i = 1,ma
         u(i) = p1(i)*pf(i,ja)*rp(i)
         w(i) = p2(i)*pf(i,ja)*rp(i)
 250  continue
      call yint(u,w,x,y,ma,h)
      do 300 i = 1,ma
         u(i) = (p1(i)*y(i)+p2(i)*x(i))/wr
         x(i) = p(i)*u(i)*rp(i) 
 300  continue
      de = 0.25d0*(anorm - 1.0d0)/rint(x,1,ma,7,h)
      wa = wa + de
      return
 901  return 1
      end

       subroutine wron(wr,p1,q1,p2,q2,m)
************************************************************************
*
*  Calculates and monitors the wronskian 'wr' of two independent 
*  solutions (p1,q1) and (p2,q2) to the radial Schrodinger equation
*          wr = q2*p1 - p2*q1
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500)
      dimension p1(NGP),q1(NGP),p2(NGP),q2(NGP)
      data eps/1.0d-4/,nn/40/

      ma = m/2
      wr = q2(ma)*p1(ma) - p2(ma)*q1(ma)
      ma = nn
      mb = m - nn
      mc = nn

      do 100  i = ma,mb,mc
         w = q2(i)*p1(i) - p2(i)*q1(i)
         dw = dabs((w-wr)/wr)
         if(dw.gt.eps) go to 200
 100  continue
      return
 200  write(6,1000)
 1000 format(/' Constancy of the wronskian :')
      do 300 i = ma,mb,mc
         w = (q2(i)*p1(i) - p2(i)*q1(i))/wr
         write(6,1200) i,w
 1200    format(i6,1pd20.7)
 300  continue 
      return
      end

      function clam(l1,l2,l)
************************************************************************
*
*    clam(l1,l2,l) = threej(l1, l2, l, 0, 0, 0)
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NF=100)
      dimension fact(0:NF)
      data init /0/
      save init, fact

*** initialize factorial array if neccessary
*   fact(i) = (i)!
*   fact(0) = 1
*
      if(init.eq.0) then
         init = 1
         fact(0) = 1.0d0
         do 100 i = 1,NF
           fact(i) = i*fact(i-1)
 100     continue
      endif
***  test triangle inequality
      if(l.gt.(l1+l2).or.l.lt.iabs(l1-l2)) then
         clam = 0.0d0
         return
      endif
***  test parity selection rules
      lg = l1+l2+l
      if(mod(lg,2).eq.0) then
*  even parity 
         jg = lg/2
         x = sqrt(fact(lg-2*l1)*fact(lg-2*l2)*fact(lg-2*l)/fact(lg+1))
         y = fact(jg)/(fact(jg-l1)*fact(jg-l2)*fact(jg-l))
         clam = (-1)**jg*x*y
         return
      else
* odd parity
         clam = 0.0d0
         return
      endif
      end

      subroutine yfun(x,y,l,m,*)
************************************************************************
*
*  this program calculates the Hartree's y-functions
*
*  x  : input function
*  y  : output hartree's y-function : y(l,r)/r =
*            r**(-l-1) * Int {0,r} [ r**l x(r) dr ] 
*               + r**l * Int{r,inf}[ r**(-l-1) x(r) ] 
*  l  : order of the y function (must be .ge. 0)
*  m  : number of tabulation points for the input function x
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      dimension x(NGP),y(NGP)
      dimension u(NGP),v(NGP),w(NGP),s(NGP)
      data ja/6/
*
*  note that ja is determined by the order of the integration method
*  used in subroutine yint
*
      v(1)=0.0
      w(1)=0.0
      if(l.lt.0) then
         write(6,*) ' Error in yint: l < 0 '
         return 1
      elseif(l.eq.0) then 
         do 100  i= 2,m
            v(i) = x(i)*rp(i)
            w(i) = x(i)*rpor(i)
 100     continue
         call yint(v,w,y,u,m,h)
         ym = y(m)
         y(1) = u(1)
         do 120 i = 2,m
            y(i) = y(i)/r(i) + u(i)
 120     continue 
         m1 = m + 1
         if(m1.le.NGP) then
           do 140 i = m1,NGP
              y(i) = ym/r(i)
 140       continue
         endif
         return
      else
         if(l.gt.ja) then
            write(6,1000) l
 1000       format('  Warning in yint:  the value of l =',
     &      i4,' is too large')
         endif
         s(1)=0.0
         do 160 i = 2,NGP
            s(i)=r(i)**l
 160     continue
         do 180 i = 2,m
            v(i) = x(i)*rp(i)*s(i)
            w(i) = x(i)*rpor(i)/s(i)
 180     continue
         call yint(v,w,y,u,m,h)
         ym = y(m)
         do 200 i = 2,m
            y(i) = y(i)/(r(i)*s(i)) + u(i)*s(i)
 200     continue
         m1 = m + 1
         if(m1.le.NGP) then
            do 210  i = m1,NGP
               y(i) = ym/(r(i)*s(i))
 210        continue
         endif
      endif
      return
      end

      subroutine yint (v,w,y,z,m,h)
************************************************************************
*  this program calculates the indefinite integrals y and z using the
*  lagrange integration formula
*
*  y(r) = integral of v from 0 to r
*  z(r) = integral of w from r to infinity
*  m is the maximum tabulation point of v and w (virtual infinity)
*  h is the step size of the radial grid
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500)
      dimension v(NGP),w(NGP),y(NGP),z(NGP)
      dimension dy(NGP),dz(NGP)
*
*     lagrange 6 point integration formula
*     ************************************
*      parameter(NO=6,NP=NO/2)
*      dimension aa(NP,NO),a(NO,NP),b(NP)
*      data  da/1440./
*      data  aa/  475.,   -27.,    11.,
*     2          1427.,   637.,   -93.,
*     3          -798.,  1022.,   802.,
*     4           482.,  -258.,   802.,
*     5          -173.,    77.,   -93.,
*     6            27.,   -11.,    11./
*      data   b/  802.,   -93.,    11./
************************************************************************
*     lagrange 8 point integration formula
*     ************************************
*      parameter(NO=8,NP=NO/2)
*      dimension aa(NP,NO),a(NO,NP),b(NP)
*      data  da/120960./
*      data  aa/  36799.,  -1375.,   351.,  -191.,
*     2          139849.,  47799., -4183.,  1879.,
*     3         -121797., 101349., 57627., -9531.,
*     4          123133., -44797., 81693., 68323.,
*     5          -88547.,  26883.,-20227., 68323.,
*     6           41499., -11547.,  7227., -9531.,
*     7          -11351.,   2999., -1719.,  1879.,
*     8            1375.,   -351.,   191.,  -191./   
*      data   b/  68323.,  -9531.,  1879.,  -191./
************************************************************************
*     lagrange 10 point integration formula
*     ************************************
      parameter(NO=10,NP=NO/2)
      dimension aa(NP,NO),a(NO,NP),b(NP)
      data da/ 7257600./
      data aa/ 2082753.,  -57281.,   10625.,   -3969.,    2497., 
     2         9449717., 2655563., -163531.,   50315.,  -28939.,
     3       -11271304., 6872072., 3133688., -342136.,  162680.,
     4        16002320.,-4397584., 5597072., 3609968., -641776.,
     5       -17283646., 3973310.,-2166334., 4763582., 4134338.,
     6        13510082.,-2848834., 1295810.,-1166146., 4134338., 
     7        -7394032., 1481072., -617584.,  462320., -641776., 
     8         2687864., -520312.,  206072., -141304.,  162680.,
     9         -583435.,  110219.,  -42187.,   27467.,  -28939.,
     a           57281.,  -10625.,    3969.,   -2497.,    2497./
      data  b/ 4134338., -641776.,  162680.,  -28939.,    2497./
************************************************************************
      data h0/0./
*
*  note that a different even order method can be used by replacing the
*  dimension and data statements in above block 
************************************************************************
      save h0, a, b, hd 
      if(h.eq.h0) go to 180
      hd = h/da
      do 150 i = 1,NP
         do 100 j = 1,NO
            a(j,i) = aa(i,j)*hd
 100     continue
         b(i) = b(i)*hd   
 150  continue
      h0 = h
 180  continue
      y(1) = 0.
      z(m) = 0.
      do 200 i = 2,NP
         k = m - i + 1
         y(i) = y(i-1)
         z(k) = z(k+1)
         ii = i - 1
         do 190 j = 1,NO
            y(i) = y(i) + a(j,ii)*v(j)
            z(k) = z(k) + a(j,ii)*w(m-j+1)
 190     continue
 200  continue
      im = NP + 1
      in = m - NP + 1

       do 240 i = im,in
          k = m - i + 1
          dy(i) = 0.0
          dz(k) = 0.0
          do 230 j = 1,NP
             dy(i) = dy(i) + b(j)*(v(i+j-1) + v(i-j))
             dz(k) = dz(k) + b(j)*(w(k-j+1) + w(k+j))
  230     continue        
  240  continue

      do 250 i = im,in
         k = m - i + 1
         y(i) = y(i-1) + dy(i)
         z(k) = z(k+1) + dz(k)
 250  continue

      in = in + 1
      do 300 i = in,m
         k = m - i + 1
         y(i) = y(i-1)
         z(k) = z(k+1)
         do 280 j = 1,NO
            y(i) = y(i) + a(j,k)*v(m-j+1)
            z(k) = z(k) + a(j,k)*w(j)
 280     continue
 300  continue
      return
      end

      subroutine master(p,q,z,eau,n,l,m,k,*)
************************************************************************
*
*   find radial eigenfunction and eigenvalue for
*   Schroedinger equation
*   p(i) = radial wave function                   (output)
*   q(i) = (dp/dr)(i)                             (output)
*   z(i) = effective charge;  V(r) = -z(r)/r      (input)
*   eau = eigenvalue:  input is a guess           (input) and
*                      output is precise value    (output)
*   n = principal quantum number                  (input)
*   l = angular momentum quantum number           (input)
*   m = practical infinity for orbital            (output)
*   k = number of trys at solution <= ntry        (output)
*   (*)  is an error exit
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      dimension p(NGP),q(NGP),z(NGP)
      dimension u(NGP)
      data del/1.0d-11/,ntry/30/,alr/800.0d0/
      
      eupper=0.0d0
      elower=0.0d0
      more=0
      less=0
***   nint = number of internal nodes = n-l-1
      nint = n - l - 1

****  k = number of trys at solution  <= ntry  (ntry = 30)
      k = 0
 100  k = k + 1
      if(k.gt.ntry) then
         k = k - 1
         return 1
      endif

****  seek practical infinity:  m
* 
*       p(r(m)) = exp(- alam r(m) )
*       alam = sqrt( -2 (eau + z(m)/r(m)) )
*       asymptotic region: alam r(m) > 40  (alr = 40*40/2)
*       ==>   -2 (eau*r(m) + z(m))*r(m) > 1600
*       ==>   (eau*r(m) + z(m))*r(m) + 800 < 0
*
      m = NGP + 1
 110  continue
         m = m - 1
      if( ( ( eau*r(m) + z(m) )*r(m) + alr ).lt.0.0d0 )  go to 110

****  seek classical turning point: mtp
*
*     classically forbidden region: eau + z(mtp)/r(mtp) < 0
*     ==>  eau*r(mtp) + z(mtp) < 0
* 
      mtp = m + 1
 120  continue
          mtp = mtp - 1
      if( ( eau*r(mtp) + z(mtp) ).lt.0.0d0 ) go to 120

**** integrate inward from practical infinity to classical turning point
      call insch(p,q,z,eau,l,m,mtp,*901)

**** save p(mtp) 
      ptp = p(mtp)
      qtp = q(mtp)
**** integrate outward from origin to classical turning point
      call outsch(p,q,z,eau,l,mtp)

**** match p(mct) to make solution continuous

      rat = p(mtp)/ptp
      qtp = rat*qtp
      
      mtp1 = mtp + 1
      do 130 i = mtp1,m
         p(i) = rat*p(i)
         q(i) = rat*q(i)
 130  continue

**** find number of zeros in p(r)
 
      nzero =0

      sp = dsign(1.0d0,p(2))

      do 140 i = 3,m
      if(dsign(1.0d0,p(i)).ne.sp) then
        sp = -sp
        nzero = nzero + 1
      endif 
 140  continue
 
**** compare nzero with nint
*
*    nzero > nint  ==> eau too high, 
*                       decrease and try again
*    nzero < nint  ==> eau too low,
*                       increase and try again
*    nzero = nint  ==> eau in correct correct range,
*                       use perturbation theory
*
      if(nzero.gt.nint) then
         more = more + 1
*** record upper bound on eau
         if((more.eq.1).or.(eau.lt.eupper)) eupper = eau
         if(less.ne.0) then
**** use average of upper and lower bounds
            eau = 0.5d0*(eupper+elower)
            go to 100
          endif
**** otherwise take 10% less than upper bound
          eau = 1.10d0*eupper
          go to 100
      elseif(nzero.lt.nint) then
          less = less + 1
          if((less.eq.1).or.(eau.gt.elower)) elower = eau
          if(more.ne.0) then
**** use average of upper and lower bounds
             eau = 0.5d0*(eupper+elower)
             go to 100
          endif
**** otherwise take 10% more than lower bound
          eau = 0.90d0*elower
          go to 100
      endif

****  calculate normalization integral
      do 150 i = 1,m
         u(i) = p(i)*p(i)*rp(i)
 150  continue
      anorm = rint(u,1,m,7,h)
**** use perturbation theory to calculate energy change
*
*    de = p(mtp) [ q(mpt-) - q(mtp+) ] / [ 2 int( p(r)**2 ) ]
*
      de = p(mtp)*(q(mtp) - qtp)/(2.0d0*anorm)
      eau = eau + de
      if((less.ne.0).and.(eau.lt.elower)) then
         eau = 0.5d0*(eau - de + elower)
         go to 100
      elseif(eau.gt.eupper) then
         eau = 0.5d0*(eau - de + eupper)
         go to 100
      elseif(dabs(de/eau).gt.del) then
         go to 100
      endif
      eau = eau - de

      an = 1.0d0/dsqrt(anorm)

      do 200 i = 1,m
         p(i) = an*p(i)
         q(i) = an*q(i)
 200  continue
      return
 901  return 1
      end

      subroutine adams(p,q,z,eau,l,na,nb)
************************************************************************
*   Subroutine to solve the radial Schroedinger equation using
*   Adams-Moulton NO-point method:  NO = 6,7,8,9
*   used in conjunction with 
*   subroutines outsch and insch
*
*   p(i) = radial wave function  
*   q(i) = (dp/dr)(i)
*   dp/dr = q(r)
*   dq/dr = -2(eau + z(r)/r - l(l+1)/r**2) p(r)
*  if nb > na then
*         p(na-NO)...p(na-1) given on input
*         q(na-NO)...q(na-1) given on input
*         p(na)...p(nb)  found by forward integration
*         q(na)...q(nb)  found by forward integration
*   else
*         p(nb+1)...p(nb+NO) given on input
*         q(nb+1)...q(nb+NO) given on input
*         p(nb)...p(na)  found by backward integration
*         q(nb)...q(na)  found by backward integration
*
*   z(i) effective charge, V(r) = - z(r)/r
*   eau = energy in a.u.
*   l = orbital angular momentum quantum number
*
************************************************************************ 
      implicit real*8(a-h,o-z)
      parameter(NGP=500,NO=8)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      dimension p(NGP),q(NGP),z(NGP)
      dimension dp(NO),dq(NO),aa(NO)
      dimension ia(NO)
******** six-point adams method ***************************************
*     parameter(NO=5)
*      data ia /  27,-173, 482,-798,1427/
*      data id /1440/,iaa/ 475/
******** seven-point adams method **************************************
*     parameter(NO=6)
*      data ia /  -863,  6312,-20211, 37504,-46461, 65112/
*      data id /60480/,  iaa/19087/
******** eight-point adams method **************************************
*     parameter(NO=7)
*      data ia /   1375, -11351,  41499, -88547, 123133,-121797, 139849/
*      data id / 120960/,   iaa/  36799/
********** nine-point adams method *************************************
*     parameter(NO=8)
      data ia/  -33953,  312874,-1291214, 3146338,-5033120, 5595358,
     1        -4604594, 4467094/         
      data id /3628800/,   iaa /1070017/
************************************************************************
      cof = h/dfloat(id)
      ang = 0.5d0*l*(l+1)

*    fill in the preliminary arrays for derivatives

      if(nb.gt.na) then
        inc = 1
        mstep = nb-na+1
      else
        inc = -1
        mstep = na-nb+1
      endif
      i = na-inc*(NO+1)
      do 100 k = 1,NO
         i = i + inc
         dp(k) = inc*rp(i)*q(i)
         dq(k) = -2*inc*( eau*rp(i) +
     &                (z(i)-ang/r(i))*rpor(i) )*p(i)
         aa(k) = cof*ia(k)
 100  continue
      a0 = cof*iaa    
      i = na - inc
      do 500 ii = 1,mstep
         i = i + inc
         dpi = inc*rp(i)
         dqi = -2*inc*( eau*rp(i)+ (z(i)-ang/r(i))*rpor(i) )
         b = a0*dpi
         c = a0*dqi
         det = 1.0d0 - b*c
         sp = p(i-inc)
         sq = q(i-inc)
         do 200 k = 1,NO
            sp = sp + aa(k)*dp(k)
            sq = sq + aa(k)*dq(k)
 200     continue
         p(i) = (sp + b*sq)/det
         q(i) = (c*sp + sq)/det
         do 300 k = 1,NO-1
             dp(k) = dp(k+1)
             dq(k) = dq(k+1)
 300     continue
         dp(NO) = dpi*q(i)
         dq(NO) = dqi*p(i)
 500  continue
      return
      end

      subroutine insch(p,q,z,eau,l,nb,na,*)
************************************************************************
*   routine to integrate the radial Schroedinger equation
*   inward from r(nb) = "practical infinity"  to r(na)
*                  
*   p(i) = radial wave function  
*   q(i) = (dp/dr)(i)
*   dp/dr = q(r)
*   dq/dr = -2(eau + z(r)/r - l(l+1)/r**2) p(r)
*   
*   z(i) = effective charge, v(r) = -z(r)/r
*   eau = energy in a.u.
*   l = orbital angular momentum quantum number
*   na : r(na) = final point for integration
*   nb : r(nb) = initial point = practical infinity 
*   nb > na
*
*   Method:
*
*   a)   Asymptotic exspansion
*        p(r) = r**sig exp(-lam r ) [ a0 + a1/r + a2/r**2 + ...
*
*           p(nb)...p(nb-NO)  found this way
*           q(nb)...q(nb-NO) found this way
*
*    b)  Continue solution using Adams-Moulton interpolation scheme
*
*           p(nb-NO-1)..p(na) found this way
*           q(nb-NO-1)..q(na) found this way
*
************************************************************************
      implicit real*8(a-h,o-z)
      parameter(NGP=500,NO=8,NX=15)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      dimension p(NGP),q(NGP),z(NGP)
      dimension ax(0:NX),bx(0:NX)
      data eps /1.0d-8/, epr/1d-3/
************************************************************************
*   a)  Asymptotic series
*
*       p(r) = r**sig exp(-lam r) [ 1 + a1/r + a2/r**2 + ... ]
*       lam = sqrt(-2 e)
*       sig = z/lam
*       a0 = 1
*       a(k+1) = [l(l+1)-(sigk)(sig-k-1)] a(k) /[2 lam (k+1) ]
*
************************************************************************
      if(eau.ge.0.0d0) then
      write(6,1010) eau
 1010 format('  Error in insch: eau = ',1pd16.6,' > 0' )
      return 1
      endif
      alam = dsqrt(-2.0d0*eau)
      sig = z(nb)/alam
      ang = l*(l+1)
      ax(0) = 1.0d0
      bx(0) = -alam
      do 100 k = 1,NX
         ax(k) = (ang - (sig-k+1)*(sig-k))*ax(k-1)/(2.0d0*k*alam)
         bx(k) = ((sig-k+1)*(sig+k) - ang)*ax(k-1)/(2.0d0*k)
 100  continue 

      do 300 i = nb,nb-NO,-1
         rfac = r(i)**sig*dexp(-alam*r(i))
         ps = ax(0)
         qs = bx(0)
         rk = 1.0d0
         do 200 k = 1,NX
             rk = rk*r(i)
             pt = ax(k)/rk
             qt = bx(k)/rk
             ps = ps + pt
             qs = qs + qt
             xe = dmax1(dabs(pt/ps),dabs(qt/qs))
             if(xe.lt.eps) go to 250
 200     continue            
         if(xe.gt.epr) write(6,1000)
 1000    format(' Warning: asymptotic series not converging ')
 250     continue
         p(i) = rfac*ps
         q(i) = rfac*qs
 300  continue
************************************************************************
*
*  b)   Continue solution inward using Adams method
*       n.b.  be sure NO is consistent in the two routines !
*
************************************************************************
      if((nb-NO-1).ge.na) then
        call adams(p,q,z,eau,l,nb-NO-1,na)
      endif
      return
      end

      subroutine outsch(p,q,z,eau,l,nb)
************************************************************************
*   routine to integrate the radial Schroedinger equation
*   outward from r(1) = 0 to r(nb)
*                  
*   p(i) = radial wave function  
*   q(i) = (dp/dr)(i)
*   dp/dr = q(r)
*   dq/dr = -2(eau + z(r)/r - l(l+1)/r**2) p(r)
*   
*   z(i) = effective charge, v(r) = -z(r)/r
*   eau = energy in a.u.
*   l = orbital angular momentum quantum number
*   nb : r(nb) = final point for integration
*
*   Method:
*
*   a)   Factor r**(l+1) and start integration using NO-order Lagrangian
*        differentiation formulas NO = 6,7,8,9 for LO loops
*
*           p(1)...p(LO*NO+1) found this way
*           q(1)...q(LO*NO+1) found this way
*
*    b)  Continue solution using Adams-Moulton interpolation scheme
*
*           p(LO*NO+2)..p(nb) found this way
*           q(LO*NO+2)..q(nb) found this way
*
************************************************************************

      implicit real*8(a-h,o-z)
      parameter(NGP=500,NO=8,LO=3)
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      dimension p(NGP),q(NGP),z(NGP)
      dimension us(NO),vs(NO),b(NO),c(NO),d(NO),s(NO)
      dimension em(NO,NO),fm(NO,NO)
      dimension ia(NO),ie(NO,NO),ipvt(NO)
      dimension work(NO),det(2)
*
*      the following coefficients are suitable for a sixth-order method *
*      use with six-point adams method, set NO = 5
*
*      data ie  / -65,  -30,   15,  -20,   75,
*     1           120,  -20,  -60,   60, -200,
*     2           -60,   60,   20, -120,  300,
*     3            20,  -15,   30,   65, -300,
*     4            -3,    2,   -3,   12,  137/
*      data ia  / -12,    3,   -2,    3,  -12/
*      data id  /  60/
*
*      the following coefficients are suitable for a seventh-order method *
*      use with seven-point adams method, set NO = 6
*
*      data ie  / -77,  -24,    9,   -8,   15,  -72,
*     1           150,  -35,  -45,   30,  -50,  225,
*     2          -100,   80,    0,  -80,  100, -400,
*     3            50,  -30,   45,   35, -150,  450,
*     4           -15,    8,   -9,   24,   77, -360,
*     5             2,   -1,    1,   -2,   10,  147/
*      data ia  / -10,    2,   -1,    1,   -2,   10/
*      data id  /  60/
*
*     the following coefficients are suitable for a eighth-order method*
*     use with eight-point adams method, set NO = 7
*
*      data ie / -609, -140,   42,  -28,    35,   -84,   490,
*     1          1260, -329, -252,  126,  -140,   315, -1764, 
*     2         -1050,  700, -105, -420,   350,  -700,  3675,
*     3           700, -350,  420,  105,  -700,  1050, -4900,
*     4          -315,  140, -126,  252,   329, -1260,  4410,
*     5            84,  -35,   28,  -42,   140,   609, -2940,
*     6           -10,    4,   -3,    4,   -10,    60,  1089/
*      data ia /  -60,   10,   -4,    3,    -4,    10,   -60/
*      data id / 420/
*
*     the following coefficients are suitable for a ninth-order method*
*      use with nine-point adams method, set NO = 8
*
      data ie / -1338,  -240,    60,   -32,    30,   -48,   140,  -960,
     1           2940,  -798,  -420,   168,  -140,   210,  -588,  3920,
     2          -2940,  1680,  -378,  -672,   420,  -560,  1470, -9408,
     3           2450, -1050,  1050,     0, -1050,  1050, -2450, 14700,
     4          -1470,   560,  -420,   672,   378, -1680,  2940,-15680,
     5            588,  -210,   140,  -168,   420,   798, -2940, 11760,
     6           -140,    48,   -30,    32,   -60,   240,  1338, -6720,
     7             15,    -5,     3,    -3,     5,   -15,   105,  2283/
      data ia /  -105,    15,    -5,     3,    -3,     5,   -15,   105/
      data id /   840/
*
      data job /01/, nd /NO/, lda /NO/

************************************************************************
*  a)    Lagrangian differentiation formulas
*
*             p(r) = r**(l+1) u(r)
*
*             du/dr = v(r)
*             dv/dr =  -2(l+1)/r v(r) -2[ e + z/r ] u(r)
*
*             u =    1     + a1 r + a2 r**2 + ...
*             v = -z/(l+1) + b1 r + b2 r**2 + ...
*
************************************************************************

       u0 = 1.0d0
       v0 = -z(1)/(l+1)

       p(1) = 0.0
       if(l.eq.0) then
          q(1) = 1.0d0
       else
          q(1) = 0.0d0
       endif 

       do 550 nl = 1,LO
         i0 = 1 + (nl-1)*NO
         do 150  i = 1,NO
             b(i) = id*h*rp(i+i0)
             c(i) = -2*id*h*(eau*rp(i+i0) + z(i+i0)*rpor(i+i0))
             d(i) = -2*id*h*(l+1)*rpor(i+i0)
             do 100 j = 1,NO
                if(j.eq.i) then
                   em(i,j) = ie(i,j) - d(i)
                else
                   em(i,j) = ie(i,j)
                endif
 100         continue
 150     continue

         call dgefa(em,lda,nd,ipvt,info)
         call dgedi(em,lda,nd,ipvt,det,work,job)

         do 250 i = 1,NO
            s(i) = -ia(i)*u0
            do 200 j = 1,NO
               fm(i,j) = ie(i,j) - b(i)*em(i,j)*c(j)
               s(i) = s(i) - b(i)*em(i,j)*ia(j)*v0
 200        continue
 250     continue

         call dgefa(fm,lda,nd,ipvt,info)
         call dgedi(fm,lda,nd,ipvt,det,work,job)

         do 350 i = 1,NO
            us(i) = 0.0d0
            do 300 j = 1,NO
               us(i) = us(i) + fm(i,j)*s(j)
 300        continue 
 350     continue

         do 450 i = 1,NO
            vs(i) = 0.0d0
            do 400 j = 1,NO
               vs(i) = vs(i) + em(i,j)*(c(j)*us(j)-ia(j)*v0)
 400        continue
 450     continue

         do 500 i = 1,NO
            p(i+i0) = r(i+i0)**(l+1)*us(i)
            q(i+i0) = r(i+i0)**l*(r(i+i0)*vs(i)+(l+1)*us(i))
 500     continue

         u0 = us(NO)
         v0 = vs(NO)

 550  continue  

***********************************************************************
*
*  b)   Continue solution outward using Adams method
*       n.b.  be sure NO is consistent in the two routines !
*
************************************************************************

         na = LO*NO + 2
         if(nb.gt.na) then
            call adams(p,q,z,eau,l,na,nb)
         endif

      return
      end

      function rint (f,na,nb,nq,h)
      implicit real*8(a-h,o-z)
c
c  this program calculates the integral of the function f from point na
c  to point nb using a nq points quadrature ( nq is any integer between
c  1 and 14 ).  h is the grid size.
c                                      written by c. c. j. roothaan
c
      dimension c(105),c1(78),c2(27),d(14),f(nb)
      equivalence (c1(1),c(1)) ,(c2(1),c(79))
      data c1/
     a 1.,
     b 2.,1.,
     c 23.,28.,9.,
     d 25.,20.,31.,8.,
     e 1413.,1586.,1104.,1902.,475.,
     f 1456.,1333.,1746.,944.,1982.,459.,
     g 119585.,130936.,89437.,177984.,54851.,176648.,36799.,
     h 122175.,111080.,156451.,46912.,220509.,29336.,185153.,35584.,
     i 7200319.,7783754.,5095890.,12489922.,-1020160.,16263486.,261166.,
     i 11532470.,2082753.,
     j 7305728 ,6767167.,9516362.,1053138.,18554050.,-7084288.,
     j 20306238.,-1471442.,11965622.,2034625.,
     k 952327935.,1021256716.,636547389.,1942518504.,-1065220914.,
     k 3897945600.,-2145575886.,3373884696.,-454944189.,1637546484.,
     k 262747265.,
     l 963053825.,896771060.,1299041091.,-196805736.,3609224754.,
     l-3398609664.,6231334350.,-3812282136.,4207237821.,-732728564.,
     l 1693103359., 257696640./
      data c2 / 5206230892907.,5551687979302.,3283609164916.,
     m 12465244770050.,-13155015007785.,39022895874876.,
     m-41078125154304.,53315213499588.,-32865015189975.,28323664941310.,
     m-5605325192308.,9535909891802.,1382741929621.,
     n 5252701747968.,4920175305323.,7268021504806.,-3009613761932.,
     n 28198302087170.,-41474518178601.,76782233435964.,
     n-78837462715392.,81634716670404.,-48598072507095.,
     n 34616887868158.,-7321658717812.,9821965479386.,1360737653653./
      data d/2.,2.,24.,24.,1440.,1440.,120960.,120960.,7257600.,
     a  7257600.,958003200.,958003200.,5230697472000.,5230697472000./

      a = 0.0
      l = na
      m = nb
      i = nq*(nq+1)/2
      do 100 j = 1,nq
         a = a + c(i)*( f(l) + f(m) )
         l = l + 1
         m = m - 1
         i = i - 1
 100  continue
      a = a/d(nq)
      do 200 n = l,m
        a = a + f(n)
 200  continue
      rint = a*h
      return
      end

      subroutine enclab(n,k,*)
c
c  this program encodes the principal and angular quantum numbers 'n'
c  and 'l' into an orbital label 'lab'
c
      implicit character*1(l)
      common/orblab/lab(4)
      dimension l1(10),l2(10),l3(10),l4(2)
      data l1/' ','1','2','3','4','5','6','7','8','9'/
      data l2/'0','1','2','3','4','5','6','7','8','9'/
      data l3/'s','p','d','f','g','h','i','j','k','l'/
      data l4/' ','*'/
      data nx/10/

      if(n.lt.1.or.n.gt.99) then
       write(6,1000)
 1000  format('Error in enclab: principal quantum number out of range')
       return 1
      endif
      n1=n/10+1
      n2=n-n1*10+11
      lab(1)=l1(n1)
      lab(2)=l2(n2)
      n3=k+1
      n4=1
      if(n3.lt.0.or.n3.gt.min0(n,nx)) then
         write(6,1010)
 1010    format('Error in enclab: angular quantum number out of range')
         return 1
      endif
      lab(3)=l3(n3)
      lab(4)=l4(n4)
      return
      end

****  contains routines dgedi - dgefa - dgesl 
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,*),det(2),work(*)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      subroutine dgesl(a,lda,n,ipvt,b,job)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
      integer lda,n,ipvt(*),job
      double precision a(lda,*),b(*)
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      double precision a(lda,*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end

      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
