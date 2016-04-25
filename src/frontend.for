      SUBROUTINE ATOM(IZ, C)
      implicit real*8(a-h,o-z)
      integer C, IZ
      parameter(NGP=500,NOR=30)
      character*4 idn,lab,label
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/labdat/idn,lab(NOR)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/charge/z(NGP)
      common/orblab/label
      data r0def /0.0005d0/, hdef/0.03d0/
      jz = IZ
      jx = 1
      jc = C
      ns (NOR) = 0 ! placeholder for ncore
      END

      SUBROUTINE ADD_ORBITAL(N, L, O)
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
      ns(jx) = N
      ls(jx) = L
      js(jx) = O
      call enclab(N, L, *901)
      lab(jx) = label
      if (jc.eq.0.or.jx.lt.jc) then
         js(jx) = 4*ls(jx) + 2
         ns(NOR) = ns(NOR) + js(jx)
         wi(jx)=-jz**2/(2.0d0*N**2)
      else
         if (jx.ge.jc) then
            js(jx) = 0
            wi(jx) = -(jz-ncore)**2/(2.0d0*N**2)
         endif
      endif
      jx = jx + 1
      return
 901  return 1
      end

      SUBROUTINE SET_GRID(R0, IH)
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

      ns (NOR) = 0 ! clean placeholder ncore
      if (jc.eq.0) then
         jc = jx - 1
      ENDIF
      jx = jx - 1
      if(r0.eq.0.0)  r0 = r0def
      if(h.eq.0.0)   h = hdef
      else h = ih
      r(1)=0d0
      rp(1)=r0
      rpor(1)=0d0
      z(1) = jz
      do 200 i=2,NGP
         rp(i) = r0*exp((i-1)*h)
         r(i) = rp(i) - r0
         rpor(i) = rp(i)/r(i)
         z(i) = jz
 200  continue
      end
