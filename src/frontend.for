      SUBROUTINE ATOM(Z, C)
      parameter(NGP=500,NOR=30)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      jz = Z
      jx = 1
      jc = C
      ns (NOR) = 0 ! placeholder for ncore
      END

      SUBROUTINE ADD_ORBITAL(N, L, O)
      parameter(NGP=500,NOR=30)
      common/intdat/jx,jz,jc,ns(NOR),ls(NOR),js(NOR),ms(NOR),ks(NOR)
      common/labdat/idn,lab(NOR)
      common/fixdat/wi(NOR),wf(NOR)
      common/orblab/label
      ns(jx) = N
      ls(jx) = L
      js(jx) = O
      call enclab(N, L, *901)
      lab(jx) = label
      if (jc.eq.0) then
         js(jx) = 4*ls(jx) + 2
         ns(NOR) = ns(NOR) + js(jx)
         wi(jx)=-jz**2/(2.0d0*N**2)
      else
         if (ja.gt.jc) then
            js(jx) = 1
            wi(jx) = -(jz-ncore)**2/(2.0d0*N**2)
         endif
      endif
      jx = jx + 1
      return
 901  return 1
      end

      SUBROUTINE SETGRID(R0, IH)
      parameter(NGP=500,NOR=30)
      character*4 idn,lab,label
      common/radial/r(NGP),rp(NGP),rpor(NGP),h
      common/charge/z(NGP)

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
