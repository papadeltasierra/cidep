$LARGE
      subroutine first(same,nflag)
c
c.....Purpose : to calc the degeneracies for a 1st order esr spectrum
c
c.....Input   : ngroup(2),quant(20,2),nnuc(20,2)
c
c.....Output  : ideg(20,20,2)
c
c.....jzcr,qn(k,j,i),nk   : counters
c
c.....N.B. highest "energy" states lie at "lowest" field.
c
      common /b/quant(20,2),a(20,2),nnuc(20,2),lines(20,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /f/ideg(20,20,2),nbuff(30),nk(10),jzcr(20),
     &          sumstm(1024,2),stmnet(2)
      common /h/qn(20,20,2),stmpol(20,20,2)
c
      character same*1
c
      i=0
1     i=i+1
      ngi=ngroup(i)
      do 11 j=1,ngi
        intq=int(2.0*quant(j,i))
        lines(j,i)=nnuc(j,i)*intq+1
        naji=lines(j,i)
        var1=-1.0-float(nnuc(j,i))*quant(j,i)
        do 6 k=1,naji
          ideg(k,j,i)=0
6         qn(k,j,i)=float(k)+var1
        do 7 ncr=1,nnuc(j,i)
7         jzcr(ncr)=1
        njz=intq+1
8       n=1
        k=1-nnuc(j,i)
        do 9 ncr=1,nnuc(j,i)
9         k=k+jzcr(ncr)
        ideg(k,j,i)=ideg(k,j,i)+1
        jzcr(n)=jzcr(n)+1
        if(jzcr(n)-njz) 8,8,10
10      jzcr(n)=1
        n=n+1
        if(n-nnuc(j,i).gt.0) goto 11
        jzcr(n)=jzcr(n)+1
        if(jzcr(n)-njz) 8,8,10
c.....  nflag=1
11    continue
      if(i.eq.2) return
      if(same.eq.'N') goto 1
c.....radicals not identical
      do 16 j=1,ngi
        lines(j,2)=lines(j,1)
        naji=lines(j,1)
        do 15 k=1,naji
          qn(k,j,2)=qn(k,j,1)
15        ideg(k,j,2)=ideg(k,j,1)
16    continue
      return
      end
c
      subroutine second(same)
c
c.....Purpose : to calc the degeneracies for a second order esr spectrum
c.....          (and st-minus contributions on each component if istm(i)=1)
c
c.....Input   : ngroup(2),qn(20,20,2),ideg(20,20,2),same
c
c.....Output  : ideg(20,20,2)
c
      common /b/quant(20,2),a(20,2),nnuc(20,2),lines(20,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /f/ideg(20,20,2),nbuff(30),nk(10),jzcr(20),
     &          sumstm(1024,2),stmnet(2)
      common /h/qn(20,20,2),stmpol(20,20,2)
c
      character same*1
      real jmax,jm,jjm
      i=0
1     i=i+1
      if(istm(i).eq.1) nsec(i)=2
      if(nsec(i).ne.2) goto 75
      ngi=ngroup(i)
c
      do 70 j=1,ngi
c.........determine j(max) and set up counter for j-values....
        jmax=abs(qn(1,j,i))
        jctr=int(jmax+1.01)
c.........set up degeneracies nbuff(jk) for each j-quantum number....
        do 40 jk=1,jctr
          if(jk-1) 20,20,10
10        nbuff(jk)=iabs(ideg(jk,j,i)-ideg(jk-1,j,i))
          goto 30
20        nbuff(jk)=1
30        qm(jk)=abs(qn(jk,j,i))
40      continue
        isos=0
c.........initialise counter for the m-quantum number and value of j(max)....
        mctr=int(2.0*jmax+1.01)
        jjm=jmax
c.........determine all the lines within the group and place in qn()..m-values,
c.....ideg()..degeneracies, qj()..j-values....
        do 60 jk=1,jctr
        jm=jjm
        do 50 mk=1,mctr
        isos=isos+1
        qn(isos,j,i)=-jm
        ideg(isos,j,i)=nbuff(jk)
        qj(isos,j,i)=qm(jk)
c
        if(istm(i).ne.1) goto 50
c.....  goto 50 if st-minus not wanted
        stmpol(isos,j,i)=(qj(isos,j,i)-qn(isos,j,i)+1)
        stmpol(isos,j,i)=stmpol(isos,j,i)*(qj(isos,j,i)+qn(isos,j,i))
        stmpol(isos,j,i)=(stmpol(isos,j,i))*a(j,i)**2
44      stmpol(isos,j,i)=stmpol(isos,j,i)*ideg(isos,j,i)
c
50      jm=jm-1
        mctr=mctr-2
60      jjm=jjm-1
70      lines(j,i)=isos
75    if(i.eq.2) return
c
      if(same.eq.'N') goto 1
      do 90 j=1,ngi
        lines(j,2)=lines(j,1)
        isos=lines(j,1)
        do 80 k=1,isos
          qn(k,j,2)=qn(k,j,1)
          qj(k,j,2)=qj(k,j,1)
80      ideg(k,j,2)=ideg(k,j,1)
90    continue
      return
      end
c
      subroutine ami(same,ef)
c
c.....Purpose : having calculated the degeneracies, from either first or
c.....          second, go on to find the total degeneracy of each line
c.....          (itotal(i,k)) and its hyperfine energy (sumami(k,i))
c.....          (and the st-minus enhancement if istm(i)=1).
c
c.....Input   : ideg(20,20,2),lines(20,2)
c
c.....Input from Cifer  : swp, freq, stfld, npt
c
c.....Output  : itotal(1024,2), sumami(1024,2), nlines(2)
c
c.....          jzcr, qn(k,j,i), nk - counters
c
      common /b/quant(20,2),a(20,2),nnuc(20,2),lines(20,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /f/ideg(20,20,2),nbuff(30),nk(10),jzcr(20),
     &          sumstm(1024,2),stmnet(2)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /h/qn(20,20,2),stmpol(20,20,2)
c
      character same*1,esc*1
      esc=char(27)
c
      print 100
      read(*,110,err=9001) swpD
      if(swpD.ne.0.0) swp=swpD
9001  print *,esc,'[2A'
      print 101,swp
      read(*,111,err=9002) stfldD
      if(stfldD.ne.0.0) stfld=stfldD
9002  print *,esc,'[2A'
      print 102,swp,stfld
      print *,'   '
      print 200
      read(*,210,err=9003) freqD
      if(freqD.ne.0.0) freq=freqD
9003  print *,esc,'[2A'
      print 201,freq
      read(*,211,err=9004) nptD
      if(nptD.ne.0) npt=nptD
9004  print *,esc,'[2A'
      print 202,freq,npt
      print *,'   '
c
      i=0
50    i=i+1
      lx=1
      nlines(i)=1
      ngi=ngroup(i)
      do 140 l=1,ngi
        nlines(i)=nlines(i)*lines(l,i)
140     nk(l)=1
150   n=1
      sumami(lx,i)=0.0
      itotal(lx,i)=1
      sumstm(lx,i)=0.0
      setfld=0.5/stfld
c
      do 160 j1=1,ngi
        nkj=nk(j1)
        sumami(lx,i)=sumami(lx,i)+a(j1,i)*qn(nkj,j1,i)
        if(nsec(i).eq.2) sumami(lx,i)=sumami(lx,i)
     &             + (a(j1,i)**2)*(qj(nkj,j1,i)*(qj(nkj,j1,i)+1.0)
     &             - qn(nkj,j1,i)**2)*setfld
        if(istm(i).eq.0) goto 160
c.....if st-minus not wanted goto 160
        sumstm(lx,i)=sumstm(lx,i)+stmpol(nkj,j1,i)
160     itotal(lx,i)=itotal(lx,i)*ideg(nkj,j1,i)
c
      lx=lx+1
c
      if(lx.gt.1024) then
        print *,' TOO MANY STATES '
        ef=1
        return
      end if
c
      nk(n)=nk(n)+1
170   if((nk(n)-lines(n,i)).le.0) goto 150
      if(n-ngi.ne.0) goto 180
      if(nk(n)-lines(ngi,i).gt.0) goto 190
180   nk(n)=1
      n=n+1
      nk(n)=nk(n)+1
      goto 170
190   continue
c
      nli=nlines(i)
      if(i.eq.2) return
      if(same.eq.'N') goto 50
c.....radicals not identical
      nlines(2)=nlines(1)
      do 220 j=1,nli
        itotal(j,2)=itotal(j,1)
220     sumami(j,2)=sumami(j,1)
      return
c
100   format(' ',' field sweep(Gauss)    ',\)
110   format(F7.0)
101   format(' ',' field sweep(Gauss)    ',F7.2,2X,
     &'; centre field(Gauss)   ',\)
102     format(' ',' field sweep(Gauss)    ',F7.2,2X,
     &'; centre field(Gauss)   ',F7.1,\)
111   format(F7.0)
200   format(' ',' freq(GHz)             ',\)
210   format(F6.0)
201   format(' ',' freq(GHz)              ',F6.4,'  ; no. of points  '
     &,10X,\)
202   format(' ',' freq(GHz)              ',F6.4,'  ; no. of points  '
     &,10X,I4,/,\)
211   format(I4)
c
      end

