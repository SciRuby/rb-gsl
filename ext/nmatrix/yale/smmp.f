c=======================================================================
c Sparse Matrix Multiplication Package
c
c Randolph E. Bank and Craig C. Douglas
c
c na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
c
c Compile this with the following command (or a similar one):
c
c     f77 -c -O smmp.f
c
c=======================================================================
        subroutine symbmm 
     *                   (n, m, l, 
     *                    ia, ja, diaga, 
     *                    ib, jb, diagb,
     *                    ic, jc, diagc,
     *                    index)
c
            integer       ia(*), ja(*), diaga,
     *                    ib(*), jb(*), diagb,
     *                    ic(*), jc(*), diagc,
     *                    index(*)
c
c       symbolic matrix multiply c=a*b
c
        maxlmn = max(l,m,n)
        do 10 i=1,maxlmn
   10       index(i)=0
        if (diagc.eq.0) then
            ic(1)=1
        else
            ic(1)=n+2
        endif
        minlm = min(l,m)
        minmn = min(m,n)
c
c    main loop
c
        do 50 i=1,n
            istart=-1
            length=0
c
c    merge row lists
c
            do 30 jj=ia(i),ia(i+1)
c    a = d + ...
                if (jj.eq.ia(i+1)) then
                    if (diaga.eq.0 .or. i.gt.minmn) goto 30
                    j = i
                else
                    j=ja(jj)
                endif
c    b = d + ...
                if (index(j).eq.0 .and. diagb.eq.1 .and. j.le.minlm)then
                    index(j)=istart
                    istart=j
                    length=length+1
                endif
                do 20 k=ib(j),ib(j+1)-1 
                    if(index(jb(k)).eq.0) then
                        index(jb(k))=istart
                        istart=jb(k)
                        length=length+1
                    endif
   20           continue
   30       continue
c
c   row i of jc
c
            if (diagc.eq.1 .and. index(i).ne.0) length = length - 1
            ic(i+1)=ic(i)+length
            do 40 j= ic(i),ic(i+1)-1
                if (diagc.eq.1 .and. istart.eq.i) then
                    istart = index(istart)
                    index(i) = 0
                endif
                jc(j)=istart
                istart=index(istart)
                index(jc(j))=0
   40       continue
            index(i) = 0
   50   continue
        return
        end
        subroutine numbmm 
     *                   (n, m, l,
     *                    ia, ja, diaga, a,
     *                    ib, jb, diagb, b,
     *                    ic, jc, diagc, c,
     *                    temp)
c
            integer       ia(*), ja(*), diaga,
     *                    ib(*), jb(*), diagb,
     *                    ic(*), jc(*), diagc 
c
            real          a(*), b(*), c(*), temp(*)
c
c       numeric matrix multiply c=a*b
c
        maxlmn = max(l,m,n)
        do 10 i = 1,maxlmn
 10         temp(i) = 0.
        minlm = min(l,m)
        minln = min(l,n)
        minmn = min(m,n)
c
c   c = a*b
c
        do 50 i = 1,n
             do 30 jj = ia(i),ia(i+1)
c    a = d + ...
                if (jj.eq.ia(i+1)) then
                    if (diaga.eq.0 .or. i.gt.minmn) goto 30
                    j = i
                    ajj = a(i)
                else
                    j=ja(jj)
                    ajj = a(jj)
                endif
c    b = d + ...
                if (diagb.eq.1 .and. j.le.minlm) 
     *              temp(j) = temp(j) + ajj * b(j)
                do 20 k = ib(j),ib(j+1)-1
 20                 temp(jb(k)) = temp(jb(k)) + ajj * b(k)
 30         continue
c    c = d + ...
            if (diagc.eq.1 .and. i.le.minln) then
                c(i) = temp(i)
                temp(i) = 0.
            endif
            do 40 j = ic(i),ic(i+1)-1
                c(j) = temp(jc(j))
 40             temp(jc(j)) = 0.
 50     continue
        return
        end
        subroutine transp
     *                   (n, m,
     *                    ia, ja, diaga, a,
     *                    ib, jb,        b,
     *                    move)
c
            integer       ia(*), ja(*), diaga,
     *                    ib(*), jb(*),
     *                    move  
c
            real          a(*), b(*)
c
c       compute b = a(transpose)
c
c       first make ib
c
        do 10 i=1,m+1
   10       ib(i)=0
        if (move.eq.1) then
            do 15 i =1,m+1
   15           b(i) = 0.
        endif
        if (diaga.eq.1) then
            ib(1)=m + 2
        else
            ib(1)=1
        endif
c
c       count indices for each column 
c
        do 30 i=1,n   
            do 20 j=ia(i),ia(i+1)-1
                ib(ja(j)+1)=ib(ja(j)+1)+1
   20       continue
   30   continue
        do 40 i=1,m
   40      ib(i+1)=ib(i)+ib(i+1)
c
c       now make jb
c
        do 60 i=1,n   
            do 50 j=ia(i),ia(i+1)-1
                index=ja(j)
                jb(ib(index))=i
                if (move.eq.1) b(ib(index)) = a(j)
                ib(index)=ib(index)+1
   50       continue
   60   continue
c
c       now fixup ib 
c
        do 70 i=m,2,-1
   70       ib(i)=ib(i-1)
        if (diaga.eq.1) then
            if (move.eq.1) then
                j = min(n,m)
                do 80 i = 1,j
   80               b(i) = a(i)
            endif
            ib(1)=m + 2
        else
            ib(1)=1
        endif
        return
        end
        subroutine ytobs
     *                  (n,
     *                   ia, ja, diaga, syma, a,
     *                   ib, jb,              b,
     *                   move)
c
            integer     ia(*), ja(*), diaga, syma,
     *                  ib(*), jb(*), move
c
            real        a(*), b(*)
c
c       create the bank-smith data structures b from the
c       corresponding yale data structures a
c
        do 10 i=1,n
   10       ib(i+1)=ia(i+1)-ia(i)
c
c       look for upper triangular entries and duplicate entries
c
        do 50 i=1,n
            do 40 jj=ia(i),ia(i+1)-1
                j=ja(jj)
                if (i.eq.j) then
                    ib(i+1)=ib(i+1)-1
                    ja(jj) = -j
                endif
                if(j.gt.i) then
                    ib(i+1)=ib(i+1)-1
                    ib(j+1)=ib(j+1)+1
c
c       check for duplicates
c
                    do 20 k=ia(j),ia(j+1)-1
                        if(ja(k).eq.i) then
                            ib(j+1)=ib(j+1)-1
                            ja(jj)=-j
                            go to 30
                        endif
   20               continue
   30               continue
                endif
   40       continue
   50   continue
c
c       compute ib
c
        ib(1)=n + 2
        do 60 i=1,n
   60       ib(i+1)=ib(i+1)+ib(i)
c
c       initialize b if move = 1
c
        if (move.eq.1) then
            lshift = 0
            if (syma.eq.0) lshift = ib(n+1) - ib(1)
            do 62 ii = 1,ib(n+1)+lshift-1
   62           b(ii) = 0.
            if (diaga.eq.1) then
                do 64 ii = 1,n
   64               b(ii) = a(ii)
            endif
        endif
c
c       compute jb
c
        do 80 i=1,n
            do 70 jj=ia(i),ia(i+1)-1
                j=ja(jj)
                if(j.gt.i) then
                    jb(ib(j))=i
                    if (move.eq.1) b(ib(j)) = a(jj)
                    ib(j)=ib(j)+1
                else
                    if(j.le.0) then
                        ja(jj)=-j
                        if (move.eq.1 .and. i.eq.-j) b(i) = a(jj)
                    else
                        jb(ib(i))=j
                        if (move.eq.1) b(ib(i)+lshift) = a(jj)
                        ib(i)=ib(i)+1
                    endif
                endif
   70       continue
   80   continue
c
c       fixup ib
c
        do 90 i=n,2,-1
   90       ib(i)=ib(i-1)
        ib(1)=n + 2
        return
        end
        subroutine bstoy
     *                  (n,
     *                   ia, ja,        syma, a,
     *                   ib, jb, diagb,       b,
     *                   move)
c
            integer      ia(*), ja(*),        syma,
     *                   ib(*), jb(*), diagb, 
     *                   move
c
            real         a(*), b(*)
c
c       create the yale data structures b from the
c       corresponding bank-smith data structures a
c
c       compute ib
c
        if (diagb.eq.1) then
            ib(1) = n + 2
            icor = 0
            if (move.eq.1) then
                lshift = 0
                if (syma.eq.0) lshift = ia(n+1) - ia(1)
                do 2 i = 1,n
    2               b(i) = a(i)
            endif
        else
            ib(1) = 1
            icor = 1
        endif
        do 10 i=1,n
   10       ib(i+1)=ia(i+1)-ia(i)+icor
        do 30 i=1,n
            do 20 j=ia(i),ia(i+1)-1
                ib(ja(j)+1)=ib(ja(j)+1)+1
   20       continue
   30   continue
c
        do 40 i=1,n
   40       ib(i+1)=ib(i+1)+ib(i)
        if (diagb.eq.0) then
            do 45 i = 1,n
                jb(ib(i)) = i
                if (move.eq.1) b(ib(i)) = a(i)
   45           ib(i) = ib(i) + 1
        endif
c
c       now compute jb
c
        do 60 i=1,n
            do 50 jj=ia(i),ia(i+1)-1
                j = ja(jj)
                jb(ib(j))=i
                jb(ib(i))=j
                if (move.eq.1) then
                    b(ib(j)) = a(jj)
                    b(ib(i)) = a(jj+lshift)
                endif
                ib(i)=ib(i)+1
                ib(j)=ib(j)+1
   50       continue
   60   continue
c
c       fixup ib
c
        do 70 i=n,2,-1
   70       ib(i)=ib(i-1)
        if (diagb.eq.1) then
            ib(1)=n+2
        else
            ib(1)=1
        endif
        return
        end
