
!!$*DECK DUMACH

function dumach ()
real*8::dumach
real*8::u,comp
u = 1.0d0
do
u = u*0.5d0
call dumsum(1.0d0,u,comp)
if (comp .eq. 1.0d0) exit
end do
dumach = u*2.0d0
return
end function dumach

!*******************************
subroutine dumsum(a,b,c)
real*8::a,b,c
c = a + b
return
end subroutine dumsum


!!$*DECK DCFODE

subroutine dcfode (meth,elco,tesco)
integer::meth
integer::i,ib,nq,nqm1,nqp1
real*8::elco(13,12),tesco(3,12)
real*8::agamq,fnq,fnqm1,pc(12),pint,ragq,&
rqfac,rq1fac,tsign,xpin

if(meth==1) then
elco(1,1) = 1.0d0
elco(2,1) = 1.0d0
tesco(1,1) = 0.0d0
tesco(2,1) = 2.0d0
tesco(1,2) = 1.0d0
tesco(3,12) = 0.0d0
pc(1) = 1.0d0
rqfac = 1.0d0
do nq = 2,12
rq1fac = rqfac
rqfac = rqfac/nq
nqm1 = nq - 1
fnqm1 = nqm1
nqp1 = nq + 1
pc(nq) = 0.0d0
do ib = 1,nqm1
i = nqp1 - ib
pc(i) = pc(i-1) + fnqm1*pc(i)
end do
pc(1) = fnqm1*pc(1)
pint = pc(1)
xpin = pc(1)/2.0d0
tsign = 1.0d0
do i = 2,nq
tsign = -tsign
pint = pint + tsign*pc(i)/i
xpin = xpin + tsign*pc(i)/(i+1)
end do
elco(1,nq) = pint*rq1fac
elco(2,nq) = 1.0d0
do i = 2,nq
elco(i+1,nq) = rq1fac*pc(i)/i
end do
agamq = rqfac*xpin
ragq = 1.0d0/agamq
tesco(2,nq) = ragq
if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/nqp1
tesco(3,nqm1) = ragq
end do
return
elseif(meth==2) then
pc(1) = 1.0d0
rq1fac = 1.0d0
do nq = 1,5
fnq = nq
nqp1 = nq + 1
pc(nqp1) = 0.0d0
do ib = 1,nq
i = nq + 2 - ib
pc(i) = pc(i-1) + fnq*pc(i)
end do
pc(1) = fnq*pc(1)
do i = 1,nqp1
elco(i,nq) = pc(i)/pc(2)
end do
elco(2,nq) = 1.0d0
tesco(1,nq) = rq1fac
tesco(2,nq) = nqp1/elco(1,nq)
tesco(3,nq) = (nq+2)/elco(1,nq)
rq1fac = rq1fac/fnq
end do
else
print *,"ERROR_90: WRONG METH",METH
stop
end if
return
end subroutine dcfode

!!***************************************
!!$*DECK DINTDY
subroutine dintdy (t,k,yh,nyh,dky,iflag)
integer::k,nyh,iflag
real*8::t,yh(nyh,*),dky(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,ic,j,jb,jb2,jj,jj1,jp1
real*8::c,r,s,tp
character*80 msg
iflag = 0
if (k .lt. 0 .or. k .gt. nq) then
msg = 'dintdy- k (=i1) illegal '
call xerrwd (msg,30,51,0,1,k,0,0,0.0d0,0.0d0)
iflag = -1
return
end if

tp = tn - hu - 100.0d0*uround*sign(abs(tn) + abs(hu),hu)

if ((t-tp)*(t-tn) .gt. 0.0d0) then
msg = 'dintdy- t (=r1) illegal '
call xerrwd (msg,30,52,0,0,0,0,1,t,0.0d0)
msg=' t not in interval tcur - hu (= r1) to tcur (=r2) '
call xerrwd (msg,60,52,0,0,0,0,2,tp,tn)
iflag = -2
return
end if
s = (t - tn)/h
ic = 1
if (k .ne. 0) then
jj1 = l - k
do jj = jj1,nq
ic = ic*jj
end do
end if
c = ic
do i = 1,n
dky(i) = c*yh(i,l)
end do
if (k .ne. nq) then
jb2 = nq - k
do jb = 1,jb2
j = nq - jb
jp1 = j + 1
ic = 1
if (k .eq. 0) then
jj1 = jp1 - k
do jj = jj1,j
ic = ic*jj
end do
end if
c = ic
do i = 1,n
dky(i) = c*yh(i,jp1) + s*dky(i)
end do
end do
if (k .eq. 0) return
end if
r = h**(-k)
do i = 1,n
dky(i) = r*dky(i)
end do
return
end subroutine dintdy
!!*************************
!!$*DECK DPREPJ
subroutine dprepj (neq,y,yh,nyh,ewt,ftem,savf,wm,iwm,&
f,jac)
external f,jac
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),ewt(*),ftem(*),savf(*),wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,i1,i2,ier,ii,j,j1,jj,lenp,&
mba,mband,meb1,meband,ml,ml3,mu,np1
real*8::con,di,fac,hl0,r,r0,srur,yi,yj,yjj,&
dvnorm
nje = nje + 1
ierpj = 0
jcur = 1
hl0 = h*el0
if(miter==1) then
lenp = n*n
do i = 1,lenp
wm(i+2) = 0.0d0
end do
call jac (neq,tn,y,0,0,wm(3),n)
con = -hl0
do i = 1,lenp
wm(i+2) = wm(i+2)*con
end do
j = 3
np1 = n + 1
do i = 1,n
wm(j) = wm(j) + 1.0d0
j = j + np1
end do
call dgefa (wm(3),n,n,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
elseif(miter==2) then
fac = dvnorm (n,savf,ewt)
r0 = 1000.0d0*abs(h)*uround*n*fac
if (r0 .eq. 0.0d0) r0 = 1.0d0
srur = wm(1)
j1 = 2
do j = 1,n
yj = y(j)
r = max(srur*abs(yj),r0/ewt(j))
y(j) = y(j) + r
fac = -hl0/r
call f (neq,tn,y,ftem)
do i = 1,n
wm(i+j1) = (ftem(i) - savf(i))*fac
end do
y(j) = yj
j1 = j1 + n
end do
nfe = nfe + n
j = 3
np1 = n + 1
do i = 1,n
wm(j) = wm(j) + 1.0d0
j = j + np1
end do
call dgefa (wm(3),n,n,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
elseif(miter==3) then
wm(2) = hl0
r = el0*0.1d0
do i = 1,n
y(i) = y(i) + r*(h*savf(i) - yh(i,2))
end do
call f (neq,tn,y,wm(3))
nfe = nfe + 1
do i = 1,n
r0 = h*savf(i) - yh(i,2)
di = 0.1d0*r0 - h*(wm(i+2) - savf(i))
wm(i+2) = 1.0d0
if (abs(r0) .lt. uround/ewt(i)) cycle
if (abs(di) .eq. 0.0d0) exit
wm(i+2) = 0.1d0*r0/di
end do
ierpj = 1
return
elseif(miter==400) then
ml = iwm(1)
mu = iwm(2)
ml3 = ml + 3
mband = ml + mu + 1
meband = mband + ml
lenp = meband*n
do i = 1,lenp
wm(i+2) = 0.0d0
end do
call jac (neq,tn,y,ml,mu,wm(ml3),meband)
con = -hl0
do i = 1,lenp
wm(i+2) = wm(i+2)*con
end do
ii = mband + 2
do i = 1,n
wm(ii) = wm(ii) + 1.0d0
ii = ii + meband
end do
call dgbfa (wm(3),meband,n,ml,mu,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
elseif(miter==5) then
ml = iwm(1)
mu = iwm(2)
mband = ml + mu + 1
mba = min(mband,n)
meband = mband + ml
meb1 = meband - 1
srur = wm(1)
fac = dvnorm (n,savf,ewt)
r0 = 1000.0d0*abs(h)*uround*n*fac
if (r0 .eq. 0.0d0) r0 = 1.0d0
do j = 1,mba
do i = j,n,mband
yi = y(i)
r = max(srur*abs(yi),r0/ewt(i))
y(i) = y(i) + r
end do
call f (neq,tn,y,ftem)
do jj = j,n,mband
y(jj) = yh(jj,1)
yjj = y(jj)
r = max(srur*abs(yjj),r0/ewt(jj))
fac = -hl0/r
i1 = max(jj-mu,1)
i2 = min(jj+ml,n)
ii = jj*meb1 - ml + 2
do i = i1,i2
wm(ii+i) = (ftem(i) - savf(i))*fac
end do
end do
end do
nfe = nfe + mba
ii = mband + 2
do i = 1,n
wm(ii) = wm(ii) + 1.0d0
ii = ii + meband
end do
call dgbfa (wm(3),meband,n,ml,mu,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
end if
end subroutine dprepj

!****************************
!!$*DECK DSOLSY

subroutine dsolsy (wm,iwm,x,tem)
integer::iwm(*)
real*8:: wm(*),x(*),tem(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,meband,ml,mu
real*8::di,hl0,phl0,r
iersl = 0
if(miter==1.or.miter==2) then
call dgesl (wm(3),n,n,iwm(21),x,0)
return
elseif(miter==3) then
phl0 = wm(2)
hl0 = h*el0
wm(2) = hl0
if (hl0 .ne. phl0) then
r = hl0/phl0
do i = 1,n
di = 1.0d0 - r*(1.0d0 - 1.0d0/wm(i+2))
if (abs(di) .eq. 0.0d0) go to 390
wm(i+2) = 1.0d0/di
end do
end if
do i = 1,n
x(i) = wm(i+2)*x(i)
end do
return
390 iersl = 1
return
elseif(miter==4) then
ml = iwm(1)
mu = iwm(2)
meband = 2*ml + mu + 1
call dgbsl (wm(3),meband,n,ml,mu,iwm(21),x,0)
return
end if
end subroutine dsolsy

!******************************
!!$*DECK DSRCOM
subroutine dsrcom (rsav,isav,job)
integer::isav(*),job
integer::ils
integer::i,lenils,lenrls
real*8::rsav(*),rls
save lenrls,lenils
common /dls001/ rls(218),ils(37)
data lenrls/218/,lenils/37/
if (job .ne. 2) then
do i = 1,lenrls
rsav(i) = rls(i)
end do
do i = 1,lenils
isav(i) = ils(i)
end do
return
end if

do i = 1,lenrls
rls(i) = rsav(i)
end do
do i = 1,lenils
ils(i) = isav(i)
end do
return
end subroutine dsrcom

!!******************************
!!$*DECK DSTODE
subroutine dstode (neq,y,yh,nyh,yh1,ewt,savf,acor,&
wm,iwm,f,jac,pjac,slvs)
external f,jac,pjac,slvs
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),&
acor(*),wm(*)
integer::iownd,ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,i1,iredo,iret,j,jb,m,ncf,newq
real*8::conit,crate,el,elco,hold,rmax,tesco,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::dcon,ddn,del,delp,dsm,dup,exdn,exsm,exup,&
r,rh,rhdn,rhsm,rhup,told,dvnorm
common /dls001/ conit,crate,el(13),elco(13,12),&
hold,rmax,tesco(3,12),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
kflag = 0
told = tn
ncf = 0
ierpj = 0
iersl = 0
jcur = 0
icf = 0
delp = 0.0d0
if (jstart .gt. 0) go to 200
if (jstart .eq. -1) go to 100
if (jstart .eq. -2) go to 160
lmax = maxord + 1
nq = 1
l = 2
ialth = 2
rmax = 10000.0d0
rc = 0.0d0
el0 = 1.0d0
crate = 0.7d0
hold = h
meo = meth
nslp = 0
ipup = miter
iret = 3
go to 140
100 ipup = miter
lmax = maxord + 1
if (ialth .eq. 1) ialth = 2
if (meth .eq. meo) go to 110
call dcfode (meth,elco,tesco)
meo = meth
if (nq .gt. maxord) go to 120
ialth = l
iret = 1
go to 150
110 if (nq .le. maxord) go to 160
120 nq = maxord
l = lmax
do i = 1,l
el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
ddn = dvnorm (n,savf,ewt)/tesco(1,l)
exdn = 1.0d0/l
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
rh = min(rhdn,1.0d0)
iredo = 3
if (h .eq. hold) go to 170
rh = min(rh,abs(h/hold))
h = hold
go to 175
140 call dcfode (meth,elco,tesco)
150 do i = 1,l
el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
go to (160,170,200),iret
160 if (h .eq. hold) go to 200
rh = h/hold
h = hold
iredo = 3
go to 175
170 rh = max(rh,hmin/abs(h))
175 rh = min(rh,rmax)
rh = rh/max(1.0d0,abs(h)*hmxi*rh)
r = 1.0d0
do j = 2,l
r = r*rh
do i = 1,n
180   yh(i,j) = yh(i,j)*r
end do
end do
h = h*rh
rc = rc*rh
ialth = l
if (iredo .eq. 0) go to 690
200 if (abs(rc-1.0d0) .gt. ccmax) ipup = miter
if (nst .ge. nslp+msbp) ipup = miter
tn = tn + h
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
210   yh1(i) = yh1(i) + yh1(i+nyh)
end do
215 continue
end do
220 m = 0
do i = 1,n
230 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
if (ipup .le. 0) go to 250
call pjac (neq,y,yh,nyh,ewt,acor,savf,wm,iwm,f,jac)
ipup = 0
rc = 1.0d0
nslp = nst
crate = 0.7d0
if (ierpj .ne. 0) go to 430
250 do i = 1,n
260 acor(i) = 0.0d0
end do
270 if (miter .ne. 0) go to 350
do i = 1,n
savf(i) = h*savf(i) - yh(i,2)
290 y(i) = savf(i) - acor(i)
end do
del = dvnorm (n,y,ewt)
do i = 1,n
y(i) = yh(i,1) + el(1)*savf(i)
300 acor(i) = savf(i)
end do
go to 400
350 do i = 1,n
360 y(i) = h*savf(i) - (yh(i,2) + acor(i))
end do
call slvs (wm,iwm,y,savf)
if (iersl .lt. 0) go to 430
if (iersl .gt. 0) go to 410
del = dvnorm (n,y,ewt)
do i = 1,n
acor(i) = acor(i) + y(i)
380 y(i) = yh(i,1) + el(1)*acor(i)
end do
400 if (m .ne. 0) crate = max(0.2d0*crate,del/delp)
dcon = del*min(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
if (dcon .le. 1.0d0) go to 450
m = m + 1
if (m .eq. maxcor) go to 410
if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
delp = del
call f (neq,tn,y,savf)
nfe = nfe + 1
go to 270
410 if (miter .eq. 0 .or. jcur .eq. 1) go to 430
icf = 1
ipup = miter
go to 220
430 icf = 2
ncf = ncf + 1
rmax = 2.0d0
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
440   yh1(i) = yh1(i) - yh1(i+nyh)
end do
445 continue
end do
if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
if (abs(h) .le. hmin*1.00001d0) go to 670
if (ncf .eq. mxncf) go to 670
rh = 0.25d0
ipup = miter
iredo = 1
go to 170
450 jcur = 0
if (m .eq. 0) dsm = del/tesco(2,nq)
if (m .gt. 0) dsm = dvnorm (n,acor,ewt)/tesco(2,nq)
if (dsm .gt. 1.0d0) go to 500
kflag = 0
iredo = 0
nst = nst + 1
hu = h
nqu = nq
do j = 1,l
do i = 1,n
470   yh(i,j) = yh(i,j) + el(j)*acor(i)
end do
end do
ialth = ialth - 1
if (ialth .eq. 0) go to 520
if (ialth .gt. 1) go to 700
if (l .eq. lmax) go to 700
do i = 1,n
490 yh(i,lmax) = acor(i)
end do
go to 700
500 kflag = kflag - 1
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
510   yh1(i) = yh1(i) - yh1(i+nyh)
end do
515 continue
end do
rmax = 2.0d0
if (abs(h) .le. hmin*1.00001d0) go to 660
if (kflag .le. -3) go to 640
iredo = 2
rhup = 0.0d0
go to 540
520 rhup = 0.0d0
if (l .eq. lmax) go to 540
do i = 1,n
530 savf(i) = acor(i) - yh(i,lmax)
end do
dup = dvnorm (n,savf,ewt)/tesco(3,nq)
exup = 1.0d0/(l+1)
rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
540 exsm = 1.0d0/l
rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rhdn = 0.0d0
if (nq .eq. 1) go to 560
ddn = dvnorm (n,yh(1,l),ewt)/tesco(1,nq)
exdn = 1.0d0/nq
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
560 if (rhsm .ge. rhup) go to 570
if (rhup .gt. rhdn) go to 590
go to 580
570 if (rhsm .lt. rhdn) go to 580
newq = nq
rh = rhsm
go to 620
580 newq = nq - 1
rh = rhdn
if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
go to 620
590 newq = l
rh = rhup
if (rh .lt. 1.1d0) go to 610
r = el(l)/l
do i = 1,n
600 yh(i,newq+1) = acor(i)*r
end do
go to 630
610 ialth = 3
go to 700
620 if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
if (kflag .le. -2) rh = min(rh,0.2d0)
if (newq .eq. nq) go to 170
630 nq = newq
l = nq + 1
iret = 2
go to 150
640 if (kflag .eq. -10) go to 660
rh = 0.1d0
rh = max(hmin/abs(h),rh)
h = h*rh
do i = 1,n
645 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
do i = 1,n
650 yh(i,2) = h*savf(i)
end do
ipup = miter
ialth = 5
if (nq .eq. 1) go to 200
nq = 1
l = 2
iret = 3
go to 150
660 kflag = -1
go to 720
670 kflag = -2
go to 720
680 kflag = -3
go to 720
690 rmax = 10.0d0
700 r = 1.0d0/tesco(2,nqu)
do i = 1,n
710 acor(i) = acor(i)*r
end do
720 hold = h
jstart = 1
return
end subroutine dstode

!!$*DECK DEWSET

subroutine dewset (n,itol,rtol,atol,ycur,ewt)
integer::n,itol
integer::i
real*8::rtol(*),atol(*),ycur(n),ewt(n)
go to (10,20,30,40),itol
10 continue
do i = 1,n
15  ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
end do
return
20 continue
do i = 1,n
25  ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
end do
return
30 continue
do i = 1,n
35  ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
end do
return
40 continue
do i = 1,n
45  ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
end do
return
end subroutine dewset

!!$*DECK DVNORM

function dvnorm (n,v,w)
real*8::dvnorm
integer::n,i
real*8::v(n),w(n),sum
sum = 0.0d0
do i = 1,n
10  sum = sum + (v(i)*w(i))**2
end do
dvnorm = sqrt(sum/n)
return
end function dvnorm

!!$*DECK DIPREP

subroutine diprep (neq,y,rwork,ia,ja,ipflag,f,jac)
external f,jac
integer::neq(*),ia(*),ja(*),ipflag
real*8::y(*),rwork(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rlss
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,imax,lewtn,lyhd,lyhn
ipflag = 0
call dprep (neq,y,rwork(lyh),rwork(lsavf),rwork(lewt),&
rwork(lacor),ia,ja,rwork(lwm),rwork(lwm),ipflag,f,jac)
lenwk = max(lreq,lwmin)
if (ipflag .lt. 0) return
lyhn = lwm + lenwk
if (lyhn .gt. lyh) return
lyhd = lyh - lyhn
if (lyhd .eq. 0) go to 20
imax = lyhn - 1 + lenyhm
do i = lyhn,imax
10  rwork(i) = rwork(i+lyhd)
end do
lyh = lyhn
20 lsavf = lyh + lenyh
lewtn = lsavf + n
lacor = lewtn + n
if (istatc .eq. 3) go to 40
if (lewtn .gt. lewt) return
do i = 1,n
30  rwork(i+lewtn-1) = rwork(i+lewt-1)
end do
40 lewt = lewtn
return
end subroutine diprep

!!$*DECK DPREP

subroutine dprep (neq,y,yh,savf,ewt,ftem,ia,ja,&
wk,iwk,ipper,f,jac)
external f,jac
integer::neq(*),ia(*),ja(*),iwk(*),ipper
real*8::y(*),yh(*),savf(*),ewt(*),ftem(*),wk(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::con0,conmin,ccmxj,psmall,rbig,seth
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ con0,conmin,ccmxj,psmall,rbig,seth,&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,ibr,ier,ipil,ipiu,iptt1,iptt2,j,jfound,k,&
knew,kmax,kmin,ldif,lenigp,liwk,maxg,np1,nzsut
real*8::dq,dyj,erwt,fac,yj
ibian = lrat*2
ipian = ibian + 1
np1 = n + 1
ipjan = ipian + np1
ibjan = ipjan - 1
liwk = lenwk*lrat
if (ipjan+n-1 .gt. liwk) go to 210
if (moss .eq. 0) go to 30
if (istatc .eq. 3) go to 20
do i = 1,n
erwt = 1.0d0/ewt(i)
fac = 1.0d0 + 1.0d0/(i + 1.0d0)
y(i) = y(i) + fac*sign(erwt,y(i))
10  continue
end do
go to (70,100),moss
20 continue
do i = 1,n
25  y(i) = yh(i)
end do
go to (70,100),moss
30 knew = ipjan
kmin = ia(1)
iwk(ipian) = 1
do j = 1,n
jfound = 0
kmax = ia(j+1) - 1
if (kmin .gt. kmax) go to 45
do k = kmin,kmax
i = ja(k)
if (i .eq. j) jfound = 1
if (knew .gt. liwk) go to 210
iwk(knew) = i
knew = knew + 1
40   continue
end do
if (jfound .eq. 1) go to 50
45  if (knew .gt. liwk) go to 210
iwk(knew) = j
knew = knew + 1
50  iwk(ipian+j) = knew + 1 - ipjan
kmin = kmax + 1
60  continue
end do
go to 140
70 continue
call f (neq,tn,y,savf)
k = ipjan
iwk(ipian) = 1
do j = 1,n
if (k .gt. liwk) go to 210
iwk(k) = j
k = k + 1
do i = 1,n
75   savf(i) = 0.0d0
end do
call jac (neq,tn,y,j,iwk(ipian),iwk(ipjan),savf)
do i = 1,n
if (abs(savf(i)) .le. seth) go to 80
if (i .eq. j) go to 80
if (k .gt. liwk) go to 210
iwk(k) = i
k = k + 1
80   continue
end do
iwk(ipian+j) = k + 1 - ipjan
90  continue
end do
go to 140
100 k = ipjan
iwk(ipian) = 1
call f (neq,tn,y,savf)
do j = 1,n
if (k .gt. liwk) go to 210
iwk(k) = j
k = k + 1
yj = y(j)
erwt = 1.0d0/ewt(j)
dyj = sign(erwt,yj)
y(j) = yj + dyj
call f (neq,tn,y,ftem)
y(j) = yj
do i = 1,n
dq = (ftem(i) - savf(i))/dyj
if (abs(dq) .le. seth) go to 110
if (i .eq. j) go to 110
if (k .gt. liwk) go to 210
iwk(k) = i
k = k + 1
110   continue
end do
iwk(ipian+j) = k + 1 - ipjan
120 continue
end do
140 continue
if (moss .eq. 0 .or. istatc .ne. 1) go to 150
do i = 1,n
145 y(i) = yh(i)
end do
150 nnz = iwk(ipian+n) - 1
lenigp = 0
ipigp = ipjan + nnz
if (miter .ne. 2) go to 160
maxg = np1
ipjgp = ipjan + nnz
ibjgp = ipjgp - 1
ipigp = ipjgp + n
iptt1 = ipigp + np1
iptt2 = iptt1 + n
lreq = iptt2 + n - 1
if (lreq .gt. liwk) go to 220
call jgroup (n,iwk(ipian),iwk(ipjan),maxg,ngp,iwk(ipigp),&
iwk(ipjgp),iwk(iptt1),iwk(iptt2),ier)
if (ier .ne. 0) go to 220
lenigp = ngp + 1
160 ipr = ipigp + lenigp
ipc = ipr
ipic = ipc + n
ipisp = ipic + n
iprsp = (ipisp - 2)/lrat + 2
iesp = lenwk + 1 - iprsp
if (iesp .lt. 0) go to 230
ibr = ipr - 1
do i = 1,n
170 iwk(ibr+i) = i
end do
nsp = liwk + 1 - ipisp
call odrv (n,iwk(ipian),iwk(ipjan),wk,iwk(ipr),iwk(ipic),&
nsp,iwk(ipisp),1,iys)
if (iys .eq. 11*n+1) go to 240
if (iys .ne. 0) go to 230
ipa = lenwk + 1 - nnz
nsp = ipa - iprsp
lreq = max(12*n/lrat,6*n/lrat+2*n+nnz) + 3
lreq = lreq + iprsp - 1 + nnz
if (lreq .gt. lenwk) go to 250
iba = ipa - 1
do i = 1,nnz
180 wk(iba+i) = 0.0d0
end do
ipisp = lrat*(iprsp - 1) + 1
call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
lreq = lenwk - iesp
if (iys .eq. 10*n+1) go to 250
if (iys .ne. 0) go to 260
ipil = ipisp
ipiu = ipil + 2*n + 1
nzu = iwk(ipil+n) - iwk(ipil)
nzl = iwk(ipiu+n) - iwk(ipiu)
if (lrat .gt. 1) go to 190
call adjlr (n,iwk(ipisp),ldif)
lreq = lreq + ldif
190 continue
if (lrat .eq. 2 .and. nnz .eq. n) lreq = lreq + 1
nsp = nsp + lreq - lenwk
ipa = lreq + 1 - nnz
iba = ipa - 1
ipper = 0
return
210 ipper = -1
lreq = 2 + (2*n + 1)/lrat
lreq = max(lenwk+1,lreq)
return
220 ipper = -2
lreq = (lreq - 1)/lrat + 1
return
230 ipper = -3
call cntnzu (n,iwk(ipian),iwk(ipjan),nzsut)
lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
return
240 ipper = -4
return
250 ipper = -5
return
260 ipper = -6
lreq = lenwk
return
end subroutine dprep

!!$*DECK JGROUP

subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)
integer::n,ia(*),ja(*),maxg,ngrp,igp(*)
integer::jgp(*),incl(*),jdone(*),ier
integer::i,j,k,kmin,kmax,ncol,ng
ier = 0
do j = 1,n
10  jdone(j) = 0
end do
ncol = 1
do ng = 1,maxg
igp(ng) = ncol
do i = 1,n
20   incl(i) = 0
end do
do j = 1,n
if (jdone(j) .eq. 1) go to 50
kmin = ia(j)
kmax = ia(j+1) - 1
do k = kmin,kmax
i = ja(k)
if (incl(i) .eq. 1) go to 50
30     continue
end do
jgp(ncol) = j
ncol = ncol + 1
jdone(j) = 1
do k = kmin,kmax
i = ja(k)
40     incl(i) = 1
end do
50   continue
end do
if (ncol .eq. igp(ng)) go to 70
60  continue
end do
if (ncol .le. n) go to 80
ng = maxg
70 ngrp = ng - 1
return
80 ier = 1
return
end subroutine jgroup

!!$*DECK ADJLR

subroutine adjlr (n,isp,ldif)
integer::n,isp(*),ldif
integer::ip,jlmax,jumax,lnfc,lsfc,nzlu
ip = 2*n + 1
jlmax = isp(ip)
jumax = isp(ip+ip)
nzlu = isp(n+1) - isp(1) + isp(ip+n+1) - isp(ip+1)
lsfc = 12*n + 3 + 2*max(jlmax,jumax)
lnfc = 9*n + 2 + jlmax + jumax + nzlu
ldif = max(0,lsfc - lnfc)
return
end subroutine adjlr

!!$*DECK CNTNZU

subroutine cntnzu (n,ia,ja,nzsut)
integer::n,ia(*),ja(*),nzsut
integer::ii,jj,j,jmin,jmax,k,kmin,kmax,num
num = 0
do ii = 1,n
jmin = ia(ii)
jmax = ia(ii+1) - 1
if (jmin .gt. jmax) go to 50
do j = jmin,jmax
if (ja(j) - ii) 10,40,30
10   jj =ja(j)
kmin = ia(jj)
kmax = ia(jj+1) - 1
if (kmin .gt. kmax) go to 30
do k = kmin,kmax
if (ja(k) .eq. ii) go to 40
20     continue
end do
30   num = num + 1
40   continue
end do
50  continue
end do
nzsut = num
return
end subroutine cntnzu

!!$*DECK DPRJS

subroutine dprjs (neq,y,yh,nyh,ewt,ftem,savf,wk,iwk,f,jac)
external f,jac
integer::neq(*),nyh,iwk(*)
real*8::y(*),yh(nyh,*),ewt(*),ftem(*),savf(*),wk(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::con0,conmin,ccmxj,psmall,rbig,seth
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ con0,conmin,ccmxj,psmall,rbig,seth,&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,imul,j,jj,jok,jmax,jmin,k,kmax,kmin,ng
real*8::con,di,fac,hl0,pij,r,r0,rcon,rcont,&
srur,dvnorm
hl0 = h*el0
con = -hl0
if (miter .eq. 3) go to 300
jok = 1
if (nst .eq. 0 .or. nst .ge. nslj+msbj) jok = 0
if (icf .eq. 1 .and. abs(rc - 1.0d0) .lt. ccmxj) jok = 0
if (icf .eq. 2) jok = 0
if (jok .eq. 1) go to 250
20 jcur = 1
nje = nje + 1
nslj = nst
iplost = 0
conmin = abs(con)
go to (100,200),miter
100 continue
kmin = iwk(ipian)
do j = 1,n
kmax = iwk(ipian+j) - 1
do i = 1,n
110   ftem(i) = 0.0d0
end do
call jac (neq,tn,y,j,iwk(ipian),iwk(ipjan),ftem)
do k = kmin,kmax
i = iwk(ibjan+k)
wk(iba+k) = ftem(i)*con
if (i .eq. j) wk(iba+k) = wk(iba+k) + 1.0d0
120   continue
end do
kmin = kmax + 1
130 continue
end do
go to 290
200 continue
fac = dvnorm(n,savf,ewt)
r0 = 1000.0d0 * abs(h) * uround * n * fac
if (r0 .eq. 0.0d0) r0 = 1.0d0
srur = wk(1)
jmin = iwk(ipigp)
do ng = 1,ngp
jmax = iwk(ipigp+ng) - 1
do j = jmin,jmax
jj = iwk(ibjgp+j)
r = max(srur*abs(y(jj)),r0/ewt(jj))
210   y(jj) = y(jj) + r
end do
call f (neq,tn,y,ftem)
do j = jmin,jmax
jj = iwk(ibjgp+j)
y(jj) = yh(jj,1)
r = max(srur*abs(y(jj)),r0/ewt(jj))
fac = -hl0/r
kmin =iwk(ibian+jj)
kmax =iwk(ibian+jj+1) - 1
do k = kmin,kmax
i = iwk(ibjan+k)
wk(iba+k) = (ftem(i) - savf(i))*fac
if (i .eq. jj) wk(iba+k) = wk(iba+k) + 1.0d0
220    continue
end do
230   continue
end do
jmin = jmax + 1
240 continue
end do
nfe = nfe + ngp
go to 290
250 jcur = 0
rcon = con/con0
rcont = abs(con)/conmin
if (rcont .gt. rbig .and. iplost .eq. 1) go to 20
kmin = iwk(ipian)
do j = 1,n
kmax = iwk(ipian+j) - 1
do k = kmin,kmax
i = iwk(ibjan+k)
pij = wk(iba+k)
if (i .ne. j) go to 260
pij = pij - 1.0d0
if (abs(pij) .ge. psmall) go to 260
iplost = 1
conmin = min(abs(con0),conmin)
260   pij = pij*rcon
if (i .eq. j) pij = pij + 1.0d0
wk(iba+k) = pij
270   continue
end do
kmin = kmax + 1
275 continue
end do
290 nlu = nlu + 1
con0 = con
ierpj = 0
do i = 1,n
295 ftem(i) = 0.0d0
end do
call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
if (iys .eq. 0) return
imul = (iys - 1)/n
ierpj = -2
if (imul .eq. 8) ierpj = 1
if (imul .eq. 10) ierpj = -1
return
300 continue
jcur = 1
nje = nje + 1
wk(2) = hl0
ierpj = 0
r = el0*0.1d0
do i = 1,n
310 y(i) = y(i) + r*(h*savf(i) - yh(i,2))
end do
call f (neq,tn,y,wk(3))
nfe = nfe + 1
do i = 1,n
r0 = h*savf(i) - yh(i,2)
di = 0.1d0*r0 - h*(wk(i+2) - savf(i))
wk(i+2) = 1.0d0
if (abs(r0) .lt. uround/ewt(i)) go to 320
if (abs(di) .eq. 0.0d0) go to 330
wk(i+2) = 0.1d0*r0/di
320 continue
end do
return
330 ierpj = 2
return
end subroutine dprjs

!!$*DECK DSOLSS

subroutine dsolss (wk,iwk,x,tem)
integer::iwk(*)
real*8::wk(*),x(*),tem(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rlss
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i
real*8::di,hl0,phl0,r
iersl = 0
go to (100,100,300),miter
100 call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
if (iersl .ne. 0) iersl = -1
return
300 phl0 = wk(2)
hl0 = h*el0
wk(2) = hl0
if (hl0 .eq. phl0) go to 330
r = hl0/phl0
do i = 1,n
di = 1.0d0 - r*(1.0d0 - 1.0d0/wk(i+2))
if (abs(di) .eq. 0.0d0) go to 390
320 wk(i+2) = 1.0d0/di
end do
330 do i = 1,n
340 x(i) = wk(i+2)*x(i)
end do
return
390 iersl = 1
return
end subroutine dsolss

!!$*DECK DSRCMS

subroutine dsrcms (rsav,isav,job)
integer::isav(*),job
integer::ils,ilss
integer::i,lenils,leniss,lenrls,lenrss
real*8::rsav(*),rls,rlss
save lenrls,lenils,lenrss,leniss
common /dls001/ rls(218),ils(37)
common /dlss01/ rlss(6),ilss(34)
data lenrls/218/,lenils/37/,lenrss/6/,leniss/34/
if (job .eq. 2) go to 100
do i = 1,lenrls
10  rsav(i) = rls(i)
end do
do i = 1,lenrss
15  rsav(lenrls+i) = rlss(i)
end do
do i = 1,lenils
20  isav(i) = ils(i)
end do
do i = 1,leniss
25  isav(lenils+i) = ilss(i)
end do
return
100 continue
do i = 1,lenrls
110 rls(i) = rsav(i)
end do
do i = 1,lenrss
115 rlss(i) = rsav(lenrls+i)
end do
do i = 1,lenils
120 ils(i) = isav(i)
end do
do i = 1,leniss
125 ilss(i) = isav(lenils+i)
end do
return
end subroutine dsrcms

!!$*DECK ODRV

subroutine odrv   (n,ia,ja,a,p,ip,nsp,isp,path,flag)
integer::ia(*),ja(*),p(*),ip(*),isp(*),path,flag, v,l,head,tmp,q
integer::max,nsp,n,next
real*8::a(*)
logical dflag
flag = 0
if (path.lt.1 .or. 5.lt.path) go to 111
if ((path-1) * (path-2) * (path-4) .ne. 0) go to 1
max = (nsp-n)/2
v = 1
l = v + max
head = l + max
next = head + n
if (max.lt.n) go to 110
call md   (n,ia,ja,max,isp(v),isp(l),isp(head),p,ip,isp(v),flag)
if (flag.ne.0) go to 100
1 if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0) go to 2
tmp = (nsp+1) - n
q = tmp - (ia(n+1)-1)
if (q.lt.1) go to 110
dflag = path.eq.4 .or. path.eq.5
call sro   (n,ip,ia,ja,a,isp(tmp),isp(q),dflag)
2 return
100 return
110 flag = 10*n + 1
return
111 flag = 11*n + 1
return
end subroutine odrv

subroutine md (n,ia,ja,max,v,l,head,last,next,mark,flag)
integer::ia(*),ja(*),v(*),l(*),head(*),last(*),next(*), mark(*),flag,tag,dmin,vk,ek,tail
integer::n,max,k
equivalence (vk,ek)
tag = 0
call mdi (n,ia,ja,max,v,l,head,last,next,mark,tag,flag)
if (flag.ne.0) return
k = 0
dmin = 1
1 if (k.ge.n) go to 4
2 if (head(dmin).gt.0) go to 3
dmin = dmin + 1
go to 2
3 vk = head(dmin)
head(dmin) = next(vk)
if (head(dmin).gt.0) last(head(dmin)) = -dmin
k = k+1
next(vk) = -k
last(ek) = dmin - 1
tag = tag + last(ek)
mark(vk) = tag
call mdm   (vk,tail,v,l,last,next,mark)
call mdp   (k,ek,tail,v,l,head,last,next,mark)
call mdu   (ek,dmin,v,l,head,last,next,mark)
go to 1
4 do k=1,n
next(k) = -next(k)
5  last(next(k)) = k
end do
return
end subroutine md

subroutine mdi   (n,ia,ja,max,v,l,head,last,next,mark,tag,flag)
integer::ia(*),ja(*),v(*),l(*),head(*),last(*),next(*), &
     mark(*),tag,flag,sfs,vi,dvi,vj
integer::n,jmin,jmax,j,lwk,kmax,k,max,nextvi,lvk
do vi=1,n
mark(vi) = 1
l(vi) = 0
1  head(vi) = 0
end do
sfs = n+1
do vi=1,n
jmin = ia(vi)
jmax = ia(vi+1) - 1
if (jmin.gt.jmax) go to 6
do j=jmin,jmax
vj = ja(j)
if (vj-vi) 2,5,4
2    lvk = vi
kmax = mark(vi) - 1
if (kmax .eq. 0) go to 4
do k=1,kmax
lvk = l(lvk)
if (v(lvk).eq.vj) go to 5
3     continue
end do
4    if (sfs.ge.max) go to 101
mark(vi) = mark(vi) + 1
v(sfs) = vj
l(sfs) = l(vi)
l(vi) = sfs
sfs = sfs+1
mark(vj) = mark(vj) + 1
v(sfs) = vi
l(sfs) = l(vj)
l(vj) = sfs
sfs = sfs+1
5    continue
end do
6  continue
end do
do vi=1,n
dvi = mark(vi)
next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
nextvi = next(vi)
if (nextvi.gt.0) last(nextvi) = vi
7  mark(vi) = tag
end do
return
101 flag = 9*n + vi
return
end subroutine mdi

subroutine mdm   (vk,tail,v,l,last,next,mark)
integer::vk,tail,v(*),l(*),last(*),next(*),mark(*), tag,s,ls,vs,es,b,lb,vb,blp,blpmax
equivalence (vs,es)
tag = mark(vk)
tail = vk
ls = l(vk)
1 s = ls
if (s.eq.0) go to 5
ls = l(s)
vs = v(s)
if (next(vs).lt.0) go to 2
mark(vs) = tag
l(tail) = s
tail = s
go to 4
2 lb = l(es)
blpmax = last(es)
do blp=1,blpmax
b = lb
lb = l(b)
vb = v(b)
if (mark(vb).ge.tag) go to 3
mark(vb) = tag
l(tail) = b
tail = b
3  continue
end do
mark(es) = tag
4 go to 1
5 l(tail) = 0
return
end subroutine mdm

subroutine mdp   (k,ek,tail,v,l,head,last,next,mark)
integer::ek,tail,v(*),l(*),head(*),last(*),next(*), mark(*),tag,free,li,vi,lvi,evi,s,ls,es,ilp,ilpmax
integer::i,k
tag = mark(ek)
li = ek
ilpmax = last(ek)
if (ilpmax.le.0) go to 12
do ilp=1,ilpmax
i = li
li = l(i)
vi = v(li)
if (last(vi).eq.0) go to 3
if (last(vi).gt.0) go to 1
head(-last(vi)) = next(vi)
go to 2
1  next(last(vi)) = next(vi)
2  if (next(vi).gt.0) last(next(vi)) = last(vi)
3  ls = vi
4  s = ls
ls = l(s)
if (ls.eq.0) go to 6
es = v(ls)
if (mark(es).lt.tag) go to 5
free = ls
l(s) = l(ls)
ls = s
5  go to 4
6  lvi = l(vi)
if (lvi.ne.0) go to 7
l(i) = l(li)
li = i
k = k+1
next(vi) = -k
last(ek) = last(ek) - 1
go to 11
7  if (l(lvi).ne.0) go to 9
evi = v(lvi)
if (next(evi).ge.0) go to 9
if (mark(evi).lt.0) go to 8
last(vi) = evi
mark(evi) = -1
l(tail) = li
tail = li
l(i) = l(li)
li = i
go to 10
8  last(vi) = 0
mark(evi) = mark(evi) - 1
go to 10
9  last(vi) = -ek
10  v(free) = ek
l(free) = l(vi)
l(vi) = free
11  continue
end do
12 l(tail) = 0
return
end subroutine mdp

subroutine mdu   (ek,dmin,v,l,head,last,next,mark)
integer::ek,dmin,v(*),l(*),head(*),last(*),next(*), mark(*),tag,vi,evi,dvi,s,vs,es,b,vb,ilp,ilpmax, blp,blpmax
integer::i
equivalence (vs,es)
tag = mark(ek) - last(ek)
i = ek
ilpmax = last(ek)
if (ilpmax.le.0) go to 11
do ilp=1,ilpmax
i = l(i)
vi = v(i)
if (last(vi)) 1,10,8
1  tag = tag + 1
dvi = last(ek)
s = l(vi)
2  s = l(s)
if (s.eq.0) go to 9
vs = v(s)
if (next(vs).lt.0) go to 3
mark(vs) = tag
dvi = dvi + 1
go to 5
3  if (mark(es).lt.0) go to 6
b = es
blpmax = last(es)
do blp=1,blpmax
b = l(b)
vb = v(b)
if (mark(vb).ge.tag) go to 4
mark(vb) = tag
dvi = dvi + 1
4    continue
end do
5  go to 2
6  last(vi) = 0
mark(es) = mark(es) - 1
7  s = l(s)
if (s.eq.0) go to 10
es = v(s)
if (mark(es).lt.0) mark(es) = mark(es) - 1
go to 7
8  evi = last(vi)
dvi = last(ek) + last(evi) + mark(evi)
mark(evi) = 0
9  next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
if (next(vi).gt.0) last(next(vi)) = vi
if (dvi.lt.dmin) dmin = dvi
10  continue
end do
11 return
end subroutine mdu

subroutine sro(n,ip,ia,ja,a,q,r,dflag)
integer::ip(*),ia(*),ja(*),q(*),r(*)
integer::i,jmin,jmax,j,k,ilast,jdummy,jak
integer::n
real*8::a(*),ak
logical dflag
do i=1,n
1  q(i) = 0
end do
do i=1,n
jmin = ia(i)
jmax = ia(i+1) - 1
if (jmin.gt.jmax) go to 3
do j=jmin,jmax
k = ja(j)
if (ip(k).lt.ip(i)) ja(j) = i
if (ip(k).ge.ip(i)) k = i
r(j) = k
2    q(k) = q(k) + 1
end do
3  continue
end do
do i=1,n
ia(i+1) = ia(i) + q(i)
4  q(i) = ia(i+1)
end do
ilast = 0
jmin = ia(1)
jmax = ia(n+1) - 1
j = jmax
do jdummy=jmin,jmax
i = r(j)
if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast) go to 5
r(j) = ia(i)
ilast = i
go to 6
5  q(i) = q(i) - 1
r(j) = q(i)
6  j = j-1
end do
do j=jmin,jmax
7  if (r(j).eq.j) go to 8
k = r(j)
r(j) = r(k)
r(k) = k
jak = ja(k)
ja(k) = ja(j)
ja(j) = jak
ak = a(k)
a(k) = a(j)
a(j) = ak
go to 7
8  continue
end do
return
end subroutine sro

!!$*DECK CDRV

subroutine cdrv   (n,r,c,ic,ia,ja,a,b,z,nsp,isp,rsp,esp,path,flag)
integer::r(*),c(*),ic(*),ia(*),ja(*),isp(*),esp,path, flag,d,u,q,row,tmp,ar,umax
integer::lratio,il,ijl,n,iu,iju,irl,jrl,jl,max,nsp,jlmax,ira,jra,irac,iru
integer::jru,jutmp,jumax,i,ju,j,l,lmax
real*8::a(*),b(*),z(*),rsp(*)
data lratio/2/
if (path.lt.1 .or. 5.lt.path) go to 111
il = 1
ijl = il + (n+1)
iu = ijl + n
iju = iu + (n+1)
irl = iju + n
jrl = irl + n
jl = jrl + n
if ((path-1) * (path-5) .ne. 0) go to 5
max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
jlmax = max/2
q = jl + jlmax
ira = q + (n+1)
jra = ira + n
irac = jra + n
iru = irac + n
jru = iru + n
jutmp = jru + n
jumax = lratio*nsp + 1 - jutmp
esp = max/lratio
if (jlmax.le.0 .or. jumax.le.0) go to 110
do i=1,n
if (c(i).ne.i) go to 2
1  continue
end do
go to 3
2 ar = nsp + 1 - n
call nroc   (n,ic,ia,ja,a,isp(il),rsp(ar),isp(iu),flag)
if (flag.ne.0) go to 100
3 call nsfc   (n,r,ic,ia,ja,  jlmax,isp(il),isp(jl),isp(ijl),&
       jumax,isp(iu),isp(jutmp),isp(iju),&
       isp(q),isp(ira),isp(jra),isp(irac),&
       isp(irl),isp(jrl),isp(iru),isp(jru),flag)
if(flag .ne. 0) go to 100
jlmax = isp(ijl+n-1)
ju = jl + jlmax
jumax = isp(iju+n-1)
if (jumax.le.0) go to 5
do j=1,jumax
4  isp(ju+j-1) = isp(jutmp+j-1)
end do
5 jlmax = isp(ijl+n-1)
ju = jl + jlmax
jumax = isp(iju+n-1)
l = (ju + jumax - 2 + lratio) / lratio + 1
lmax = isp(il+n) - 1
d = l + lmax
u = d + n
row = nsp + 1 - n
tmp = row - n
umax = tmp - u
esp = umax - (isp(iu+n) - 1)
if ((path-1) * (path-2) .ne. 0) go to 6
if (umax.lt.0) go to 110
call nnfc   (n,r,c,ic,ia,ja,a,z,b,&
     lmax,isp(il),isp(jl),isp(ijl),rsp(l),rsp(d),&
     umax,isp(iu),isp(ju),isp(iju),rsp(u),&
     rsp(row),rsp(tmp),isp(irl),isp(jrl),flag)
if(flag .ne. 0) go to 100
6 if ((path-3) .ne. 0) go to 7
call nnsc   (n,r,c,isp(il),isp(jl),isp(ijl),rsp(l),&
     rsp(d), isp(iu),isp(ju),isp(iju),rsp(u),  z,b,rsp(tmp))
7 if ((path-4) .ne. 0) go to 8
call nntc   (n,r,c,isp(il),isp(jl),isp(ijl),rsp(l),&
     rsp(d), isp(iu),isp(ju),isp(iju),rsp(u),  z,b,rsp(tmp))
8 return
100 return
110 flag = 10*n + 1
return
111 flag = 11*n + 1
return
end subroutine cdrv

subroutine nroc (n,ic,ia,ja,a,jar,ar,p,flag)
integer::ic(*),ia(*),ja(*),jar(*),p(*),flag
integer::k,n,jmin,jmax,j,newj,i
real*8::a(*),ar(*)
do k=1,n
jmin = ia(k)
jmax = ia(k+1) - 1
if(jmin .gt. jmax) go to 5
p(n+1) = n + 1
do j=jmin,jmax
newj = ic(ja(j))
i = n + 1
1    if(p(i) .ge. newj) go to 2
i = p(i)
go to 1
2    if(p(i) .eq. newj) go to 102
p(newj) = p(i)
p(i) = newj
jar(newj) = ja(j)
ar(newj) = a(j)
3    continue
end do
i = n + 1
do j=jmin,jmax
i = p(i)
ja(j) = jar(i)
4    a(j) = ar(i)
end do
5  continue
end do
flag = 0
return
102 flag = n + k
return
end subroutine nroc

subroutine nsfc(n,r,ic,ia,ja,jlmax,il,jl,ijl,&
  jumax,iu,ju,iju,q,ira,jra,irac,irl,jrl,iru,jru,flag)
integer::cend,qm,rend,rk,vj
integer::ia(*),ja(*),ira(*),jra(*),il(*),jl(*),ijl(*)
integer::iu(*),ju(*),iju(*),irl(*),jrl(*),iru(*),jru(*)
integer::r(*),ic(*),q(*),irac(*),flag
integer::np1,n,jlmin,jlptr,jumin,juptr,k
integer::iak,jaiak,luk,m,lastid,lasti,i
integer::jmin,jmax,long,jtmp,j,irll
integer::jlmax,i1,irul,jumax,irai,jairai
np1 = n + 1
jlmin = 1
jlptr = 0
il(1) = 1
jumin = 1
juptr = 0
iu(1) = 1
do k=1,n
irac(k) = 0
jra(k) = 0
jrl(k) = 0
1  jru(k) = 0
end do
do k=1,n
rk = r(k)
iak = ia(rk)
if (iak .ge. ia(rk+1)) go to 101
jaiak = ic(ja(iak))
if (jaiak .gt. k) go to 105
jra(k) = irac(jaiak)
irac(jaiak) = k
2  ira(k) = iak
end do
do k=1,n
q(np1) = np1
luk = -1
vj = irac(k)
if (vj .eq. 0) go to 5
3  qm = np1
4  m = qm
qm = q(m)
if (qm .lt. vj) go to 4
if (qm .eq. vj) go to 102
luk = luk + 1
q(m) = vj
q(vj) = qm
vj = jra(vj)
if (vj .ne. 0) go to 3
5  lastid = 0
lasti = 0
ijl(k) = jlptr
i = k
6  i = jru(i)
if (i .eq. 0) go to 10
qm = np1
jmin = irl(i)
jmax = ijl(i) + il(i+1) - il(i) - 1
long = jmax - jmin
if (long .lt. 0) go to 6
jtmp = jl(jmin)
if (jtmp .ne. k) long = long + 1
if (jtmp .eq. k) r(i) = -r(i)
if (lastid .ge. long) go to 7
lasti = i
lastid = long
7  do j=jmin,jmax
vj = jl(j)
8    m = qm
qm = q(m)
if (qm .lt. vj) go to 8
if (qm .eq. vj) go to 9
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
9    continue
end do
go to 6
10  qm = q(np1)
if (qm .ne. k) go to 105
if (luk .eq. 0) go to 17
if (lastid .ne. luk) go to 11
irll = irl(lasti)
ijl(k) = irll + 1
if (jl(irll) .ne. k) ijl(k) = ijl(k) - 1
go to 17
11  if (jlmin .gt. jlptr) go to 15
qm = q(qm)
do j=jlmin,jlptr
if (jl(j) - qm) 12,13,15
12   continue
end do
go to 15
13  ijl(k) = j
do i=j,jlptr
if (jl(i) .ne. qm) go to 15
qm = q(qm)
if (qm .gt. n) go to 17
14   continue
end do
jlptr = j - 1
15  jlmin = jlptr + 1
ijl(k) = jlmin
if (luk .eq. 0) go to 17
jlptr = jlptr + luk
if (jlptr .gt. jlmax) go to 103
qm = q(np1)
do j=jlmin,jlptr
qm = q(qm)
16   jl(j) = qm
end do
17  irl(k) = ijl(k)
il(k+1) = il(k) + luk
q(np1) = np1
luk = -1
rk = r(k)
jmin = ira(k)
jmax = ia(rk+1) - 1
if (jmin .gt. jmax) go to 20
do j=jmin,jmax
vj = ic(ja(j))
qm = np1
18   m = qm
qm = q(m)
if (qm .lt. vj) go to 18
if (qm .eq. vj) go to 102
luk = luk + 1
q(m) = vj
q(vj) = qm
19   continue
end do
20  lastid = 0
lasti = 0
iju(k) = juptr
i = k
i1 = jrl(k)
21  i = i1
if (i .eq. 0) go to 26
i1 = jrl(i)
qm = np1
jmin = iru(i)
jmax = iju(i) + iu(i+1) - iu(i) - 1
long = jmax - jmin
if (long .lt. 0) go to 21
jtmp = ju(jmin)
if (jtmp .eq. k) go to 22
long = long + 1
cend = ijl(i) + il(i+1) - il(i)
irl(i) = irl(i) + 1
if (irl(i) .ge. cend) go to 22
j = jl(irl(i))
jrl(i) = jrl(j)
jrl(j) = i
22  if (lastid .ge. long) go to 23
lasti = i
lastid = long
23  do j=jmin,jmax
vj = ju(j)
24   m = qm
qm = q(m)
if (qm .lt. vj) go to 24
if (qm .eq. vj) go to 25
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
25   continue
end do
go to 21
26  if (il(k+1) .le. il(k)) go to 27
j = jl(irl(k))
jrl(k) = jrl(j)
jrl(j) = k
27  qm = q(np1)
if (qm .ne. k) go to 105
if (luk .eq. 0) go to 34
if (lastid .ne. luk) go to 28
irul = iru(lasti)
iju(k) = irul + 1
if (ju(irul) .ne. k) iju(k) = iju(k) - 1
go to 34
28  if (jumin .gt. juptr) go to 32
qm = q(qm)
do j=jumin,juptr
if (ju(j) - qm) 29,30,32
29   continue
end do
go to 32
30  iju(k) = j
do i=j,juptr
if (ju(i) .ne. qm) go to 32
qm = q(qm)
if (qm .gt. n) go to 34
31   continue
end do
juptr = j - 1
32  jumin = juptr + 1
iju(k) = jumin
if (luk .eq. 0) go to 34
juptr = juptr + luk
if (juptr .gt. jumax) go to 106
qm = q(np1)
do j=jumin,juptr
qm = q(qm)
33   ju(j) = qm
end do
34  iru(k) = iju(k)
iu(k+1) = iu(k) + luk
i = k
35  i1 = jru(i)
if (r(i) .lt. 0) go to 36
rend = iju(i) + iu(i+1) - iu(i)
if (iru(i) .ge. rend) go to 37
j = ju(iru(i))
jru(i) = jru(j)
jru(j) = i
go to 37
36  r(i) = -r(i)
37  i = i1
if (i .eq. 0) go to 38
iru(i) = iru(i) + 1
go to 35
38  i = irac(k)
if (i .eq. 0) go to 41
39  i1 = jra(i)
ira(i) = ira(i) + 1
if (ira(i) .ge. ia(r(i)+1)) go to 40
irai = ira(i)
jairai = ic(ja(irai))
if (jairai .gt. i) go to 40
jra(i) = irac(jairai)
irac(jairai) = i
40  i = i1
if (i .ne. 0) go to 39
41  continue
end do
ijl(n) = jlptr
iju(n) = juptr
flag = 0
return
101 flag = n + rk
return
102 flag = 2*n + rk
return
103 flag = 3*n + k
return
105 flag = 5*n + k
return
106 flag = 6*n + k
return
end subroutine nsfc

subroutine nnfc(n,r,c,ic,ia,ja,a,z,b,lmax,il,jl,ijl,l,d,umax,iu,ju,iju,u,row,tmp,irl,jrl,flag)
integer::rk,umax
integer::r(*),c(*),ic(*),ia(*),ja(*),il(*),jl(*),ijl(*)
integer::iu(*),ju(*),iju(*),irl(*),jrl(*),flag
integer::n,lmax,k,i1,i,i2
integer::jmin,jmax,j,mu,ijlb
real*8::a(*),l(*),d(*),u(*),z(*),b(*),row(*)
real*8::tmp(*),lki,sum,dk
if(il(n+1)-1 .gt. lmax) go to 104
if(iu(n+1)-1 .gt. umax) go to 107
do k=1,n
irl(k) = il(k)
jrl(k) = 0
1  continue
end do
do k=1,n
row(k) = 0
i1 = 0
if (jrl(k) .eq. 0) go to 3
i = jrl(k)
2  i2 = jrl(i)
jrl(i) = i1
i1 = i
row(i) = 0
i = i2
if (i .ne. 0) go to 2
3  jmin = iju(k)
jmax = jmin + iu(k+1) - iu(k) - 1
if (jmin .gt. jmax) go to 5
do j=jmin,jmax
4    row(ju(j)) = 0
end do
5  rk = r(k)
jmin = ia(rk)
jmax = ia(rk+1) - 1
do j=jmin,jmax
row(ic(ja(j))) = a(j)
6    continue
end do
sum = b(rk)
i = i1
if (i .eq. 0) go to 10
7  lki = -row(i)
l(irl(i)) = -lki
sum = sum + lki * tmp(i)
jmin = iu(i)
jmax = iu(i+1) - 1
if (jmin .gt. jmax) go to 9
mu = iju(i) - jmin
do j=jmin,jmax
8    row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
end do
9  i = jrl(i)
if (i .ne. 0) go to 7
10  if (row(k) .eq. 0.0d0) go to 108
dk = 1.0d0 / row(k)
d(k) = dk
tmp(k) = sum * dk
if (k .eq. n) go to 19
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax) go to 12
mu = iju(k) - jmin
do j=jmin,jmax
11   u(j) = row(ju(mu+j)) * dk
end do
12  continue
i = i1
if (i .eq. 0) go to 18
14  irl(i) = irl(i) + 1
i1 = jrl(i)
if (irl(i) .ge. il(i+1)) go to 17
ijlb = irl(i) - il(i) + ijl(i)
j = jl(ijlb)
15  if (i .gt. jrl(j)) go to 16
j = jrl(j)
go to 15
16  jrl(i) = jrl(j)
jrl(j) = i
17  i = i1
if (i .ne. 0) go to 14
18  if (irl(k) .ge. il(k+1)) go to 19
j = jl(ijl(k))
jrl(k) = jrl(j)
jrl(j) = k
19  continue
end do
k = n
do i=1,n
sum = tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax) go to 21
mu = iju(k) - jmin
do j=jmin,jmax
20   sum = sum - u(j) * tmp(ju(mu+j))
end do
21  tmp(k) = sum
z(c(k)) = sum
22  k = k-1
end do
flag = 0
return
104 flag = 4*n + 1
return
107 flag = 7*n + 1
return
108 flag = 8*n + k
return
end subroutine nnfc

subroutine nnsc(n,r,c,il,jl,ijl,l,d,iu,ju,iju,u,z,b,tmp)
integer::r(*),c(*),il(*),jl(*),ijl(*),iu(*),ju(*),iju(*)
integer::k,n,jmin,jmax,ml,j,i,mu
real*8::l(*),d(*),u(*),b(*),z(*),tmp(*),tmpk,sum
do k=1,n
1  tmp(k) = b(r(k))
end do
do k=1,n
jmin = il(k)
jmax = il(k+1) - 1
tmpk = -d(k) * tmp(k)
tmp(k) = -tmpk
if (jmin .gt. jmax) go to 3
ml = ijl(k) - jmin
do j=jmin,jmax
2    tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
end do
3  continue
end do
k = n
do i=1,n
sum = -tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax) go to 5
mu = iju(k) - jmin
do j=jmin,jmax
4    sum = sum + u(j) * tmp(ju(mu+j))
end do
5  tmp(k) = -sum
z(c(k)) = -sum
k = k - 1
6  continue
end do
return
end subroutine nnsc

subroutine nntc(n,r,c,il,jl,ijl,l,d,iu,ju,iju,u,z,b,tmp)
integer::r(*),c(*),il(*),jl(*),ijl(*),iu(*),ju(*),iju(*)
integer::k,n,jmin,jmax,mu,j,i,ml
real*8::l(*),d(*),u(*),b(*),z(*),tmp(*),tmpk,sum
do k=1,n
1  tmp(k) = b(c(k))
end do
do k=1,n
jmin = iu(k)
jmax = iu(k+1) - 1
tmpk = -tmp(k)
if (jmin .gt. jmax) go to 3
mu = iju(k) - jmin
do j=jmin,jmax
2    tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
end do
3  continue
end do
k = n
do i=1,n
sum = -tmp(k)
jmin = il(k)
jmax = il(k+1) - 1
if (jmin .gt. jmax) go to 5
ml = ijl(k) - jmin
do j=jmin,jmax
4    sum = sum + l(j) * tmp(jl(ml+j))
end do
5  tmp(k) = -sum * d(k)
z(r(k)) = tmp(k)
k = k - 1
6  continue
end do
return
end subroutine nntc

!!$*DECK DSTODA

subroutine dstoda (neq,y,yh,nyh,yh1,ewt,savf,acor,&
wm,iwm,f,jac,pjac,slvs)
external f,jac,pjac,slvs
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),&
acor(*),wm(*)
integer::iownd,ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iownd2,icount,irflag,jtyp,mused,mxordn,mxords
real*8::conit,crate,el,elco,hold,rmax,tesco,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rownd2,cm1,cm2,pdest,pdlast,ratio,&
pdnorm
common /dls001/ conit,crate,el(13),elco(13,12),&
hold,rmax,tesco(3,12),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlsa01/ rownd2,cm1(12),cm2(5),pdest,pdlast,ratio,&
pdnorm,&
iownd2(3),icount,irflag,jtyp,mused,mxordn,mxords
integer::i,i1,iredo,iret,j,jb,m,ncf,newq
integer::lm1,lm1p1,lm2,lm2p1,nqm1,nqm2
real*8::dcon,ddn,del,delp,dsm,dup,exdn,exsm,exup,&
r,rh,rhdn,rhsm,rhup,told,dmnorm
real*8::alpha,dm1,dm2,exm1,exm2,&
pdh,pnorm,rate,rh1,rh1it,rh2,rm,sm1(12)
save sm1
data sm1/0.5d0,0.575d0,0.55d0,0.45d0,0.35d0,0.25d0,&
0.20d0,0.15d0,0.10d0,0.075d0,0.050d0,0.025d0/
kflag = 0
told = tn
ncf = 0
ierpj = 0
iersl = 0
jcur = 0
icf = 0
delp = 0.0d0
if (jstart .gt. 0) go to 200
if (jstart .eq. -1) go to 100
if (jstart .eq. -2) go to 160
lmax = maxord + 1
nq = 1
l = 2
ialth = 2
rmax = 10000.0d0
rc = 0.0d0
el0 = 1.0d0
crate = 0.7d0
hold = h
nslp = 0
ipup = miter
iret = 3
icount = 20
irflag = 0
pdest = 0.0d0
pdlast = 0.0d0
ratio = 5.0d0
call dcfode (2,elco,tesco)
do i = 1,5
10  cm2(i) = tesco(2,i)*elco(i+1,i)
end do
call dcfode (1,elco,tesco)
do i = 1,12
20  cm1(i) = tesco(2,i)*elco(i+1,i)
end do
go to 150
100 ipup = miter
lmax = maxord + 1
if (ialth .eq. 1) ialth = 2
if (meth .eq. mused) go to 160
call dcfode (meth,elco,tesco)
ialth = l
iret = 1
150 do i = 1,l
155 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
go to (160,170,200),iret
160 if (h .eq. hold) go to 200
rh = h/hold
h = hold
iredo = 3
go to 175
170 rh = max(rh,hmin/abs(h))
175 rh = min(rh,rmax)
rh = rh/max(1.0d0,abs(h)*hmxi*rh)
if (meth .eq. 2) go to 178
irflag = 0
pdh = max(abs(h)*pdlast,0.000001d0)
if (rh*pdh*1.00001d0 .lt. sm1(nq)) go to 178
rh = sm1(nq)/pdh
irflag = 1
178 continue
r = 1.0d0
do j = 2,l
r = r*rh
do i = 1,n
180   yh(i,j) = yh(i,j)*r
end do
end do
h = h*rh
rc = rc*rh
ialth = l
if (iredo .eq. 0) go to 690
200 if (abs(rc-1.0d0) .gt. ccmax) ipup = miter
if (nst .ge. nslp+msbp) ipup = miter
tn = tn + h
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
210   yh1(i) = yh1(i) + yh1(i+nyh)
end do
215 continue
end do
pnorm = dmnorm (n,yh1,ewt)
220 m = 0
rate = 0.0d0
del = 0.0d0
do i = 1,n
230 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
if (ipup .le. 0) go to 250
call pjac (neq,y,yh,nyh,ewt,acor,savf,wm,iwm,f,jac)
ipup = 0
rc = 1.0d0
nslp = nst
crate = 0.7d0
if (ierpj .ne. 0) go to 430
250 do i = 1,n
260 acor(i) = 0.0d0
end do
270 if (miter .ne. 0) go to 350
do i = 1,n
savf(i) = h*savf(i) - yh(i,2)
290 y(i) = savf(i) - acor(i)
end do
del = dmnorm (n,y,ewt)
do i = 1,n
y(i) = yh(i,1) + el(1)*savf(i)
300 acor(i) = savf(i)
end do
go to 400
350 do i = 1,n
360 y(i) = h*savf(i) - (yh(i,2) + acor(i))
end do
call slvs (wm,iwm,y,savf)
if (iersl .lt. 0) go to 430
if (iersl .gt. 0) go to 410
del = dmnorm (n,y,ewt)
do i = 1,n
acor(i) = acor(i) + y(i)
380 y(i) = yh(i,1) + el(1)*acor(i)
end do
400 continue
if (del .le. 100.0d0*pnorm*uround) go to 450
if (m .eq. 0 .and. meth .eq. 1) go to 405
if (m .eq. 0) go to 402
rm = 1024.0d0
if (del .le. 1024.0d0*delp) rm = del/delp
rate = max(rate,rm)
crate = max(0.2d0*crate,rm)
402 dcon = del*min(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
if (dcon .gt. 1.0d0) go to 405
pdest = max(pdest,rate/abs(h*el(1)))
if (pdest .ne. 0.0d0) pdlast = pdest
go to 450
405 continue
m = m + 1
if (m .eq. maxcor) go to 410
if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
delp = del
call f (neq,tn,y,savf)
nfe = nfe + 1
go to 270
410 if (miter .eq. 0 .or. jcur .eq. 1) go to 430
icf = 1
ipup = miter
go to 220
430 icf = 2
ncf = ncf + 1
rmax = 2.0d0
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
440   yh1(i) = yh1(i) - yh1(i+nyh)
end do
445 continue
end do
if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
if (abs(h) .le. hmin*1.00001d0) go to 670
if (ncf .eq. mxncf) go to 670
rh = 0.25d0
ipup = miter
iredo = 1
go to 170
450 jcur = 0
if (m .eq. 0) dsm = del/tesco(2,nq)
if (m .gt. 0) dsm = dmnorm (n,acor,ewt)/tesco(2,nq)
if (dsm .gt. 1.0d0) go to 500
kflag = 0
iredo = 0
nst = nst + 1
hu = h
nqu = nq
mused = meth
do j = 1,l
do i = 1,n
460   yh(i,j) = yh(i,j) + el(j)*acor(i)
end do
end do
icount = icount - 1
if (icount .ge. 0) go to 488
if (meth .eq. 2) go to 480
if (nq .gt. 5) go to 488
if (dsm .gt. 100.0d0*pnorm*uround .and. pdest .ne. 0.0d0)&
go to 470
if (irflag .eq. 0) go to 488
rh2 = 2.0d0
nqm2 = min(nq,mxords)
go to 478
470 continue
exsm = 1.0d0/l
rh1 = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rh1it = 2.0d0*rh1
pdh = pdlast*abs(h)
if (pdh*rh1 .gt. 0.00001d0) rh1it = sm1(nq)/pdh
rh1 = min(rh1,rh1it)
if (nq .le. mxords) go to 474
nqm2 = mxords
lm2 = mxords + 1
exm2 = 1.0d0/lm2
lm2p1 = lm2 + 1
dm2 = dmnorm (n,yh(1,lm2p1),ewt)/cm2(mxords)
rh2 = 1.0d0/(1.2d0*dm2**exm2 + 0.0000012d0)
go to 476
474 dm2 = dsm*(cm1(nq)/cm2(nq))
rh2 = 1.0d0/(1.2d0*dm2**exsm + 0.0000012d0)
nqm2 = nq
476 continue
if (rh2 .lt. ratio*rh1) go to 488
478 rh = rh2
icount = 20
meth = 2
miter = jtyp
pdlast = 0.0d0
nq = nqm2
l = nq + 1
go to 170
480 continue
exsm = 1.0d0/l
if (mxordn .ge. nq) go to 484
nqm1 = mxordn
lm1 = mxordn + 1
exm1 = 1.0d0/lm1
lm1p1 = lm1 + 1
dm1 = dmnorm (n,yh(1,lm1p1),ewt)/cm1(mxordn)
rh1 = 1.0d0/(1.2d0*dm1**exm1 + 0.0000012d0)
go to 486
484 dm1 = dsm*(cm2(nq)/cm1(nq))
rh1 = 1.0d0/(1.2d0*dm1**exsm + 0.0000012d0)
nqm1 = nq
exm1 = exsm
486 rh1it = 2.0d0*rh1
pdh = pdnorm*abs(h)
if (pdh*rh1 .gt. 0.00001d0) rh1it = sm1(nqm1)/pdh
rh1 = min(rh1,rh1it)
rh2 = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
if (rh1*ratio .lt. 5.0d0*rh2) go to 488
alpha = max(0.001d0,rh1)
dm1 = (alpha**exm1)*dm1
if (dm1 .le. 1000.0d0*uround*pnorm) go to 488
rh = rh1
icount = 20
meth = 1
miter = 0
pdlast = 0.0d0
nq = nqm1
l = nq + 1
go to 170
488 continue
ialth = ialth - 1
if (ialth .eq. 0) go to 520
if (ialth .gt. 1) go to 700
if (l .eq. lmax) go to 700
do i = 1,n
490 yh(i,lmax) = acor(i)
end do
go to 700
500 kflag = kflag - 1
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
510   yh1(i) = yh1(i) - yh1(i+nyh)
end do
515 continue
end do
rmax = 2.0d0
if (abs(h) .le. hmin*1.00001d0) go to 660
if (kflag .le. -3) go to 640
iredo = 2
rhup = 0.0d0
go to 540
520 rhup = 0.0d0
if (l .eq. lmax) go to 540
do i = 1,n
530 savf(i) = acor(i) - yh(i,lmax)
end do
dup = dmnorm (n,savf,ewt)/tesco(3,nq)
exup = 1.0d0/(l+1)
rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
540 exsm = 1.0d0/l
rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rhdn = 0.0d0
if (nq .eq. 1) go to 550
ddn = dmnorm (n,yh(1,l),ewt)/tesco(1,nq)
exdn = 1.0d0/nq
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
550 if (meth .eq. 2) go to 560
pdh = max(abs(h)*pdlast,0.000001d0)
if (l .lt. lmax) rhup = min(rhup,sm1(l)/pdh)
rhsm = min(rhsm,sm1(nq)/pdh)
if (nq .gt. 1) rhdn = min(rhdn,sm1(nq-1)/pdh)
pdest = 0.0d0
560 if (rhsm .ge. rhup) go to 570
if (rhup .gt. rhdn) go to 590
go to 580
570 if (rhsm .lt. rhdn) go to 580
newq = nq
rh = rhsm
go to 620
580 newq = nq - 1
rh = rhdn
if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
go to 620
590 newq = l
rh = rhup
if (rh .lt. 1.1d0) go to 610
r = el(l)/l
do i = 1,n
600 yh(i,newq+1) = acor(i)*r
end do
go to 630
610 ialth = 3
go to 700
620 if (meth .eq. 2) go to 622
if (rh*pdh*1.00001d0 .ge. sm1(newq)) go to 625
622 if (kflag .eq. 0 .and. rh .lt. 1.1d0) go to 610
625 if (kflag .le. -2) rh = min(rh,0.2d0)
if (newq .eq. nq) go to 170
630 nq = newq
l = nq + 1
iret = 2
go to 150
640 if (kflag .eq. -10) go to 660
rh = 0.1d0
rh = max(hmin/abs(h),rh)
h = h*rh
do i = 1,n
645 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
do i = 1,n
650 yh(i,2) = h*savf(i)
end do
ipup = miter
ialth = 5
if (nq .eq. 1) go to 200
nq = 1
l = 2
iret = 3
go to 150
660 kflag = -1
go to 720
670 kflag = -2
go to 720
680 kflag = -3
go to 720
690 rmax = 10.0d0
700 r = 1.0d0/tesco(2,nqu)
do i = 1,n
710 acor(i) = acor(i)*r
end do
720 hold = h
jstart = 1
return
end subroutine dstoda

!!$*DECK DPRJA

subroutine dprja (neq,y,yh,nyh,ewt,ftem,savf,wm,iwm,&
f,jac)
external f,jac
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),ewt(*),ftem(*),savf(*),&
wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iownd2,iowns2,jtyp,mused,mxordn,mxords
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rownd2,rowns2,pdnorm
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlsa01/ rownd2,rowns2(20),pdnorm,&
iownd2(3),iowns2(2),jtyp,mused,mxordn,mxords
integer::i,i1,i2,ier,ii,j,j1,jj,lenp,&
mba,mband,meb1,meband,ml,ml3,mu,np1
real*8::con,fac,hl0,r,r0,srur,yi,yj,yjj,&
dmnorm,dfnorm,dbnorm
nje = nje + 1
ierpj = 0
jcur = 1
hl0 = h*el0
go to (100,200,300,400,500),miter
100 lenp = n*n
do i = 1,lenp
110 wm(i+2) = 0.0d0
end do
call jac (neq,tn,y,0,0,wm(3),n)
con = -hl0
do i = 1,lenp
120 wm(i+2) = wm(i+2)*con
end do
go to 240
200 fac = dmnorm (n,savf,ewt)
r0 = 1000.0d0*abs(h)*uround*n*fac
if (r0 .eq. 0.0d0) r0 = 1.0d0
srur = wm(1)
j1 = 2
do j = 1,n
yj = y(j)
r = max(srur*abs(yj),r0/ewt(j))
y(j) = y(j) + r
fac = -hl0/r
call f (neq,tn,y,ftem)
do i = 1,n
220   wm(i+j1) = (ftem(i) - savf(i))*fac
end do
y(j) = yj
j1 = j1 + n
230 continue
end do
nfe = nfe + n
240 continue
pdnorm = dfnorm (n,wm(3),ewt)/abs(hl0)
j = 3
np1 = n + 1
do i = 1,n
wm(j) = wm(j) + 1.0d0
250 j = j + np1
end do
call dgefa (wm(3),n,n,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
300 return
400 ml = iwm(1)
mu = iwm(2)
ml3 = ml + 3
mband = ml + mu + 1
meband = mband + ml
lenp = meband*n
do i = 1,lenp
410 wm(i+2) = 0.0d0
end do
call jac (neq,tn,y,ml,mu,wm(ml3),meband)
con = -hl0
do i = 1,lenp
420 wm(i+2) = wm(i+2)*con
end do
go to 570
500 ml = iwm(1)
mu = iwm(2)
mband = ml + mu + 1
mba = min(mband,n)
meband = mband + ml
meb1 = meband - 1
srur = wm(1)
fac = dmnorm (n,savf,ewt)
r0 = 1000.0d0*abs(h)*uround*n*fac
if (r0 .eq. 0.0d0) r0 = 1.0d0
do j = 1,mba
do i = j,n,mband
yi = y(i)
r = max(srur*abs(yi),r0/ewt(i))
530   y(i) = y(i) + r
end do
call f (neq,tn,y,ftem)
do jj = j,n,mband
y(jj) = yh(jj,1)
yjj = y(jj)
r = max(srur*abs(yjj),r0/ewt(jj))
fac = -hl0/r
i1 = max(jj-mu,1)
i2 = min(jj+ml,n)
ii = jj*meb1 - ml + 2
do i = i1,i2
540    wm(ii+i) = (ftem(i) - savf(i))*fac
end do
550   continue
end do
560 continue
end do
nfe = nfe + mba
570 continue
pdnorm = dbnorm (n,wm(ml+3),meband,ml,mu,ewt)/abs(hl0)
ii = mband + 2
do i = 1,n
wm(ii) = wm(ii) + 1.0d0
580 ii = ii + meband
end do
call dgbfa (wm(3),meband,n,ml,mu,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
end subroutine dprja

!!$*DECK DMNORM

function dmnorm (n,v,w)
real*8::dmnorm
integer::n,i
real*8::v(n),w(n),vm
vm = 0.0d0
do i = 1,n
10  vm = max(vm,abs(v(i))*w(i))
end do
dmnorm = vm
return
end function dmnorm

!!$*DECK DFNORM

function dfnorm (n,a,w)
real*8::dfnorm
integer::n,i,j
real*8::a(n,n),w(n),an,sum
an = 0.0d0
do i = 1,n
sum = 0.0d0
do j = 1,n
10   sum = sum + abs(a(i,j))/w(j)
end do
an = max(an,sum*w(i))
20  continue
end do
dfnorm = an
return
end function dfnorm

!!$*DECK DBNORM

function dbnorm (n,a,nra,ml,mu,w)
real*8::dbnorm
integer::n,nra,ml,mu
integer::i,i1,jlo,jhi,j
real*8::an,w(n)
real*8::a(nra,n),sum
an = 0.0d0
do i = 1,n
sum = 0.0d0
i1 = i + mu + 1
jlo = max(i-ml,1)
jhi = min(i+mu,n)
do j = jlo,jhi
10   sum = sum + abs(a(i1-j,j))/w(j)
end do
an = max(an,sum*w(i))
20  continue
end do
dbnorm = an
return
end function dbnorm

!!$*DECK DSRCMA

subroutine dsrcma (rsav,isav,job)
integer::isav(*),job
integer::ils,ilsa
integer::i,lenrls,lenils,lenrla,lenila
real*8::rsav(*)
real*8::rls,rlsa
save lenrls,lenils,lenrla,lenila
common /dls001/ rls(218),ils(37)
common /dlsa01/ rlsa(22),ilsa(9)
data lenrls/218/,lenils/37/,lenrla/22/,lenila/9/
if (job .eq. 2) go to 100
do i = 1,lenrls
10  rsav(i) = rls(i)
end do
do i = 1,lenrla
15  rsav(lenrls+i) = rlsa(i)
end do
do i = 1,lenils
20  isav(i) = ils(i)
end do
do i = 1,lenila
25  isav(lenils+i) = ilsa(i)
end do
return
100 continue
do i = 1,lenrls
110 rls(i) = rsav(i)
end do
do i = 1,lenrla
115 rlsa(i) = rsav(lenrls+i)
end do
do i = 1,lenils
120 ils(i) = isav(i)
end do
do i = 1,lenila
125 ilsa(i) = isav(lenils+i)
end do
return
end subroutine dsrcma

!!$*DECK DRCHEK

subroutine drchek (job,g,neq,y,yh,nyh,g0,g1,gx,jroot,irt)
external g
integer::job,neq(*),nyh,jroot(*),irt
real*8::y(*),yh(nyh,*),g0(*),g1(*),gx(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iownd3,iownr3,irfnd,itaskc,ngc,nge
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rownr3,t0,tlast,toutc
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlsr01/ rownr3(2),t0,tlast,toutc,&
iownd3(3),iownr3(2),irfnd,itaskc,ngc,nge
integer::i,iflag,jflag
real*8::hming,t1,temp1,temp2,x
logical zroot
irt = 0
do i = 1,ngc
10  jroot(i) = 0
end do
hming = (abs(tn) + abs(h))*uround*100.0d0
go to (100,200,300),job
100 continue
t0 = tn
call g (neq,t0,y,ngc,g0)
nge = 1
zroot = .false.
do i = 1,ngc
110 if (abs(g0(i)) .le. 0.0d0) zroot = .true.
end do
if (.not. zroot) go to 190
temp2 = max(hming/abs(h),0.1d0)
temp1 = temp2*h
t0 = t0 + temp1
do i = 1,n
120 y(i) = y(i) + temp2*yh(i,2)
end do
call g (neq,t0,y,ngc,g0)
nge = nge + 1
zroot = .false.
do i = 1,ngc
130 if (abs(g0(i)) .le. 0.0d0) zroot = .true.
end do
if (.not. zroot) go to 190
irt = -1
return
190 continue
return
200 continue
if (irfnd .eq. 0) go to 260
call dintdy (t0,0,yh,nyh,y,iflag)
call g (neq,t0,y,ngc,g0)
nge = nge + 1
zroot = .false.
do i = 1,ngc
210 if (abs(g0(i)) .le. 0.0d0) zroot = .true.
end do
if (.not. zroot) go to 260
temp1 = sign(hming,h)
t0 = t0 + temp1
if ((t0 - tn)*h .lt. 0.0d0) go to 230
temp2 = temp1/h
do i = 1,n
220 y(i) = y(i) + temp2*yh(i,2)
end do
go to 240
230 call dintdy (t0,0,yh,nyh,y,iflag)
240 call g (neq,t0,y,ngc,g0)
nge = nge + 1
zroot = .false.
do i = 1,ngc
if (abs(g0(i)) .gt. 0.0d0) go to 250
jroot(i) = 1
zroot = .true.
250 continue
end do
if (.not. zroot) go to 260
irt = 1
return
260 if (tn .eq. tlast) go to 390
300 continue
if (itaskc.eq.2 .or. itaskc.eq.3 .or. itaskc.eq.5) go to 310
if ((toutc - tn)*h .ge. 0.0d0) go to 310
t1 = toutc
if ((t1 - t0)*h .le. 0.0d0) go to 390
call dintdy (t1,0,yh,nyh,y,iflag)
go to 330
310 t1 = tn
do i = 1,n
320 y(i) = yh(i,1)
end do
330 call g (neq,t1,y,ngc,g1)
nge = nge + 1
jflag = 0
350 continue
call droots (ngc,hming,jflag,t0,t1,g0,g1,gx,x,jroot)
if (jflag .gt. 1) go to 360
call dintdy (x,0,yh,nyh,y,iflag)
call g (neq,x,y,ngc,gx)
nge = nge + 1
go to 350
360 t0 = x
call dcopy (ngc,gx,1,g0,1)
if (jflag .eq. 4) go to 390
call dintdy (x,0,yh,nyh,y,iflag)
irt = 1
return
390 continue
return
end subroutine drchek

!!$*DECK DROOTS

subroutine droots (ng,hmin,jflag,x0,x1,g0,g1,gx,x,jroot)
integer::ng,jflag,jroot(ng)
real*8::hmin,x0,x1,g0(ng),g1(ng),gx(ng),x
integer::iownd3,imax,last,idum3
real*8::alpha,x2,rdum3
common /dlsr01/ alpha,x2,rdum3(3),&
iownd3(3),imax,last,idum3(4)
integer::i,imxold,nxlast
real*8::t2,tmax,fracint,fracsub,zero,half,tenth,five
logical zroot,sgnchg,xroot
save zero,half,tenth,five
data zero/0.0d0/,half/0.5d0/,tenth/0.1d0/,five/5.0d0/
if (jflag .eq. 1) go to 200
imax = 0
tmax = zero
zroot = .false.
do i = 1,ng
if (abs(g1(i)) .gt. zero) go to 110
zroot = .true.
go to 120
110 if (sign(1.0d0,g0(i)) .eq. sign(1.0d0,g1(i))) go to 120
t2 = abs(g1(i)/(g1(i)-g0(i)))
if (t2 .le. tmax) go to 120
tmax = t2
imax = i
120 continue
end do
if (imax .gt. 0) go to 130
sgnchg = .false.
go to 140
130 sgnchg = .true.
140 if (.not. sgnchg) go to 400
xroot = .false.
nxlast = 0
last = 1
150 continue
if (xroot) go to 300
if (nxlast .eq. last) go to 160
alpha = 1.0d0
go to 180
160 if (last .eq. 0) go to 170
alpha = 0.5d0*alpha
go to 180
170 alpha = 2.0d0*alpha
180 x2 = x1 - (x1 - x0)*g1(imax) / (g1(imax) - alpha*g0(imax))
if (abs(x2 - x0) < half*hmin) then
fracint = abs(x1 - x0)/hmin
fracsub = tenth
if (fracint .le. five) fracsub = half/fracint
x2 = x0 + fracsub*(x1 - x0)
endif
if (abs(x1 - x2) < half*hmin) then
fracint = abs(x1 - x0)/hmin
fracsub = tenth
if (fracint .le. five) fracsub = half/fracint
x2 = x1 - fracsub*(x1 - x0)
endif
jflag = 1
x = x2
return
200 imxold = imax
imax = 0
tmax = zero
zroot = .false.
do i = 1,ng
if (abs(gx(i)) .gt. zero) go to 210
zroot = .true.
go to 220
210 if (sign(1.0d0,g0(i)) .eq. sign(1.0d0,gx(i))) go to 220
t2 = abs(gx(i)/(gx(i) - g0(i)))
if (t2 .le. tmax) go to 220
tmax = t2
imax = i
220 continue
end do
if (imax .gt. 0) go to 230
sgnchg = .false.
imax = imxold
go to 240
230 sgnchg = .true.
240 nxlast = last
if (.not. sgnchg) go to 250
x1 = x2
call dcopy (ng,gx,1,g1,1)
last = 1
xroot = .false.
go to 270
250 if (.not. zroot) go to 260
x1 = x2
call dcopy (ng,gx,1,g1,1)
xroot = .true.
go to 270
260 continue
call dcopy (ng,gx,1,g0,1)
x0 = x2
last = 0
xroot = .false.
270 if (abs(x1-x0) .le. hmin) xroot = .true.
go to 150
300 jflag = 2
x = x1
call dcopy (ng,g1,1,gx,1)
do i = 1,ng
jroot(i) = 0
if (abs(g1(i)) .gt. zero) go to 310
jroot(i) = 1
go to 320
310 if (sign(1.0d0,g0(i)) .ne. sign(1.0d0,g1(i))) jroot(i) = 1
320 continue
end do
return
400 if (.not. zroot) go to 420
x = x1
call dcopy (ng,g1,1,gx,1)
do i = 1,ng
jroot(i) = 0
if (abs(g1(i)) .le. zero) jroot (i) = 1
410 continue
end do
jflag = 3
return
420 call dcopy (ng,g1,1,gx,1)
x = x1
jflag = 4
return
end subroutine droots

!!$*DECK DSRCAR

subroutine dsrcar (rsav,isav,job)
integer::isav(*),job
integer::ils,ilsa,ilsr
integer::i,ioff,lenrls,lenils,lenrla,lenila,lenrlr,lenilr
real*8::rsav(*)
real*8::rls,rlsa,rlsr
save lenrls,lenils,lenrla,lenila,lenrlr,lenilr
common /dls001/ rls(218),ils(37)
common /dlsa01/ rlsa(22),ilsa(9)
common /dlsr01/ rlsr(5),ilsr(9)
data lenrls/218/,lenils/37/,lenrla/22/,lenila/9/
data lenrlr/5/,lenilr/9/
if (job .eq. 2) go to 100
do i = 1,lenrls
10  rsav(i) = rls(i)
end do
do i = 1,lenrla
15  rsav(lenrls+i) = rlsa(i)
end do
ioff = lenrls + lenrla
do i = 1,lenrlr
20  rsav(ioff+i) = rlsr(i)
end do
do i = 1,lenils
30  isav(i) = ils(i)
end do
do i = 1,lenila
35  isav(lenils+i) = ilsa(i)
end do
ioff = lenils + lenila
do i = 1,lenilr
40  isav(ioff+i) = ilsr(i)
end do
return
100 continue
do i = 1,lenrls
110 rls(i) = rsav(i)
end do
do i = 1,lenrla
115 rlsa(i) = rsav(lenrls+i)
end do
ioff = lenrls + lenrla
do i = 1,lenrlr
120 rlsr(i) = rsav(ioff+i)
end do
do i = 1,lenils
130 ils(i) = isav(i)
end do
do i = 1,lenila
135 ilsa(i) = isav(lenils+i)
end do
ioff = lenils + lenila
do i = 1,lenilr
140 ilsr(i) = isav(ioff+i)
end do
return
end subroutine dsrcar

!!$*DECK DSTODPK

subroutine dstodpk (neq,y,yh,nyh,yh1,ewt,savf,savx,acor,&
wm,iwm,f,jac,psol)
external f,jac,psol
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),&
savx(*),acor(*),wm(*)
integer::iownd,ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
real*8::conit,crate,el,elco,hold,rmax,tesco,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::delt,epcon,sqrtn,rsqrtn
common /dls001/ conit,crate,el(13),elco(13,12),&
hold,rmax,tesco(3,12),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlpk01/ delt,epcon,sqrtn,rsqrtn,&
jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
integer::i,i1,iredo,iret,j,jb,m,ncf,newq
real*8::dcon,ddn,del,delp,dsm,dup,exdn,exsm,exup,&
r,rh,rhdn,rhsm,rhup,told,dvnorm
kflag = 0
told = tn
ncf = 0
ierpj = 0
iersl = 0
jcur = 0
icf = 0
delp = 0.0d0
if (jstart .gt. 0) go to 200
if (jstart .eq. -1) go to 100
if (jstart .eq. -2) go to 160
lmax = maxord + 1
nq = 1
l = 2
ialth = 2
rmax = 10000.0d0
rc = 0.0d0
el0 = 1.0d0
crate = 0.7d0
hold = h
meo = meth
nslp = 0
ipup = miter
iret = 3
go to 140
100 ipup = miter
lmax = maxord + 1
if (ialth .eq. 1) ialth = 2
if (meth .eq. meo) go to 110
call dcfode (meth,elco,tesco)
meo = meth
if (nq .gt. maxord) go to 120
ialth = l
iret = 1
go to 150
110 if (nq .le. maxord) go to 160
120 nq = maxord
l = lmax
do i = 1,l
125 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
epcon = conit*tesco(2,nq)
ddn = dvnorm (n,savf,ewt)/tesco(1,l)
exdn = 1.0d0/l
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
rh = min(rhdn,1.0d0)
iredo = 3
if (h .eq. hold) go to 170
rh = min(rh,abs(h/hold))
h = hold
go to 175
140 call dcfode (meth,elco,tesco)
150 do i = 1,l
155 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
epcon = conit*tesco(2,nq)
go to (160,170,200),iret
160 if (h .eq. hold) go to 200
rh = h/hold
h = hold
iredo = 3
go to 175
170 rh = max(rh,hmin/abs(h))
175 rh = min(rh,rmax)
rh = rh/max(1.0d0,abs(h)*hmxi*rh)
r = 1.0d0
do j = 2,l
r = r*rh
do i = 1,n
180   yh(i,j) = yh(i,j)*r
end do
end do
h = h*rh
rc = rc*rh
ialth = l
if (iredo .eq. 0) go to 690
200 if (jacflg .ne. 0) go to 202
ipup = 0
crate = 0.7d0
go to 205
202 if (abs(rc-1.0d0) .gt. ccmax) ipup = miter
if (nst .ge. nslp+msbp) ipup = miter
205 tn = tn + h
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
210   yh1(i) = yh1(i) + yh1(i+nyh)
end do
215 continue
end do
220 m = 0
mnewt = 0
do i = 1,n
230 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
if (ipup .le. 0) go to 250
call dpkset (neq,y,yh1,ewt,acor,savf,wm,iwm,f,jac)
ipup = 0
rc = 1.0d0
nslp = nst
crate = 0.7d0
if (ierpj .ne. 0) go to 430
250 do i = 1,n
260 acor(i) = 0.0d0
end do
270 if (miter .ne. 0) go to 350
do i = 1,n
savf(i) = h*savf(i) - yh(i,2)
290 y(i) = savf(i) - acor(i)
end do
del = dvnorm (n,y,ewt)
do i = 1,n
y(i) = yh(i,1) + el(1)*savf(i)
300 acor(i) = savf(i)
end do
go to 400
350 do i = 1,n
360 savx(i) = h*savf(i) - (yh(i,2) + acor(i))
end do
call dsolpk (neq,y,savf,savx,ewt,wm,iwm,f,psol)
if (iersl .lt. 0) go to 430
if (iersl .gt. 0) go to 410
del = dvnorm (n,savx,ewt)
do i = 1,n
acor(i) = acor(i) + savx(i)
380 y(i) = yh(i,1) + el(1)*acor(i)
end do
400 if (m .ne. 0) crate = max(0.2d0*crate,del/delp)
dcon = del*min(1.0d0,1.5d0*crate)/epcon
if (dcon .le. 1.0d0) go to 450
m = m + 1
if (m .eq. maxcor) go to 410
if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
mnewt = m
delp = del
call f (neq,tn,y,savf)
nfe = nfe + 1
go to 270
410 if (miter.eq.0 .or. jcur.eq.1 .or. jacflg.eq.0) go to 430
icf = 1
ipup = miter
go to 220
430 icf = 2
ncf = ncf + 1
ncfn = ncfn + 1
rmax = 2.0d0
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
440   yh1(i) = yh1(i) - yh1(i+nyh)
end do
445 continue
end do
if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
if (abs(h) .le. hmin*1.00001d0) go to 670
if (ncf .eq. mxncf) go to 670
rh = 0.5d0
ipup = miter
iredo = 1
go to 170
450 jcur = 0
if (m .eq. 0) dsm = del/tesco(2,nq)
if (m .gt. 0) dsm = dvnorm (n,acor,ewt)/tesco(2,nq)
if (dsm .gt. 1.0d0) go to 500
kflag = 0
iredo = 0
nst = nst + 1
hu = h
nqu = nq
do j = 1,l
do i = 1,n
470   yh(i,j) = yh(i,j) + el(j)*acor(i)
end do
end do
ialth = ialth - 1
if (ialth .eq. 0) go to 520
if (ialth .gt. 1) go to 700
if (l .eq. lmax) go to 700
do i = 1,n
490 yh(i,lmax) = acor(i)
end do
go to 700
500 kflag = kflag - 1
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
510   yh1(i) = yh1(i) - yh1(i+nyh)
end do
515 continue
end do
rmax = 2.0d0
if (abs(h) .le. hmin*1.00001d0) go to 660
if (kflag .le. -3) go to 640
iredo = 2
rhup = 0.0d0
go to 540
520 rhup = 0.0d0
if (l .eq. lmax) go to 540
do i = 1,n
530 savf(i) = acor(i) - yh(i,lmax)
end do
dup = dvnorm (n,savf,ewt)/tesco(3,nq)
exup = 1.0d0/(l+1)
rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
540 exsm = 1.0d0/l
rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rhdn = 0.0d0
if (nq .eq. 1) go to 560
ddn = dvnorm (n,yh(1,l),ewt)/tesco(1,nq)
exdn = 1.0d0/nq
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
560 if (rhsm .ge. rhup) go to 570
if (rhup .gt. rhdn) go to 590
go to 580
570 if (rhsm .lt. rhdn) go to 580
newq = nq
rh = rhsm
go to 620
580 newq = nq - 1
rh = rhdn
if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
go to 620
590 newq = l
rh = rhup
if (rh .lt. 1.1d0) go to 610
r = el(l)/l
do i = 1,n
600 yh(i,newq+1) = acor(i)*r
end do
go to 630
610 ialth = 3
go to 700
620 if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
if (kflag .le. -2) rh = min(rh,0.2d0)
if (newq .eq. nq) go to 170
630 nq = newq
l = nq + 1
iret = 2
go to 150
640 if (kflag .eq. -10) go to 660
rh = 0.1d0
rh = max(hmin/abs(h),rh)
h = h*rh
do i = 1,n
645 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
do i = 1,n
650 yh(i,2) = h*savf(i)
end do
ipup = miter
ialth = 5
if (nq .eq. 1) go to 200
nq = 1
l = 2
iret = 3
go to 150
660 kflag = -1
go to 720
670 kflag = -2
go to 720
680 kflag = -3
go to 720
690 rmax = 10.0d0
700 r = 1.0d0/tesco(2,nqu)
do i = 1,n
710 acor(i) = acor(i)*r
end do
720 hold = h
jstart = 1
return
end subroutine dstodpk

!!$*DECK DPKSET

subroutine dpkset (neq,y,ysv,ewt,ftem,savf,wm,iwm,f,jac)
external f,jac
integer::neq(*),iwm(*)
real*8::y(*),ysv(*),ewt(*),ftem(*),savf(*),&
wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::delt,epcon,sqrtn,rsqrtn
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlpk01/ delt,epcon,sqrtn,rsqrtn,&
jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
integer::ier
real*8::hl0
ierpj = 0
jcur = 1
hl0 = el0*h
call jac (f,neq,tn,y,ysv,ewt,savf,ftem,hl0,&
wm(locwp),iwm(lociwp),ier)
nje = nje + 1
if (ier .eq. 0) return
ierpj = 1
return
end subroutine dpkset

!!$*DECK DSOLPK

subroutine dsolpk (neq,y,savf,x,ewt,wm,iwm,f,psol)
external f,psol
integer::neq(*),iwm(*)
real*8::y(*),savf(*),x(*),ewt(*),wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::delt,epcon,sqrtn,rsqrtn
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlpk01/ delt,epcon,sqrtn,rsqrtn,&
jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
integer::iflag,lb,ldl,lhes,liom,lgmr,lpcg,lp,lq,lr,&
lv,lw,lwk,lz,maxlp1,npsl
real*8::delta,hl0
iersl = 0
hl0 = h*el0
delta = delt*epcon
go to (100,200,300,400,900,900,900,900,900),miter
100 continue
lv = 1
lb = lv + n*maxl
lhes = lb + n
lwk = lhes + maxl*maxl
call dcopy (n,x,1,wm(lb),1)
call dscal (n,rsqrtn,ewt,1)
call dspiom (neq,tn,y,savf,wm(lb),ewt,n,maxl,kmp,delta,&
hl0,jpre,mnewt,f,psol,npsl,x,wm(lv),wm(lhes),iwm,&
liom,wm(locwp),iwm(lociwp),wm(lwk),iflag)
nni = nni + 1
nli = nli + liom
nps = nps + npsl
call dscal (n,sqrtn,ewt,1)
if (iflag .ne. 0) ncfl = ncfl + 1
if (iflag .ge. 2) iersl = 1
if (iflag .lt. 0) iersl = -1
return
200 continue
maxlp1 = maxl + 1
lv = 1
lb = lv + n*maxl
lhes = lb + n + 1
lq = lhes + maxl*maxlp1
lwk = lq + 2*maxl
ldl = lwk + min(1,maxl-kmp)*n
call dcopy (n,x,1,wm(lb),1)
call dscal (n,rsqrtn,ewt,1)
call dspigmr (neq,tn,y,savf,wm(lb),ewt,n,maxl,maxlp1,kmp,&
delta,hl0,jpre,mnewt,f,psol,npsl,x,wm(lv),wm(lhes),&
wm(lq),lgmr,wm(locwp),iwm(lociwp),wm(lwk),wm(ldl),iflag)
nni = nni + 1
nli = nli + lgmr
nps = nps + npsl
call dscal (n,sqrtn,ewt,1)
if (iflag .ne. 0) ncfl = ncfl + 1
if (iflag .ge. 2) iersl = 1
if (iflag .lt. 0) iersl = -1
return
300 continue
lr = 1
lp = lr + n
lw = lp + n
lz = lw + n
lwk = lz + n
call dcopy (n,x,1,wm(lr),1)
call dpcg (neq,tn,y,savf,wm(lr),ewt,n,maxl,delta,hl0,&
jpre,mnewt,f,psol,npsl,x,wm(lp),wm(lw),wm(lz),&
lpcg,wm(locwp),iwm(lociwp),wm(lwk),iflag)
nni = nni + 1
nli = nli + lpcg
nps = nps + npsl
if (iflag .ne. 0) ncfl = ncfl + 1
if (iflag .ge. 2) iersl = 1
if (iflag .lt. 0) iersl = -1
return
400 continue
lr = 1
lp = lr + n
lw = lp + n
lz = lw + n
lwk = lz + n
call dcopy (n,x,1,wm(lr),1)
call dpcgs (neq,tn,y,savf,wm(lr),ewt,n,maxl,delta,hl0,&
jpre,mnewt,f,psol,npsl,x,wm(lp),wm(lw),wm(lz),&
lpcg,wm(locwp),iwm(lociwp),wm(lwk),iflag)
nni = nni + 1
nli = nli + lpcg
nps = nps + npsl
if (iflag .ne. 0) ncfl = ncfl + 1
if (iflag .ge. 2) iersl = 1
if (iflag .lt. 0) iersl = -1
return
900 continue
lb = 1
lwk = lb + n
call dcopy (n,x,1,wm(lb),1)
call dusol (neq,tn,y,savf,wm(lb),ewt,n,delta,hl0,mnewt,&
psol,npsl,x,wm(locwp),iwm(lociwp),wm(lwk),iflag)
nni = nni + 1
nps = nps + npsl
if (iflag .ne. 0) ncfl = ncfl + 1
if (iflag .eq. 3) iersl = 1
if (iflag .lt. 0) iersl = -1
return
end subroutine dsolpk

!!$*DECK DSPIOM

subroutine dspiom (neq,tn,y,savf,b,wght,n,maxl,kmp,delta,&
hl0,jpre,mnewt,f,psol,npsl,x,v,hes,ipvt,&
liom,wp,iwp,wk,iflag)
external f,psol
integer::neq(*),n,maxl,kmp,jpre,mnewt,npsl,ipvt(*),liom,iwp(*),iflag
real*8::tn,y(*),savf(*),b(*),wght(*),delta,hl0,x(*),v(n,*),hes(maxl,maxl),wp(*),wk(*)
integer::i,ier,info,j,k,ll,lm1
real*8::bnrm,bnrm0,prod,rho,snormw,dnrm2,tem
iflag = 0
liom = 0
npsl = 0
do i = 1,n
10  v(i,1) = b(i)*wght(i)
end do
bnrm0 = dnrm2 (n,v,1)
bnrm = bnrm0
if (bnrm0 .gt. delta) go to 30
if (mnewt .gt. 0) go to 20
call dcopy (n,b,1,x,1)
return
20 do i = 1,n
25  x(i) = 0.0d0
end do
return
30 continue
ier = 0
if (jpre .eq. 0 .or. jpre .eq. 2) go to 55
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,b,1,ier)
npsl = 1
if (ier .ne. 0) go to 300
do i = 1,n
50  v(i,1) = b(i)*wght(i)
end do
bnrm = dnrm2(n,v,1)
delta = delta*(bnrm/bnrm0)
55 tem = 1.0d0/bnrm
call dscal (n,tem,v(1,1),1)
do j = 1,maxl
do i = 1,maxl
60   hes(i,j) = 0.0d0
end do
65  continue
end do
prod = 1.0d0
do ll = 1,maxl
liom = ll
call datv (neq,y,savf,v(1,ll),wght,x,f,psol,v(1,ll+1),&
wk,wp,iwp,hl0,jpre,ier,npsl)
if (ier .ne. 0) go to 300
call dorthog (v(1,ll+1),v,hes,n,ll,maxl,kmp,snormw)
call dhefa (hes,maxl,ll,ipvt,info,ll)
lm1 = ll - 1
if (ll .gt. 1 .and. ipvt(lm1) .eq. lm1) prod = prod*hes(ll,lm1)
if (info .ne. ll) go to 70
if (snormw .eq. 0.0d0) go to 120
if (ll .eq. maxl) go to 120
go to 80
70  continue
rho = bnrm*snormw*abs(prod/hes(ll,ll))
if (rho .le. delta) go to 200
if (ll .eq. maxl) go to 100
80  continue
hes(ll+1,ll) = snormw
tem = 1.0d0/snormw
call dscal (n,tem,v(1,ll+1),1)
90  continue
end do
100 continue
if (rho .le. 1.0d0) go to 150
if (rho .le. bnrm .and. mnewt .eq. 0) go to 150
120 continue
iflag = 2
return
150 iflag = 1
200 continue
ll = liom
do k = 1,ll
210 b(k) = 0.0d0
end do
b(1) = bnrm
call dhesl (hes,maxl,ll,ipvt,b)
do k = 1,n
220 x(k) = 0.0d0
end do
do i = 1,ll
call daxpy (n,b(i),v(1,i),1,x,1)
230 continue
end do
do i = 1,n
240 x(i) = x(i)/wght(i)
end do
if (jpre .le. 1) return
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,x,2,ier)
npsl = npsl + 1
if (ier .ne. 0) go to 300
return
300 continue
if (ier .lt. 0) iflag = -1
if (ier .gt. 0) iflag = 3
return
end subroutine dspiom

!!$*DECK DATV

subroutine datv (neq,y,savf,v,wght,ftem,f,psol,z,vtem,&
wp,iwp,hl0,jpre,ier,npsl)
external f,psol
integer::neq(*),iwp(*),jpre,ier,npsl
real*8::hl0
real*8::y(*),savf(*),v(*),wght(*),ftem(*),z(*),&
vtem(*),wp(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i
real*8::fac,rnorm,dnrm2,tempn
do i = 1,n
10  vtem(i) = v(i)/wght(i)
end do
ier = 0
if (jpre .ge. 2) go to 30
call dcopy (n,y,1,z,1)
do i = 1,n
20  y(i) = z(i) + vtem(i)
end do
fac = hl0
go to 60
30 continue
call psol (neq,tn,y,savf,ftem,hl0,wp,iwp,vtem,2,ier)
npsl = npsl + 1
if (ier .ne. 0) return
do i = 1,n
40  z(i) = vtem(i)*wght(i)
end do
tempn = dnrm2 (n,z,1)
rnorm = 1.0d0/tempn
call dcopy (n,y,1,z,1)
do i = 1,n
50  y(i) = z(i) + vtem(i)*rnorm
end do
fac = hl0*tempn
60 continue
call f (neq,tn,y,ftem)
nfe = nfe + 1
call dcopy (n,z,1,y,1)
do i = 1,n
70  z(i) = ftem(i) - savf(i)
end do
do i = 1,n
80  z(i) = vtem(i) - fac*z(i)
end do
if (jpre .eq. 0 .or. jpre .eq. 2) go to 85
call psol (neq,tn,y,savf,ftem,hl0,wp,iwp,z,1,ier)
npsl = npsl + 1
if (ier .ne. 0) return
85 continue
do i = 1,n
90  z(i) = z(i)*wght(i)
end do
return
end subroutine datv

!!$*DECK DORTHOG

subroutine dorthog (vnew,v,hes,n,ll,ldhes,kmp,snormw)
integer::n,ll,ldhes,kmp
real*8::vnew(*),v(n,*),hes(ldhes,*),snormw
integer::i,i0
real*8::arg,ddot,dnrm2,sumdsq,tem,vnrm
vnrm = dnrm2 (n,vnew,1)
i0 = max(1,ll-kmp+1)
do i = i0,ll
hes(i,ll) = ddot (n,v(1,i),1,vnew,1)
tem = -hes(i,ll)
call daxpy (n,tem,v(1,i),1,vnew,1)
10  continue
end do
snormw = dnrm2 (n,vnew,1)
if (vnrm + 0.001d0*snormw .ne. vnrm) return
sumdsq = 0.0d0
do i = i0,ll
tem = -ddot (n,v(1,i),1,vnew,1)
if (hes(i,ll) + 0.001d0*tem .eq. hes(i,ll)) go to 30
hes(i,ll) = hes(i,ll) - tem
call daxpy (n,tem,v(1,i),1,vnew,1)
sumdsq = sumdsq + tem**2
30  continue
end do
if (sumdsq .eq. 0.0d0) return
arg = max(0.0d0,snormw**2 - sumdsq)
snormw = sqrt(arg)
return
end subroutine dorthog

!!$*DECK DSPIGMR

subroutine dspigmr (neq,tn,y,savf,b,wght,n,maxl,maxlp1,&
kmp,delta,hl0,jpre,mnewt,f,psol,npsl,x,v,hes,q,&
lgmr,wp,iwp,wk,dl,iflag)
external f,psol
integer::neq(*),n,maxl,maxlp1,kmp,jpre,mnewt,npsl,lgmr,iwp(*),iflag
real*8::tn,y(*),savf(*),b(*),wght(*),delta,hl0,x(*),v(n,*),hes(maxlp1,*)
real*8::q(*),wp(*),wk(*),dl(*)
integer::i,ier,info,ip1,i2,j,k,ll,llp1
real*8::bnrm,bnrm0,c,dlnrm,prod,rho,s,snormw,dnrm2,tem
iflag = 0
lgmr = 0
npsl = 0
do i = 1,n
10  v(i,1) = b(i)*wght(i)
end do
bnrm0 = dnrm2 (n,v,1)
bnrm = bnrm0
if (bnrm0 .gt. delta) go to 30
if (mnewt .gt. 0) go to 20
call dcopy (n,b,1,x,1)
return
20 do i = 1,n
25  x(i) = 0.0d0
end do
return
30 continue
ier = 0
if (jpre .eq. 0 .or. jpre .eq. 2) go to 55
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,b,1,ier)
npsl = 1
if (ier .ne. 0) go to 300
do i = 1,n
50  v(i,1) = b(i)*wght(i)
end do
bnrm = dnrm2 (n,v,1)
delta = delta*(bnrm/bnrm0)
55 tem = 1.0d0/bnrm
call dscal (n,tem,v(1,1),1)
do j = 1,maxl
do i = 1,maxlp1
60   hes(i,j) = 0.0d0
end do
65  continue
end do
prod = 1.0d0
do ll = 1,maxl
lgmr = ll
call datv (neq,y,savf,v(1,ll),wght,x,f,psol,v(1,ll+1),&
wk,wp,iwp,hl0,jpre,ier,npsl)
if (ier .ne. 0) go to 300
call dorthog (v(1,ll+1),v,hes,n,ll,maxlp1,kmp,snormw)
hes(ll+1,ll) = snormw
call dheqr (hes,maxlp1,ll,q,info,ll)
if (info .eq. ll) go to 120
prod = prod*q(2*ll)
rho = abs(prod*bnrm)
if ((ll.gt.kmp) .and. (kmp.lt.maxl)) then
if (ll .eq. kmp+1) then
call dcopy (n,v(1,1),1,dl,1)
do i = 1,kmp
ip1 = i + 1
i2 = i*2
s = q(i2)
c = q(i2-1)
do k = 1,n
70        dl(k) = s*dl(k) + c*v(k,ip1)
end do
75      continue
end do
endif
s = q(2*ll)
c = q(2*ll-1)/snormw
llp1 = ll + 1
do k = 1,n
80     dl(k) = s*dl(k) + c*v(k,llp1)
end do
dlnrm = dnrm2 (n,dl,1)
rho = rho*dlnrm
endif
if (rho .le. delta) go to 200
if (ll .eq. maxl) go to 100
tem = 1.0d0/snormw
call dscal (n,tem,v(1,ll+1),1)
90  continue
end do
100 continue
if (rho .le. 1.0d0) go to 150
if (rho .le. bnrm .and. mnewt .eq. 0) go to 150
120 continue
iflag = 2
return
150 iflag = 1
200 continue
ll = lgmr
llp1 = ll + 1
do k = 1,llp1
210 b(k) = 0.0d0
end do
b(1) = bnrm
call dhels (hes,maxlp1,ll,q,b)
do k = 1,n
220 x(k) = 0.0d0
end do
do i = 1,ll
call daxpy (n,b(i),v(1,i),1,x,1)
230 continue
end do
do i = 1,n
240 x(i) = x(i)/wght(i)
end do
if (jpre .le. 1) return
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,x,2,ier)
npsl = npsl + 1
if (ier .ne. 0) go to 300
return
300 continue
if (ier .lt. 0) iflag = -1
if (ier .gt. 0) iflag = 3
return
end subroutine dspigmr

!!$*DECK DPCG

subroutine dpcg (neq,tn,y,savf,r,wght,n,maxl,delta,hl0,&
jpre,mnewt,f,psol,npsl,x,p,w,z,lpcg,wp,iwp,wk,iflag)
external f,psol
integer::neq(*),n,maxl,jpre,mnewt,npsl,lpcg,iwp(*),iflag
real*8::tn,delta,hl0
real*8::y(*),savf(*),r(*),wght(*),x(*),p(*),w(*),&
z(*),wp(*),wk(*)
integer::i,ier
real*8::alpha,beta,bnrm,ptw,rnrm,ddot,dvnorm,ztr,ztr0
iflag = 0
npsl = 0
lpcg = 0
do i = 1,n
10  x(i) = 0.0d0
end do
bnrm = dvnorm (n,r,wght)
if (bnrm .gt. delta) go to 20
if (mnewt .gt. 0) return
call dcopy (n,r,1,x,1)
return
20 ztr = 0.0d0
30 continue
lpcg = lpcg + 1
call dcopy (n,r,1,z,1)
ier = 0
if (jpre .eq. 0) go to 40
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,z,3,ier)
npsl = npsl + 1
if (ier .ne. 0) go to 100
40 continue
ztr0 = ztr
ztr = ddot (n,z,1,r,1)
if (lpcg .ne. 1) go to 50
call dcopy (n,z,1,p,1)
go to 70
50 continue
if (ztr0 .eq. 0.0d0) go to 200
beta = ztr/ztr0
do i = 1,n
60  p(i) = z(i) + beta*p(i)
end do
70 continue
call datp (neq,y,savf,p,wght,hl0,wk,f,w)
ptw = ddot (n,p,1,w,1)
if (ptw .eq. 0.0d0) go to 200
alpha = ztr/ptw
call daxpy (n,alpha,p,1,x,1)
alpha = -alpha
call daxpy (n,alpha,w,1,r,1)
rnrm = dvnorm (n,r,wght)
if (rnrm .le. delta) return
if (lpcg .lt. maxl) go to 30
iflag = 2
if (rnrm .le. 1.0d0) iflag = 1
if (rnrm .le. bnrm .and. mnewt .eq. 0) iflag = 1
return
100 continue
if (ier .lt. 0) iflag = -1
if (ier .gt. 0) iflag = 3
return
200 continue
iflag = 4
return
end subroutine dpcg

!!$*DECK DPCGS

subroutine dpcgs (neq,tn,y,savf,r,wght,n,maxl,delta,hl0,&
jpre,mnewt,f,psol,npsl,x,p,w,z,lpcg,wp,iwp,wk,iflag)
external f,psol
integer::neq(*),n,maxl,jpre,mnewt,npsl,lpcg,iwp(*),iflag
real*8::tn,delta,hl0
real*8::y(*),savf(*),r(*),wght(*),x(*),p(*),w(*),&
z(*),wp(*),wk(*)
integer::i,ier
real*8::alpha,beta,bnrm,ptw,rnrm,dvnorm,ztr,ztr0
iflag = 0
npsl = 0
lpcg = 0
do i = 1,n
10  x(i) = 0.0d0
end do
bnrm = dvnorm (n,r,wght)
if (bnrm .gt. delta) go to 20
if (mnewt .gt. 0) return
call dcopy (n,r,1,x,1)
return
20 ztr = 0.0d0
30 continue
lpcg = lpcg + 1
call dcopy (n,r,1,z,1)
ier = 0
if (jpre .eq. 0) go to 40
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,z,3,ier)
npsl = npsl + 1
if (ier .ne. 0) go to 100
40 continue
ztr0 = ztr
ztr = 0.0d0
do i = 1,n
45  ztr = ztr + z(i)*r(i)*wght(i)**2
end do
if (lpcg .ne. 1) go to 50
call dcopy (n,z,1,p,1)
go to 70
50 continue
if (ztr0 .eq. 0.0d0) go to 200
beta = ztr/ztr0
do i = 1,n
60  p(i) = z(i) + beta*p(i)
end do
70 continue
call datp (neq,y,savf,p,wght,hl0,wk,f,w)
ptw = 0.0d0
do i = 1,n
80  ptw = ptw + p(i)*w(i)*wght(i)**2
end do
if (ptw .eq. 0.0d0) go to 200
alpha = ztr/ptw
call daxpy (n,alpha,p,1,x,1)
alpha = -alpha
call daxpy (n,alpha,w,1,r,1)
rnrm = dvnorm (n,r,wght)
if (rnrm .le. delta) return
if (lpcg .lt. maxl) go to 30
iflag = 2
if (rnrm .le. 1.0d0) iflag = 1
if (rnrm .le. bnrm .and. mnewt .eq. 0) iflag = 1
return
100 continue
if (ier .lt. 0) iflag = -1
if (ier .gt. 0) iflag = 3
return
200 continue
iflag = 4
return
end subroutine dpcgs

!!$*DECK DATP

subroutine datp (neq,y,savf,p,wght,hl0,wk,f,w)
external f
integer::neq(*)
real*8::hl0
real*8::y(*),savf(*),p(*),wght(*),wk(*),w(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i
real*8::fac,pnrm,rpnrm,dvnorm
pnrm = dvnorm (n,p,wght)
rpnrm = 1.0d0/pnrm
call dcopy (n,y,1,w,1)
do i = 1,n
20  y(i) = w(i) + p(i)*rpnrm
end do
call f (neq,tn,y,wk)
nfe = nfe + 1
call dcopy (n,w,1,y,1)
fac = hl0*pnrm
do i = 1,n
40  w(i) = p(i) - fac*(wk(i) - savf(i))
end do
return
end subroutine datp

!!$*DECK DUSOL

subroutine dusol (neq,tn,y,savf,b,wght,n,delta,hl0,mnewt,&
psol,npsl,x,wp,iwp,wk,iflag)
external psol
integer::neq(*),n,mnewt,npsl,iwp(*),iflag
real*8::tn,delta,hl0
real*8::y(*),savf(*),b(*),wght(*),x(*),&
wp(*),wk(*)
integer::i,ier
real*8::bnrm,dvnorm
iflag = 0
npsl = 0
bnrm = dvnorm (n,b,wght)
if (bnrm .gt. delta) go to 30
if (mnewt .gt. 0) go to 10
call dcopy (n,b,1,x,1)
return
10 do i = 1,n
20  x(i) = 0.0d0
end do
return
30 ier = 0
call psol (neq,tn,y,savf,wk,hl0,wp,iwp,b,0,ier)
npsl = 1
if (ier .ne. 0) go to 100
call dcopy (n,b,1,x,1)
return
100 continue
if (ier .lt. 0) iflag = -1
if (ier .gt. 0) iflag = 3
return
end subroutine dusol

!!$*DECK DSRCPK

subroutine dsrcpk (rsav,isav,job)
integer::isav(*),job
integer::ils,ilsp
integer::i,lenilp,lenrlp,lenils,lenrls
real*8::rsav(*),rls,rlsp
save lenrls,lenils,lenrlp,lenilp
common /dls001/ rls(218),ils(37)
common /dlpk01/ rlsp(4),ilsp(13)
data lenrls/218/,lenils/37/,lenrlp/4/,lenilp/13/
if (job .eq. 2) go to 100
call dcopy (lenrls,rls,1,rsav,1)
call dcopy (lenrlp,rlsp,1,rsav(lenrls+1),1)
do i = 1,lenils
20  isav(i) = ils(i)
end do
do i = 1,lenilp
40  isav(lenils+i) = ilsp(i)
end do
return
100 continue
call dcopy (lenrls,rsav,1,rls,1)
call dcopy (lenrlp,rsav(lenrls+1),1,rlsp,1)
do i = 1,lenils
120 ils(i) = isav(i)
end do
do i = 1,lenilp
140 ilsp(i) = isav(lenils+i)
end do
return
end subroutine dsrcpk

!!$*DECK DHEFA

subroutine dhefa (a,lda,n,ipvt,info,job)
integer::lda,n,ipvt(*),info,job
real*8::a(lda,*)
integer::idamax,j,k,km1,kp1,l,nm1
real*8::t
if (job .gt. 1) go to 80
info = 0
nm1 = n - 1
if (nm1 .lt. 1) go to 70
do k = 1,nm1
kp1 = k + 1
l = idamax (2,a(k,k),1) + k - 1
ipvt(k) = l
if (a(l,k) .eq. 0.0d0) go to 40
if (l .eq. k) go to 10
t = a(l,k)
a(l,k) = a(k,k)
a(k,k) = t
10  continue
t = -1.0d0/a(k,k)
a(k+1,k) = a(k+1,k)*t
do j = kp1,n
t = a(l,j)
if (l .eq. k) go to 20
a(l,j) = a(k,j)
a(k,j) = t
20   continue
call daxpy (n-k,t,a(k+1,k),1,a(k+1,j),1)
30   continue
end do
go to 50
40  continue
info = k
50  continue
60  continue
end do
70 continue
ipvt(n) = n
if (a(n,n) .eq. 0.0d0) info = n
return
80 continue
nm1 = n - 1
if (nm1 .le. 1) go to 105
do k = 2,nm1
km1 = k - 1
l = ipvt(km1)
t = a(l,n)
if (l .eq. km1) go to 90
a(l,n) = a(km1,n)
a(km1,n) = t
90  continue
a(k,n) = a(k,n) + a(k,km1)*t
100 continue
end do
105 continue
info = 0
l = idamax (2,a(nm1,nm1),1) + nm1 - 1
ipvt(nm1) = l
if (a(l,nm1) .eq. 0.0d0) go to 140
if (l .eq. nm1) go to 110
t = a(l,nm1)
a(l,nm1) = a(nm1,nm1)
a(nm1,nm1) = t
110 continue
t = -1.0d0/a(nm1,nm1)
a(n,nm1) = a(n,nm1)*t
t = a(l,n)
if (l .eq. nm1) go to 120
a(l,n) = a(nm1,n)
a(nm1,n) = t
120 continue
a(n,n) = a(n,n) + t*a(n,nm1)
go to 150
140 continue
info = nm1
150 continue
ipvt(n) = n
if (a(n,n) .eq. 0.0d0) info = n
return
end subroutine dhefa

!!$*DECK DHESL

subroutine dhesl (a,lda,n,ipvt,b)
integer::lda,n,ipvt(*)
real*8::a(lda,*),b(*)
integer::k,kb,l,nm1
real*8::t
nm1 = n - 1
if (nm1 .lt. 1) go to 30
do k = 1,nm1
l = ipvt(k)
t = b(l)
if (l .eq. k) go to 10
b(l) = b(k)
b(k) = t
10  continue
b(k+1) = b(k+1) + t*a(k+1,k)
20  continue
end do
30 continue
do kb = 1,n
k = n + 1 - kb
b(k) = b(k)/a(k,k)
t = -b(k)
call daxpy (k-1,t,a(1,k),1,b(1),1)
40  continue
end do
return
end subroutine dhesl

!!$*DECK DHEQR

subroutine dheqr (a,lda,n,q,info,ijob)
integer::lda,n,info,ijob
real*8::a(lda,*),q(*)
integer::i,iq,j,k,km1,kp1,nm1
real*8::c,s,t,t1,t2
if (ijob .gt. 1) go to 70
info = 0
do k = 1,n
km1 = k - 1
kp1 = k + 1
if (km1 .lt. 1) go to 20
do j = 1,km1
i = 2*(j-1) + 1
t1 = a(j,k)
t2 = a(j+1,k)
c = q(i)
s = q(i+1)
a(j,k) = c*t1 - s*t2
a(j+1,k) = s*t1 + c*t2
10   continue
end do
20  continue
iq = 2*km1 + 1
t1 = a(k,k)
t2 = a(kp1,k)
if (t2 .ne. 0.0d0) go to 30
c = 1.0d0
s = 0.0d0
go to 50
30  continue
if (abs(t2) .lt. abs(t1)) go to 40
t = t1/t2
s = -1.0d0/sqrt(1.0d0+t*t)
c = -s*t
go to 50
40  continue
t = t2/t1
c = 1.0d0/sqrt(1.0d0+t*t)
s = -c*t
50  continue
q(iq) = c
q(iq+1) = s
a(k,k) = c*t1 - s*t2
if (a(k,k) .eq. 0.0d0) info = k
60  continue
end do
return
70 continue
nm1 = n - 1
do k = 1,nm1
i = 2*(k-1) + 1
t1 = a(k,n)
t2 = a(k+1,n)
c = q(i)
s = q(i+1)
a(k,n) = c*t1 - s*t2
a(k+1,n) = s*t1 + c*t2
100 continue
end do
info = 0
t1 = a(n,n)
t2 = a(n+1,n)
if (t2 .ne. 0.0d0) go to 110
c = 1.0d0
s = 0.0d0
go to 130
110 continue
if (abs(t2) .lt. abs(t1)) go to 120
t = t1/t2
s = -1.0d0/sqrt(1.0d0+t*t)
c = -s*t
go to 130
120 continue
t = t2/t1
c = 1.0d0/sqrt(1.0d0+t*t)
s = -c*t
130 continue
iq = 2*n - 1
q(iq) = c
q(iq+1) = s
a(n,n) = c*t1 - s*t2
if (a(n,n) .eq. 0.0d0) info = n
return
end subroutine dheqr

!!$*DECK DHELS

subroutine dhels (a,lda,n,q,b)
integer::lda,n
real*8::a(lda,*),b(*),q(*)
integer::iq,k,kb,kp1
real*8::c,s,t,t1,t2
do k = 1,n
kp1 = k + 1
iq = 2*(k-1) + 1
c = q(iq)
s = q(iq+1)
t1 = b(k)
t2 = b(kp1)
b(k) = c*t1 - s*t2
b(kp1) = s*t1 + c*t2
20  continue
end do
do kb = 1,n
k = n + 1 - kb
b(k) = b(k)/a(k,k)
t = -b(k)
call daxpy (k-1,t,a(1,k),1,b(1),1)
40  continue
end do
return
end subroutine dhels

!!$*DECK DLHIN

subroutine dlhin (neq,n,t0,y0,ydot,f,tout,uround,&
ewt,itol,atol,y,temp,h0,niter,ier)
external f
real*8::t0,tout,uround,h0
integer::neq(*),n,itol,niter,ier
real*8::y0(*),ydot(*),ewt(*),atol(*),y(*),temp(*)
real*8::afi,atoli,delyi,half,hg,hlb,hnew,hrat,&
hub,hun,pt1,t1,tdist,tround,two,dvnorm,yddnrm
integer::i,iter
save half,hun,pt1,two
data half /0.5d0/,hun /100.0d0/,pt1 /0.1d0/,two /2.0d0/
niter = 0
tdist = abs(tout - t0)
tround = uround*max(abs(t0),abs(tout))
if (tdist .lt. two*tround) go to 100
hlb = hun*tround
hub = pt1*tdist
atoli = atol(1)
do i = 1,n
if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
delyi = pt1*abs(y0(i)) + atoli
afi = abs(ydot(i))
if (afi*hub .gt. delyi) hub = delyi/afi
10  continue
end do
iter = 0
hg = sqrt(hlb*hub)
if (hub .lt. hlb) then
h0 = hg
go to 90
endif
50 continue
t1 = t0 + hg
do i = 1,n
60  y(i) = y0(i) + hg*ydot(i)
end do
call f (neq,t1,y,temp)
do i = 1,n
70  temp(i) = (temp(i) - ydot(i))/hg
end do
yddnrm = dvnorm (n,temp,ewt)
if (yddnrm*hub*hub .gt. two) then
hnew = sqrt(two/yddnrm)
else
hnew = sqrt(hg*hub)
endif
iter = iter + 1
if (iter .ge. 4) go to 80
hrat = hnew/hg
if ( (hrat .gt. half) .and. (hrat .lt. two) ) go to 80
if ( (iter .ge. 2) .and. (hnew .gt. two*hg) ) then
hnew = hg
go to 80
endif
hg = hnew
go to 50
80 h0 = hnew*half
if (h0 .lt. hlb) h0 = hlb
if (h0 .gt. hub) h0 = hub
90 h0 = sign(h0,tout - t0)
call dcopy (n,y0,1,y,1)
niter = iter
ier = 0
return
100 ier = -1
return
end subroutine dlhin

!!$*DECK DSTOKA

subroutine dstoka (neq,y,yh,nyh,yh1,ewt,savf,savx,acor,&
wm,iwm,f,jac,psol)
external f,jac,psol
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),&
savx(*),acor(*),wm(*)
integer::iownd,ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::newt,nsfi,nslj,njev
integer::jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
real*8::conit,crate,el,elco,hold,rmax,tesco,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::stifr
real*8::delt,epcon,sqrtn,rsqrtn
common /dls001/ conit,crate,el(13),elco(13,12),&
hold,rmax,tesco(3,12),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dls002/ stifr,newt,nsfi,nslj,njev
common /dlpk01/ delt,epcon,sqrtn,rsqrtn,&
jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
integer::i,i1,iredo,iret,j,jb,jok,m,ncf,newq,nslow
real*8::dcon,ddn,del,delp,drc,dsm,dup,exdn,exsm,&
exup,dfnorm,r,rh,rhdn,rhsm,rhup,roc,stiff,told,dvnorm
kflag = 0
told = tn
ncf = 0
ierpj = 0
iersl = 0
jcur = 0
icf = 0
delp = 0.0d0
if (jstart .gt. 0) go to 200
if (jstart .eq. -1) go to 100
if (jstart .eq. -2) go to 160
lmax = maxord + 1
nq = 1
l = 2
ialth = 2
rmax = 10000.0d0
rc = 0.0d0
el0 = 1.0d0
crate = 0.7d0
hold = h
meo = meth
nslp = 0
nslj = 0
ipup = 0
iret = 3
newt = 0
stifr = 0.0d0
go to 140
100 ipup = miter
lmax = maxord + 1
if (ialth .eq. 1) ialth = 2
if (meth .eq. meo) go to 110
call dcfode (meth,elco,tesco)
meo = meth
if (nq .gt. maxord) go to 120
ialth = l
iret = 1
go to 150
110 if (nq .le. maxord) go to 160
120 nq = maxord
l = lmax
do i = 1,l
125 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
epcon = conit*tesco(2,nq)
ddn = dvnorm (n,savf,ewt)/tesco(1,l)
exdn = 1.0d0/l
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
rh = min(rhdn,1.0d0)
iredo = 3
if (h .eq. hold) go to 170
rh = min(rh,abs(h/hold))
h = hold
go to 175
140 call dcfode (meth,elco,tesco)
150 do i = 1,l
155 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
epcon = conit*tesco(2,nq)
go to (160,170,200),iret
160 if (h .eq. hold) go to 200
rh = h/hold
h = hold
iredo = 3
go to 175
170 rh = max(rh,hmin/abs(h))
175 rh = min(rh,rmax)
rh = rh/max(1.0d0,abs(h)*hmxi*rh)
r = 1.0d0
do j = 2,l
r = r*rh
do i = 1,n
180   yh(i,j) = yh(i,j)*r
end do
end do
h = h*rh
rc = rc*rh
ialth = l
if (iredo .eq. 0) go to 690
200 if (newt .eq. 0 .or. jacflg .eq. 0) then
drc = 0.0d0
ipup = 0
crate = 0.7d0
else
drc = abs(rc - 1.0d0)
if (drc .gt. ccmax) ipup = miter
if (nst .ge. nslp+msbp) ipup = miter
endif
tn = tn + h
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
210   yh1(i) = yh1(i) + yh1(i+nyh)
end do
215 continue
end do
220 m = 0
mnewt = 0
stiff = 0.0d0
roc = 0.05d0
nslow = 0
do i = 1,n
230 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
if (newt .eq. 0 .or. ipup .le. 0) go to 250
jok = 1
if (nst .eq. 0 .or. nst .gt. nslj+50) jok = -1
if (icf .eq. 1 .and. drc .lt. 0.2d0) jok = -1
if (icf .eq. 2) jok = -1
if (jok .eq. -1) then
nslj = nst
njev = njev + 1
endif
call dsetpk (neq,y,yh1,ewt,acor,savf,jok,wm,iwm,f,jac)
ipup = 0
rc = 1.0d0
drc = 0.0d0
nslp = nst
crate = 0.7d0
if (ierpj .ne. 0) go to 430
250 do i = 1,n
260 acor(i) = 0.0d0
end do
270 if (newt .ne. 0) go to 350
do i = 1,n
savf(i) = h*savf(i) - yh(i,2)
290 y(i) = savf(i) - acor(i)
end do
del = dvnorm (n,y,ewt)
do i = 1,n
y(i) = yh(i,1) + el(1)*savf(i)
300 acor(i) = savf(i)
end do
stiff = 1.0d0
go to 400
350 do i = 1,n
360 savx(i) = h*savf(i) - (yh(i,2) + acor(i))
end do
dfnorm = dvnorm (n,savx,ewt)
call dsolpk (neq,y,savf,savx,ewt,wm,iwm,f,psol)
if (iersl .lt. 0) go to 430
if (iersl .gt. 0) go to 410
del = dvnorm (n,savx,ewt)
if (del .gt. 1.0d-8) stiff = max(stiff,dfnorm/del)
do i = 1,n
acor(i) = acor(i) + savx(i)
380 y(i) = yh(i,1) + el(1)*acor(i)
end do
400 if (m .ne. 0) then
roc = max(0.05d0,del/delp)
crate = max(0.2d0*crate,roc)
endif
dcon = del*min(1.0d0,1.5d0*crate)/epcon
if (dcon .le. 1.0d0) go to 450
m = m + 1
if (m .eq. maxcor) go to 410
if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
if (roc .gt. 10.0d0) go to 410
if (roc .gt. 0.8d0) nslow = nslow + 1
if (nslow .ge. 2) go to 410
mnewt = m
delp = del
call f (neq,tn,y,savf)
nfe = nfe + 1
go to 270
410 icf = 1
if (newt .eq. 0) then
if (nst .eq. 0) go to 430
if (miter .eq. 0) go to 430
newt = miter
stifr = 1023.0d0
ipup = miter
go to 220
endif
if (jcur.eq.1 .or. jacflg.eq.0) go to 430
ipup = miter
go to 220
430 icf = 2
ncf = ncf + 1
ncfn = ncfn + 1
rmax = 2.0d0
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
440   yh1(i) = yh1(i) - yh1(i+nyh)
end do
445 continue
end do
if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
if (abs(h) .le. hmin*1.00001d0) go to 670
if (ncf .eq. mxncf) go to 670
rh = 0.5d0
ipup = miter
iredo = 1
go to 170
450 jcur = 0
if (newt .gt. 0) stifr = 0.5d0*(stifr + stiff)
if (m .eq. 0) dsm = del/tesco(2,nq)
if (m .gt. 0) dsm = dvnorm (n,acor,ewt)/tesco(2,nq)
if (dsm .gt. 1.0d0) go to 500
kflag = 0
iredo = 0
nst = nst + 1
if (newt .eq. 0) nsfi = nsfi + 1
if (newt .gt. 0 .and. stifr .lt. 1.5d0) newt = 0
hu = h
nqu = nq
do j = 1,l
do i = 1,n
470   yh(i,j) = yh(i,j) + el(j)*acor(i)
end do
end do
ialth = ialth - 1
if (ialth .eq. 0) go to 520
if (ialth .gt. 1) go to 700
if (l .eq. lmax) go to 700
do i = 1,n
490 yh(i,lmax) = acor(i)
end do
go to 700
500 kflag = kflag - 1
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
510   yh1(i) = yh1(i) - yh1(i+nyh)
end do
515 continue
end do
rmax = 2.0d0
if (abs(h) .le. hmin*1.00001d0) go to 660
if (kflag .le. -3) go to 640
iredo = 2
rhup = 0.0d0
go to 540
520 rhup = 0.0d0
if (l .eq. lmax) go to 540
do i = 1,n
530 savf(i) = acor(i) - yh(i,lmax)
end do
dup = dvnorm (n,savf,ewt)/tesco(3,nq)
exup = 1.0d0/(l+1)
rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
540 exsm = 1.0d0/l
rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rhdn = 0.0d0
if (nq .eq. 1) go to 560
ddn = dvnorm (n,yh(1,l),ewt)/tesco(1,nq)
exdn = 1.0d0/nq
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
560 if (rhsm .ge. rhup) go to 570
if (rhup .gt. rhdn) go to 590
go to 580
570 if (rhsm .lt. rhdn) go to 580
newq = nq
rh = rhsm
go to 620
580 newq = nq - 1
rh = rhdn
if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
go to 620
590 newq = l
rh = rhup
if (rh .lt. 1.1d0) go to 610
r = el(l)/l
do i = 1,n
600 yh(i,newq+1) = acor(i)*r
end do
go to 630
610 ialth = 3
go to 700
620 if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
if (kflag .le. -2) rh = min(rh,0.2d0)
if (newq .eq. nq) go to 170
630 nq = newq
l = nq + 1
iret = 2
go to 150
640 if (kflag .eq. -10) go to 660
rh = 0.1d0
rh = max(hmin/abs(h),rh)
h = h*rh
do i = 1,n
645 y(i) = yh(i,1)
end do
call f (neq,tn,y,savf)
nfe = nfe + 1
do i = 1,n
650 yh(i,2) = h*savf(i)
end do
ipup = miter
ialth = 5
if (nq .eq. 1) go to 200
nq = 1
l = 2
iret = 3
go to 150
660 kflag = -1
go to 720
670 kflag = -2
go to 720
680 kflag = -3
go to 720
690 rmax = 10.0d0
700 r = 1.0d0/tesco(2,nqu)
do i = 1,n
710 acor(i) = acor(i)*r
end do
720 hold = h
jstart = 1
return
end subroutine dstoka

!!$*DECK DSETPK

subroutine dsetpk (neq,y,ysv,ewt,ftem,savf,jok,wm,iwm,&
f,jac)
external f,jac
integer::neq(*),jok,iwm(*)
real*8::y(*),ysv(*),ewt(*),ftem(*),savf(*),wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::delt,epcon,sqrtn,rsqrtn
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlpk01/ delt,epcon,sqrtn,rsqrtn,&
jpre,jacflg,locwp,lociwp,lsavx,kmp,maxl,mnewt,&
nni,nli,nps,ncfn,ncfl
integer::ier
real*8::hl0
ierpj = 0
jcur = 0
if (jok .eq. -1) jcur = 1
hl0 = el0*h
call jac (f,neq,tn,y,ysv,ewt,savf,ftem,hl0,jok,&
wm(locwp),iwm(lociwp),ier)
nje = nje + 1
if (ier .eq. 0) return
ierpj = 1
return
end subroutine dsetpk

!!$*DECK DSRCKR

subroutine dsrckr (rsav,isav,job)
integer::isav(*),job
integer::ils,ils2,ilsr,ilsp
integer::i,ioff,lenilp,lenrlp,lenils,lenrls,lenilr,lenrlr
real*8::rsav(*),rls,rls2,rlsr,rlsp
save lenrls,lenils,lenrlp,lenilp,lenrlr,lenilr
common /dls001/ rls(218),ils(37)
common /dls002/ rls2,ils2(4)
common /dlsr01/ rlsr(5),ilsr(9)
common /dlpk01/ rlsp(4),ilsp(13)
data lenrls/218/,lenils/37/,lenrlp/4/,lenilp/13/
data lenrlr/5/,lenilr/9/
if (job .eq. 2) go to 100
call dcopy (lenrls,rls,1,rsav,1)
rsav(lenrls+1) = rls2
call dcopy (lenrlr,rlsr,1,rsav(lenrls+2),1)
call dcopy (lenrlp,rlsp,1,rsav(lenrls+lenrlr+2),1)
do i = 1,lenils
20  isav(i) = ils(i)
end do
isav(lenils+1) = ils2(1)
isav(lenils+2) = ils2(2)
isav(lenils+3) = ils2(3)
isav(lenils+4) = ils2(4)
ioff = lenils + 2
do i = 1,lenilr
30  isav(ioff+i) = ilsr(i)
end do
ioff = ioff + lenilr
do i = 1,lenilp
40  isav(ioff+i) = ilsp(i)
end do
return
100 continue
call dcopy (lenrls,rsav,1,rls,1)
rls2 = rsav(lenrls+1)
call dcopy (lenrlr,rsav(lenrls+2),1,rlsr,1)
call dcopy (lenrlp,rsav(lenrls+lenrlr+2),1,rlsp,1)
do i = 1,lenils
120 ils(i) = isav(i)
end do
ils2(1) = isav(lenils+1)
ils2(2) = isav(lenils+2)
ils2(3) = isav(lenils+3)
ils2(4) = isav(lenils+4)
ioff = lenils + 2
do i = 1,lenilr
130 ilsr(i) = isav(ioff+i)
end do
ioff = ioff + lenilr
do i = 1,lenilp
140 ilsp(i) = isav(ioff+i)
end do
return
end subroutine dsrckr

!!$*DECK DAINVG

subroutine dainvg (res,adda,neq,t,y,ydot,miter,&
ml,mu,pw,ipvt,ier )
external res,adda
integer::neq,miter,ml,mu,ipvt(*),ier
integer::i,lenpw,mlp1,nrowpw
real*8::t,y(*),ydot(*),pw(*)
if (miter .ge. 4) go to 100
lenpw = neq*neq
do i = 1,lenpw
10  pw(i) = 0.0d0
end do
ier = 1
call res ( neq,t,y,pw,ydot,ier )
if (ier .gt. 1) return
call adda ( neq,t,y,0,0,pw,neq )
call dgefa ( pw,neq,neq,ipvt,ier )
if (ier .eq. 0) go to 20
ier = -ier
return
20 call dgesl ( pw,neq,neq,ipvt,ydot,0 )
return
100 continue
nrowpw = 2*ml + mu + 1
lenpw = neq * nrowpw
do i = 1,lenpw
110 pw(i) = 0.0d0
end do
ier = 1
call res ( neq,t,y,pw,ydot,ier )
if (ier .gt. 1) return
mlp1 = ml + 1
call adda ( neq,t,y,ml,mu,pw(mlp1),nrowpw )
call dgbfa ( pw,nrowpw,neq,ml,mu,ipvt,ier )
if (ier .eq. 0) go to 120
ier = -ier
return
120 call dgbsl ( pw,nrowpw,neq,ml,mu,ipvt,ydot,0 )
return
end subroutine dainvg

!!$*DECK DSTODI

subroutine dstodi (neq,y,yh,nyh,yh1,ewt,savf,savr,&
acor,wm,iwm,res,adda,jac,pjac,slvs )
external res,adda,jac,pjac,slvs
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),&
savr(*),acor(*),wm(*)
integer::iownd,ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::conit,crate,el,elco,hold,rmax,tesco,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ conit,crate,el(13),elco(13,12),&
hold,rmax,tesco(3,12),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),ialth,ipup,lmax,meo,nqnyh,nslp,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,i1,iredo,ires,iret,j,jb,kgo,m,ncf,newq
real*8::dcon,ddn,del,delp,dsm,dup,&
eljh,el1h,exdn,exsm,exup,&
r,rh,rhdn,rhsm,rhup,told,dvnorm
kflag = 0
told = tn
ncf = 0
ierpj = 0
iersl = 0
jcur = 0
icf = 0
delp = 0.0d0
if (jstart .gt. 0) go to 200
if (jstart .eq. -1) go to 100
if (jstart .eq. -2) go to 160
lmax = maxord + 1
nq = 1
l = 2
ialth = 2
rmax = 10000.0d0
rc = 0.0d0
el0 = 1.0d0
crate = 0.7d0
hold = h
meo = meth
nslp = 0
ipup = miter
iret = 3
go to 140
100 ipup = miter
lmax = maxord + 1
if (ialth .eq. 1) ialth = 2
if (meth .eq. meo) go to 110
call dcfode (meth,elco,tesco)
meo = meth
if (nq .gt. maxord) go to 120
ialth = l
iret = 1
go to 150
110 if (nq .le. maxord) go to 160
120 nq = maxord
l = lmax
do i = 1,l
125 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
ddn = dvnorm (n,savf,ewt)/tesco(1,l)
exdn = 1.0d0/l
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
rh = min(rhdn,1.0d0)
iredo = 3
if (h .eq. hold) go to 170
rh = min(rh,abs(h/hold))
h = hold
go to 175
140 call dcfode (meth,elco,tesco)
150 do i = 1,l
155 el(i) = elco(i,nq)
end do
nqnyh = nq*nyh
rc = rc*el(1)/el0
el0 = el(1)
conit = 0.5d0/(nq+2)
go to (160,170,200),iret
160 if (h .eq. hold) go to 200
rh = h/hold
h = hold
iredo = 3
go to 175
170 rh = max(rh,hmin/abs(h))
175 rh = min(rh,rmax)
rh = rh/max(1.0d0,abs(h)*hmxi*rh)
r = 1.0d0
do j = 2,l
r = r*rh
do i = 1,n
180   yh(i,j) = yh(i,j)*r
end do
end do
h = h*rh
rc = rc*rh
ialth = l
if (iredo .eq. 0) go to 690
200 if (abs(rc-1.0d0) .gt. ccmax) ipup = miter
if (nst .ge. nslp+msbp) ipup = miter
tn = tn + h
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
210   yh1(i) = yh1(i) + yh1(i+nyh)
end do
215 continue
end do
220 m = 0
do i = 1,n
savf(i) = yh(i,2) / h
230 y(i) = yh(i,1)
end do
if (ipup .le. 0) go to 240
call pjac (neq,y,yh,nyh,ewt,acor,savr,savf,wm,iwm,&
res,jac,adda )
ipup = 0
rc = 1.0d0
nslp = nst
crate = 0.7d0
if (ierpj .eq. 0) go to 250
if (ierpj .lt. 0) go to 435
ires = ierpj
go to (430,435,430),ires
240 ires = 1
call res ( neq,tn,y,savf,savr,ires )
nfe = nfe + 1
kgo = abs(ires)
go to ( 250,435,430 ) ,kgo
250 do i = 1,n
260 acor(i) = 0.0d0
end do
270 continue
call slvs (wm,iwm,savr,savf)
if (iersl .lt. 0) go to 430
if (iersl .gt. 0) go to 410
el1h = el(1) * h
del = dvnorm (n,savr,ewt) * abs(h)
do i = 1,n
acor(i) = acor(i) + savr(i)
savf(i) = acor(i) + yh(i,2)/h
380 y(i) = yh(i,1) + el1h*acor(i)
end do
if (m .ne. 0) crate = max(0.2d0*crate,del/delp)
dcon = del*min(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
if (dcon .le. 1.0d0) go to 460
m = m + 1
if (m .eq. maxcor) go to 410
if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
delp = del
ires = 1
call res ( neq,tn,y,savf,savr,ires )
nfe = nfe + 1
kgo = abs(ires)
go to ( 270,435,410 ) ,kgo
410 icf = 1
if (jcur .eq. 1) go to 430
ipup = miter
go to 220
430 icf = 2
ncf = ncf + 1
rmax = 2.0d0
435 tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
440   yh1(i) = yh1(i) - yh1(i+nyh)
end do
445 continue
end do
if (ires .eq. 2) go to 680
if (ierpj .lt. 0 .or. iersl .lt. 0) go to 685
if (abs(h) .le. hmin*1.00001d0) go to 450
if (ncf .eq. mxncf) go to 450
rh = 0.25d0
ipup = miter
iredo = 1
go to 170
450 if (ires .eq. 3) go to 680
go to 670
460 jcur = 0
if (m .eq. 0) dsm = del/tesco(2,nq)
if (m .gt. 0) dsm = abs(h) * dvnorm (n,acor,ewt)/tesco(2,nq)
if (dsm .gt. 1.0d0) go to 500
kflag = 0
iredo = 0
nst = nst + 1
hu = h
nqu = nq
do j = 1,l
eljh = el(j)*h
do i = 1,n
470   yh(i,j) = yh(i,j) + eljh*acor(i)
end do
end do
ialth = ialth - 1
if (ialth .eq. 0) go to 520
if (ialth .gt. 1) go to 700
if (l .eq. lmax) go to 700
do i = 1,n
490 yh(i,lmax) = acor(i)
end do
go to 700
500 kflag = kflag - 1
tn = told
i1 = nqnyh + 1
do jb = 1,nq
i1 = i1 - nyh
do i = i1,nqnyh
510   yh1(i) = yh1(i) - yh1(i+nyh)
end do
515 continue
end do
rmax = 2.0d0
if (abs(h) .le. hmin*1.00001d0) go to 660
if (kflag .le. -7) go to 660
iredo = 2
rhup = 0.0d0
go to 540
520 rhup = 0.0d0
if (l .eq. lmax) go to 540
do i = 1,n
530 savf(i) = acor(i) - yh(i,lmax)
end do
dup = abs(h) * dvnorm (n,savf,ewt)/tesco(3,nq)
exup = 1.0d0/(l+1)
rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
540 exsm = 1.0d0/l
rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
rhdn = 0.0d0
if (nq .eq. 1) go to 560
ddn = dvnorm (n,yh(1,l),ewt)/tesco(1,nq)
exdn = 1.0d0/nq
rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
560 if (rhsm .ge. rhup) go to 570
if (rhup .gt. rhdn) go to 590
go to 580
570 if (rhsm .lt. rhdn) go to 580
newq = nq
rh = rhsm
go to 620
580 newq = nq - 1
rh = rhdn
if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
go to 620
590 newq = l
rh = rhup
if (rh .lt. 1.1d0) go to 610
r = h*el(l)/l
do i = 1,n
600 yh(i,newq+1) = acor(i)*r
end do
go to 630
610 ialth = 3
go to 700
620 if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
if (kflag .le. -2) rh = min(rh,0.1d0)
if (newq .eq. nq) go to 170
630 nq = newq
l = nq + 1
iret = 2
go to 150
660 kflag = -1
go to 720
670 kflag = -2
go to 720
680 kflag = -1 - ires
go to 720
685 kflag = -5
go to 720
690 rmax = 10.0d0
700 r = h/tesco(2,nqu)
do i = 1,n
710 acor(i) = acor(i)*r
end do
720 hold = h
jstart = 1
return
end subroutine dstodi

!!$*DECK DPREPJI

subroutine dprepji (neq,y,yh,nyh,ewt,rtem,savr,s,wm,iwm,&
res,jac,adda)
external res,jac,adda
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),ewt(*),rtem(*),&
s(*),savr(*),wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,i1,i2,ier,ii,ires,j,j1,jj,lenp,&
mba,mband,meb1,meband,ml,ml3,mu
real*8::con,fac,hl0,r,srur,yi,yj,yjj
nje = nje + 1
hl0 = h*el0
ierpj = 0
jcur = 1
go to (100,200,300,400,500),miter
100 ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
lenp = n*n
do i = 1,lenp
110 wm(i+2) = 0.0d0
end do
call jac ( neq,tn,y,s,0,0,wm(3),n )
con = -hl0
do i = 1,lenp
120 wm(i+2) = wm(i+2)*con
end do
go to 240
200 continue
ires = -1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
srur = wm(1)
j1 = 2
do j = 1,n
yj = y(j)
r = max(srur*abs(yj),0.01d0/ewt(j))
y(j) = y(j) + r
fac = -hl0/r
call res ( neq,tn,y,s,rtem,ires )
nfe = nfe + 1
if (ires .gt. 1) go to 600
do i = 1,n
220   wm(i+j1) = (rtem(i) - savr(i))*fac
end do
y(j) = yj
j1 = j1 + n
230 continue
end do
ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
240 continue
call adda(neq,tn,y,0,0,wm(3),n)
call dgefa (wm(3),n,n,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
300 return
400 ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
ml = iwm(1)
mu = iwm(2)
ml3 = ml + 3
mband = ml + mu + 1
meband = mband + ml
lenp = meband*n
do i = 1,lenp
410 wm(i+2) = 0.0d0
end do
call jac ( neq,tn,y,s,ml,mu,wm(ml3),meband)
con = -hl0
do i = 1,lenp
420 wm(i+2) = wm(i+2)*con
end do
go to 570
500 continue
ires = -1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
ml = iwm(1)
mu = iwm(2)
ml3 = ml + 3
mband = ml + mu + 1
mba = min(mband,n)
meband = mband + ml
meb1 = meband - 1
srur = wm(1)
do j = 1,mba
do i = j,n,mband
yi = y(i)
r = max(srur*abs(yi),0.01d0/ewt(i))
530   y(i) = y(i) + r
end do
call res ( neq,tn,y,s,rtem,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
do jj = j,n,mband
y(jj) = yh(jj,1)
yjj = y(jj)
r = max(srur*abs(yjj),0.01d0/ewt(jj))
fac = -hl0/r
i1 = max(jj-mu,1)
i2 = min(jj+ml,n)
ii = jj*meb1 - ml + 2
do i = i1,i2
540    wm(ii+i) = (rtem(i) - savr(i))*fac
end do
550   continue
end do
560 continue
end do
ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
570 continue
call adda(neq,tn,y,ml,mu,wm(ml3),meband)
call dgbfa (wm(3),meband,n,ml,mu,iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
600 ierpj = ires
return
end subroutine dprepji

!!$*DECK DAIGBT

subroutine daigbt (res,adda,neq,t,y,ydot,&
mb,nb,pw,ipvt,ier )
external res,adda
integer::neq(*),mb,nb,ipvt(*),ier
integer::i,lenpw,lblox,lpb,lpc
real*8::t,y(*),ydot(*),pw(*)
lblox = mb*mb*nb
lpb = 1 + lblox
lpc = lpb + lblox
lenpw = 3*lblox
do i = 1,lenpw
10  pw(i) = 0.0d0
end do
ier = 1
call res (neq,t,y,pw,ydot,ier)
if (ier .gt. 1) return
call adda (neq,t,y,mb,nb,pw(1),pw(lpb),pw(lpc) )
call ddecbt (mb,nb,pw,pw(lpb),pw(lpc),ipvt,ier)
if (ier .eq. 0) go to 20
ier = -ier
return
20 call dsolbt (mb,nb,pw,pw(lpb),pw(lpc),ydot,ipvt)
return
end subroutine daigbt

!!$*DECK DPJIBT

subroutine dpjibt (neq,y,yh,nyh,ewt,rtem,savr,s,wm,iwm,&
res,jac,adda)
external res,jac,adda
integer::neq(*),nyh,iwm(*)
real*8::y(*),yh(nyh,*),ewt(*),rtem(*),&
s(*),savr(*),wm(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::i,ier,iia,iib,iic,ipa,ipb,ipc,ires,j,j1,j2,&
k,k1,lenp,lblox,lpb,lpc,mb,mbsq,mwid,nb
real*8::con,fac,hl0,r,srur
nje = nje + 1
hl0 = h*el0
ierpj = 0
jcur = 1
mb = iwm(1)
nb = iwm(2)
mbsq = mb*mb
lblox = mbsq*nb
lpb = 3 + lblox
lpc = lpb + lblox
lenp = 3*lblox
go to (100,200),miter
100 ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
do i = 1,lenp
110 wm(i+2) = 0.0d0
end do
call jac (neq,tn,y,s,mb,nb,wm(3),wm(lpb),wm(lpc))
con = -hl0
do i = 1,lenp
120 wm(i+2) = wm(i+2)*con
end do
go to 260
200 continue
ires = -1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
mwid = 3*mb
srur = wm(1)
do i = 1,lenp
205 wm(2+i) = 0.0d0
end do
do k = 1,3
do j = 1,mb
j1 = j+(k-1)*mb
do i = j1,n,mwid
r = max(srur*abs(y(i)),0.01d0/ewt(i))
y(i) = y(i) + r
210    continue
end do
call res (neq,tn,y,s,rtem,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
do i = 1,n
215    rtem(i) = rtem(i) - savr(i)
end do
k1 = k
do i = j1,n,mwid
y(i) = yh(i,1)
r = max(srur*abs(y(i)),0.01d0/ewt(i))
fac = -hl0/r
iia = i - j
ipa = 2 + (j-1)*mb + (k1-1)*mbsq
do j2 = 1,mb
221      wm(ipa+j2) = rtem(iia+j2)*fac
end do
if (k1 .le. 1) go to 223
iib = iia - mb
ipb = ipa + lblox - mbsq
do j2 = 1,mb
222      wm(ipb+j2) = rtem(iib+j2)*fac
end do
223    continue
if (k1 .ge. nb) go to 225
iic = iia + mb
ipc = ipa + 2*lblox + mbsq
do j2 = 1,mb
224      wm(ipc+j2) = rtem(iic+j2)*fac
end do
225    continue
if (k1 .ne. 3) go to 227
ipc = ipa - 2*mbsq + 2*lblox
do j2 = 1,mb
226      wm(ipc+j2) = rtem(j2)*fac
end do
227    continue
if (k1 .ne. nb-2) go to 229
iib = n - mb
ipb = ipa + 2*mbsq + lblox
do j2 = 1,mb
228      wm(ipb+j2) = rtem(iib+j2)*fac
end do
229    k1 = k1 + 3
230    continue
end do
240   continue
end do
250 continue
end do
ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
260 continue
call adda (neq,tn,y,mb,nb,wm(3),wm(lpb),wm(lpc))
call ddecbt (mb,nb,wm(3),wm(lpb),wm(lpc),iwm(21),ier)
if (ier .ne. 0) ierpj = 1
return
600 ierpj = ires
return
end subroutine dpjibt

!!$*DECK DSLSBT

subroutine dslsbt (wm,iwm,x,tem)
integer::iwm(*)
integer::lblox,lpb,lpc,mb,nb
real*8::wm(*),x(*),tem(*)
mb = iwm(1)
nb = iwm(2)
lblox = mb*mb*nb
lpb = 3 + lblox
lpc = lpb + lblox
call dsolbt (mb,nb,wm(3),wm(lpb),wm(lpc),x,iwm(21))
return
end subroutine dslsbt

!!$*DECK DDECBT

subroutine ddecbt (m,n,a,b,c,ip,ier)
integer::m,n,ip(m,n),ier
real*8::a(m,m,n),b(m,m,n),c(m,m,n)
integer::nm1,nm2,km1,i,j,k
real*8::dp,ddot
if (m .lt. 1 .or. n .lt. 4) go to 210
nm1 = n - 1
nm2 = n - 2
call dgefa (a,m,m,ip,ier)
k = 1
if (ier .ne. 0) go to 200
do j = 1,m
call dgesl (a,m,m,ip,b(1,j,1),0)
call dgesl (a,m,m,ip,c(1,j,1),0)
10  continue
end do
do j = 1,m
do i = 1,m
dp = ddot (m,c(i,1,2),m,c(1,j,1),1)
b(i,j,2) = b(i,j,2) - dp
30   continue
end do
40  continue
end do
do k = 2,nm1
km1 = k - 1
do j = 1,m
do i = 1,m
dp = ddot (m,c(i,1,k),m,b(1,j,km1),1)
a(i,j,k) = a(i,j,k) - dp
60     continue
end do
70   continue
end do
call dgefa (a(1,1,k),m,m,ip(1,k),ier)
if (ier .ne. 0) go to 200
do j = 1,m
80   call dgesl (a(1,1,k),m,m,ip(1,k),b(1,j,k),0)
end do
100 continue
end do
do j = 1,m
do i = 1,m
dp = ddot (m,b(i,1,n),m,b(1,j,nm2),1)
c(i,j,n) = c(i,j,n) - dp
120   continue
end do
130 continue
end do
do j = 1,m
do i = 1,m
dp = ddot (m,c(i,1,n),m,b(1,j,nm1),1)
a(i,j,n) = a(i,j,n) - dp
150   continue
end do
160 continue
end do
call dgefa (a(1,1,n),m,m,ip(1,n),ier)
k = n
if (ier .ne. 0) go to 200
return
200 ier = k
return
210 ier = -1
return
end subroutine ddecbt

!!$*DECK DSOLBT

subroutine dsolbt (m,n,a,b,c,y,ip)
integer::m,n,ip(m,n)
real*8::a(m,m,n),b(m,m,n),c(m,m,n),y(m,n)
integer::nm1,nm2,i,k,kb,km1,kp1
real*8::dp,ddot
nm1 = n - 1
nm2 = n - 2
call dgesl (a,m,m,ip,y,0)
do k = 2,nm1
km1 = k - 1
do i = 1,m
dp = ddot (m,c(i,1,k),m,y(1,km1),1)
y(i,k) = y(i,k) - dp
20   continue
end do
call dgesl (a(1,1,k),m,m,ip(1,k),y(1,k),0)
30  continue
end do
do i = 1,m
dp = ddot (m,c(i,1,n),m,y(1,nm1),1)&
+ ddot (m,b(i,1,n),m,y(1,nm2),1)
y(i,n) = y(i,n) - dp
50  continue
end do
call dgesl (a(1,1,n),m,m,ip(1,n),y(1,n),0)
do kb = 1,nm1
k = n - kb
kp1 = k + 1
do i = 1,m
dp = ddot (m,b(i,1,k),m,y(1,kp1),1)
y(i,k) = y(i,k) - dp
70   continue
end do
80  continue
end do
do i = 1,m
dp = ddot (m,c(i,1,1),m,y(1,3),1)
y(i,1) = y(i,1) - dp
100 continue
end do
return
end subroutine dsolbt

!!$*DECK DIPREPI

subroutine diprepi (neq,y,s,rwork,ia,ja,ic,jc,ipflag,&
res,jac,adda)
external res,jac,adda
integer::neq(*),ia(*),ja(*),ic(*),jc(*),ipflag
real*8::y(*),s(*),rwork(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rlss
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,imax,lewtn,lyhd,lyhn
ipflag = 0
call dprepi (neq,y,s,rwork(lyh),rwork(lsavf),rwork(lewt),&
rwork(lacor),ia,ja,ic,jc,rwork(lwm),rwork(lwm),ipflag,&
res,jac,adda)
lenwk = max(lreq,lwmin)
if (ipflag .lt. 0) return
lyhn = lwm + lenwk
if (lyhn .gt. lyh) return
lyhd = lyh - lyhn
if (lyhd .eq. 0) go to 20
imax = lyhn - 1 + lenyhm
do i=lyhn,imax
10  rwork(i) = rwork(i+lyhd)
end do
lyh = lyhn
20 lsavf = lyh + lenyh
lewtn = lsavf + n
lacor = lewtn + n
if (istatc .eq. 3) go to 40
if (lewtn .gt. lewt) return
do i=1,n
30  rwork(i+lewtn-1) = rwork(i+lewt-1)
end do
40 lewt = lewtn
return
end subroutine diprepi

!!$*DECK DPREPI

subroutine dprepi (neq,y,s,yh,savr,ewt,rtem,ia,ja,ic,jc,&
wk,iwk,ipper,res,jac,adda)
external res,jac,adda
integer::neq(*),ia(*),ja(*),ic(*),jc(*),iwk(*),ipper
real*8::y(*),s(*),yh(*),savr(*),ewt(*),rtem(*),wk(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rlss
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,ibr,ier,ipil,ipiu,iptt1,iptt2,j,k,knew,kamax,&
kamin,kcmax,kcmin,ldif,lenigp,lenwk1,liwk,ljfo,maxg,&
np1,nzsut
real*8::erwt,fac,yj
ibian = lrat*2
ipian = ibian + 1
np1 = n + 1
ipjan = ipian + np1
ibjan = ipjan - 1
lenwk1 = lenwk - n
liwk = lenwk*lrat
if (moss .eq. 0) liwk = liwk - n
if (moss .eq. 1 .or. moss .eq. 2) liwk = lenwk1*lrat
if (ipjan+n-1 .gt. liwk) go to 310
if (moss .eq. 0) go to 30
if (istatc .eq. 3) go to 20
do i=1,n
erwt = 1.0d0/ewt(i)
fac = 1.0d0 + 1.0d0/(i + 1.0d0)
y(i) = y(i) + fac*sign(erwt,y(i))
s(i) = 1.0d0 + fac*erwt
10  continue
end do
go to (70,100,150,200),moss
20 continue
do i = 1,n
y(i) = yh(i)
25  s(i) = yh(n+i)
end do
go to (70,100,150,200),moss
30 knew = ipjan
kamin = ia(1)
kcmin = ic(1)
iwk(ipian) = 1
do j = 1,n
do i = 1,n
35   iwk(liwk+i) = 0
end do
kamax = ia(j+1) - 1
if (kamin .gt. kamax) go to 45
do k = kamin,kamax
i = ja(k)
iwk(liwk+i) = 1
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
40   continue
end do
45  kamin = kamax + 1
kcmax = ic(j+1) - 1
if (kcmin .gt. kcmax) go to 55
do k = kcmin,kcmax
i = jc(k)
if (iwk(liwk+i) .ne. 0) go to 50
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
50   continue
end do
55  iwk(ipian+j) = knew + 1 - ipjan
kcmin = kcmax + 1
60  continue
end do
go to 240
70 continue
ier = 1
call res (neq,tn,y,s,savr,ier)
if (ier .gt. 1) go to 370
do i = 1,n
savr(i) = 0.0d0
75  wk(lenwk1+i) = 0.0d0
end do
k = ipjan
iwk(ipian) = 1
do j = 1,n
call adda (neq,tn,y,j,iwk(ipian),iwk(ipjan),wk(lenwk1+1))
call jac (neq,tn,y,s,j,iwk(ipian),iwk(ipjan),savr)
do i = 1,n
ljfo = lenwk1 + i
if (wk(ljfo) .eq. 0.0d0) go to 80
wk(ljfo) = 0.0d0
savr(i) = 0.0d0
go to 85
80   if (savr(i) .eq. 0.0d0) go to 90
savr(i) = 0.0d0
85   if (k .gt. liwk) go to 310
iwk(k) = i
k = k+1
90   continue
end do
iwk(ipian+j) = k + 1 - ipjan
95  continue
end do
go to 240
100 do i = 1,n
105 wk(lenwk1+i) = 0.0d0
end do
k = ipjan
iwk(ipian) = 1
ier = -1
if (miter .eq. 1) ier = 1
call res (neq,tn,y,s,savr,ier)
if (ier .gt. 1) go to 370
do j = 1,n
call adda (neq,tn,y,j,iwk(ipian),iwk(ipjan),wk(lenwk1+1))
yj = y(j)
erwt = 1.0d0/ewt(j)
y(j) = yj + sign(erwt,yj)
call res (neq,tn,y,s,rtem,ier)
if (ier .gt. 1) return
y(j) = yj
do i = 1,n
ljfo = lenwk1 + i
if (wk(ljfo) .eq. 0.0d0) go to 110
wk(ljfo) = 0.0d0
go to 115
110   if (rtem(i) .eq. savr(i)) go to 120
115   if (k .gt. liwk) go to 310
iwk(k) = i
k = k + 1
120   continue
end do
iwk(ipian+j) = k + 1 - ipjan
130 continue
end do
go to 240
150 continue
ier = 1
call res (neq,tn,y,s,savr,ier)
if (ier .gt. 1) go to 370
do i = 1,n
155 savr(i) = 0.0d0
end do
knew = ipjan
kamin = ia(1)
iwk(ipian) = 1
do j = 1,n
call jac (neq,tn,y,s,j,iwk(ipian),iwk(ipjan),savr)
kamax = ia(j+1) - 1
if (kamin .gt. kamax) go to 170
do k = kamin,kamax
i = ja(k)
savr(i) = 0.0d0
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
160   continue
end do
170 kamin = kamax + 1
do i = 1,n
if (savr(i) .eq. 0.0d0) go to 180
savr(i) = 0.0d0
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
180   continue
end do
iwk(ipian+j) = knew + 1 - ipjan
190 continue
end do
go to 240
200 knew = ipjan
kamin = ia(1)
iwk(ipian) = 1
ier = -1
if (miter .eq. 1) ier = 1
call res (neq,tn,y,s,savr,ier)
if (ier .gt. 1) go to 370
do j = 1,n
yj = y(j)
erwt = 1.0d0/ewt(j)
y(j) = yj + sign(erwt,yj)
call res (neq,tn,y,s,rtem,ier)
if (ier .gt. 1) return
y(j) = yj
kamax = ia(j+1) - 1
if (kamin .gt. kamax) go to 225
do k = kamin,kamax
i = ja(k)
rtem(i) = savr(i)
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
220   continue
end do
225 kamin = kamax + 1
do i = 1,n
if (rtem(i) .eq. savr(i)) go to 230
if (knew .gt. liwk) go to 310
iwk(knew) = i
knew = knew + 1
230   continue
end do
iwk(ipian+j) = knew + 1 - ipjan
235 continue
end do
240 continue
if (moss .eq. 0 .or. istatc .eq. 3) go to 250
do i = 1,n
245 y(i) = yh(i)
end do
250 nnz = iwk(ipian+n) - 1
ipper = 0
ngp = 0
lenigp = 0
ipigp = ipjan + nnz
if (miter .ne. 2) go to 260
maxg = np1
ipjgp = ipjan + nnz
ibjgp = ipjgp - 1
ipigp = ipjgp + n
iptt1 = ipigp + np1
iptt2 = iptt1 + n
lreq = iptt2 + n - 1
if (lreq .gt. liwk) go to 320
call jgroup (n,iwk(ipian),iwk(ipjan),maxg,ngp,iwk(ipigp),&
iwk(ipjgp),iwk(iptt1),iwk(iptt2),ier)
if (ier .ne. 0) go to 320
lenigp = ngp + 1
260 ipr = ipigp + lenigp
ipc = ipr
ipic = ipc + n
ipisp = ipic + n
iprsp = (ipisp-2)/lrat + 2
iesp = lenwk + 1 - iprsp
if (iesp .lt. 0) go to 330
ibr = ipr - 1
do i = 1,n
270 iwk(ibr+i) = i
end do
nsp = liwk + 1 - ipisp
call odrv(n,iwk(ipian),iwk(ipjan),wk,iwk(ipr),iwk(ipic),nsp,&
iwk(ipisp),1,iys)
if (iys .eq. 11*n+1) go to 340
if (iys .ne. 0) go to 330
ipa = lenwk + 1 - nnz
nsp = ipa - iprsp
lreq = max(12*n/lrat,6*n/lrat+2*n+nnz) + 3
lreq = lreq + iprsp - 1 + nnz
if (lreq .gt. lenwk) go to 350
iba = ipa - 1
do i = 1,nnz
280 wk(iba+i) = 0.0d0
end do
ipisp = lrat*(iprsp - 1) + 1
call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
lreq = lenwk - iesp
if (iys .eq. 10*n+1) go to 350
if (iys .ne. 0) go to 360
ipil = ipisp
ipiu = ipil + 2*n + 1
nzu = iwk(ipil+n) - iwk(ipil)
nzl = iwk(ipiu+n) - iwk(ipiu)
if (lrat .gt. 1) go to 290
call adjlr (n,iwk(ipisp),ldif)
lreq = lreq + ldif
290 continue
if (lrat .eq. 2 .and. nnz .eq. n) lreq = lreq + 1
nsp = nsp + lreq - lenwk
ipa = lreq + 1 - nnz
iba = ipa - 1
ipper = 0
return
310 ipper = -1
lreq = 2 + (2*n + 1)/lrat
lreq = max(lenwk+1,lreq)
return
320 ipper = -2
lreq = (lreq - 1)/lrat + 1
return
330 ipper = -3
call cntnzu (n,iwk(ipian),iwk(ipjan),nzsut)
lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
return
340 ipper = -4
return
350 ipper = -5
return
360 ipper = -6
lreq = lenwk
return
370 ipper = -ier - 5
lreq = 2 + (2*n + 1)/lrat
return
end subroutine dprepi

!!$*DECK DAINVGS

subroutine dainvgs (neq,t,y,wk,iwk,tem,ydot,ier,res,adda)
external res,adda
integer::neq,iwk(*),ier
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,imul,j,k,kmin,kmax
real*8::t,rlss
real*8::y(*),wk(*),tem(*),ydot(*)
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
do i = 1,nnz
10  wk(iba+i) = 0.0d0
end do
ier = 1
call res (neq,t,y,wk(ipa),ydot,ier)
if (ier .gt. 1) return
kmin = iwk(ipian)
do j = 1,neq
kmax = iwk(ipian+j) - 1
do k = kmin,kmax
i = iwk(ibjan+k)
15   tem(i) = 0.0d0
end do
call adda (neq,t,y,j,iwk(ipian),iwk(ipjan),tem)
do k = kmin,kmax
i = iwk(ibjan+k)
20   wk(iba+k) = tem(i)
end do
kmin = kmax + 1
30  continue
end do
nlu = nlu + 1
ier = 0
do i = 1,neq
40  tem(i) = 0.0d0
end do
call cdrv (neq,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),tem,tem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
if (iys .eq. 0) go to 50
imul = (iys - 1)/neq
ier = 5
if (imul .eq. 8) ier = 1
if (imul .eq. 10) ier = 4
return
50 call cdrv (neq,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),ydot,ydot,nsp,iwk(ipisp),wk(iprsp),iesp,4,iys)
if (iys .ne. 0) ier = 5
return
end subroutine dainvgs

!!$*DECK DPRJIS

subroutine dprjis (neq,y,yh,nyh,ewt,rtem,savr,s,wk,iwk,&
res,jac,adda)
external res,jac,adda
integer::neq(*),nyh,iwk(*)
real*8::y(*),yh(nyh,*),ewt(*),rtem(*),&
s(*),savr(*),wk(*)
integer::iownd,iowns,&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
integer::iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
real*8::rowns,&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround
real*8::rlss
common /dls001/ rowns(209),&
ccmax,el0,h,hmin,hmxi,hu,rc,tn,uround,&
iownd(6),iowns(6),&
icf,ierpj,iersl,jcur,jstart,kflag,l,&
lyh,lewt,lacor,lsavf,lwm,liwm,meth,miter,&
maxord,maxcor,msbp,mxncf,n,nq,nst,nfe,nje,nqu
common /dlss01/ rlss(6),&
iplost,iesp,istatc,iys,iba,ibian,ibjan,ibjgp,&
ipian,ipjan,ipjgp,ipigp,ipr,ipc,ipic,ipisp,iprsp,ipa,&
lenyh,lenyhm,lenwk,lreq,lrat,lrest,lwmin,moss,msbj,&
nslj,ngp,nlu,nnz,nsp,nzl,nzu
integer::i,imul,ires,j,jj,jmax,jmin,k,kmax,kmin,ng
real*8::con,fac,hl0,r,srur
hl0 = h*el0
con = -hl0
jcur = 1
nje = nje + 1
go to (100,200),miter
100 ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
kmin = iwk(ipian)
do j = 1,n
kmax = iwk(ipian+j)-1
do i = 1,n
110   rtem(i) = 0.0d0
end do
call jac (neq,tn,y,s,j,iwk(ipian),iwk(ipjan),rtem)
do i = 1,n
120   rtem(i) = rtem(i)*con
end do
call adda (neq,tn,y,j,iwk(ipian),iwk(ipjan),rtem)
do k = kmin,kmax
i = iwk(ibjan+k)
wk(iba+k) = rtem(i)
125   continue
end do
kmin = kmax + 1
130 continue
end do
go to 290
200 continue
ires = -1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
srur = wk(1)
jmin = iwk(ipigp)
do ng = 1,ngp
jmax = iwk(ipigp+ng) - 1
do j = jmin,jmax
jj = iwk(ibjgp+j)
r = max(srur*abs(y(jj)),0.01d0/ewt(jj))
210   y(jj) = y(jj) + r
end do
call res (neq,tn,y,s,rtem,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
do j = jmin,jmax
jj = iwk(ibjgp+j)
y(jj) = yh(jj,1)
r = max(srur*abs(y(jj)),0.01d0/ewt(jj))
fac = -hl0/r
kmin = iwk(ibian+jj)
kmax = iwk(ibian+jj+1) - 1
do k = kmin,kmax
i = iwk(ibjan+k)
rtem(i) = (rtem(i) - savr(i))*fac
220    continue
end do
call adda (neq,tn,y,jj,iwk(ipian),iwk(ipjan),rtem)
do k = kmin,kmax
i = iwk(ibjan+k)
wk(iba+k) = rtem(i)
225    continue
end do
230   continue
end do
jmin = jmax + 1
240 continue
end do
ires = 1
call res (neq,tn,y,s,savr,ires)
nfe = nfe + 1
if (ires .gt. 1) go to 600
290 nlu = nlu + 1
ierpj = 0
do i = 1,n
295 rtem(i) = 0.0d0
end do
call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),&
wk(ipa),rtem,rtem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
if (iys .eq. 0) return
imul = (iys - 1)/n
ierpj = -2
if (imul .eq. 8) ierpj = 1
if (imul .eq. 10) ierpj = -1
return
600 ierpj = ires
return
end subroutine dprjis
