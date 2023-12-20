
Subroutine TamNumero

use vars_RDC

implicit none

integer negativo,nzero,pp
real varaux2,limsup


if (varaux<0) then
	negativo=1
else
	negativo=0
end if
varaux2=abs(varaux)

do pp=0,10
	limsup=10.0**pp
	if (varaux2<limsup) then
		nzero=pp
		exit
	end if
end do

if (num==1) then !eh inteiro
	if (nzero==0) then
		tamnum=1+negativo
	else
		tamnum=nzero+negativo
	end if
else !eh real
	if (nzero==0) then
		tamnum=8+1+negativo
	else
		tamnum=8+nzero+negativo
	end if
end if






return

end