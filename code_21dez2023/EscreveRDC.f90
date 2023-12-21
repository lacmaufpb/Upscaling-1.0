
Subroutine EscreveRDC

use Vars_RDC

implicit none

integer PosExt


PosExt=index(nomeRST,'.rst')
nomeRDC=nomeRST(1:PosExt-1)//'.rdc'
open(87,file=nomeRDC)
write(87,'(A14,A17)') 'file format : ','IDRISI Raster A.1'
write(87,'(A14,A1)') 'file title  : ',''
if (tipodado==1) then
	write(87,'(A14,A7)') 'data type   : ','integer'  
else
	if (tipodado==2) then
		write(87,'(A14,A4)') 'data type   : ','real'
	end if
end if
write(87,'(A14,A6)')    'file type   : ','binary'
varaux=ncol3
num=1
textoaux='columns     : '
call TamNumero
call AuxRDC
varaux=nlin3
num=1
textoaux='rows        : '
call TamNumero
call AuxRDC
write(87,'(A14,A7)')    'ref. system : ',sistemaref
if (metrordc==1) then
    write(87,'(A14,A1)')    'ref. units  : ','m'
else
    write(87,'(A14,A3)')    'ref. units  : ','deg'
end if
write(87,'(A14,F9.7)')  'unit dist.  : ',1.0
varaux=xmin3
num=2
textoaux='min. X      : '
call TamNumero
call AuxRDC
varaux=xmax3
num=2
textoaux='max. X      : '
call TamNumero
call AuxRDC
varaux=ymin3
num=2
textoaux='min. Y      : '
call TamNumero
call AuxRDC
varaux=ymax3
num=2
textoaux='max. Y      : '
call TamNumero
call AuxRDC
write(87,'(A14,A7)')    "pos'n error : ",'unknown'
varaux=dx3
num=2
textoaux='resolution  : '
call TamNumero
call AuxRDC
varaux=VarMin
num=2 
textoaux='min. value  : '
call TamNumero
call AuxRDC
varaux=VarMax
num=2 
textoaux='max. value  : '
call TamNumero
call AuxRDC
varaux=VarMin
num=2 
textoaux='display min : '
call TamNumero
call AuxRDC
varaux=VarMax
num=2 
textoaux='display max : '
call TamNumero
call AuxRDC
write(87,'(A14,A11)')   'value units : ','unspecified'
write(87,'(A14,A7)')    'value error : ','unknown'
write(87,'(A14,I1)') 'flag value  : ',0
write(87,'(A14,A4)')   "flag def'n  : ",'none'
write(87,'(A14,I1)') 'legend cats : ', 0
write(87,'(A14,A100)') 'lineage     : ','This file was created automatically by an ARP FORTRAN routine'

close(87)



return


end
