
Subroutine EscreveRDC

use Vars_RDC

implicit none

integer PosExt


!identifica posicao da extensao no arquivo .rst
PosExt=index(nomeRST,'.rst')
!define nome do arquivo rdc correspondente igual ao nome do arquivo rst
nomeRDC=nomeRST(1:PosExt-1)//'.rdc'
open(87,file=nomeRDC)
!escreve linha com formato do arquivo
write(87,'(A14,A17)') 'file format : ','IDRISI Raster A.1'
!escreve linha com titulo
write(87,'(A14,A1)') 'file title  : ',''
!escreve linha com tipo de dado
if (tipodado==1) then
	write(87,'(A14,A7)') 'data type   : ','integer'  !estava A14,A10 por que?
else
	if (tipodado==2) then
		write(87,'(A14,A4)') 'data type   : ','real'
	end if
end if
!escreve linha com tipo de arquivo
write(87,'(A14,A6)')    'file type   : ','binary'
!escreve linha com numero de colunas
varaux=ncol3
num=1
textoaux='columns     : '
call TamNumero
call AuxRDC
!escreve linha com numero de linhas
varaux=nlin3
num=1
textoaux='rows        : '
call TamNumero
call AuxRDC
!escreve linha com sistema de referencia
write(87,'(A14,A7)')    'ref. system : ',sistemaref
!write(87,'(A14,A5)')    'ref. system : ','plane'
!write(87,'(A14,A7)')    'ref. system : ','latlong'
!escreve linha com unidade da referencia
if (metrordc==1) then
    !escreve linha com unidade da referencia
    write(87,'(A14,A1)')    'ref. units  : ','m'
else
    write(87,'(A14,A3)')    'ref. units  : ','deg'
end if
!write(87,'(A14,A3)')    'ref. units  : ','deg'
!write(87,'(A14,A1)')    'ref. units  : ','m'
!escreve linha com distancia unitaria da referencia
write(87,'(A14,F9.7)')  'unit dist.  : ',1.0
!escreve linha com coord xmin
varaux=xmin3
num=2
textoaux='min. X      : '
call TamNumero
call AuxRDC
!escreve linha com coord xmax
varaux=xmax3
num=2
textoaux='max. X      : '
call TamNumero
call AuxRDC
!escreve linha com coord ymin
varaux=ymin3
num=2
textoaux='min. Y      : '
call TamNumero
call AuxRDC
!escreve linha com coord ymax
varaux=ymax3
num=2
textoaux='max. Y      : '
call TamNumero
call AuxRDC
!escreve linha com erro de posicao
write(87,'(A14,A7)')    "pos'n error : ",'unknown'
!escreve linha com resolucao
varaux=dx3
num=2
textoaux='resolution  : '
call TamNumero
call AuxRDC
!escreve linha com valor minimo dos dados
varaux=VarMin
num=2 !integer=1; real=2
textoaux='min. value  : '
call TamNumero
call AuxRDC
!escreve linha com valor maximo dos dados
varaux=VarMax
num=2 !integer=1; real=2
textoaux='max. value  : '
call TamNumero
call AuxRDC
!escreve linha com valor minimo para exibicao
varaux=VarMin
num=2 !integer=1; real=2
textoaux='display min : '
call TamNumero
call AuxRDC
!escreve linha com valor maximo para exibicao
varaux=VarMax
num=2 !integer=1; real=2
textoaux='display max : '
call TamNumero
call AuxRDC
!escreve linha com unidade dos dados
write(87,'(A14,A11)')   'value units : ','unspecified'
!escreve linha com valor do erro dos dados
write(87,'(A14,A7)')    'value error : ','unknown'
!escreve linha com sinalizador
write(87,'(A14,I1)') 'flag value  : ',0
!escreve linha com definicao do sinalizador
write(87,'(A14,A4)')   "flag def'n  : ",'none'
!escreve linha com numero de categorias da legenda
write(87,'(A14,I1)') 'legend cats : ', 0
!escreve linha sobre criacao da imagem
write(87,'(A14,A100)') 'lineage     : ','This file was created automatically by an ARP FORTRAN routine'

close(87)



return


end
