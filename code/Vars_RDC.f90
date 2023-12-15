
Module Vars_RDC


integer i,narq
character(len=40),allocatable :: nomearq(:)
character(len=1) :: iauxchar1
character(len=2) :: iauxchar2
character(len=3) :: ichar
integer lin,col,nlin,ncol
real,allocatable :: Var2D(:,:)

real,allocatable :: xx(:),yy(:),P2lc(:,:)
integer arg,KC2


character(len=20) cabecalho

character(len=14) :: texto10,texto20    !texto contido no cabecalho do arquivo do ArcView
character(len=14) :: texto30,texto40    !texto contido no cabecalho do arquivo do ArcView
character(len=14) :: texto50           !texto contido no cabecalho do arquivo do ArcView
character(len=14) :: texto60          !texto contido no cabecalho do arquivo do ArcView
integer :: pixnul    !num colunas,linhas,valor pixel nulo
real :: coordxie,coordyie

real parax,paray,delta

real Xmin3,Xmax3,Ymin3,Ymax3,dx3
integer lin3,col3,nlin3,ncol3
real VarMin,VarMax
integer tipodado

integer tamnum,num

real varaux
character(len=14):: textoaux

real,allocatable :: VarMM2(:,:),VarMM3(:,:,:)
integer tipoMM

integer TipoAnaliseFaixa

character(len=33) :: nomeRST,nomeRDC
integer*2,allocatable :: cell(:,:)
integer i3

integer lin2,col2,cellaux

integer inivalueraster

integer :: linstac,linst
integer :: gerastack

character (len=7) :: sistemaref
character (len=3) :: unidaderef3

integer :: metrordc

end