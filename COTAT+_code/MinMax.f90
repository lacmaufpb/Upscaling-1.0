
Subroutine MinMax

use vars_RDC

implicit none


VarMax=-10000000.0
VarMin=10000000.0

if (tipoMM==2) then
	do col3=1,ncol3
		do lin3=1,nlin3
				if (VarMM2(lin3,col3)>VarMax) then
					VarMax=VarMM2(lin3,col3)
				end if
				if (VarMM2(lin3,col3)<VarMin) then
					VarMin=VarMM2(lin3,col3)
				end if
		end do
	end do
else
	do col3=1,ncol3
		do lin3=1,nlin3
				if (VarMM3(lin3,col3,i3)>VarMax) then
					VarMax=VarMM3(lin3,col3,i3)
				end if
					if (VarMM3(lin3,col3,i3)<VarMin) then
						VarMin=VarMM3(lin3,col3,i3)
					end if
		end do
	end do
end if


return

end