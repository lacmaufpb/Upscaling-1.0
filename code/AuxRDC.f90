
Subroutine AuxRDC

use vars_RDC

implicit none

integer varauxint


varauxint=int(varaux)

select case(TamNum)
case(1)
	write(87,'(A14,I1)')   textoaux,varauxint
case(2)
	write(87,'(A14,I2)')   textoaux,varauxint
case(3)
	write(87,'(A14,I3)')   textoaux,varauxint
case(4)
	write(87,'(A14,I4)')   textoaux,varauxint
case(5)
	write(87,'(A14,I5)')   textoaux,varauxint
case(6)
	write(87,'(A14,I6)')   textoaux,varauxint
case(7)
	write(87,'(A14,I7)')   textoaux,varauxint
case(8)
	write(87,'(A14,I8)')   textoaux,varauxint
case(9)
	write(87,'(A14,F9.7)')   textoaux,varaux
case(10)
	write(87,'(A14,F10.7)')   textoaux,varaux
case(11)
	write(87,'(A14,F11.7)')   textoaux,varaux
case(12)
	write(87,'(A14,F12.7)')   textoaux,varaux
case(13)
	write(87,'(A14,F13.7)')   textoaux,varaux
case(14)
	write(87,'(A14,F14.7)')   textoaux,varaux
case(15)
	write(87,'(A14,F15.7)')   textoaux,varaux
case(16)
	write(87,'(A14,F16.7)')   textoaux,varaux
case(17)
	write(87,'(A14,F17.7)')   textoaux,varaux
case(18)
	write(87,'(A14,F18.7)')   textoaux,varaux
case(19)
	write(87,'(A14,F19.7)')   textoaux,varaux
end select




return

end