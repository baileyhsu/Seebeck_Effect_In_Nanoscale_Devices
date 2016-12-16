      FUNCTION FERMIR(E)
        IMPLICIT NONE
        REAL*8 FERMIR,E
        real*8 VL,VR,KT1,KT2
        common /VL/VL/VR/VR/KT1/KT1/KT2/KT2
 
        IF(KT2.EQ.0.0d0) then 
        IF(E.LT.VR) then 
        FERMIR=1.0d0 
        else 
        FERMIR=0.0d0
        end if
        else if (E.LT.VR) then 
        FERMIR=1.0d0/(1.0d0+dexp((E-VR)/KT2))
        else 
        FERMIR=dexp(-(E-VR)/KT2)/(1.0d0+dexp(-(E-VR)/KT2))
        ENDIF 

        RETURN
        
        END
