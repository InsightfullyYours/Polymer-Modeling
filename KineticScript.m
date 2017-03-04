structPS0= KineticFits2(OldLampBis,WholeAreaPS0');
title('PS0')
structPSbis0p5 = KineticFits2(OldLampBis([1:7,9]),(WholeAreaPSbis0p5(:,[1:7,9]))');
title('PSbis0p5')
structPSbis1p2 = KineticFits2(OldLampBis,WholeAreaPSbis1p2');
title('PSbis1p2')
structPSbis4 = KineticFits2(OldLampBis([1:4,6:end-1]),(WholeAreaPSbis4(:,[1:4,6:end]))');
title('PSbis4')
structPSbis16 = KineticFits2(OldLampBis,WholeAreaPSbis16');
title('PSbis16')

structPSmono1 = KineticFits2(OldLampMono,WholeAreaPSmono1');
title('PSmono1')
structPSmono2p4 = KineticFits2(OldLampMono,WholeAreaPSmono2p4');
title('PSmono2p4')
structPSmono8 = KineticFits2(OldLampMono,WholeAreaPSmono8');
title('PSmono8')
structPSmono32 = KineticFits2(OldLampMono,WholeAreaPSmono32');
title('PSmono32')

structPSRIDmono = KineticFits2(NewLampMono,WholeAreaPSRIDmono');
title('PSRIDmono')
structPSRIDbis = KineticFits2(NewLampMono,WholeAreaPSRIDbis');
title('PSRIDbis')
structPSRIDctr = KineticFits2(NewLampMono,WholeAreaPSRIDctr');
title('PSRIDctr')
structPSUVmono = KineticFits2(NewLampMono,(WholeAreaPSUVmono)');
title('PSUVmono')
structPSUVbis = KineticFits2(NewLampMono,(WholeAreaPSUVbis)');
title('PSUVbis')
structPSUVctr  = KineticFits2(NewLampMono,(WholeAreaPSUVctr)');
title('PSUVctr')