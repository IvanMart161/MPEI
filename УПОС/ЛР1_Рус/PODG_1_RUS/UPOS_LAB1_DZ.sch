<Qucs Schematic 25.2.0>
<Properties>
  <View=-40,-297,3662,1465,1.24854,1356,342>
  <Grid=10,10,1>
  <DataSet=UPOS_LAB1_DZ.dat>
  <DataDisplay=UPOS_LAB1_DZ.dpl>
  <OpenDisplay=0>
  <Script=UPOS_LAB1_DZ.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Название>
  <FrameText1=Чертил:>
  <FrameText2=Дата:>
  <FrameText3=Версия:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <GND * 1 660 370 0 0 0 0>
  <GND * 1 1180 360 0 0 0 0>
  <src_eqndef B1 1 660 340 18 -26 0 1 "A * cos(2*pi*f0*time + m1 * sin(2*pi*Fm1*time))" 1 "" 0 "" 0 "" 0 "" 0>
  <src_eqndef B2 1 1180 330 18 -26 0 1 "A * cos(2*pi*f0*time + m2* sin(2*pi*Fm2*time))" 1 "" 0 "" 0 "" 0 "" 0>
  <.FFT FFT1 1 660 70 0 56 0 0 "1.00023 GHz" 1 "100 k" 1 "none" 1 "2" 0 "0" 0 "yes" 0>
  <Eqn Eqn2 1 1110 100 -37 17 0 0 "m2=2" 1 "delta_f2=500e3" 1 "Fm2=delta_f2 / m2" 1 "yes" 0>
  <Eqn Eqn1 1 910 100 -37 17 0 0 "A=1" 1 "f0=0.83e9" 1 "m1=10" 1 "delta_f1=500e3" 1 "Fm1=delta_f1 / m1" 1 "yes" 0>
  <GND * 1 1860 360 0 0 0 0>
  <Vac V1 1 1860 330 18 -26 0 1 "1 V" 1 "f0" 1 "0" 0 "0" 0 "0" 0 "0" 0>
</Components>
<Wires>
  <660 310 660 310 "Vout1" 670 260 0 "">
  <1180 300 1180 300 "Vout2" 1190 240 0 "">
  <1860 300 1860 300 "Vout3" 1890 270 0 "">
</Wires>
<Diagrams>
  <Rect 620 597 345 167 3 #c0c0c0 1 00 0 8.27e+08 1e+06 8.33e+08 1 -0.0111312 0.05 0.122684 1 -1 1 1 315 0 225 1 0 0 "f, Гц" "U, В" "">
	<"ngspice/ac.v(vout1)" #0000ff 2 3 0 0 0>
  </Rect>
  <Rect 1160 597 345 167 3 #c0c0c0 1 00 0 8.27e+08 1e+06 8.33e+08 1 -0.0162635 0.1 0.2 1 -1 1 1 315 0 225 1 0 0 "f, Гц" "U, В" "">
	<"ngspice/ac.v(vout2)" #ff0000 2 3 0 0 0>
  </Rect>
  <Rect 1690 587 345 167 3 #c0c0c0 1 00 0 8.27e+08 1e+06 8.33e+08 1 -0.0162635 0.1 0.2 1 -1 1 1 315 0 225 1 0 0 "f, Гц" "U, В" "">
	<"ngspice/ac.v(vout3)" #ff0000 2 3 0 0 0>
	  <Mkr 8.30066e+08 206 -213 3 1 0>
	  <Mkr 8.30005e+08 -27 -197 3 1 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
