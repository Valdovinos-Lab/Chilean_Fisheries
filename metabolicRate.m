function mbr = metabolicRate(avgl,a,b)

wcc = a.*avgl.^b*0.2*0.53*1e6;
mbr = 0.88*(6.4e-5./wcc).^0.11;