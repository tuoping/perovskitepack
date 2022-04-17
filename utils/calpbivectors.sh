for d in 1-2ns-T{170..310..20} 1-2ns-T350 0-2ns-T195 0-2ns-T200 0-2ns-T205 0-2ns-T215 0-2ns-T220 0-2ns-T225 0-2ns-T{295..305..5} 0-3ns-T{255..265..5} 0-3ns-T{275..285..5} 0-3ns-T320 0-3ns-T340;
do
    echo $d;
    cd $d;
    # tail -6153 md4.dump > 4020000.lammpstrj;
    # python ../perovskitelattpack.py *.lammpstrj;
    python ../readlogs.py
    # python ../drawhist.py molecule_longaxis.dat
    # python ../drawhist.py molecule_polaraxis.dat
    cd ../;
done
