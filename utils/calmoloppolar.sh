for d in annealT{52..54} annealT{56..59} annealT{61..63} annealT{66..69} annealT{71..74} annealT{76..79} annealT{81..84} annealT{86..89} annealT{91..94} annealT{96..99} annealT{101..104} annealT{136..139} annealT{141..144} annealT{147..149} annealT{152..154};
do
    echo $d;
    cd $d;
    python ../annealT50/getorderparam-mol-CHa.py md4.dump 1 > molrotop-polar.out;
    cd ../;
done
