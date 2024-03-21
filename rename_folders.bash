for i in output/Sprtn-SH2_Gd/*
do
  mv -- "$i" "${i/SH2_Gd-/SH2_Gd-SPRT}"
  #echo -- "$i" "${i/SH2_Se-/SH2_Se-SPRT}"
done