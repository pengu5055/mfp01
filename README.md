# DN 1. Izračun Airyjevih funkcij (Numerična aproksimacija)

Prva naloga mafijskega praktikuma se mi zdi res perfektna v svojem namenu. Cilj naloge je spoznavanje numeričnega reševanja in dejstva, da imajo številke v računalniku intrinzično omejeno natančnost. To je zelo pomembno dejstvo, ki se ga moramo zavedati, ko se lotimo numeričnega reševanja kakšnega problema.


## Navodila
Naloga želi, da z uporabi Maclaurinove vrste in asimptotskega razvoja poiščeš čim učinkovitejši postopek za izračun vrednosti Ariyjevih funkcij na realni osi. Željena za natančnost je, da spravimo absolutno napako pod $10^{-10}$. Radi bi storili to tudi z relativno napako. Poglej če je to mogoče.

## Napotki
1. Ali so referenče funkcije res popolnoma pravilne?
2. Za boljšo natančnost si lahko pomagaš z uporabo `decimal` ali pa meni še ljubše `mpmath` knjižnice.
3. Med drugim je cilj mafijskega praktikuma tudi ta, da se naučiš delati res hot grafe. Glej da bodo osi označene, da bo legenda, da bodo barvne kombinacije dobre. Hidden weapon za lepe barve je `cmasher` knjižnica, ki vsebuje dodatne colormap-e za `matplotlib`.