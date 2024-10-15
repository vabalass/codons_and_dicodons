# codons_and_dicodons
Universtity lab works

##aprašykite, kaip skaičiavote atstumo funkciją;

1. Paimi dviejų virusų kodonų arba dikodonų dažnių vektorius.
2. Apskaičiuoji skirtumus tarp tų dviejų vektorių elementų ir pakeli juos kvadratu.
3. Gautas reikšmes susumavus, gaunam atstumą tarp dviejų virusų.

##kokie medžiai gavosi su kodonais ir dikodonais

kodonų atstumų matrica:
bacterial1 0.000 0.019 0.016 0.010 0.023 0.028 0.032 0.031
bacterial2 0.019 0.000 0.024 0.017 0.021 0.016 0.025 0.027
bacterial3 0.016 0.024 0.000 0.021 0.027 0.029 0.039 0.031
bacterial4 0.010 0.017 0.021 0.000 0.022 0.025 0.025 0.036
mamalian1  0.023 0.021 0.027 0.022 0.000 0.023 0.032 0.030
mamalian2  0.028 0.016 0.029 0.025 0.023 0.000 0.037 0.015
mamalian3  0.032 0.025 0.039 0.025 0.032 0.037 0.000 0.049
mamalian4  0.031 0.027 0.031 0.036 0.030 0.015 0.049 0.000

dikodonų atstumų matrica:
bacterial1 0.000 0.003 0.003 0.002 0.003 0.004 0.004 0.004
bacterial2 0.003 0.000 0.003 0.002 0.002 0.002 0.003 0.003
bacterial3 0.003 0.003 0.000 0.003 0.003 0.004 0.004 0.004
bacterial4 0.002 0.002 0.003 0.000 0.003 0.003 0.003 0.004
mamalian1  0.003 0.002 0.003 0.003 0.000 0.003 0.003 0.003
mamalian2  0.004 0.002 0.004 0.003 0.003 0.000 0.004 0.002
mamalian3  0.004 0.003 0.004 0.003 0.003 0.004 0.000 0.005
mamalian4  0.004 0.003 0.004 0.004 0.003 0.002 0.005 0.000

##ar skiriasi kodonų ir dikodonų dažnis tarp žinduolių ir bakterijų virusų, kaip klasterizuojasi virusai. Gal kažkuris virusas labai išsiskyrė? Kokie kodonai/dikodonai labiausiai varijuoja?

Virusai aiškiai skiriasi pagal tipą (žinduoliai vs. bakterijos). 
Mamalian3 virusas turi didžiausią skirtumo reikšmę (0.049), kas rodo, jog jis labiausiai skiriasi nuo kitų virusų. Mamalian3 ir Mamalian4 virusai išsiskyrė labiausiai iš žinduolių, ypač pagal kodonų dažnius.
Kodonų dažniai labiausiai skiriasi tarp virusų grupių, o dikodonų vienodesni, bet vis tiek naudingi.
