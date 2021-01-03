"""
Test all geodesic functions here.

This file compares results of geodesic functions with the
Black Hole Perturbation Toolkit written in Mathematica.
"""
import pytest
import numpy as np
from mpmath import mpf, mp, almosteq
from functions import calc_constants

digits = 100  # accuracy requested
eps = 10 ** (-digits)  # number of digits which may be different
mp.dps = digits  # set precision


# -----------------------------------------------------------------------------
#   Tests of calc_constants function
# -----------------------------------------------------------------------------


def test_constants_sc_polar_circular():
    aa = 0
    slr = 6
    ecc = 0
    x = 0
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.94280904158206336586779248280646538571311458358463204878"
        + "4453158660488318974738025900258356218427715156676"
    )
    Lz = mpf("0")
    Q = mpf("12")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_circular():
    aa = 0
    slr = 6
    ecc = 0
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.94280904158206336586779248280646538571311458358463204878"
        + "4453158660488318974738025900258356218427715156676"
    )
    Lz = mpf(
        "1.73205080756887729352744634150587236694280525381038062"
        + "80558069794519330169088000370811461867572485756756261414"
    )
    Q = mpf("9")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_circular_equatorial():
    aa = 0
    slr = 6
    ecc = 0
    x = 1
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.9428090415820633658677924828064653857131145835846320487"
        + "844531586604883189747380259002583562184277151566758975"
    )
    Lz = mpf(
        "3.4641016151377545870548926830117447338856105076207612561"
        + "11613958903866033817600074162292373514497151351252283"
    )
    Q = mpf("0")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_polar_circular():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 0
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.93949949445667018404632400714789613593552722936335835298"
        + "18075269414504771831119358463712091768276399621636372"
    )
    Lz = mpf("0")
    Q = mpf(
        "11.49069917358251688551533489363432042697268698929664496570"
        + "628785810940303968313887646713785433448129671465657"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_circular():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.92926208518755385940023990416190759483286146511474553693"
        + "4033239507956850773101835323740984327261143567860022"
    )
    Lz = mpf(
        "1.51864386037595847509551522883727088566379643896240197669"
        + "357808405028886043522766119946171907171961531232087"
    )
    Q = mpf(
        "7.001744250020255211247949830215200508121452736181614498048"
        + "1905103921213089404453120895980950824845440005903"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_circular_equatorial():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 1
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.92259962627962169982949337843300709772200871590220343780"
        + "88735453946811240350585108621787789067290006977623192"
    )
    Lz = mpf(
        "2.79427836148321036281043368969640997944808168260479731853"
        + "0053829729254630064949889877752074500576020878348642"
    )
    Q = mpf("0")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_polar_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 0
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.94320311016419542510334707036486778180528031689323908838"
        + "78061618241015265193332910795870838688378743323652801967135"
    )
    Lz = mpf("0")
    Q = mpf(
        "12.0401337792642140468227424749163879598662207357"
        + "8595317725752508361204013377926421404682274247491"
        + "63879598662207357864"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_inclined_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.943203110164195425103347070364867781805280316893239088"
        + "38780616182410152651933329107958708386883787433236528019"
        + "67135"
    )
    Lz = mpf(
        "1.734944795898720666757003868204620133806227380741254"
        + "89151777193947131753531745174945379491474382700412679"
        + "59746716139"
    )
    Q = mpf(
        "9.030100334448160535117056856187290969899665551839"
        + "46488294314381270903010033444816053511705685618729"
        + "09698996655518395"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_equatorial_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 1
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.943203110164195425103347070364867781805280316893239088"
        + "38780616182410152651933329107958708386883787433236528019"
        + "67135"
    )
    Lz = mpf(
        "3.469889591797441333514007736409240267612454761482509"
        + "78303554387894263507063490349890758982948765400825359"
        + "19493432278"
    )
    Q = mpf("0")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_equatorial_medium_ecc():
    aa = 0
    slr = 6
    ecc = 0.5
    x = 1
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.953462589245592315446775921527215998613883506983185440"
        + "79614160325455084136463699219944117175744947073087697070"
        + "36711"
    )
    Lz = mpf(
        "3.618136134933163471761744803640749109738642049733840"
        + "28770038052303616506466441146913692007177631707371674"
        + "09721994845"
    )
    Q = mpf("0")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_polar_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 0
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.9999833356941505011042626308284425667213449655947521"
        + "046595466454286167504254217295586344501915823285618514122"
    )
    Lz = mpf("0")
    Q = mpf(
        "17.9982002699640049493250922374017210149071143630992618"
        + "8930737870866602694063607109607074840551319070849510403"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.999983335694150501104262630828442566721344965594752104"
        + "6595466454286167504254217295586344501915823285618514122"
    )
    Lz = mpf(
        "2.121214290799258507892773741617459872076805628175392680"
        + "3202618929948843840934521693383074613688504765958787920"
    )
    Q = mpf(
        "13.49865020247300371199381917805129076118033577232444641"
        + "6980534031499520205477053322053061304134893031371328024"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.999983335694150501104262630828442566721344965594752104"
        + "6595466454286167504254217295586344501915823285618514122"
    )
    Lz = mpf(
        "2.121214290799258507892773741617459872076805628175392680"
        + "32026189299488438409345216933830746136885047659587879199"
        + "7635786628786523303"
    )
    Q = mpf(
        "13.49865020247300371199381917805129076118033577232444641"
        + "6980534031499520205477053322053061304134893031371328024"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_sc_equatorial_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 1
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.999983335694150501104262630828442566721344965594752104"
        + "6595466454286167504254217295586344501915823285618514122"
    )
    Lz = mpf(
        "4.242428581598517015785547483234919744153611256350785360"
        + "640523785989768768186904338676614922737700953191757584"
    )
    Q = mpf("0")
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_polar_high_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.9999
    x = 0
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.999983335339928721196580291224660022139747002421683704"
        + "3269892342102857819986614440069285559744480834240033140"
    )
    Lz = mpf("0")
    Q = mpf(
        "15.447760940128635640183107341523789504572370741756299777"
        + "04359223520794709622630983426009037882933633201136405"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_high_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.9999
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.999983334762009390624963061217615163403107253700220379"
        + "80574894182978682677699041122633009873901778668809666999"
        + "80678119314950"
    )
    Lz = mpf(
        "1.679773151701191377352370006632168332785658441313245996"
        + "64649971146314487918503303276461332270840559456097663872111982"
    )
    Q = mpf(
        "8.4649337716238986206511841687108769598324549247957851881"
        + "295551604892270945362640145013991357228739478154330788373464"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)


def test_constants_kerr_even_higher_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.999999
    x = 0.5
    En_ch, Lz_ch, Q_ch = calc_constants(aa, slr, ecc, x, digits)
    En = mpf(
        "0.9999998333334762089978428690357821465253272401340593719"
        + "46059752529066513167011257218869930744732693058"
    )
    Lz = mpf(
        "1.6798101590762567517524755694236949833471071190540190363"
        + "7527258462637697164838885327026795514787789551830938510"
    )
    Q = mpf(
        "8.46528671410720657183089558985044315752321378660758373670"
        + "77178699937459442334893278595278583665315183011560589"
    )
    assert almosteq(En_ch, En, eps)
    assert almosteq(Lz_ch, Lz, eps)
    assert almosteq(Q_ch, Q, eps)
