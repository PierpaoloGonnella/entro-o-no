import math
from datetime import datetime, timedelta, timezone

import numpy as np
import pytest


def test_cf_ittc_known_value(tar):
    re = 1.0 * tar.HULL["L"] / tar.NU
    assert tar.cf_ittc(1.0) == pytest.approx(0.075 / (math.log10(re) - 2) ** 2)


def test_cf_ittc_decreases_with_speed(tar):
    assert tar.cf_ittc(0.5) > tar.cf_ittc(1.0) > tar.cf_ittc(2.0)


def test_cf_ittc_clamped_at_low_reynolds(tar):
    assert tar.cf_ittc(0.0) == pytest.approx(0.075 / 4)


def test_r_hull_zero_and_monotonic(tar):
    assert tar.r_hull(0) == 0
    vals = [tar.r_hull(v) for v in np.linspace(0.1, 3, 30)]
    assert all(b > a for a, b in zip(vals, vals[1:]))


def test_v_acqua_calm_equals_cruise(tar):
    assert tar.v_acqua(1.13, 4e-4, 0.0) == pytest.approx(1.13, rel=1e-6)


def test_v_acqua_headwind_slower_tailwind_faster(tar):
    v0, b = 1.13, 4.18e-4
    assert tar.v_acqua(v0, b, 5.0) < v0 < tar.v_acqua(v0, b, -5.0)


def test_v_acqua_stronger_headwind_slower(tar):
    v0, b = 1.13, 4.18e-4
    assert tar.v_acqua(v0, b, 10.0) < tar.v_acqua(v0, b, 5.0)


@pytest.mark.parametrize("da,rotta,expected", [
    (0, 0, 1.0), (180, 0, -1.0), (90, 0, 0.0), (45, 90, math.cos(math.radians(45)))])
def test_componente_vento(tar, da, rotta, expected):
    assert tar.componente_vento(10, da, rotta) == pytest.approx(10 * tar.KN * expected, abs=1e-9)


def test_haversine_one_degree_latitude(tar):
    assert tar._haversine((0, 0), (1, 0)) == pytest.approx(111195, rel=1e-3)


@pytest.mark.parametrize("b,exp", [((1, 0), 0), ((0, 1), 90), ((-1, 0), 180), ((0, -1), 270)])
def test_rilevamento_cardinal(tar, b, exp):
    assert tar._rilevamento((0, 0), b) == pytest.approx(exp, abs=1e-6)


def test_prepara(tar):
    righe = [
        dict(traccia="B", tipo="piedi", distanza_m="200", durata_s="100", rotta_deg="0",
             vento_kn="10", vento_da_deg="0", onda_m="0.3"),
        dict(traccia="A", tipo="deriva", distanza_m="60", durata_s="120", rotta_deg="180",
             vento_kn="10", vento_da_deg="0", onda_m=""),
    ]
    d, tracce = tar.prepara(righe)
    assert tracce == ["A", "B"]
    assert list(d["itr"]) == [1, 0]
    assert d["vg"] == pytest.approx([2.0, 0.5])
    assert d["Uh"][0] == pytest.approx(10 * tar.KN)
    assert d["Uh"][1] == pytest.approx(-10 * tar.KN)
    assert np.isnan(d["onda"][1])
    assert d["peso"][0] == pytest.approx(10.0)


def test_ora_vicina(tar):
    t0 = datetime(2025, 7, 1, 10, tzinfo=timezone.utc)
    times = [(t0 + timedelta(hours=h)).timestamp() for h in range(5)]
    assert tar._ora_vicina(times, t0 + timedelta(hours=2, minutes=20)) == 2
    assert tar._ora_vicina(times, t0 + timedelta(hours=2, minutes=40)) == 3


def _track(n, brg, step_m=5.0, start=(40.5, 17.0), t0=None):
    t0 = t0 or datetime(2025, 7, 1, 10, tzinfo=timezone.utc)
    pts, (lat, lon) = [], start
    for i in range(n):
        pts.append((t0 + timedelta(seconds=5 * i), lat, lon))
        lat += step_m * math.cos(math.radians(brg)) / 111195
        lon += step_m * math.sin(math.radians(brg)) / (111195 * math.cos(math.radians(lat)))
    return pts


def test_estrai_tratte_too_short(tar):
    assert tar.estrai_tratte(_track(5, 90)) == []


def test_estrai_tratte_single_straight_leg(tar):
    tr = tar.estrai_tratte(_track(60, 90))
    assert len(tr) == 1
    assert tr[0]["rotta_deg"] == pytest.approx(90, abs=1)
    assert tr[0]["durata_s"] == pytest.approx(295)
    assert tr[0]["distanza_m"] == pytest.approx(295, rel=0.02)


def test_estrai_tratte_splits_on_course_change(tar):
    a = _track(40, 0)
    last = a[-1]
    b = _track(40, 90, start=(last[1], last[2]), t0=last[0] + timedelta(seconds=5))
    tr = tar.estrai_tratte(a + b)
    assert len(tr) == 2
    rotte = sorted(t["rotta_deg"] for t in tr)
    assert rotte[0] == pytest.approx(0, abs=5)
    assert rotte[1] == pytest.approx(90, abs=5)


def test_leggi_gpx_sorts_and_skips_untimed(tar, tmp_path):
    gpx = tmp_path / "t.gpx"
    gpx.write_text(
        '<?xml version="1.0"?>\n'
        '<gpx xmlns="http://www.topografix.com/GPX/1/1"><trk><trkseg>\n'
        '<trkpt lat="40.6" lon="17.1"><time>2025-07-01T10:00:10Z</time></trkpt>\n'
        '<trkpt lat="40.5" lon="17.0"><time>2025-07-01T10:00:00Z</time></trkpt>\n'
        '<trkpt lat="40.7" lon="17.2"></trkpt>\n'
        '</trkseg></trk></gpx>', encoding="utf-8")
    pts = tar.leggi_gpx(str(gpx))
    assert len(pts) == 2
    assert pts[0][1:] == (40.5, 17.0)
    assert pts[0][0].tzinfo is not None


def test_previsione_drift_plus_current(tar):
    d = {"vg": np.zeros(2), "rotta": np.array([0.0, 90.0]), "itr": np.array([0, 0]),
         "tipo": np.array(["deriva", "deriva"]), "Uh": np.array([-5.0, -5.0])}
    alpha = 0.06
    par = np.concatenate([np.log([1.13, 4e-4, 2e-4, alpha]), [0.1, 0.0]])  # 0.1 m/s east
    out = tar.previsione(par, d, 1)
    assert out[0] == pytest.approx(alpha * 5.0)
    assert out[1] == pytest.approx(alpha * 5.0 + 0.1)


def test_previsione_paddling_calm_equals_cruise(tar):
    d = {"vg": np.zeros(1), "rotta": np.array([0.0]), "itr": np.array([0]),
         "tipo": np.array(["piedi"]), "Uh": np.array([0.0])}
    par = np.concatenate([np.log([1.13, 4e-4, 2e-4, 0.06]), [0.0, 0.0]])
    assert tar.previsione(par, d, 1)[0] == pytest.approx(1.13, rel=1e-6)
