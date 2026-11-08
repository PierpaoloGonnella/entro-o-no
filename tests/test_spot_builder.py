import json
import math

import pytest


def angdiff(a, b):
    return abs(((a - b + 180) % 360) - 180)


def straight_coast(brg):
    def f(a, r):
        return -10.0 if math.cos(math.radians(a - brg)) * r > 0.3 else 20.0
    return f


def grid(sb, f):
    az = [i * 360.0 / sb.N_AZ for i in range(sb.N_AZ)]
    return [[f(a, r) for r in sb.RAGGI_KM] for a in az], az


# ---------------- sposta ----------------

@pytest.mark.parametrize("brg", [0, 90, 180, 270, 37.5])
def test_sposta_distance_and_bearing_preserved(sb, tar, brg):
    lat2, lon2 = sb.sposta(40.5, 17.0, brg, 5.0)
    d = tar._haversine((40.5, 17.0), (lat2, lon2)) / 1000
    assert d == pytest.approx(5.0, rel=1e-3)
    assert angdiff(tar._rilevamento((40.5, 17.0), (lat2, lon2)), brg) < 0.1


def test_sposta_north_increases_latitude(sb):
    lat2, lon2 = sb.sposta(0.0, 0.0, 0, 111.195)
    assert lat2 == pytest.approx(1.0, abs=1e-3)
    assert lon2 == pytest.approx(0.0, abs=1e-9)


def test_sposta_wraps_longitude(sb):
    _, lon2 = sb.sposta(0.0, 179.99, 90, 10)
    assert -180 <= lon2 < -179


# ---------------- media_circolare ----------------

def test_media_circolare_across_north(sb):
    ang, conc = sb.media_circolare([(350, 1), (10, 1)])
    assert angdiff(ang, 0) < 1e-6
    assert conc == pytest.approx(math.cos(math.radians(10)))


def test_media_circolare_opposite_directions_zero_concentration(sb):
    _, conc = sb.media_circolare([(90, 1), (270, 1)])
    assert conc == pytest.approx(0, abs=1e-9)


def test_media_circolare_weights(sb):
    ang, _ = sb.media_circolare([(0, 3), (90, 1)])
    assert ang == pytest.approx(math.degrees(math.atan2(1, 3)))


def test_media_circolare_empty(sb):
    _, conc = sb.media_circolare([])
    assert conc == 0


# ---------------- analizza ----------------

@pytest.mark.parametrize("brg", [0, 45, 158, 270, 330])
def test_analizza_straight_coast(sb, brg):
    q, az = grid(sb, straight_coast(brg))
    r = sb.analizza(q, az, sb.RAGGI_KM)
    assert angdiff(r["coast_normal_deg"], brg) <= 6
    assert r["concentrazione"] > 0.7
    assert len(r["open_water_km"]) == sb.SETTORI
    assert max(r["open_water_km"]) == max(sb.RAGGI_KM)
    assert min(r["open_water_km"]) == 0


def test_analizza_all_land_returns_none(sb):
    q, az = grid(sb, lambda a, r: 50.0)
    assert sb.analizza(q, az) is None


def test_analizza_missing_elevations_treated_as_land(sb):
    q, az = grid(sb, lambda a, r: None)
    assert sb.analizza(q, az) is None


def test_analizza_bay_open_water_limited(sb):
    def baia(a, r):
        if angdiff(a, 200) > 45:
            return 30.0
        return -8.0 if r <= 4 else 15.0
    q, az = grid(sb, baia)
    r = sb.analizza(q, az)
    assert angdiff(r["coast_normal_deg"], 200) <= 6
    assert max(r["open_water_km"]) == 4


# ---------------- punto_query ----------------

def test_punto_query_picks_radius_closest_to_3km(sb):
    q = [-1] * len(sb.RAGGI_KM)
    p, r = sb.punto_query(40.0, 17.0, 90, q)
    assert r == 3
    assert p[0] == pytest.approx(40.0, abs=1e-3)
    assert p[1] > 17.0


def test_punto_query_skips_land(sb):
    q = [5 if r < 6 else -5 for r in sb.RAGGI_KM]
    _, r = sb.punto_query(40.0, 17.0, 90, q)
    assert r == 6


def test_punto_query_no_water(sb):
    assert sb.punto_query(40, 17, 0, [10] * len(sb.RAGGI_KM)) == (None, None)


# ---------------- avvisi ----------------

def base(**kw):
    r = {"open_water_km": [35] * 8 + [0] * 8, "frazione_mare": 0.5,
         "concentrazione": 0.9}
    r.update(kw)
    return r


def test_avvisi_clean_coast(sb):
    assert sb.avvisi(base(), "x") == []


def test_avvisi_inland_short_circuits(sb):
    out = sb.avvisi(base(open_water_km=[0.5] * 16, concentrazione=0.1), "x")
    assert len(out) == 1 and "entroterra" in out[0]


def test_avvisi_bay_and_uncertain_normal(sb):
    out = sb.avvisi(base(open_water_km=[4] * 4 + [0] * 12, concentrazione=0.3), "x")
    assert any("baia" in a for a in out)
    assert any("normale incerta" in a for a in out)


def test_avvisi_island(sb):
    out = sb.avvisi(base(open_water_km=[35] * 16, frazione_mare=0.95), "x")
    assert any("isola" in a for a in out)


# ---------------- cmd_check / selftest ----------------

def make_spot(sb, normal, dist_km):
    qlat, qlon = sb.sposta(40.5, 17.0, normal, dist_km)
    return {"name": "t", "coast_normal_deg": normal,
            "shore": {"lat": 40.5, "lon": 17.0},
            "query": {"lat": qlat, "lon": qlon},
            "open_water_km": [10] * 16}


def test_cmd_check_consistent(sb, tmp_path, capsys):
    p = tmp_path / "spots.json"
    p.write_text(json.dumps({"spots": [make_spot(sb, 158, 3)]}), encoding="utf-8")
    assert sb.cmd_check(str(p)) == 0
    assert "COERENZA OK" in capsys.readouterr().out


@pytest.mark.parametrize("normal,query_brg,dist", [(158, 338, 3), (90, 90, 20)])
def test_cmd_check_detects_inconsistency(sb, tmp_path, normal, query_brg, dist):
    s = make_spot(sb, query_brg, dist)
    s["coast_normal_deg"] = normal
    p = tmp_path / "spots.json"
    p.write_text(json.dumps({"spots": [s]}), encoding="utf-8")
    assert sb.cmd_check(str(p)) == 1


def test_builtin_selftest_passes(sb):
    assert sb.selftest() == 0
