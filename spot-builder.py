#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spot-builder.py — genera spots.json dalle sole coordinate di partenza.

PERCHE' ESISTE
Il campo coast_normal_deg (il rilevamento che punta dalla riva verso il largo)
per i primi tre spot l'ho ricavato a mano dalla geometria della costa. Non
scala, e soprattutto e' l'unico campo che se sbagliato INVERTE il verdetto: un
errore di 180 gradi trasforma "il vento ti spinge al largo" in "il vento ti
riporta a riva". Qui viene derivato dall'elevazione del terreno e poi
verificato, e i casi dubbi vengono segnalati invece di essere indovinati.

COME
Si campiona un anello di punti attorno al punto di partenza interrogando
/v1/elevation di Open-Meteo (Copernicus DEM GLO-90, 90 m, nessuna chiave, fino
a 100 coordinate per richiesta). Il mare e' dove l'elevazione e' <= 0.
Da qui si ricavano tre cose:

  coast_normal_deg  media circolare delle direzioni di mare, pesata sulla
                    frazione di mare in quella direzione.
  query             punto spostato al largo lungo la normale, verificato in
                    acqua: e' quello che si interroga per il meteo, perche' la
                    cella di terra sovrastima le raffiche e non ha onda.
  open_water_km     per ciascuno dei 16 settori, quanti km di mare aperto ci
                    sono prima di ritrovare terra. Serve a distinguere una
                    deriva verso cento chilometri di Golfo da una deriva verso
                    la sponda opposta di una baia: stesso vento, conseguenza
                    diversa.

USO
  python spot-builder.py --selftest
      Verifica il metodo su geometrie sintetiche di normale nota. Da eseguire
      prima di fidarsi del resto: non richiede rete.

  python spot-builder.py --input partenze.csv --out spots.json
      partenze.csv: id,name,sub,area,lat,lon,hazards
      Interroga l'API e scrive spots.json. Richiede rete.

  python spot-builder.py --check spots.json
      Rilegge un spots.json e ne verifica la coerenza interna.
"""

import argparse, csv, json, math, sys, time, urllib.parse, urllib.request

R_TERRA = 6371.0                      # km
RAGGI_KM = [0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 25, 35]
N_AZ = 36                             # un raggio ogni 10 gradi
SETTORI = 16                          # per open_water_km
LOTTO = 100                           # coordinate per richiesta


# ==========================================================================
# geometria
# ==========================================================================

def sposta(lat, lon, brg_deg, dist_km):
    """Punto a dist_km lungo il rilevamento brg_deg, su sfera."""
    b, la = math.radians(brg_deg), math.radians(lat)
    d = dist_km / R_TERRA
    la2 = math.asin(math.sin(la)*math.cos(d) + math.cos(la)*math.sin(d)*math.cos(b))
    lo2 = math.radians(lon) + math.atan2(math.sin(b)*math.sin(d)*math.cos(la),
                                         math.cos(d) - math.sin(la)*math.sin(la2))
    return math.degrees(la2), (math.degrees(lo2) + 540) % 360 - 180


def media_circolare(dirs_pesi):
    """Media circolare di (direzione, peso). Restituisce (dir, concentrazione).

    La concentrazione va da 0 a 1: vicino a 1 le direzioni di mare sono
    concordi (costa rettilinea), vicino a 0 sono sparse in ogni verso (punta,
    istmo, laguna) e la normale non e' definita in modo affidabile.
    """
    sx = sum(p*math.sin(math.radians(d)) for d, p in dirs_pesi)
    sy = sum(p*math.cos(math.radians(d)) for d, p in dirs_pesi)
    tot = sum(p for _, p in dirs_pesi) or 1.0
    ang = (math.degrees(math.atan2(sx, sy)) + 360) % 360
    return ang, math.hypot(sx, sy)/tot


# ==========================================================================
# analisi di un anello
# ==========================================================================

def analizza(quote, azimut=None, raggi=None):
    """Da una griglia di quote (azimut x raggi) ricava normale e mare aperto.

    quote[i][j] = elevazione in metri all'azimut azimut[i], raggio raggi[j].
    Nessuna rete: la funzione e' pura, ed e' quella che l'autotest verifica.
    """
    az = azimut if azimut is not None else [i*360.0/N_AZ for i in range(N_AZ)]
    rg = raggi if raggi is not None else RAGGI_KM

    frazione, aperto = [], []
    for i, a in enumerate(az):
        riga = quote[i]
        mare = [q is not None and q <= 0 for q in riga]
        frazione.append(sum(mare)/len(mare))
        # mare aperto: fin dove il mare e' continuo partendo dal raggio minore
        km = 0.0
        for j, m in enumerate(mare):
            if not m:
                break
            km = rg[j]
        aperto.append(km)

    # la normale pesa ogni direzione sulla sua frazione di mare, elevata al
    # quadrato: cosi' una direzione mezza-terra conta molto meno di una tutta
    # mare, e le direzioni radenti alla costa non trascinano la media.
    pesi = [(a, f**2) for a, f in zip(az, frazione) if f > 0.35]
    if not pesi:
        return None
    normale, conc = media_circolare(pesi)

    # open_water_km aggregato in 16 settori, prendendo il minimo del settore:
    # se in un settore c'e' un promontorio, e' quello che conta.
    settori = []
    for s in range(SETTORI):
        c = s*360.0/SETTORI
        vals = [aperto[i] for i, a in enumerate(az)
                if min(abs(((a-c+180) % 360)-180), 360) <= 360/SETTORI/2]
        settori.append(round(min(vals), 1) if vals else 0.0)

    return {
        "coast_normal_deg": round(normale, 1),
        "concentrazione": round(conc, 3),
        "frazione_mare": round(sum(frazione)/len(frazione), 3),
        "open_water_km": settori,
    }


def punto_query(lat, lon, normale, quote_lungo_normale, raggi=RAGGI_KM, voluto=3.0):
    """Sceglie il punto al largo: il raggio piu' vicino a 3 km che sia in acqua."""
    in_acqua = [r for r, q in zip(raggi, quote_lungo_normale)
                if q is not None and q <= 0]
    if not in_acqua:
        return None, None
    r = min(in_acqua, key=lambda x: abs(x - voluto))
    qlat, qlon = sposta(lat, lon, normale, r)
    return (round(qlat, 4), round(qlon, 4)), r


# ==========================================================================
# rete
# ==========================================================================

def quote_api(punti, pausa=1.1):
    """Elevazione per una lista di (lat, lon), a lotti da 100."""
    out = []
    for k in range(0, len(punti), LOTTO):
        blocco = punti[k:k+LOTTO]
        q = urllib.parse.urlencode({
            "latitude":  ",".join(f"{a:.5f}" for a, _ in blocco),
            "longitude": ",".join(f"{b:.5f}" for _, b in blocco),
        })
        url = "https://api.open-meteo.com/v1/elevation?" + q
        with urllib.request.urlopen(url, timeout=45) as r:
            out.extend(json.loads(r.read().decode())["elevation"])
        if k + LOTTO < len(punti):
            time.sleep(pausa)          # servizio gratuito: non lo si martella
    return out


def costruisci(lat, lon):
    az = [i*360.0/N_AZ for i in range(N_AZ)]
    punti = [sposta(lat, lon, a, r) for a in az for r in RAGGI_KM]
    piatte = quote_api(punti)
    n = len(RAGGI_KM)
    quote = [piatte[i*n:(i+1)*n] for i in range(len(az))]
    ris = analizza(quote, az, RAGGI_KM)
    if ris is None:
        return None, quote, az
    # per il punto di query servono le quote LUNGO la normale, non su un raggio
    # dell'anello: interpolo prendendo l'azimut campionato piu' vicino.
    i = min(range(len(az)), key=lambda k: abs(((az[k]-ris["coast_normal_deg"]+180) % 360)-180))
    ris["query_lungo"] = quote[i]
    return ris, quote, az


# ==========================================================================
# avvisi
# ==========================================================================

def avvisi(r, nome):
    """Segnala i casi in cui la derivazione automatica non e' affidabile.

    L'ordine conta: una baia ha poca frazione di mare come l'entroterra, ma le
    due cose si distinguono guardando ANCHE la concentrazione e il mare aperto.
    Una baia ha mare concentrato in un settore stretto e corto; l'entroterra non
    ha mare in nessun settore. Prima distinguevo solo sulla frazione, e la baia
    dell'autotest si prendeva l'avviso "coordinate nell'entroterra".
    """
    out = []
    apmax, apmin = max(r["open_water_km"]), min(r["open_water_km"])
    frz, conc = r["frazione_mare"], r["concentrazione"]

    if apmax < 0.6:
        out.append("nessuna direzione con mare continuo: coordinate quasi "
                   "certamente nell'entroterra.")
        return out

    if apmax < 8:
        out.append("mare aperto al massimo %.1f km: bacino chiuso o baia. Il "
                   "modello d'onda di Open-Meteo qui puo' non dare dati, e una "
                   "deriva finisce sulla sponda opposta invece che al largo — "
                   "meno grave, ma da sapere." % apmax)

    if conc < 0.55:
        out.append("normale incerta (concentrazione %.2f): mare su troppi lati, "
                   "probabile punta o istmo. Verificala a mano su una mappa: "
                   "sbagliarla di 180 gradi capovolge il verdetto." % conc)

    if frz > 0.85 and apmin > 3:
        out.append("mare in quasi tutte le direzioni e nessun settore chiuso: "
                   "il punto sembra su un'isola o al largo, non su una spiaggia.")

    return out


# ==========================================================================
# comandi
# ==========================================================================

def cmd_input(path, out):
    with open(path, newline="", encoding="utf-8") as f:
        righe = list(csv.DictReader(f))
    spots, problemi = [], 0
    for r in righe:
        lat, lon = float(r["lat"]), float(r["lon"])
        print(f"\n{r['name']} ({lat}, {lon})")
        try:
            ris, quote, az = costruisci(lat, lon)
        except Exception as e:
            print(f"  ERRORE di rete: {e}")
            problemi += 1
            continue
        if ris is None:
            print("  nessuna direzione di mare trovata: coordinate sbagliate?")
            problemi += 1
            continue
        q, raggio = punto_query(lat, lon, ris["coast_normal_deg"], ris["query_lungo"])
        if q is None:
            print("  nessun punto in acqua lungo la normale: salto.")
            problemi += 1
            continue
        print(f"  normale {ris['coast_normal_deg']}°  "
              f"concentrazione {ris['concentrazione']}  "
              f"query {q} a {raggio} km")
        for a in avvisi(ris, r["name"]):
            print("  ! " + a)
            problemi += 1
        spots.append({
            "id": r["id"], "name": r["name"], "sub": r.get("sub", ""),
            "area": r.get("area", ""),
            "shore": {"lat": lat, "lon": lon},
            "query": {"lat": q[0], "lon": q[1]},
            "coast_normal_deg": ris["coast_normal_deg"],
            "open_water_km": ris["open_water_km"],
            "derivazione": {"concentrazione": ris["concentrazione"],
                            "frazione_mare": ris["frazione_mare"],
                            "query_km": raggio,
                            "fonte": "Copernicus DEM GLO-90 via Open-Meteo /v1/elevation"},
            "hazards": r.get("hazards", ""),
        })
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"spots": spots}, f, ensure_ascii=False, indent=1)
    print(f"\nScritti {len(spots)} spot in {out}. Avvisi: {problemi}.")
    print("Controlla a mano ogni spot con un avviso: la normale sbagliata di 180° "
          "capovolge il verdetto, ed e' l'errore piu' pericoloso possibile qui.")
    return 0


def cmd_check(path):
    d = json.load(open(path, encoding="utf-8"))
    ok = True
    for s in d["spots"]:
        n = s["coast_normal_deg"]
        # il punto di query deve stare dal lato mare rispetto alla riva
        dlat = s["query"]["lat"] - s["shore"]["lat"]
        dlon = (s["query"]["lon"] - s["shore"]["lon"]) * math.cos(math.radians(s["shore"]["lat"]))
        brg = (math.degrees(math.atan2(dlon, dlat)) + 360) % 360
        scarto = abs(((brg - n + 180) % 360) - 180)
        km = math.hypot(dlat*111.32, dlon*111.32)
        stato = "ok" if scarto < 5 and 1.0 < km < 8.0 else "INCOERENTE"
        if stato != "ok":
            ok = False
        print(f"  {s['name']:22s} normale {n:5.1f}°  query a {brg:5.1f}° "
              f"({km:4.1f} km)  scarto {scarto:4.1f}°  {stato}")
        if min(s.get("open_water_km", [99])) < 1:
            print("      settore con meno di 1 km di mare: fondo chiuso, "
                  "attenzione ai dati d'onda.")
    print("\nCOERENZA", "OK" if ok else "DA CORREGGERE")
    return 0 if ok else 1


# ==========================================================================
# autotest: geometrie sintetiche di normale nota, nessuna rete
# ==========================================================================

def _quote_finte(f, az, raggi):
    return [[f(a, r) for r in raggi] for a in az]


def selftest():
    az = [i*360.0/N_AZ for i in range(N_AZ)]
    casi = []

    # 1. costa rettilinea: mare in un semipiano, normale attesa = brg
    for brg in (0, 45, 158, 270, 330):
        def f(a, r, b=brg):
            comp = math.cos(math.radians(a - b))    # +1 verso il largo
            return -10.0 if comp*r > 0.3 else 20.0
        casi.append((f"costa rettilinea, normale {brg}°", f, brg, "conc>0.7"))

    # 2. baia: mare solo in un settore di 90° centrato su 200°
    def baia(a, r):
        d = abs(((a - 200 + 180) % 360) - 180)
        if d > 45:
            return 30.0
        return -8.0 if r <= 4 else 15.0            # sponda opposta a 4 km
    casi.append(("baia chiusa aperta a 200°", baia, 200, "aperto<=4"))

    # 3. punta: mare su tre lati, normale non definita
    def punta(a, r):
        d = abs(((a - 90 + 180) % 360) - 180)
        return 25.0 if d > 130 else -10.0
    casi.append(("punta, mare su 260°", punta, 90, "conc bassa"))

    print("AUTOTEST — geometrie sintetiche, nessuna rete\n")
    ok = True
    for nome, f, atteso, nota in casi:
        q = _quote_finte(f, az, RAGGI_KM)
        r = analizza(q, az, RAGGI_KM)
        if r is None:
            print(f"  {nome:34s} nessuna direzione di mare  FALLITO")
            ok = False
            continue
        err = abs(((r["coast_normal_deg"] - atteso + 180) % 360) - 180)
        buono = err <= 6
        if not buono:
            ok = False
        print(f"  {nome:34s} normale {r['coast_normal_deg']:6.1f}° "
              f"(attesa {atteso:5.1f}°, errore {err:4.1f}°)  "
              f"conc {r['concentrazione']:.2f}  "
              f"mare aperto min/max {min(r['open_water_km']):4.1f}/"
              f"{max(r['open_water_km']):4.1f} km  "
              f"{'ok' if buono else 'FALLITO'}")
        for a in avvisi(r, nome):
            print(f"      ! {a}")

    # il punto di query deve cadere in acqua
    def retta158(a, r):
        return -10.0 if math.cos(math.radians(a - 158))*r > 0.3 else 20.0
    q = _quote_finte(retta158, az, RAGGI_KM)
    r = analizza(q, az, RAGGI_KM)
    i = min(range(len(az)), key=lambda k: abs(((az[k]-r["coast_normal_deg"]+180) % 360)-180))
    p, raggio = punto_query(40.5148, 17.0553, r["coast_normal_deg"], q[i])
    print(f"\n  punto di query: {p} a {raggio} km lungo {r['coast_normal_deg']}°")
    if raggio is None or not (2.0 <= raggio <= 4.0):
        print("      FALLITO: raggio fuori dall'intervallo voluto")
        ok = False

    print("\nAUTOTEST " + ("SUPERATO" if ok else "FALLITO"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--input", metavar="CSV")
    ap.add_argument("--out", default="spots.json")
    ap.add_argument("--check", metavar="JSON")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.input:
        return cmd_input(a.input, a.out)
    if a.check:
        return cmd_check(a.check)
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
