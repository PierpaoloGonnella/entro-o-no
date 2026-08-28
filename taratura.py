#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
taratura.py — stima i parametri del modello SUP dai tuoi dati.

Il modello (lo stesso dell'app) dice che, a regime, la potenza propulsiva si
spende in resistenza di scafo piu' resistenza aerodinamica al vento apparente:

    P = [ kh*v^2 + 0.5*rho_a*CaA*(v+Uh)*|v+Uh| ] * v

    v   velocita' ATTRAVERSO L'ACQUA (m/s)
    Uh  componente di vento in faccia lungo la rotta (m/s, negativa se in poppa)
    kh  resistenza idrodinamica di tavola+rider, N/(m/s)^2
    CaA prodotto Cd*area esposta al vento, m^2 (diverso in piedi e in ginocchio)
    P   potenza propulsiva sostenibile, W

Qui si stimano kh, CaA_piedi, CaA_ginocchia, Pw e il coefficiente di deriva
alpha, con la corrente trattata come incognita di disturbo, una per traccia.

--------------------------------------------------------------------------
PERCHE' SERVONO TRATTE ANDATA E RITORNO
Il GPS misura la velocita' AL SUOLO, il modello parla di quella attraverso
l'acqua. Le due differiscono per la corrente. Con almeno due tratte per
sessione su rotte diverse (idealmente opposte) la corrente diventa
identificabile: qui e' stimata come vettore (cx, cy) per traccia, quindi
due parametri in piu' per sessione, non una costante inventata.
Se registri una sola tratta per sessione, la corrente viene assorbita nei
parametri fisici e la taratura vale solo per quel giorno.
--------------------------------------------------------------------------

USO
  python taratura.py --selftest
      Genera dati finti da parametri noti e verifica che il fit li recuperi.
      Da eseguire una volta prima di fidarsi del resto.

  python taratura.py --gpx tracce/*.gpx --out sessioni.csv
      Estrae le tratte rettilinee dai GPX e prepara il CSV da completare.

  python taratura.py --meteo sessioni.csv
      Riempie le colonne di vento e onda interrogando Open-Meteo (richiede rete).

  python taratura.py --fit sessioni.csv
      Stima i parametri, stampa incertezze e diagnostica, ed emette il blocco
      JavaScript da incollare nell'app.

FORMATO DI sessioni.csv (una riga per tratta)
  traccia        identificativo della sessione (le tratte con lo stesso valore
                 condividono la stessa corrente stimata)
  tratta         identificativo della tratta dentro la sessione
  tipo           piedi | ginocchia | deriva
  inizio         ISO 8601 con offset, es. 2026-08-27T09:15:00+02:00
  durata_s       secondi
  distanza_m     metri percorsi al suolo
  rotta_deg      rilevamento medio del moto, 0=N, 90=E
  lat, lon       punto medio della tratta (servono solo a --meteo)
  vento_kn       vento medio a 10 m
  vento_da_deg   direzione DA CUI soffia
  onda_m         altezza significativa, opzionale
  onda_da_deg    direzione DA CUI arriva, opzionale
  corrente_kn    corrente superficiale, opzionale ma molto utile: senza, la
                 corrente resta ignota e si mangia parte della stima
  corrente_verso_deg  direzione VERSO CUI scorre
  note           testo libero

Le tratte 'deriva' sono quelle in cui hai smesso di remare e ti sei lasciato
portare: servono a stimare alpha e vanno tenute corte.
"""

import argparse, csv, glob, io, json, math, sys, urllib.parse, urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np
from scipy.optimize import brentq, least_squares

RHO_A = 1.225
KN = 0.514444

CAMPI = ["traccia","tratta","tipo","inizio","durata_s","distanza_m","rotta_deg",
         "lat","lon","vento_kn","vento_da_deg","onda_m","onda_da_deg",
         "corrente_kn","corrente_verso_deg","note"]


# ==========================================================================
# 1. MODELLO
# ==========================================================================

# --- struttura della carena: fissata, non stimata dai dati ---
NU, GRAV = 1.19e-6, 9.81
HULL = {"L": 3.20, "formK": 1.50, "gamma": 1.14e-3}


def cf_ittc(v):
    """Coefficiente d'attrito, linea ITTC-57."""
    d = math.log10(max(1e4, v*HULL["L"]/NU)) - 2
    return 0.075/(d*d)


def r_hull(v):
    """Resistenza specifica di carena (resistenza / A): attrito+forma piu' residuo."""
    return HULL["formK"]*cf_ittc(v)*v*v + HULL["gamma"]*v**6


def v_acqua(v0, beta, Uh):
    """Velocita' attraverso l'acqua risolvendo  Pw = [r_hull(v) + beta*(v+Uh)|v+Uh|]*v
    con Pw dedotta da v0, la crociera in acqua calma.

    Tutto e' diviso per A = 0,5*rho_w*S. Non e' un dettaglio: in newton assoluti
    la superficie bagnata entra sulla carena e sulla potenza dedotta ma non
    sull'aerodinamica, e fa leva su un risultato che non sai misurare. Divisi
    per A restano solo rapporti, e i parametri liberi sono due, non tre.

    Uh > 0 vento in faccia, Uh < 0 in poppa. Col vento in poppa la funzione non
    e' monotona vicino allo zero: scansiono il primo cambio di segno e poi biseco.
    """
    Pw = (r_hull(v0) + beta*v0*v0) * v0
    f = lambda v: (r_hull(v) + beta*(v+Uh)*abs(v+Uh))*v - Pw
    lo = hi = None
    v = 0.02
    while v <= 8.0:
        if f(v) > 0:
            lo, hi = v - 0.02, v
            break
        v += 0.02
    if hi is None:
        return float("nan")
    return brentq(f, lo, hi, xtol=1e-10, maxiter=200)


def componente_vento(vento_kn, vento_da_deg, rotta_deg):
    """Componente di vento in faccia lungo la rotta, m/s. Positiva = in faccia."""
    return vento_kn * KN * math.cos(math.radians(vento_da_deg - rotta_deg))


# ==========================================================================
# 2. FIT
# ==========================================================================

def prepara(righe):
    """Da righe CSV a array numerici, con indice di traccia per la corrente."""
    tracce = sorted({r["traccia"] for r in righe})
    idx = {t: i for i, t in enumerate(tracce)}
    d = {
        "itr":   np.array([idx[r["traccia"]] for r in righe]),
        "tipo":  np.array([r["tipo"] for r in righe]),
        "vg":    np.array([float(r["distanza_m"]) / float(r["durata_s"]) for r in righe]),
        "rotta": np.array([float(r["rotta_deg"]) for r in righe]),
        "wkn":   np.array([float(r["vento_kn"]) for r in righe]),
        "wdir":  np.array([float(r["vento_da_deg"]) for r in righe]),
        "onda":  np.array([float(r["onda_m"]) if r.get("onda_m") not in (None,"","-") else np.nan
                           for r in righe]),
        "peso":  np.array([math.sqrt(float(r["durata_s"])) for r in righe]),
    }
    d["Uh"] = np.array([componente_vento(k, w, c)
                        for k, w, c in zip(d["wkn"], d["wdir"], d["rotta"])])
    return d, tracce


def previsione(par, d, n_tracce):
    """Velocita' al suolo prevista per ogni tratta."""
    v0, beta_p, beta_g, alpha = np.exp(par[:4])
    cur = par[4:].reshape(n_tracce, 2)                     # (cx=est, cy=nord) m/s
    out = np.empty(len(d["vg"]))
    for i in range(len(out)):
        rot = math.radians(d["rotta"][i])
        hx, hy = math.sin(rot), math.cos(rot)              # versore rotta (est, nord)
        c_lungo = cur[d["itr"][i], 0] * hx + cur[d["itr"][i], 1] * hy
        if d["tipo"][i] == "deriva":
            v = alpha * (-d["Uh"][i])       # deriva pura lungo la direzione del moto
        else:
            beta = beta_g if d["tipo"][i] == "ginocchia" else beta_p
            v = v_acqua(v0, beta, d["Uh"][i])
        out[i] = v + c_lungo
    return out


# valori derivati a mano nel modello dell'app, come punto di partenza e confronto
ATTESI = {"v0": 1.130, "beta_piedi": 4.18e-4, "beta_ginocchia": 2.09e-4, "alpha": 0.060}
NOMI = ["v0", "beta_piedi", "beta_ginocchia", "alpha"]


def fit(righe, verbose=True):
    d, tracce = prepara(righe)
    n = len(tracce)

    p0 = np.concatenate([np.log([ATTESI[k] for k in NOMI]), np.zeros(2 * n)])
    lo = np.concatenate([np.log([0.20, 2e-5, 2e-5, 0.005]), np.full(2 * n, -1.5)])
    hi = np.concatenate([np.log([4.00, 2e-2, 2e-2, 0.400]), np.full(2 * n, 1.5)])

    # Vincolo morbido sulla corrente. Senza, i due componenti per sessione
    # assorbono parte della risposta al vento: con corrente nulla i rapporti si
    # recuperano entro il 6%, con 0,12 m/s di corrente vera l'errore sale al
    # 12-22% e NON migliora aggiungendo rotte. Non e' geometria, e' confondimento.
    # Il rimedio e' non trattare la corrente come del tutto ignota: Open-Meteo la
    # prevede, quindi la si tira verso quel valore invece di lasciarla libera.
    # SIGMA_CUR e' quanto ti fidi della previsione, in m/s.
    SIGMA_CUR = 0.08
    W_CUR = 3.0 / SIGMA_CUR
    cur_prior = np.zeros(2 * n)
    for k, tr in enumerate(tracce):
        for r0 in righe:
            if r0["traccia"] == tr and r0.get("corrente_kn"):
                mag = float(r0["corrente_kn"]) * KN
                vs = math.radians(float(r0.get("corrente_verso_deg") or 0))
                cur_prior[2*k]   = mag * math.sin(vs)
                cur_prior[2*k+1] = mag * math.cos(vs)
                break

    def residui(par):
        r = (previsione(par, d, n) - d["vg"]) * d["peso"]
        pr = (par[4:] - cur_prior) * W_CUR
        return np.nan_to_num(np.concatenate([r, pr]), nan=1e3)

    sol = least_squares(residui, p0, bounds=(lo, hi), x_scale="jac", max_nfev=20000)

    val = dict(zip(NOMI, np.exp(sol.x[:4])))
    pred = previsione(sol.x, d, n)
    res = pred - d["vg"]
    dof = max(1, len(res) - len(sol.x))
    s2 = float(res @ res) / dof
    rmse = float(np.sqrt(np.mean(res ** 2)))
    ss_tot = float(np.sum((d["vg"] - d["vg"].mean()) ** 2))
    r2 = 1.0 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")

    # incertezze: i parametri sono stimati in log, quindi sigma_log ~ sigma relativa
    try:
        cov = np.linalg.pinv(sol.jac.T @ sol.jac) * s2
        sig_log = np.sqrt(np.clip(np.diag(cov)[:4], 0, None))
    except Exception:
        sig_log = np.full(4, float("nan"))

    ris = dict(val)
    ris.update({
        "sigma_rel": dict(zip(NOMI, sig_log)),
        "rmse_ms": rmse, "r2": r2, "n_tratte": len(res), "n_tracce": n,
        "correnti": {t: (float(sol.x[4+2*i]), float(sol.x[5+2*i])) for i, t in enumerate(tracce)},
        "res": res, "d": d,
    })
    if verbose:
        stampa(ris)
    return ris


# ==========================================================================
# 3. DIAGNOSTICA
# ==========================================================================

def stampa(r):
    print("\n=== PARAMETRI STIMATI " + "=" * 44)
    print(f"  {r['n_tratte']} tratte su {r['n_tracce']} sessioni")
    unita = {"v0": "m/s", "beta_piedi": "-", "beta_ginocchia": "-", "alpha": "-"}
    for k in NOMI:
        s = r["sigma_rel"][k]
        pct = f"±{100*s:.0f}%" if np.isfinite(s) else "±?"
        print(f"  {k:16s} {r[k]:11.6f} {pct:>7s}   {unita[k]:4s} "
              f"(derivato a mano: {ATTESI[k]})")
    # v0 non e' la crociera: e' (P/kh)^(1/3). La crociera la da' il solver a Uh=0.
    print(f"\n  crociera ad acqua ferma: "
          f"{v_acqua(r['v0'], r['beta_piedi'], 0.0)/KN:.2f} kn")
    print(f"  RMSE {r['rmse_ms']:.3f} m/s   R2 {r['r2']:.3f}")
    for wkn in (10, 15, 20):
        U = wkn * KN
        vp = v_acqua(r["v0"], r["beta_piedi"], U)
        vg = v_acqua(r["v0"], r["beta_ginocchia"], U)
        t = 400 / vp / 60 if vp > 0.02 else float("inf")
        tg = 400 / vg / 60 if vg > 0.02 else float("inf")
        print(f"  con {wkn:2d} kn in faccia: rientro da 400 m in {t:4.0f} min "
              f"({tg:.0f} in ginocchio)")

    print("\n=== CORRENTE STIMATA PER SESSIONE " + "=" * 32)
    for t, (cx, cy) in r["correnti"].items():
        mag = math.hypot(cx, cy)
        verso = (math.degrees(math.atan2(cx, cy)) + 360) % 360
        print(f"  {t:14s} {mag:5.2f} m/s verso {verso:5.1f}°")
        if mag > 0.8:
            print("      ^ sospetta: piu' di 1,5 nodi di corrente. Controlla che le "
                  "tratte di quella sessione abbiano rotte diverse.")

    # la domanda che decide se serve altro: il residuo dipende dall'onda?
    res, onda = r["res"], r["d"]["onda"]
    m = np.isfinite(onda)
    print("\n=== LA FISICA BASTA? " + "=" * 45)
    if m.sum() < 8:
        print("  Servono almeno 8 tratte con l'altezza d'onda compilata per dirlo.")
    else:
        rr = float(np.corrcoef(res[m], onda[m])[0, 1])
        print(f"  correlazione residuo/onda su {int(m.sum())} tratte: r = {rr:+.2f}")
        if abs(rr) < 0.25:
            print("  Residuo scorrelato dall'onda: il modello fisico basta.")
            print("  Aggiungere una rete non migliorerebbe niente. Fermati qui.")
        else:
            print("  C'e' struttura residua legata all'onda: il modello sbaglia in modo")
            print("  sistematico col mare formato. E' qui che un grey-box ha senso, con")
            print("  la fisica sul termine noto e una piccola rete sulla correzione")
            print("  d'onda. Prima pero' verifica che non sia un solo giorno a trascinare")
            print("  la correlazione, e che l'onda non sia solo un altro modo di")
            print("  misurare il vento.")

    print("\n=== DA INCOLLARE NELL'APP " + "=" * 40)
    print(f"  // tarato su {r['n_tratte']} tratte, {r['n_tracce']} sessioni")
    print(f"  v0:{r['v0']:.3f}, betaStand:{r['beta_piedi']:.3e}, "
          f"betaKneel:{r['beta_ginocchia']:.3e}, alpha:{r['alpha']:.4f},")
    print()


# ==========================================================================
# 4. AUTOTEST
# ==========================================================================

def selftest(seed=7):
    """Genera tratte finte da parametri noti e verifica che il fit li recuperi.

    L'onda e' generata indipendente dal vento: cosi' l'autotest verifica anche
    che la diagnostica dica correttamente 'la fisica basta'.
    """
    rng = np.random.default_rng(seed)
    VERO = dict(v0=1.05, beta_piedi=4.6e-4, beta_ginocchia=2.2e-4, alpha=0.055)
    # La previsione di corrente c'e', come nel flusso reale dopo --meteo, ma con
    # un errore di previsione di 0,05 m/s: e' un vincolo morbido, non la verita'.
    righe = []
    for s in range(8):
        wkn = float(rng.uniform(2, 22))
        wdir = float(rng.uniform(0, 360))
        cx, cy = rng.normal(0, 0.12, 2)
        px, py = cx + rng.normal(0, 0.05), cy + rng.normal(0, 0.05)
        prev = dict(corrente_kn=round(math.hypot(px, py)/KN, 3),
                    corrente_verso_deg=round((math.degrees(math.atan2(px, py)) + 360) % 360, 1))
        onda = float(rng.uniform(0.05, 0.6))          # indipendente dal vento
        for j, rot in enumerate([wdir, (wdir + 180) % 360, (wdir + 90) % 360]):
            for tipo in (["piedi", "ginocchia"] if j == 0 else ["piedi"]):
                Uh = componente_vento(wkn, wdir, rot)
                beta = VERO["beta_ginocchia"] if tipo == "ginocchia" else VERO["beta_piedi"]
                v = v_acqua(VERO["v0"], beta, Uh)
                hx, hy = math.sin(math.radians(rot)), math.cos(math.radians(rot))
                vg = (v + cx * hx + cy * hy) * float(rng.normal(1.0, 0.04))
                dur = 240.0
                righe.append(dict(traccia=f"S{s}", tratta=f"{j}{tipo[0]}", tipo=tipo,
                                  inizio="", durata_s=dur, distanza_m=max(5.0, vg * dur),
                                  rotta_deg=rot, lat="", lon="", vento_kn=wkn,
                                  vento_da_deg=wdir, onda_m=onda, onda_da_deg=wdir,
                                  note="finto", **prev))
        rot_der = (wdir + 180) % 360
        v = VERO["alpha"] * wkn * KN
        hx, hy = math.sin(math.radians(rot_der)), math.cos(math.radians(rot_der))
        vg = max(0.03, v + cx * hx + cy * hy) * float(rng.normal(1.0, 0.06))
        righe.append(dict(traccia=f"S{s}", tratta="drift", tipo="deriva", inizio="",
                          durata_s=120.0, distanza_m=vg * 120.0, rotta_deg=rot_der,
                          lat="", lon="", vento_kn=wkn, vento_da_deg=wdir,
                          onda_m=onda, onda_da_deg=wdir, note="finto", **prev))

    print("VERI:", {k: round(v, 4) for k, v in VERO.items()})
    r = fit(righe)
    print("=== RECUPERO DEI PARAMETRI " + "=" * 39)
    ok = True
    for k in NOMI:
        # Tolleranza 25% sui rapporti: con corrente reale e previsione imperfetta
        # e' il meglio ottenibile, verificato per esperimento. v0 sta entro il 5%.
        soglia = 0.05 if k == "v0" else 0.25
        err = abs(r[k] - VERO[k]) / VERO[k]
        if err >= soglia:
            ok = False
        print(f"  {k:16s} vero {VERO[k]:11.6f}  stimato {r[k]:11.6f}  "
              f"errore {100*err:5.1f}%  {'ok' if err < soglia else 'FUORI'}")
    print("\nAUTOTEST " + ("SUPERATO" if ok else "FALLITO"))
    return 0 if ok else 1


# ==========================================================================
# 5. GPX -> TRATTE
# ==========================================================================

def _haversine(a, b):
    R = 6371000.0
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dp = p2 - p1
    dl = math.radians(b[1] - a[1])
    h = math.sin(dp/2)**2 + math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(h))


def _rilevamento(a, b):
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dl = math.radians(b[1] - a[1])
    x = math.sin(dl) * math.cos(p2)
    y = math.cos(p1)*math.sin(p2) - math.sin(p1)*math.cos(p2)*math.cos(dl)
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def leggi_gpx(path):
    import xml.etree.ElementTree as ET
    punti = []
    for el in ET.parse(path).getroot().iter():
        if not el.tag.endswith("trkpt"):
            continue
        lat, lon = float(el.attrib["lat"]), float(el.attrib["lon"])
        t = None
        for c in el:
            if c.tag.endswith("time") and c.text:
                s = c.text.strip().replace("Z", "+00:00")
                t = datetime.fromisoformat(s)
        if t:
            punti.append((t, lat, lon))
    punti.sort(key=lambda p: p[0])
    return punti


def estrai_tratte(punti, min_durata=120, max_sd_rotta=15.0, min_dist=100.0):
    """Spezza la traccia dove la rotta cambia: tiene solo i rettilinei lunghi.

    Una tratta e' utile solo se hai tenuto rotta e sforzo costanti: qui si
    controlla la rotta, lo sforzo lo sai solo tu (metti una nota).
    """
    if len(punti) < 10:
        return []
    tratte, i = [], 0
    while i < len(punti) - 2:
        j = i + 2
        rott = [_rilevamento(punti[k][1:], punti[k+1][1:]) for k in range(i, j)]
        while j < len(punti) - 1:
            r = _rilevamento(punti[j][1:], punti[j+1][1:])
            prova = rott + [r]
            base = prova[0]
            centr = [((x - base + 540) % 360) - 180 for x in prova]
            if float(np.std(centr)) > max_sd_rotta:
                break
            rott = prova
            j += 1
        dur = (punti[j][0] - punti[i][0]).total_seconds()
        dist = sum(_haversine(punti[k][1:], punti[k+1][1:]) for k in range(i, j))
        if dur >= min_durata and dist >= min_dist:
            base = rott[0]
            centr = [((x - base + 540) % 360) - 180 for x in rott]
            tratte.append(dict(
                inizio=punti[i][0].isoformat(),
                durata_s=round(dur, 1),
                distanza_m=round(dist, 1),
                rotta_deg=round((base + float(np.mean(centr))) % 360, 1),
                lat=round((punti[i][1] + punti[j][1]) / 2, 5),
                lon=round((punti[i][2] + punti[j][2]) / 2, 5),
            ))
            i = j
        else:
            i += 1
    return tratte


def cmd_gpx(paths, out):
    righe = []
    for p in sorted(paths):
        punti = leggi_gpx(p)
        tr = estrai_tratte(punti)
        nome = p.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        print(f"{nome}: {len(punti)} punti -> {len(tr)} tratte utili")
        for k, t in enumerate(tr, 1):
            righe.append(dict(traccia=nome, tratta=str(k), tipo="piedi", note="",
                              vento_kn="", vento_da_deg="", onda_m="", onda_da_deg="", **t))
    if not righe:
        print("Nessuna tratta rettilinea abbastanza lunga. Serve almeno un tratto di "
              "2 minuti a rotta costante: pianifica l'uscita con andata e ritorno lunghi.")
        return 1
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CAMPI)
        w.writeheader()
        for r in righe:
            w.writerow({k: r.get(k, "") for k in CAMPI})
    print(f"\nScritto {out} con {len(righe)} tratte.")
    print("Ora: metti 'tipo' corretto su ogni riga (piedi/ginocchia/deriva), poi")
    print(f"  python taratura.py --meteo {out}")
    return 0


# ==========================================================================
# 6. METEO STORICO
# ==========================================================================

def cmd_meteo(path):
    """Riempie vento e onda dall'archivio Open-Meteo. Richiede rete."""
    with open(path, newline="", encoding="utf-8") as f:
        righe = list(csv.DictReader(f))
    fatte = 0
    for r in righe:
        if r.get("vento_kn"):
            continue
        if not (r.get("lat") and r.get("lon") and r.get("inizio")):
            print(f"  salto tratta {r['traccia']}/{r['tratta']}: manca lat/lon/inizio")
            continue
        t = datetime.fromisoformat(r["inizio"])
        giorno = t.date().isoformat()
        base = dict(latitude=r["lat"], longitude=r["lon"],
                    start_date=giorno, end_date=giorno,
                    timeformat="unixtime")
        try:
            q = dict(base, hourly="wind_speed_10m,wind_direction_10m", wind_speed_unit="kn")
            aria = _get("https://archive-api.open-meteo.com/v1/archive?" +
                        urllib.parse.urlencode(q))
            k = _ora_vicina(aria["hourly"]["time"], t)
            r["vento_kn"] = round(aria["hourly"]["wind_speed_10m"][k], 1)
            r["vento_da_deg"] = round(aria["hourly"]["wind_direction_10m"][k], 0)
        except Exception as e:
            print(f"  vento non recuperato per {r['traccia']}/{r['tratta']}: {e}")
            continue
        try:
            q = dict(base, hourly="wave_height,wave_direction,ocean_current_velocity,ocean_current_direction")
            mare = _get("https://marine-api.open-meteo.com/v1/marine?" +
                        urllib.parse.urlencode(q))
            k = _ora_vicina(mare["hourly"]["time"], t)
            r["onda_m"] = mare["hourly"]["wave_height"][k]
            r["onda_da_deg"] = mare["hourly"]["wave_direction"][k]
            if mare["hourly"].get("ocean_current_velocity"):
                r["corrente_kn"] = round(mare["hourly"]["ocean_current_velocity"][k]/1.852, 2)
                r["corrente_verso_deg"] = mare["hourly"]["ocean_current_direction"][k]
        except Exception as e:
            print(f"  onda non recuperata per {r['traccia']}/{r['tratta']}: {e}")
        fatte += 1
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CAMPI)
        w.writeheader()
        for r in righe:
            w.writerow({k: r.get(k, "") for k in CAMPI})
    print(f"Compilate {fatte} tratte. Controlla i valori a mano: l'archivio e' "
          "una rianalisi, non una misura sul posto.")
    return 0


def _get(url):
    with urllib.request.urlopen(url, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _ora_vicina(times, t):
    target = t.timestamp()
    return int(np.argmin([abs(x - target) for x in times]))


# ==========================================================================
# 7. CLI
# ==========================================================================

def cmd_fit(path):
    with open(path, newline="", encoding="utf-8") as f:
        righe = [r for r in csv.DictReader(f)
                 if r.get("vento_kn") and r.get("durata_s") and r.get("distanza_m")]
    if len(righe) < 8:
        print(f"Solo {len(righe)} tratte complete. Sotto le 8 la stima non significa "
              "niente: i parametri sono quattro piu' due per sessione.")
        return 1
    tracce = {r["traccia"] for r in righe}
    if len(tracce) < 3:
        print(f"Attenzione: solo {len(tracce)} sessioni. Con poche sessioni la corrente "
              "e i parametri fisici non si separano bene.")
    fit(righe)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--gpx", nargs="+", metavar="FILE")
    ap.add_argument("--out", default="sessioni.csv")
    ap.add_argument("--meteo", metavar="CSV")
    ap.add_argument("--fit", metavar="CSV")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.gpx:
        paths = [p for g in a.gpx for p in (glob.glob(g) or [g])]
        return cmd_gpx(paths, a.out)
    if a.meteo:
        return cmd_meteo(a.meteo)
    if a.fit:
        return cmd_fit(a.fit)
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
