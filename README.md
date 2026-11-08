# Entro o no?

[![tests](https://github.com/PierpaoloGonnella/entro-o-no/actions/workflows/tests.yml/badge.svg)](https://github.com/PierpaoloGonnella/entro-o-no/actions/workflows/tests.yml)

Semaforo per il SUP sulla costa ionica, e su qualunque altro punto di costa.

Non risponde alla domanda "quanto soffia", che sanno già tutte le app meteo.
Risponde a **"quanto ci metto a rientrare da 400 metri con il vento che c'è"**,
che è la domanda che decide se torni a riva.

Un solo file HTML, nessun server, nessuna chiave API, nessun tracciamento.

---

## Pubblicarla

1. Crea un repository pubblico su GitHub, per esempio `entro-o-no`.
2. Metti dentro `entro-o-no.html` **rinominato `index.html`**. È l'unico file
   necessario: gli altri sono strumenti di lavoro, non servono al sito.
3. Repository → **Settings → Pages**. In *Source* scegli *Deploy from a branch*,
   ramo `main`, cartella `/ (root)`. Salva.
4. Dopo un minuto il sito è su `https://<tuo-utente>.github.io/entro-o-no/`.

Se vuoi un indirizzo tuo, in *Settings → Pages → Custom domain* metti il
dominio e aggiungi un record `CNAME` che punta a `<tuo-utente>.github.io`.

### Perché pubblicarla e non passare il file

Non è solo comodità. **La geolocalizzazione richiede HTTPS**: il pulsante
"Vicino a me" non funziona se il file viene aperto localmente o servito in
chiaro. GitHub Pages dà HTTPS gratis, e questo è il motivo tecnico per cui
vale la pena mettercela.

### Quanto costa in traffico verso Open-Meteo

Niente, e per una ragione che vale conoscere: **le chiamate partono dal browser
di chi visita**, non da un tuo server. Ogni visitatore usa la propria quota
giornaliera di Open-Meteo, quindi il sito regge un numero qualunque di persone
senza che tu paghi o configuri nulla.

In cambio hai due obblighi, entrambi già rispettati nel piede della pagina:
i dati Open-Meteo sono sotto licenza **CC BY 4.0** e vanno attribuiti, e il
piano gratuito è per **uso non commerciale**. Se un giorno ci metti pubblicità
o un abbonamento, quel presupposto cade e servono i termini a pagamento.

Non aggiungere Analytics o simili: la pagina oggi non manda niente a nessuno
tranne le richieste meteo, ed è una proprietà che si perde con una riga.

---

## Cosa c'è nel repository

| File | A cosa serve |
|---|---|
| `index.html` | L'applicazione. È autosufficiente: quattro spot già dentro, e la capacità di derivarne di nuovi. |
| `spot-builder.py` | Deriva l'orientamento della costa e il mare aperto dall'elevazione del terreno, per preparare spot verificati con calma. La stessa logica è dentro l'HTML, ma qui l'anello è più fitto. |
| `taratura.py` | Stima i parametri del modello dalle tue tracce GPS. Senza questo, l'app usa i miei numeri derivati a tavolino. |
| `spots.json` | Non più necessario all'app. Resta come formato di scambio per una lista curata. |

Entrambi gli script hanno `--selftest`: eseguilo prima di fidarti dei risultati.

---

## Link condivisibili

- `?spot=lizzano` apre uno spot già catalogato.
- `?lat=40.3261&lon=17.4518` analizza un punto qualunque e lo apre.
- `?pos=40.33,17.45` finge di essere lì, per provare il comportamento altrove.

Il primo e il secondo sono utili da mandare in chat. Il terzo mostra sempre una
fascia gialla, perché una posizione simulata dimenticata è un modo di guardare
le condizioni della spiaggia sbagliata.

---

## Cosa devi dire a chi la usa

Questa è la parte che conta più del codice, perché lo strumento parla a persone
che stanno decidendo se entrare in acqua.

**I parametri del pagaiatore non sono i loro.** `v0` (velocità di crociera in
acqua calma) e `beta` (quanto crolli col vento) nell'app sono valori che ho
derivato da un bilancio di forze, non misurati. L'incertezza propagata sui
minuti di rientro è del 50% a vento debole e del 75% a 25 nodi: per questo il
rientro è mostrato come **intervallo** e non come numero, e il colore lo decide
l'estremo pessimistico. Chi vuole numeri suoi usa `taratura.py` e sostituisce
quattro valori in cima allo script.

**L'orientamento della costa dei punti derivati non è verificato da nessuno.**
È l'unico dato che, sbagliato di 180°, capovolge il verdetto. Per questo i punti
automatici lo dichiarano in pagina e invitano a controllare il quadrante: la
metà scura deve stare dove c'è il mare. Nei casi ambigui — punte, istmi, baie
strette — l'app dice che il valore è poco affidabile invece di tacere.

**Le insidie locali non sono derivabili.** Correnti alle foci, scogli,
rientranze: per i quattro spot catalogati c'è una riga scritta a mano, per i
punti automatici c'è scritto che nessuno le ha controllate.

**E il verde non è un permesso.** È scritto nel piede della pagina e va lasciato
lì: è una stima da modelli meteo su un punto a 3 km dalla riva, le condizioni
cambiano in minuti e da un punto all'altro della stessa spiaggia.

---

## Contribuire uno spot

Il modo utile non è mandare coordinate: quelle l'app le deriva da sola. È
mandare **le due cose che non si possono calcolare**: la conferma che
l'orientamento derivato sia giusto, e la riga di insidie locali di chi conosce
quella riva.

E se pagai regolarmente, la cosa più preziosa sono le tracce GPS con andata e
ritorno sulla stessa linea: sono quelle che trasformano un modello fondato in un
modello tarato.
