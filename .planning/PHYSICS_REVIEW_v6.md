---
status: complete
reviewed: 2026-07-24
scope: independent physics review of ETNA (~6 weeks after PHYSICS_AUDIT_v5) — synthesis + next-steps plan
reviewer: Fable (synthesis), 5x Sonnet sub-agents (evidence gathering), Opus advisor (verification)
---

# Physics Review v6 — ETNA Simulator

Requested as "una review fable di tutta la fisica e i risultati fisici" with next-steps
planning. Process: 5 Sonnet sub-agents independently gathered file:line-cited evidence
across five physics domains (electrostatics/C-V, CCE/dark-current,
radiation-damage/FLASH, microdosimetry/MC-coupling, 2D-structures/honesty-audit) →
a Fable agent synthesized an independent physics judgment and a phased plan → Opus
advisor verified the synthesis and flagged one unverified claim, which was then
reproduced live (see Addendum). This is a deliberate follow-up to
`PHYSICS_AUDIT_v4.md`/`v5.md` (2026-06-13/14): it does not repeat those audits, it
assesses what changed since (~6 weeks) and what's newly found.

---

## Valutazione fisica indipendente (Fable)

Il quadro generale è quello di un simulatore la cui **fisica di base è solida** —
elettrostatica, drift-diffusion, raccolta di carica, danno da radiazione strutturale —
ma con **due categorie di problemi che non sono equivalenti in gravità** e che gli
audit precedenti tendono a trattare come fossero sullo stesso piano: (1) fisica
mancante per limiti di dati esterni (kappa, NIEL, Ramo), che è un blocco onesto e già
ben documentato nel codice; (2) **rappresentazione documentale che diverge dalla
fisica sottostante** — e qui trovo i due elementi nuovi più seri dell'intero ciclo di
audit.

**Dominio 1-2 (elettrostatica, DD, raccolta di carica): pronto per pubblicazione con
una correzione di etichettatura.** `README.md:19-20` e
`DESIGN-1-dosimetry-pn.md:118-119` dichiarano "R²=0.998 vs C-V misurato" come fatto
validato, senza alcun cenno alla circolarità (gli stessi tre punti sperimentali
calibrano il profilo di doping E fungono da target di validazione) né alla verifica
indipendente reale (Bruzzi & Verroi 2023, 2% di scarto) che invece esisterebbe ed è la
prova migliore disponibile. Questo non è un errore di fisica — la DD e la Poisson
equation sono corrette contro Sze & Ng, 54/54 test passano — è un **errore di standard
epistemico interno incoerente**: lo stesso README copre correttamente la corrente di
buio ("calibrazione, non predizione") ma non il C-V. Mj-4 (Ramo/displacement current
assente) è confermato aperto ma limitato a I(t)/tempi di raccolta — non tocca Q
raccolta né CCE, quindi non inficia i due design deliverable.

**Dominio 3 (danno da radiazione, FLASH): onestamente disclosurato, in buono stato.**
La correzione EH6/7 (sigma_n) è ben fatta e ben commentata come "trascrizione fedele
ma fisicamente implausibile". Il FLASH claim è stato de-escalato coerentemente in
tutti i file controllati. L'unico residuo è cosmetico: la frase nell'intro del
notebook 04 non allineata col resto del notebook.

**Dominio 4 (microdosimetria): qui il giudizio deve essere più severo di quanto
indicato da "ancora bloccato dai dati".** Il fatto che kappa fittizio (~0.58, invertito
rispetto al valore fisico atteso ~1.13-1.24) sia usato come _default silenzioso_ di
`tissue_equivalence_correction()` era già noto. Nuovo e più grave: il notebook 18,
cella 13, contiene un output eseguito e committato che stampa una "Microdosimetric
Quantities Summary" con `y_F Tissue-Eq. = 12.041 keV/um` e `y_D` corrispondente, a
quattro cifre significative, **senza una sola parola "placeholder" o "fabricated" in
tutto il notebook**. Chiunque apra questo notebook — un revisore di paper, un
collaboratore, la stessa Giada tra sei mesi — vede numeri di aspetto pubblicabile che
sono fisicamente sbagliati e con il segno invertito. Non è un bug di codice; è un
artefatto di presentazione che ha ereditato silenziosamente la fabbricazione. La UI
Streamlit è pulita (non chiama questa funzione), quindi il rischio è isolato a uso
diretto della libreria e ai notebook — ma i notebook sono esattamente il tipo di file
che finisce allegato a un paper o mostrato a un collaboratore.

> **Correzione 2026-10-21 (testo sopra lasciato com'era).** Il kappa da 0.58 era
> fabbricato ma **non invertito**. I dati reali NIST PSTAR ora nel repo
> (`data/nist_star/README.md`, `docs/USAGE.md`) danno il kappa per l'energia lineale
> (rapporto LET SiC -> acqua) pari a 0.426 a 1 MeV, in calo fino a 0.375 a 150 MeV. Il
> valore "~1.13-1.24" citato qui è un rapporto di stopping power **massico** (conversione
> di dose, per g/cm²), non il fattore per y; e anche come rapporto massico i dati danno
> 1.37 (1 MeV) -> 1.20 (150 MeV). La tabella piatta 0.58 era quindi fabbricata e troppo
> alta (circa 1.4-1.5x), ma dalla stessa parte di 1 del rapporto lineale reale. Stato
> aggiornato degli item in fondo al documento.

**Dominio 5 (2D, notebook, honesty audit): le fix strutturali (Mj-1/2/3) sono
solide** — un test si chiama letteralmente `test_deceptive_undepleted_ranks_below_valid`.
Ma i notebook 19/20 non sono stati toccati dalla campagna di fix fisici, solo dal
rename meccanico petringa→etna, e contengono numeri di full-depletion voltage (50V)
che il codice live smentisce ora (~10.5V a 10um epi). Meno grave del caso kappa
(errore "troppo conservativo" non "invertito"), ma comunque un disallineamento fisico
visibile.

---

## Findings classificati per priorità (severità × probabilità di fuorviare qualcuno)

1. **Notebook 18/19: kappa fabbricato/invertito presentato come output pubblicabile,
   zero disclaimer.** Il caso peggiore del repo: dato fisicamente sbagliato (segno
   invertito), formattato come tabella di risultati finale. Aggravante:
   `scripts/create_notebook_18.py:315-341` genera prosa che afferma come fatto la
   fisica invertita. Priorità massima, azione immediata, basso costo.
   _Correzione 2026-10-21:_ fabbricato sì, invertito no (vedi la nota al Dominio 4);
   risolto il 2026-10-20 con dati NIST reali (C-1, f6d7ade, merge 15b736f).

2. **README.md:19-20 e DESIGN-1 §6: claim di validazione C-V circolare presentato
   senza hedge**, mentre lo stesso documento hedgia esplicitamente altrove (dark
   current). Rischio: è la vetrina principale del progetto e il documento che va nelle
   mani del fab. Fix a costo quasi zero: una frase, più citare Bruzzi & Verroi come
   prova indipendente reale (che esiste ed è positiva).

3. **Legacy kappa come default silenzioso** (`microdosimetry.py:392-407`,
   `kappa_constant=0.58`). Root cause del punto 1: finché il default restituisce il
   numero fabbricato senza errore, ogni futuro notebook/script ripete l'errore. Fix di
   codice quasi banale, non bloccato dai dati esterni.

4. **Notebook 20: numeri di full-depletion voltage obsoleti** (50V vs ~10.5V reale) e
   raccomandazione guard-ring basata su numeri inseriti a mano, non da simulazione live
   nella cella. Severità media (errore conservativo, non ingannevole in senso
   opposto).

5. **Mj-4 (Ramo current assente), derivata Newton hardcoded a zero in
   `dark_current.py:341-348`, NIEL placeholder.** Correttamente minori/data-blocked.
   Non toccano Q raccolta, CCE, o convergenza (solo velocità/I(t) transitorio).

**Dissenso rispetto ai sub-agent:** il Dominio 4 dice giustamente "presentation-blocking,
peggio di data-blocked" — concordo e lo rimarco: non è "manca un dato", è "un dato
sbagliato è già stato pubblicato con l'aspetto di un dato buono". Questo item va sopra
qualunque altro nel lavoro immediato, sopra anche Mj-4, che potrebbe sembrare
prioritario solo perché è "un bug di fisica vero" — non lo è nel senso che conta (non
tocca CCE/Q).

---

## Piano dei prossimi passi fisici

**(a) Fix documentali/honesty — eseguibili oggi da un agente, priorità immediata:**

- Notebook 18/19: warning esplicito prima/dopo l'uso di `compute_kappa_table`/
  `tissue_equivalence_correction` senza `source="bragg"`; correggere la prosa di
  `scripts/create_notebook_18.py:315-341` che afferma l'inverso come fatto.
- README.md:19-20 e DESIGN-1 §6: aggiungere la clausola di circolarità e citare Bruzzi
  & Verroi 2023 (2% scarto) come prova indipendente reale.
- Notebook 04: correggere la frase intro (riga 14) per coerenza col caveat già presente
  altrove nello stesso notebook.
- Notebook 20: aggiornare i numeri di full-depletion voltage stantii (50V→~10.5V) e
  propagare il caveat "guard-ring non modellato / numeri inseriti a mano".
- `poisson.py:317`: correggere il docstring (dice 1%, il codice usa 50%).

**(b) Lavoro codice/fisica reale, fattibile ora, nessun blocco dati:**

- Rendere `kappa_constant` default non silenzioso (richiedere `source="bragg"`
  esplicito con eccezione, o warning non ignorabile).
- SRV Newton derivative in `dark_current.py:341-348`: implementare la derivata
  analitica al posto dello zero hardcoded (qualità numerica, basso rischio).
- Mj-4 (Ramo/displacement current): fattibile ma non urgente — da programmare solo se
  emerge un caso d'uso concreto su timing/pulse shape.
- **Nuovo, minore (vedi Addendum):** `optimization.py:148` — `center_cce`/`edge_cce`
  possono essere `None` (non solo NaN) per configurazioni non depletabili a bassa
  tensione, causando `TypeError` in `center_cce >= CCE_FLOOR` invece di un fallimento
  gestito. Il gate di validità (`is_valid=False`) è comunque fisicamente corretto
  quando la scan fallisce del tutto (via l'except-block), quindi non è un bug di
  fisica — è un buco di robustezza minore da chiudere con un `None`-check esplicito.

**(c) Bloccato su dato esterno reale — non schedulabile come lavoro-agente:**

- Stopping power PSTAR/SRIM (acqua/Si/C) → sblocca kappa reale.
- SR-NIEL hardness factors → sblocca fattori di danno da radiazione realistici.
- Fermi da oltre un mese per mancanza di accesso rete in questo ambiente: un umano
  deve scaricare i CSV una volta e depositarli in `data/srim/`; l'integrazione
  successiva è agente-eseguibile in un ciclo singolo.

**(d) Da NON fare ora, e perché:**

- Non prioritizzare Mj-4 (Ramo current) — non altera Q/CCE, distoglierebbe tempo da
  fix a più alto rischio-di-fraintendimento.
- Non rincorrere fisica FLASH "vera" (screening, trasporto ambipolare) in questo
  ciclo — correttamente derisked come esplorativo, nessun consumer chiaro identificato.
- Non estendere l'honesty-audit a tutti i 22 notebook esaustivamente subito — 18/19/20/04
  sono i portatori concreti noti; uno scan completo ha basso ritorno marginale ora.

**In sintesi:** i due deliverable fab-ready (DESIGN-1/2) restano fisicamente
difendibili — non si revoca "done" — ma **DESIGN-1 è "done dopo il fix di onestà C-V
in un rigo"**, non "done senza riserve": va reso esplicito a chiunque lo legga come
pronto per il fab senza ulteriori controlli. Prima di mostrare il repo a un revisore
esterno o allegare notebook a un paper, i due item in cima alla lista (kappa
fabbricato nei notebook, claim C-V circolare) vanno sistemati: costo bassissimo,
rischio reputazionale sproporzionato se scoperti da altri.

---

## Addendum — verifica Opus advisor

L'advisor (Opus) ha approvato il flusso Sonnet→Fable→Opus e la sostanza della review,
ma ha bloccato la finalizzazione su un punto non verificato: nel run del notebook 20
allegato ai findings del dominio 5, 34/36 configurazioni dello sweep 2D di
microdosimetria fallivano con `DEVSIM FATAL: Solver "custom" specified, but
"solver_callback" not set`, ogni volta subito dopo un "Resetting DEVSIM". L'ipotesi
dell'advisor: potrebbe essere un vero bug di regressione nel percorso di reset di
`reset_devsim_fully` (usato dallo sweep di ottimizzazione 2D dietro DESIGN-2), mascherato
dal fatto che i test di `optimization.py` verificano solo il ranking analitico
(`full_depletion_voltage_graded`/`_rank_sweep_results`), non il DD solve stesso.

**Verifica eseguita:** rieseguito `microdosimetric_sweep()` in un processo Python
pulito su 4 configurazioni (stesso spazio parametri del notebook, sottoinsieme). Le 4
configurazioni hanno completato tutti i solve DD (equilibrio + lateral scan a 5 punti)
con convergenza pulita — **zero occorrenze** di `solver_callback not set`. Una
configurazione ha invece incontrato un errore diverso e reale (§b sopra, punto
`None`-check), non correlato al problema segnalato dall'advisor.

**Conclusione:** il fallimento `solver_callback not set` osservato nel notebook 20 **non
si riproduce** sul codice corrente in un processo pulito — non è una regressione live
di `reset_devsim_fully`. È coerente con uno stato di sessione ipickernel/devsim
persistente, specifico del notebook (probabilmente legato al fatto che il notebook
esegue anche celle precedenti con costrutti di device diversi/cilindrici prima dello
sweep, in un unico kernel di lunga durata), non con un bug nella libreria. La
classificazione di questo item resta quella di Fable: parte del disallineamento
notebook-20 già coperto in (a), nessuna azione di codice richiesta in (b) per questo
punto specifico.

---

## Aggiornamento di stato (2026-10-21)

Aggiunto a posteriori; il testo sopra resta com'era. Fonti: `data/nist_star/README.md`,
`data/niel/`, `docs/USAGE.md`, `.issues/TECH-DEBT-2026-10-18.md`,
`.issues/audit-v7/issues/{A,B,D}-INDEX.md`. Dettaglio in `PHYSICS_AUDIT_v5.md` §10.

- **Correzione kappa tessuto-equivalente:** fabbricato, non invertito. Rapporto LET
  SiC -> acqua (fattore per y) 0.426 @1 MeV, 0.391 @10, 0.378 @62, 0.375 @150 MeV;
  rapporto massico 1.368 -> 1.255 -> 1.215 -> 1.204. Il "~1.13-1.24" era un rapporto
  massico non verificato, che non coincide con nessuno dei due.
- **C-1** (dati kappa): risolto (fix 2026-10-20, merge 2026-10-21), NIST PSTAR/ASTAR
  reali, SiC per additività di Bragg (f6d7ade, merge 15b736f). Nessun default silenzioso
  di `kappa_constant` già dalla fix di onestà v6.
- **C-2** (kappa indicizzato all'energia depositata): risolto 2026-10-20, lookup
  all'energia cinetica, eventi MC con energia cinetica primaria (5ba7c78).
- **Mj-5** (profili di generazione): alpha = curva di Bragg NIST ASTAR di default;
  protoni: dipendenza dall'energia PSTAR opt-in (B-10), modo di default invariato
  (8acb0ec, merge 15b736f).
- **Mj-4** (Ramo/displacement current): risolto 2026-10-20, corrente ai terminali =
  conduzione + spostamento; CCE e Q invariati (eba0801, merge d958408).
- **C-5** (NIEL): parzialmente risolto (fix 2026-10-20, merge 2026-10-21), dati SR-NIEL-7 per protoni
  1-1000 MeV nella convenzione Si 1-MeV-neq (kappa(62 MeV) = 1.942; B-07, merge
  1142e9c). **Aperto:** normalizzazione ai neutroni nel SiC (serve un account SR-NIEL
  registrato o NJOY).
- **FLASH:** fisica di plasma / alta iniezione ancora non implementata (bound di
  sensibilità). Termine Auger: era inattivo nello sweep stazionario, corretto
  2026-10-20 (afbf055), e nello sweep transiente, corretto 2026-10-21 (49cb003, merge
  bf6e970).
- **audit-v7:** B-07, B-08, B-10, B-11 fixed; D-15 closed con follow-up fixed; A-02
  mitigato **opt-in** (`DeviceConfig(incomplete_ionization=True)`, default OFF,
  calibrazione con il modello attivo ancora aperta).
