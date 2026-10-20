"""Fetch proton NIEL in SiC from the SR-NIEL web calculator (provenance script).

Usage:  uv run python data/niel/fetch_sr_niel.py <Ed_Si_eV> <Ed_C_eV> <tag>
        uv run python data/niel/fetch_sr_niel.py --si <Ed_Si_eV> <tag>   (pure Si target)

POSTs the public "Protons & Ions NIEL Calculator" form
(https://www.sr-niel.org/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator)
with: incident particle 1:H, model Hadron + Coulomb, hadronic contribution No Scaling,
target Compound (IV) SiC (Si:C = 1:1), the given displacement thresholds, energy range
0.1-1000 MeV plus the extra energies below. Saves the raw result page
(sr_niel_proton_sic_<tag>.html) and the parsed table (sr_niel_proton_sic_<tag>.csv).
Needs network; nothing in the package calls this at run time.
"""

import csv
import re
import sys
import urllib.parse
import urllib.request
from http.cookiejar import CookieJar
from pathlib import Path

FORM = (
    "https://www.sr-niel.org/index.php/sr-niel-web-calculators/"
    "niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator"
)
EXTRA = [
    1,
    2,
    3,
    4,
    5,
    5.6,
    6,
    8,
    10,
    12,
    15,
    20,
    25,
    30,
    40,
    50,
    62,
    70,
    80,
    100,
    120,
    150,
    200,
    250,
]
ROW = re.compile(
    r"(\d\.\d{4}e[+-]\d\d)\s+(\d\.\d{4}e[+-]\d\d)\s+\d\.\d{4}e[+-]\d\d"
    r"\s+\d\.\d{4}e[+-]\d\d"
)


def main(ed_si, ed_c, tag, material="sic"):
    here = Path(__file__).resolve().parent
    opener = urllib.request.build_opener(
        urllib.request.HTTPCookieProcessor(CookieJar())
    )
    opener.open(FORM, timeout=60).read()
    data = {
        "Zi": "1",
        "Ionmodel": "1",
        "Scale": "0",
        **(
            # form order for the SiC compound: element 1 = 14:Si, element 2 = 6:C
            {"TARGET": "1", "COMP": "1", "Th1": str(ed_si), "Th2": str(ed_c)}
            if material == "sic"
            # user-defined single-element target: Si, stoichiometric index 1
            else {"TARGET": "0", "NELEM": "1", "Z1": "14", "Stoich1": "1",
                  "Th1": str(ed_si)}
        ),
        "Emin": "0.1",
        "Emax": "1000",
        "Fluence": "1",
        "AddEnergy": "\n".join(str(e) for e in EXTRA),
        "formSubmit": "CALCULATE",
    }
    body = urllib.parse.urlencode(data).encode()
    page = opener.open(FORM, body, timeout=120).read().decode()
    m = re.search(r"window\.open\('(http[^']+niel_p\.html)'", page)
    if not m:
        raise SystemExit("no result link in the response")
    url = m.group(1).replace("http://", "https://")
    html = opener.open(url, timeout=60).read().decode()
    (here / f"sr_niel_proton_{material}_{tag}.html").write_text(html)
    text = re.sub(r"<[^>]+>", " ", html)
    header = re.search(r"Ions&Protons:[^\n]*", text)
    thr = re.search(r"Displacement Threshold Energies \[eV\]:\s*([^\n]*)", text)
    rows = ROW.findall(text)
    with open(here / f"sr_niel_proton_{material}_{tag}.csv", "w", newline="") as f:
        f.write(
            f"# SR-NIEL proton NIEL in {material} (Hadron+Coulomb, no scaling); "
            f"{header.group(0).strip() if header else ''}\n"
        )
        f.write(
            "# Displacement thresholds Si; C [eV]: "
            f"{thr.group(1).strip() if thr else ''}\n"
        )
        f.write(f"# source result page: {url}\n")
        w = csv.writer(f)
        w.writerow(["energy_MeV", "niel_MeV_cm2_per_g"])
        for e, n in rows:
            w.writerow([e, n])
    print(url, len(rows), "rows")


if __name__ == "__main__":
    if sys.argv[1] == "--si":
        main(float(sys.argv[2]), None, sys.argv[3], material="si")
    else:
        main(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3])
