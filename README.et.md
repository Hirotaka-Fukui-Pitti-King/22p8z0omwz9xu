# Teadusartiklite reprodutseerimine

**Keel:** [한국어](README.md) | Eesti

See repositoorium koondab erinevate teadusartiklite empiiriliste analüüside, teoreetiliste mudelite ja arvutuslike tulemuste reprodutseerimiseks kasutatavat koodi ning sellega seotud materjale.

Eesmärk on mõista teadusartiklites kasutatud meetodeid, sõltumatult kontrollida nende peamisi tulemusi ning edendada reprodutseeritavat ja läbipaistvat teadustööd.

## Eesmärgid

- **Tulemuste reprodutseerimine (Replication):** Teadusartiklite peamiste tulemuste reprodutseerimine algsetes uuringutes kasutatud meetodite ja andmete põhjal.
- **Metoodika mõistmine (Methodological Understanding):** Ökonomeetriliste hindamismeetodite, struktuursete mudelite ja arvutuslike algoritmide rakendamise õppimine.
- **Tulemuste kontrollimine (Verification):** Reprodutseeritud tulemuste võrdlemine avaldatud tulemustega ning võimalike erinevuste põhjuste dokumenteerimine.
- **Uuringute laiendamine (Extension):** Alternatiivsete eelduste, hindamismeetodite ja täiendavate eksperimentide uurimine.

## Repositooriumi struktuur

Iga teadusartiklit käsitletakse eraldi kataloogis. Võimaluse korral kasutatakse järgmist struktuuri:

```text
.
├── README.md
├── README.et.md
├── papers/
│   ├── paper-01/
│   │   ├── README.md
│   │   ├── code/
│   │   ├── data/
│   │   ├── results/
│   │   └── figures/
│   ├── paper-02/
│   └── ...
└── LICENSE
```

Iga artikli kataloog sisaldab võimaluse korral bibliograafilist teavet, reprodutseerimise eesmärke, kasutatud tarkvarakeskkonna kirjeldust, käivitamisjuhiseid ja tulemuste võrdlust.

Kui algandmete kasutamine on piiratud, ei avaldata neid repositooriumis. Selle asemel esitatakse viited ametlikele andmeallikatele ja juurdepääsutingimustele.

## Reprodutseerimise staatus

Iga artikli reprodutseerimise edenemist dokumenteeritakse vastava kataloogi failis `README.md`.

| Staatus | Kirjeldus |
|---|---|
| Not Started | Reprodutseerimist ei ole veel alustatud |
| In Progress | Koodi rakendamine või tulemuste kontrollimine on pooleli |
| Replicated | Peamised tulemused on edukalt reprodutseeritud |
| Partially Replicated | Reprodutseeritud on ainult osa tulemustest |
| Not Replicated | Peamisi tulemusi ei õnnestunud reprodutseerida |

Reprodutseerimise edukust hinnatakse algses artiklis avaldatud ja tegelikult saadud tulemuste võrdlemise põhjal. Kontrollimata tulemusi ei märgita edukalt reprodutseerituks.

## Kasutatavad töövahendid

Sõltuvalt uuritavast artiklist võidakse kasutada järgmisi programmeerimiskeeli ja analüüsivahendeid:

- Julia
- Python
- R
- MATLAB
- Stata

Tarkvara versioonid, vajalikud paketid ja käivitamisjuhised dokumenteeritakse iga projekti juures eraldi.

## Olulised märkused

See repositoorium on mõeldud akadeemiliseks õppe- ja uurimistööks ning teadustulemuste sõltumatuks kontrollimiseks.

Algsete teadusartiklite autoriõigused kuuluvad nende autoritele ja vastavatele kirjastajatele. Algse koodi ja andmete kasutamisel järgitakse kohaldatavaid litsentse ning kasutustingimusi.

Kui ei ole märgitud teisiti, ei ole siin avaldatud reprodutseerimiskood algsete autorite koostatud ega ametlikult kinnitatud. Reprodutseeritud tulemuste tõlgendamise eest vastutab selle repositooriumi haldaja.

## Viitamine

Selle repositooriumi materjalide kasutamisel tuleks eelkõige viidata vastavatele algsetele teadusartiklitele.

Üksikasjalik bibliograafiline teave on esitatud iga projekti kataloogis.
