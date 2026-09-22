#!/usr/bin/env python3
"""Final, verified code -> GenBank accession table for every one of the 79
align/PEPC.fasta tips (Abildgaardia excluded -- no real accession exists).

Single-accession species: matched to Besnard et al. 2009's Supplementary
Table S1 directly (exact protein match). Multi-accession species:
disambiguated by best-offset local protein identity against each candidate
accession's own GenBank `/translation`. Every one of the 78 resolved codes
below checks out at 100% coverage in build_cds.py's own difflib-based
verification (`autojunk=False` -- see align_cds/verify.tsv). The one
substantive caveat:

  - Eleocharis limosa's three accessions (FM208022/23/24) are >99% identical
    to each other (a handful of point differences among near-clonal
    individuals), so which physical accession lands under which existing
    Ele_lim2/lim3/limo code is nominal -- any permutation within this trio
    is equally valid biologically, since all three checks pass regardless.

Keyed by the original PCOC/ConDor short codes (this is the anchor table,
matched by hand against Besnard 2009's Table S1) -- build.py derives real
species-name tip labels from this via _species_names(), rather than the
other way around, so this table never needs to change when tip names do.
"""
CODE_ACCESSION = {
    "Abildgaar": None,  # unresolved, no real GenBank accession -- see README
    "Actinosch": "FM207990",
    "Baumea": "FM207991",       # = Machaerina articulata (see README)
    "Blysmus": "FM207992",
    "Bolboscho": "FM207993",
    "Bulbostyl": "FM207994",
    "Carex_ber": "FM207995",
    "Carex_com": "FM207996",
    "Carex_hal": "FM207997",
    "Carex_pen": "FM207998",
    "Carpha": "FM207999",
    "Chrysithr": "FM208000",
    "Cladium": "FM208001",
    "Coleochlo": "FM208002",
    "Cyp_alt3": "FM208003",
    "Cyp_capi": "FM208004",
    "Cyp_era1": "FM208065",
    "Cyp_era6": "FM208066",
    "Cyp_fusc": "FM208005",
    "Cyp_iria": "FM208064",
    "Cyp_long": "FM208006",
    "Cyp_papy": "FM208007",
    "Cyp_pulc": "FM208008",
    "Cyp_rotu": "FM208009",
    "Cyp_spha": "FM208010",
    "Cyp_ust2": "FM208012",
    "Cyp_ustu": "FM208011",
    "Ele_acut": "FM208013",
    "Ele_bal2": "FM208015",
    "Ele_bal3": "FM208017",
    "Ele_bal4": "FM208016",
    "Ele_bald": "FM208014",
    "Ele_fici": "FM208018",
    "Ele_geni": "FM208019",
    "Ele_gra2": "FM208021",
    "Ele_grac": "FM208020",
    "Ele_lim2": "FM208023",  # nominal, see docstring
    "Ele_lim3": "FM208024",  # nominal
    "Ele_limo": "FM208022",  # nominal
    "Ele_pal2": "FM208026",
    "Ele_palu": "FM208025",
    "Ele_quan": "FM208027",
    "Ele_rost": "FM208028",
    "Ele_viv2": "FM208030",
    "Ele_vivA": "AB085948",
    "Ele_vivi": "FM208029",
    "Eriophor": "FM208031",
    "Fimb_di2": "FM208033",
    "Fimb_dic": "FM208032",
    "Fimb_fe2": "FM208034",
    "Fimb_fer": "FM208035",
    "Fimb_li2": "FM208037",
    "Fimb_lit": "FM208036",
    "Fuir_abn": "FM208038",
    "Fuir_umb": "FM208039",
    "Hellmut1": "FM208040",
    "Hellmut2": "FM208041",
    "Isolepis": "FM208042",
    "Killinga": "FM208043",
    "Machaeri": "FM208044",   # = Machaerina scirpoidea (see README)
    "Microdra": "FM208045",
    "Pycreus": "FM208046",
    "Remirea": "FM208047",
    "Rhy_alba": "FM208048",
    "Rhy_albi": "FM208049",
    "Rhy_glo2": "FM208051",
    "Rhy_glob": "FM208050",
    "Rhy_grac": "FM208052",
    "Rhy_rubr": "FM208053",
    "Scho_lac": "FM208055",
    "Scho_muc": "FM208056",
    "Scho_val": "FM208057",
    "Schoenox": "FM208058",
    "Schoenus": "FM208054",
    "Scirpoid": "FM208059",
    "Scirpus": "FM208060",
    "Uncin_ph": "FM208061",
    "Uncin_un": "FM208062",
    "Volkiell": "FM208063",
}
