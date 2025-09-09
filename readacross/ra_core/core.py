#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import subprocess
import os
import re
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdmolops, inchi
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
import sygma
import requests
from urllib.parse import quote
import time
from rdkit.DataStructs import TanimotoSimilarity
from typing import Iterable, Dict, List, Set, Union
import xlsxwriter



# In[ ]:





# In[2]:


# --- 1. SETUP ---
RDLogger.DisableLog('rdApp.*')

from pathlib import Path
import os, shutil

# --- Resolve repo root: readacross/ (ra_core is inside this folder)
APP_ROOT = Path(__file__).resolve().parents[1]   # .../readacross
TOOLS_DIR = APP_ROOT / "tools"

# JVMs (env overrides allowed)
JAVA_BIN_DEFAULT  = shutil.which("java") or "/usr/bin/java"
JAVA_BIN_TOXTREE  = os.environ.get("JAVA_BIN_TOXTREE", "/opt/java/temurin-11/bin/java")
JAVA_BIN_BT       = os.environ.get("JAVA_BIN_BIOTRANSFORMER", JAVA_BIN_DEFAULT)

# Heaps (MB)
TT_HEAP_MB = int(os.environ.get("TT_HEAP_MB", "512"))   # Toxtree
BT_HEAP_MB = int(os.environ.get("BT_HEAP_MB", "1024"))  # BioTransformer

def _find_first(root: Path, pattern: str) -> Path:
    """Find first file that matches pattern anywhere under root."""
    hits = sorted(root.rglob(pattern))
    if not hits:
        # Log tree to help debugging if missing
        print(f"[JAR RESOLVE] No match for '{pattern}' under {root}")
        try:
            for p in sorted(root.glob("**/*"))[:200]:
                pass  # just force traversal for build logs if needed
        except Exception as e:
            print("[JAR RESOLVE] listing failed:", e)
        raise FileNotFoundError(f"No file matching '{pattern}' under {root}")
    return hits[0]

def _ensure_exists(p: Path, label: str):
    print(f"[JAR RESOLVE] Using {label}: {p}")
    if not p.exists():
        raise FileNotFoundError(f"{label} not found at {p}")

# --- BioTransformer: find JAR and working directory (must contain config.json)
BT_JAR = _find_first(TOOLS_DIR / "biotransformer", "BioTransformer*.jar")
_ensure_exists(BT_JAR, "BioTransformer JAR")

# Prefer the JAR's folder; if it doesn't have config.json, try parent
if (BT_JAR.parent / "config.json").exists():
    BT_DIR = BT_JAR.parent
elif (BT_JAR.parent.parent / "config.json").exists():
    BT_DIR = BT_JAR.parent.parent
else:
    raise FileNotFoundError(
        f"BioTransformer config.json not found near {BT_JAR}. "
        f"Checked: {BT_JAR.parent} and {BT_JAR.parent.parent}"
    )
print(f"[JAR RESOLVE] BioTransformer working dir: {BT_DIR}")

# --- Toxtree: find JAR and working directory (must contain ext/index.properties)
TX_JAR = _find_first(TOOLS_DIR / "toxtree", "Toxtree*.jar")
_ensure_exists(TX_JAR, "Toxtree JAR")

# Prefer the JAR's folder; if it doesn't have ext/index.properties, try parent
if (TX_JAR.parent / "ext" / "index.properties").exists():
    TX_DIR = TX_JAR.parent
elif (TX_JAR.parent.parent / "ext" / "index.properties").exists():
    TX_DIR = TX_JAR.parent.parent
else:
    raise FileNotFoundError(
        f"Toxtree ext/index.properties not found near {TX_JAR}. "
        f"Checked: {TX_JAR.parent} and {TX_JAR.parent.parent}"
    )
print(f"[JAR RESOLVE] Toxtree working dir: {TX_DIR}")

dart_alerts_dictionary = {
    "Category 1: Inorganics and Derivatives": {
        "1.a.1": {
            "tier": 1,
            "name": "Metals and Organometallics",
            "description": "Identifies chemicals containing specific metals and metalloids (As, B, Pb, Hg, Sn, V, etc.) cited as having DART activity.",
            "smarts": "[As,B,Mn,Cr,Zn,Te,Al,Cd,Cu,Ni,Pb,Hg,Sn,V,Co,Ga,Li]"
        },
        "1.b.1a": {
            "tier": 1,
            "name": "Organophosphorus (Core)",
            "description": "Identifies the core phosphoryl (P=O) or thiophosphoryl (P=S) bond common to organophosphorus compounds.",
            "smarts": "[P;v5]=O"
        },
        "1.b.1b": {
            "tier": 1,
            "name": "Organophosphorus (Core)",
            "description": "Identifies the core phosphoryl (P=O) or thiophosphoryl (P=S) bond common to organophosphorus compounds.",
            "smarts": "[P;v5]=S"
        },
        "1.b.2": {
            "tier": 2,
            "name": "Organophosphorus - Cyclophosphamide-like (Key Feature)",
            "description": "Identifies the highly potent sub-class of cyclophosphamide-like mustards containing a bis(2-chloroethyl)amino group.",
            "smarts": "[P;v5](=O)N(CC[Cl])CC[Cl]"
        },
        "1.c.1": {
            "tier": 1,
            "name": "Organosiloxane - Phenyl Siloxane (Core)",
            "description": "Identifies the core feature of a phenyl group directly attached to a silicon atom that is part of a siloxane (-Si-O-) backbone.",
            "smarts": "[Si](c1ccccc1)O"
        }
    },
    "Category 2: Estrogen and Androgen Receptor Binders": {
        "2.a.1": {
            "tier": 1,
            "name": "Steroid Nucleus (Core)",
            "description": "Identifies the general fused tetracyclic steroid nucleus.",
            "smarts": "[#6]12[#6]3[#6]4[#6][#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6][#6]4"
        },
        "2.a.2": {
            "tier": 2,
            "name": "Steroid - Estradiol-Like (Key Feature)",
            "description": "Identifies the key feature of estrogens: an aromatic A-ring with oxygens at C-3 and C-17.",
            "smarts": "c1ccc2c(c1[O;H1,R0])CC[C@H]3[C@H]2CC[C@H]4[C@@]3([#6])[#6][#6][C@@H]([O;H1,R0])C4"
        },
        "2.a.3": {
            "tier": 2,
            "name": "Steroid - Glucocorticoid (Key Feature)",
            "description": "Identifies the key features of glucocorticoids: a C-3 keto, a C-11 oxygen (keto or hydroxyl), and a ketol/chloro-ketone side chain at C-17.",
            "smarts": "[#6]12[C@H]3[C@H]([#6;!a]([#6;!a]=O)[#6;!a][C@]4([#6])[C@H]1[C@H](C(=O)C([#8,#17]))C[C@H]4[#6])[C@@H](O)[#6;!a][#6;!a]23"
        },
        "2.a.4": {
            "tier": 2,
            "name": "Steroid - Progestin-Like (Key Feature)",
            "description": "Identifies the key features of progestins: a C-3 keto and an acetyl side chain at C-17.",
            "smarts": "[#6]12[C@H]3[C@H]([#6;!a]([#6;!a]=O)C[C@]4([#6])[C@H]1[C@H](C(=O)C)CC4)[#6;!a][#6;!a][#6;!a]23"
        },
        "2.a.5": {
            "tier": 2,
            "name": "Steroid - Androgen-Like (Key Feature)",
            "description": "Identifies the key features of androgens: a C-3 keto and a hydroxyl group at C-17.",
            "smarts": "[#6]12[C@H]3[C@H]([#6;!a]([#6;!a]=O)C[C@]4([#6])[C@H]1[C@@H](O)CC4)[#6;!a][#6;!a][#6;!a]23"
        },
        "2.b.1": {
            "tier": 1,
            "name": "Flavonoids (Core)",
            "description": "Identifies the phenylchromone scaffold of flavonoids and isoflavones.",
            "smarts": "c1oc(c2ccccc12)c3ccccc3"
        },
        "2.b.2": {
            "tier": 1,
            "name": "Mycoestrogens (Core)",
            "description": "Identifies the resorcylic acid macrocyclic lactone scaffold of mycoestrogens.",
            "smarts": "c1c(O)cc(O)c2c1C(=O)O[C;R1][C;R1][C;R1][C;R1][C;R1][C;R1][C;R1][C;R1][C;R1]2"
        },
        "2.b.3": {
            "tier": 1,
            "name": "Diaryl Endocrine Disruptors (Core)",
            "description": "Identifies the general scaffolds of diaryl endocrine disruptors like diphenylalkanes and stilbenes.",
            "smarts": "c1ccccc1[#6]~[#6]c2ccccc2"
        },
        "2.b.4": {
            "tier": 2,
            "name": "DDT-like (Key Feature)",
            "description": "Identifies the key feature of a 4,4'-dichloro diphenylalkane pattern.",
            "smarts": "c1cc(Cl)ccc1[#6][#6]c2ccc(Cl)cc2"
        },
        "2.b.5": {
            "tier": 2,
            "name": "BPA/DES-like (Key Feature)",
            "description": "Identifies the key feature of para-hydroxyl or para-amino groups on a diphenylalkane or stilbene scaffold.",
            "smarts": "c1cc([O,N])ccc1[#6]~[#6]c2ccc([O,N])cc2"
        },
        "2.b.6": {
            "tier": 1,
            "name": "4-Alkylphenols (Core)",
            "description": "Identifies the general 4-alkylphenol scaffold.",
            "smarts": "[OH]c1ccc(C)cc1"
        },
        "2.b.7": {
            "tier": 2,
            "name": "Long-Chain 4-Alkylphenols (Key Feature)",
            "description": "Identifies the key feature of a long (C7-C12) alkyl chain at the para position of a phenol.",
            "smarts": "[OH]c1ccc(CCCCCCC)cc1"
        },
        "2.b.8": {
            "tier": 1,
            "name": "N-Aryl Amides/Ureas (Core)",
            "description": "Identifies the general N-Aryl amide, urea, or carbamate scaffold.",
            "smarts": "[c]N([H,C])C(=O)"
        },
        "2.b.9": {
            "tier": 2,
            "name": "Non-Steroid - Polysubstituted Aryl (Key Feature)",
            "description": "Identifies aryl rings with multiple activating/deactivating groups. This alert finds an aryl ring with at least one halogen AND at least one nitrile.",
            "smarts": "c1(C#N)cccc(c1)[F,Cl]"
        }
    },
    "Category 3: RAR/AhR Binders and Prostaglandin Agonists": {
        "3.a.1": {
            "tier": 1,
            "name": "Retinoid - Ring and Polar Group (Core)",
            "description": "Identifies molecules with a retinoid-like ring system (aromatic or non-aromatic) and a required polar functional group.",
            "smarts": "[#6]1([#6])[#6][#6][#6]([#6])[#6]1.[$([#6]=O),$(C[O;H1]),$(S(=O)=O)]"
        },
        "3.a.2": {
            "tier": 2,
            "name": "Retinoid - Long Conjugated Chain (Key Feature)",
            "description": "Identifies the key feature of a long (>3 units) conjugated polyene chain attached to a polar group, associated with high potency.",
            "smarts": "[#6]=[#6]/[#6]=[#6]/[#6]=[#6]/C(=O)"
        },
        "3.b.1": {
            "tier": 2,
            "name": "TCDD-like - Halogenated Dibenzofuran",
            "description": "Identifies the dibenzofuran ring system with at least one chlorine or bromine atom.",
            "smarts": "c12oc3ccccc3c1c(c(c(c2)))[Cl,Br]"
        },
        "3.b.2": {
            "tier": 2,
            "name": "TCDD-like - Halogenated Dibenzo-p-dioxin",
            "description": "Identifies the dibenzo-p-dioxin ring system with at least one chlorine or bromine atom.",
            "smarts": "c12oc3oc(ccc3)c1c(c(c(c2)))[Cl,Br]"
        },
        "3.b.3": {
            "tier": 1,
            "name": "HAH - Polychlorinated/brominated Biphenyl (Core)",
            "description": "Identifies a biphenyl scaffold with four or more chlorine or bromine atoms.",
            "smarts": "c1c([Cl,Br])c([Cl,Br])c(c([Cl,Br])c1[Cl,Br])-c2ccccc2"
        },
        "3.b.4": {
            "tier": 1,
            "name": "HAH - Monochlorinated/brominated Biphenyl (Core)",
            "description": "Identifies a biphenyl scaffold with at least one chlorine or bromine atom.",
            "smarts": "c1c([Cl,Br])cccc1-c2ccccc2"
        },
        "3.b.5": {
            "tier": 1,
            "name": "PAH - 3-5 Fused Rings",
            "description": "Identifies fused aromatic systems of 3 to 5 rings.",
            "smarts": "a1a2a3a(a4aa1a2a34)a"
        },
        "3.b.6": {
            "tier": 1,
            "name": "Indole-Related Chemicals",
            "description": "Identifies the core indole ring system.",
            "smarts": "c1cccc2c1c[nH]c2"
        },
        "3.c.1": {
            "tier": 2,
            "name": "Prostaglandin Receptor Agonists",
            "description": "Identifies the core prostaglandin scaffold: a 5-membered ring with adjacent carboxylic acid and alcohol side chains.",
            "smarts": "[#6;R]1[#6;R]([#6]C(=O)O)[#6;R]([#6][#6](O))[#6;R][#6;R]1"
        }
    },
    "Category 4: nAChR Binders and AChE Inhibitors": {
        "4.a.1": {
            "tier": 2,
            "name": "nAChR Binder - Atropine-like (Key Feature)",
            "description": "Identifies the rigid azabicyclo[3.2.1]octane (tropane) ester scaffold.",
            "smarts": "[#6]1[#6][#6]2[#6][#7](C)[#6]1[#6]C2OC(=O)C(c1ccccc1)CO"
        },
        "4.a.2": {
            "tier": 2,
            "name": "nAChR Binder - Diphenylhydramine-like (Key Feature)",
            "description": "Identifies the classic diaryl pharmacophore connected to an amine via a short linker.",
            "smarts": "[#7]CC[O,C](c1ccccc1)c2ccccc2"
        },
        "4.a.3": {
            "tier": 1,
            "name": "nAChR Binder - Piperidine Alkaloid (Core)",
            "description": "Identifies the general piperidine ring.",
            "smarts": "C1CCNCC1"
        },
        "4.a.4": {
            "tier": 2,
            "name": "nAChR Binder - Pyridyl Piperidine (Key Feature)",
            "description": "Identifies the key feature of a 2-(pyridin-2-yl)piperidine scaffold found in potent piperidine alkaloids.",
            "smarts": "c1cnccc1C2CCCCN2"
        },
        "4.b.1": {
            "tier": 1,
            "name": "AChE Inhibitor - O-Aryl Carbamate (Core)",
            "description": "Identifies the core O-Aryl carbamate functional group.",
            "smarts": "[#7]C(=O)Oc1ccccc1"
        },
        "4.b.2": {
            "tier": 2,
            "name": "AChE Inhibitor - Naphthyl Carbamate (Key Feature)",
            "description": "Identifies the key feature of a carbamate of a naphthol (naphthalene ring).",
            "smarts": "[#7]C(=O)Oc1c2ccccc2ccc1"
        },
        "4.b.3": {
            "tier": 2,
            "name": "AChE Inhibitor - Benzofuran Carbamate (Key Feature)",
            "description": "Identifies the key feature of a carbamate of a benzofuranol.",
            "smarts": "[#7]C(=O)Oc1c2Occcc2c(C)c1C"
        }
    },
    "Category 5: Ion Channel Modulators, Enzyme Inhibitors, and Signaling Pathway Inhibitors": {
        "5.a.1": {
            "tier": 2,
            "name": "Antiarrhythmic - IKr Blocker (Key Feature)",
            "description": "Identifies the core pharmacophore of Class III antiarrhythmic IKr blockers.",
            "smarts": "c1ccc(cc1)C(O)CCN"
        },
        "5.a.2": {
            "tier": 2,
            "name": "Antiarrhythmic - Amiloride-like (Key Feature)",
            "description": "Identifies the pyrazine carboxamide scaffold of amiloride-like sodium channel inhibitors.",
            "smarts": "c1(N)nc(N)c(C(=O)N)nc1Cl"
        },
        "5.b.1": {
            "tier": 1,
            "name": "Beta-Blocker - Aryloxypropanolamine (Core)",
            "description": "Identifies the general aryloxypropanolamine scaffold of beta-blockers.",
            "smarts": "[a]O[#6]C(O)[#6][#7]"
        },
        "5.b.2": {
            "tier": 2,
            "name": "Beta-Blocker - Fused Aryloxypropanolamine (Key Feature)",
            "description": "Identifies the key feature where the aryl group of an aryloxypropanolamine is a fused ring system (e.g., indole, naphthalene).",
            "smarts": "[a;R2]OCCC(O)CN"
        },
        "5.c.1": {
            "tier": 2,
            "name": "ACE Inhibitor - Acyl-Proline (Key Feature)",
            "description": "Identifies the N-acyl-proline moiety, a key feature of ACE inhibitors.",
            "smarts": "[#6]C(=O)N1[#6](C(=O)O)CCC1"
        },
        "5.c.2": {
            "tier": 2,
            "name": "ARB - Biphenyl-Tetrazole (Key Feature)",
            "description": "Identifies the biphenyl-tetrazole scaffold, the hallmark of the 'sartan' class of ARBs.",
            "smarts": "c1ccc(c2ccccc2-c3nnnn3)cc1"
        },
        "5.d.1": {
            "tier": 2,
            "name": "Shh Inhibitor - Cyclopamine-like (Key Feature)",
            "description": "Identifies the unique fused furanopiperidine system of teratogenic Veratrum alkaloids.",
            "smarts": "[#6]12[#6][#6][#6](O1)[#6][#7][#6][#6]2"
        },
        "5.d.2": {
            "tier": 2,
            "name": "Shh Inhibitor - AY-9944-like (Key Feature)",
            "description": "Identifies the trans-1,4-bis(2-chlorobenzylaminomethyl)cyclohexane scaffold of AY-9944.",
            "smarts": "c1(Cl)ccccc1CNCC2CCC(NCc1ccccc1Cl)CC2"
        },
        "5.d.3": {
            "tier": 2,
            "name": "Cholesterol Synth Inhibitor - Triparanol-like (Key Feature)",
            "description": "Identifies the triarylethanol scaffold of cholesterol synthesis inhibitors like triparanol.",
            "smarts": "[#6](O)(c1ccccc1)(c2ccccc2)c3ccccc3"
        }
    },
    "Category 6: Opioid Receptor Binders and Tubulin Interactors": {
        "6.a.1": {
            "tier": 2,
            "name": "Opioid - Morphinan Scaffold (Key Feature)",
            "description": "Identifies the rigid, pentacyclic 4,5a-epoxymorphinan ring system of morphine.",
            "smarts": "[#6]12[#6][#6][#7]C[#6]c3c4O[#6]([#6]1c3)C=C[#6]24"
        },
        "6.a.2": {
            "tier": 2,
            "name": "Opioid - Meperidine-like Scaffold (Key Feature)",
            "description": "Identifies the 4-phenyl-4-carbalkoxypiperidine scaffold of meperidine.",
            "smarts": "C1(c2ccccc2)(C(=O)O)CCNCC1"
        },
        "6.b.1": {
            "tier": 2,
            "name": "Tubulin - Benzimidazole Carbamate (Key Feature)",
            "description": "Identifies the 2-(methoxycarbonylamino)benzimidazole scaffold.",
            "smarts": "c12[nH]c(nc1cccc2)NC(=O)OC"
        },
        "6.b.2": {
            "tier": 2,
            "name": "Tubulin - Podophyllotoxin (Key Feature)",
            "description": "Identifies the complex polycyclic lignan scaffold of podophyllotoxin.",
            "smarts": "c1(OC)c(OC)c(OC)cc(c1)[#6]1[#6]c2c(c(O[#6]O2)cc)C(=O)O1"
        },
        "6.b.3": {
            "tier": 2,
            "name": "Tubulin - Colchicine (Key Feature)",
            "description": "Identifies the unique tropolone-containing scaffold of colchicine.",
            "smarts": "c1(OC)c(OC)c(OC)c2c(c1)C(=O)C=CC=C[#6]2"
        },
        "6.b.4": {
            "tier": 2,
            "name": "Tubulin - Vinca Alkaloid Scaffold (Key Feature)",
            "description": "Identifies the complex, fused indole-containing scaffold of Vinca alkaloids like vinblastine.",
            "smarts": "C1=CC=C2C(=C1)C3(C4=C(CCN3C)C5=C(C=CC=C5)N4)C6(C(C(C7(C(C(=O)OC)C6O)C=C)CC)N(C2)C7)O"
        },
        "6.b.5": {
            "tier": 2,
            "name": "Tubulin - Taxane (Key Feature)",
            "description": "Identifies the characteristic four-membered oxetane ring fused to the taxane core.",
            "smarts": "[#6;r4]1[#6;r4][#6;r4]O1"
        },
        "6.b.6": {
            "tier": 2,
            "name": "Tubulin - Epothilone (Key Feature)",
            "description": "Identifies the 16-membered macrocyclic lactone core of the epothilones.",
            "smarts": "C1C/C=C/C(=O)OCCCC/C=C/C(C)C(c2csc(n2)C)C1(C)O"
        }
    },
    "Category 7: Nucleotide and Nucleobase Derivatives": {
        "7.a.1": {
            "tier": 1,
            "name": "Pyrimidine Nucleoside (Core)",
            "description": "Identifies the general pyrimidine nucleoside scaffold.",
            "smarts": "c1ncc(n(c1)[#6;R]2[#6;R][#6;R][#6;R]O2)"
        },
        "7.a.2": {
            "tier": 2,
            "name": "Pyrimidine Nucleoside - C5-Halogenated (Key Feature)",
            "description": "Identifies the key feature of a C-5 halogen on a pyrimidine nucleoside.",
            "smarts": "c1c([F,Cl,Br,I])cn(C2OC(CO)C(O)C2O)c(=O)[nH]1"
        },
        "7.a.3": {
            "tier": 2,
            "name": "C-5 Methylated Pyrimidine Nucleoside (Key Feature)",
            "description": "Identifies pyrimidine nucleosides with a methyl group at the C-5 position.",
            "smarts": "c1(C)cn(C2OC(CO)C(O)C2O)c(=O)n(c1=O)"
        },
        "7.a.4": {
            "tier": 2,
            "name": "5-Aza Pyrimidine Nucleoside (Key Feature)",
            "description": "Identifies nucleosides with a 5-aza (triazine) base.",
            "smarts": "n1cn(C2OC(CO)C(O)C2O)c(=O)n(c1=O)"
        },
        "7.b.1": {
            "tier": 1,
            "name": "Purine (Core)",
            "description": "Identifies the general purine scaffold.",
            "smarts": "c1c2ncnc2[nH]c1"
        },
        "7.b.2": {
            "tier": 2,
            "name": "Acyclovir-like Purine (Key Feature)",
            "description": "Identifies the key 2-hydroxyethoxymethyl side chain on a purine base.",
            "smarts": "c1c(N)nc2n(c1=O)cn(CCOCCO)c2"
        },
        "7.b.3": {
            "tier": 2,
            "name": "5-Fluorouracil (Key Feature)",
            "description": "Identifies the potent antimetabolite 5-fluorouracil.",
            "smarts": "c1(F)c(NC(=O)NC1=O)"
        },
        "7.b.4": {
            "tier": 2,
            "name": "Thiouracil (Key Feature)",
            "description": "Identifies the propylthiouracil scaffold, a known thyroid toxicant.",
            "smarts": "c1(CCC)c(S)nc(=O)[nH]c1"
        }
    },
    "Category 8: Aromatic Compounds with Alkyl, Multi-Halogen, and Nitro Groups": {
        "8.a.1": {
            "tier": 1,
            "name": "Alkylbenzene (Core)",
            "description": "Identifies benzene rings substituted with a short alkyl chain (approx. C1-C4).",
            "smarts": "c1ccccc1[#6;!a;!R;!D3]"
        },
        "8.b.1": {
            "tier": 1,
            "name": "Nitroaromatic (Core)",
            "description": "Identifies aromatic rings containing one or more nitro groups.",
            "smarts": "c[N+](=O)[O-]"
        },
        "8.b.2": {
            "tier": 2,
            "name": "Nitrotoluene (Key Feature)",
            "description": "Identifies nitrotoluenes, a combination associated with a specific male reproductive toxicity profile.",
            "smarts": "c1c(C)cccc1[N+](=O)[O-]"
        },
        "8.c.1": {
            "tier": 1,
            "name": "Di-chlorinated Benzene (Core)",
            "description": "Identifies benzene rings having two chlorine substituents.",
            "smarts": "c1(Cl)c(Cl)cccc1"
        },
        "8.c.2": {
            "tier": 2,
            "name": "Poly-chlorinated Benzene (Key Feature)",
            "description": "Identifies benzene rings with three or more chlorine atoms.",
            "smarts": "c1(Cl)c(Cl)c(Cl)ccc1"
        },
        "8.d.1": {
            "tier": 1,
            "name": "Substituted Diphenyl Ether (Core)",
            "description": "Identifies a diphenyl ether scaffold containing at least one halogen or nitro group.",
            "smarts": "[$(c1(Oc2ccccc2)c([Cl,Br,I,F])cccc1),$(c1(Oc2ccccc2)c([N+](=O)[O-])cccc1)]"
        },
        "8.d.2": {
            "tier": 2,
            "name": "Nitrofen-like (Key Feature)",
            "description": "Identifies the key feature of a para-nitro group on one ring and halogenation on the other.",
            "smarts": "c1cc([N+](=O)[O-])ccc1Oc2ccc(Cl)c(Cl)c2"
        },
        "8.e.1": {
            "tier": 2,
            "name": "Substituted Phenol - 2,6-Dihalo (Key Feature)",
            "description": "Identifies the key feature of a phenol with halogen atoms at the 2 and 6 positions.",
            "smarts": "c1(O)c([Cl,Br])ccc([Cl,Br])c1"
        },
        "8.e.2": {
            "tier": 2,
            "name": "2,4-Dinitro Phenol (Key Feature)",
            "description": "Identifies the key feature of 2,4-dinitro substitution on a phenol or phenyl ester.",
            "smarts": "c1c(O)c([N+](=O)[O-])cc([N+](=O)[O-])c1"
        }
    },
    "Category 9: Aromatic Compounds with Oxygenated Side Chains": {
        "9.a.1": {
            "tier": 2,
            "name": "BMHCA-like (Key Feature)",
            "description": "Identifies phenylpropanals with a sterically hindered para-alkyl group.",
            "smarts": "c1cc(C(C)(C)C)ccc1CC(C)C=O"
        },
        "9.b.1": {
            "tier": 1,
            "name": "Aryl-Alkyl Acid (Core)",
            "description": "Identifies the general aryl-alkyl-acid scaffold of NSAIDs.",
            "smarts": "[a]C(C)C(=O)O"
        },
        "9.b.2": {
            "tier": 2,
            "name": "Chlorambucil (Key Feature)",
            "description": "Identifies the potent alkylating bis(2-chloroethyl)amino moiety of Chlorambucil.",
            "smarts": "c1cc(N(CCCl)CCCl)ccc1CCCC(=O)O"
        },
        "9.c.1": {
            "tier": 1,
            "name": "Aryloxy-Aliphatic Acid (Core)",
            "description": "Identifies the general aryloxy-aliphatic acid core of phenoxy herbicides.",
            "smarts": "[a]O[#6;!a]C(=O)O"
        },
        "9.c.2": {
            "tier": 2,
            "name": "Aryloxy-Aliphatic Acid - Polychlorinated (Key Feature)",
            "description": "Identifies the key feature of a polychlorinated (2,4,5-trichloro) phenoxy acetic acid.",
            "smarts": "c1(Cl)c(Cl)cc(O)c(c1)Cl.CC(=O)O"
        }
    },
    "Category 10: Aromatic Compounds with Sulfonamide and Urea Moieties, Phenytoins": {
        "10.a.1": {
            "tier": 1,
            "name": "Arylsulfonamide (Core)",
            "description": "Identifies the general arylsulfonamide scaffold.",
            "smarts": "[a]S(=O)(=O)N"
        },
        "10.a.2": {
            "tier": 2,
            "name": "Sulfa Drug (Key Feature)",
            "description": "Identifies the 4-aminophenylsulfonamide core of sulfa drugs.",
            "smarts": "c1cc(N)ccc1S(=O)(=O)N"
        },
        "10.a.3": {
            "tier": 2,
            "name": "Arylsulfonylurea (Key Feature)",
            "description": "Identifies the arylsulfonylurea moiety.",
            "smarts": "[a]S(=O)(=O)NC(=O)N"
        },
        "10.b.1": {
            "tier": 1,
            "name": "Hydantoin (Core)",
            "description": "Identifies the general hydantoin ring system.",
            "smarts": "C1C(=O)NC(=O)N1"
        },
        "10.b.2": {
            "tier": 2,
            "name": "Phenytoin (Key Feature)",
            "description": "Identifies the 5,5-diphenylhydantoin scaffold of the teratogen phenytoin.",
            "smarts": "C1(C(=O)NC(=O)N1)(c2ccccc2)c3ccccc3"
        }
    },
    "Category 11: Aromatic Compounds with Aliphatic Amine Moieties": {
        "11.a.1": {
            "tier": 1,
            "name": "Arylethane Amine (Core)",
            "description": "Identifies the general arylethane amine scaffold (Aryl-C-C-N).",
            "smarts": "[a]CC[N]"
        },
        "11.a.2": {
            "tier": 2,
            "name": "Phenyl Nitrogen Mustard (Key Feature)",
            "description": "Identifies the potent phenyl nitrogen mustard moiety of melphalan.",
            "smarts": "c1ccc(cc1)N(CCCl)CCCl"
        },
        "11.b.1": {
            "tier": 2,
            "name": "Diarylmethylpiperazine (Key Feature)",
            "description": "Identifies the diarylmethylpiperazine scaffold of antihistamines like cyclizine.",
            "smarts": "[#6](c1ccccc1)(c2ccccc2)N3CCNCC3"
        }
    },
    "Category 12: Aromatic Diamines, Azo/Diazo Compounds, and Dyes": {
        "12.a.1": {
            "tier": 1,
            "name": "Diaminopyrimidine & Analogs (Core)",
            "description": "Identifies the 2,4-diaminopyrimidine scaffold and its fused analogs (pteridine, quinazoline).",
            "smarts": "c1(N)ncc(N)nc1"
        },
        "12.a.2": {
            "tier": 2,
            "name": "5-Phenyl-Diaminopyrimidine (Key Feature)",
            "description": "Identifies the key feature of a phenyl substituent at the C-5 position of a diaminopyrimidine.",
            "smarts": "c1(N)c(c2ccccc2)c(N)nc(n1)"
        },
        "12.b.1": {
            "tier": 2,
            "name": "Benzidine-Azo Dye (Key Feature)",
            "description": "Identifies the specific scaffold of azo dyes derived from benzidine.",
            "smarts": "c1cc(N=Nc2ccc(c(c2)c3ccc(cc3)N=N)C)ccc1"
        },
        "12.b.2": {
            "tier": 2,
            "name": "Dialkylaminoazobenzene (Key Feature)",
            "description": "Identifies the specific scaffold of dialkylaminoazobenzene dyes.",
            "smarts": "c1cc(N(C)C)ccc1N=Nc2ccccc2"
        },
        "12.c.1": {
            "tier": 1,
            "name": "Triarylmethane (Core)",
            "description": "Identifies the general triarylmethane scaffold.",
            "smarts": "[#6](c1ccccc1)(c2ccccc2)c3ccccc3"
        },
        "12.c.2": {
            "tier": 2,
            "name": "Triarylmethane Dye (Key Feature)",
            "description": "Identifies the key feature of para-amino substitution on at least two of the phenyl rings.",
            "smarts": "[#6](c1ccc(N)cc1)(c2ccc(N)cc2)c3ccccc3"
        },
        "12.d.1": {
            "tier": 2,
            "name": "Aryl Triazene (Key Feature)",
            "description": "Identifies the potent aryl triazene functional group.",
            "smarts": "[a]N=NN"
        }
    },
    "Category 13: Imidazoles, Nitrofurans, and Triazoles": {
        "13.a.1": {
            "tier": 2,
            "name": "Antifungal - Phenethyl Imidazole (Key Feature)",
            "description": "Identifies the phenethyl imidazole scaffold of certain antifungal agents.",
            "smarts": "c1ccc(cc1)CC(n2ccnc2)"
        },
        "13.a.2": {
            "tier": 2,
            "name": "Antifungal - Phenoxymethyl Imidazole (Key Feature)",
            "description": "Identifies the phenoxymethyl imidazole scaffold of certain antifungal agents.",
            "smarts": "c1ccccc1OC[#6](n2ccnc2)"
        },
        "13.b.1": {
            "tier": 2,
            "name": "Nitroimidazole (Key Feature)",
            "description": "Identifies the potent nitroimidazole toxicophore.",
            "smarts": "c1c([N+](=O)[O-])[nH]cn1"
        },
        "13.b.2": {
            "tier": 2,
            "name": "Nitrofuran-CH=N- (Key Feature)",
            "description": "Identifies the nitrofurfurylideneamino moiety.",
            "smarts": "c1c(oc1[N+](=O)[O-])C=N"
        },
        "13.c.1": {
            "tier": 1,
            "name": "1,2,4-Triazole (Core)",
            "description": "Identifies the general 1,2,4-triazole ring.",
            "smarts": "c1ncn[nH]1"
        },
        "13.c.2": {
            "tier": 2,
            "name": "Antifungal - Conazole (Key Feature)",
            "description": "Identifies the key feature of conazole antifungals: a triazole ring attached to a tertiary carbon bearing an aryl group.",
            "smarts": "c1cc(ccc1)[#6](Cn2cncn2)"
        }
    },
    "Category 14: Fused Aromatic and Heterocyclic Derivatives": {
        "14.a.1": {
            "tier": 1,
            "name": "Coumarin (Core)",
            "description": "Identifies the general benzopyran-2-one (coumarin) scaffold.",
            "smarts": "c1ccc2c(c1)ccc(=O)o2"
        },
        "14.a.2": {
            "tier": 2,
            "name": "Hydroxy-Coumarin (Key Feature)",
            "description": "Identifies the key feature of a hydroxyl group on the coumarin ring, associated with teratogenicity.",
            "smarts": "c1c(O)cc2c(c1)ccc(=O)o2"
        },
        "14.b.1": {
            "tier": 2,
            "name": "Thalidomide (Key Feature)",
            "description": "Identifies the specific phthalimide-glutarimide fused structure of the potent teratogen thalidomide.",
            "smarts": "O=C1NC(=O)c2ccccc2N1C3CCC(=O)NC3=O"
        },
        "14.b.2": {
            "tier": 1,
            "name": "Quinolone Acid (Core)",
            "description": "Identifies the core 4-oxo-quinoline-3-carboxylic acid scaffold.",
            "smarts": "c12c(cccc1)c(C(=O)O)cn(c2=O)"
        },
        "14.b.3": {
            "tier": 2,
            "name": "Quinolone Acid with C7-Amine-Heterocycle (Key Feature)",
            "description": "Identifies the key feature of an amine-containing heterocycle at the C-7 position of the quinolone core.",
            "smarts": "c1(N2CCNCC2)cc2c(c(c1)F)c(C(=O)O)cn(c2=O)"
        },
        "14.b.4": {
            "tier": 2,
            "name": "Benzodiazepine (Key Feature)",
            "description": "Identifies the 1,4-benzodiazepine-2-one scaffold of anxiolytic drugs.",
            "smarts": "c1ccc2c(c1)N=C(c3ccccc3)C[#7]C2=O"
        },
        "14.c.1": {
            "tier": 2,
            "name": "Tricyclic - Phenothiazine (Key Feature)",
            "description": "Identifies the tricyclic phenothiazine ring system of antipsychotic drugs.",
            "smarts": "c1ccc2c(c1)Nc3ccccc3S2"
        },
        "14.c.2": {
            "tier": 2,
            "name": "Tricyclic - Dibenzazepine (Key Feature)",
            "description": "Identifies the tricyclic dibenzazepine ring system of antidepressant drugs.",
            "smarts": "c1ccc2c(c1)Nccc3ccccc23"
        },
        "14.d.1": {
            "tier": 2,
            "name": "Tetracycline Scaffold (Key Feature)",
            "description": "Identifies the unique, linearly fused, and highly oxygenated naphthacene carboxamide core of tetracycline antibiotics.",
            "smarts": "[#6]12C(=O)c3c(O)c4c(c(O)ccc3)[C@@](C)(O)[C@H]1[C@@H](N(C)C)C(=O)C2=C(O)C4=O"
        },
        "14.d.2": {
            "tier": 2,
            "name": "Anthracycline (Key Feature)",
            "description": "Identifies the planar, tetracyclic quinone core of anthracycline chemotherapeutics.",
            "smarts": "c1c(O)c2C(=O)c3c(C)c(O)ccc3C(=O)c2cc1"
        }
    },
    "Category 15: Miscellaneous Aromatic Chemicals": {
        "15.a.1": {
            "tier": 2,
            "name": "3-Aminopyrazine-2-carboxamide (Scaffold)",
            "description": "Alert for 3-aminopyrazine-2-carboxamide.",
            "smarts": "c1(N)c(C(=O)N)nccn1"
        },
        "15.b.1": {
            "tier": 2,
            "name": "Isoniazid (Scaffold)",
            "description": "Alert for the isonicotinohydrazide scaffold of Isoniazid.",
            "smarts": "c1cncc(c1)C(=O)NN"
        },
        "15.c.1": {
            "tier": 2,
            "name": "Aminonicotinamide (Scaffold)",
            "description": "Alert for 6-Aminonicotinamide.",
            "smarts": "c1c(N)ccnc1C(=O)N"
        },
        "15.d.1": {
            "tier": 2,
            "name": "Phenelzine (Scaffold)",
            "description": "Alert for the phenylethylhydrazine scaffold of Phenelzine.",
            "smarts": "c1ccccc1CCNN"
        },
        "15.e.1": {
            "tier": 2,
            "name": "Phencyclidine (PCP) (Scaffold)",
            "description": "Alert for the 1-(1-phenylcyclohexyl)piperidine scaffold of Phencyclidine.",
            "smarts": "C1(CCCCC1)(c2ccccc2)N3CCCCC3"
        },
        "15.f.1": {
            "tier": 2,
            "name": "Aminoglutethimide (Scaffold)",
            "description": "Alert for the aminophenyl-substituted glutarimide scaffold of Aminoglutethimide.",
            "smarts": "C1(CC(=O)NC1=O)c2ccc(N)cc2"
        },
        "15.g.1": {
            "tier": 2,
            "name": "Ketamine (Scaffold)",
            "description": "Alert for the 2-(2-chlorophenyl)-2-(methylamino)cyclohexan-1-one scaffold of Ketamine.",
            "smarts": "C1(=O)CCCCC1(NC)c2ccccc2Cl"
        },
        "15.h.1": {
            "tier": 2,
            "name": "Phenindione (Scaffold)",
            "description": "Alert for the 2-phenyl-1,3-indandione scaffold of Phenindione.",
            "smarts": "c1ccccc1C2C(=O)c3ccccc3C2=O"
        },
        "15.i.1": {
            "tier": 2,
            "name": "Haloperidol (Scaffold)",
            "description": "Alert for the butyrophenone-piperidine core of Haloperidol.",
            "smarts": "c1cc(F)ccc1C(=O)CCCCN2CCC(O)(c3ccc(Cl)cc3)CC2"
        },
        "15.j.1": {
            "tier": 2,
            "name": "Pyrethroid (Scaffold)",
            "description": "Alert for the alpha-cyano-3-phenoxybenzyl ester of a dichlorovinylcyclopropane carboxylic acid, characteristic of Type II pyrethroids.",
            "smarts": "ClC(=CC1C(C1(C)C)C(=O)OC(C#N)c2cccc(c2)Oc3ccccc3)Cl"
        }
    },
    "Category 16: Non-aromatic Cyclic Compounds and Derivatives": {
        "16.a.1": {
            "tier": 2,
            "name": "Vitamin D-like Secosteroid (Key Feature)",
            "description": "Identifies the secosteroid scaffold of Vitamin D, characterized by a broken steroid B-ring and a conjugated triene system.",
            "smarts": "C1(CCC/C1=C/C=C)=C"
        },
        "16.b.1": {
            "tier": 2,
            "name": "Long-Chain N-Alkylmorpholine (Key Feature)",
            "description": "Identifies morpholine rings substituted with a long alkyl chain (C11-C14).",
            "smarts": "C1COCCN1CCCCCCCCCCC"
        },
        "16.c.1": {
            "tier": 2,
            "name": "Oxirane (Epoxide)",
            "description": "Identifies the reactive oxirane ring.",
            "smarts": "[#6]1[#6]O1"
        },
        "16.d.1": {
            "tier": 2,
            "name": "Aminoglycoside - 2-Deoxystreptamine (Key Feature)",
            "description": "Identifies the highly functionalized 2-deoxystreptamine aminocyclitol core common to this class of antibiotics.",
            "smarts": "C1(N)C(O)C(O)C(N)C(O)C1O"
        },
        "16.e.1": {
            "tier": 2,
            "name": "Macrocyclic Lactone - Avermectin (Key Feature)",
            "description": "Identifies the complex spiroketal-containing macrocycle of Avermectins.",
            "smarts": "[#6]1O[#6]2[#6](O1)[#6]C(O[#6]3C[#6](O[#6]4CC(O)C(C)O4)C(C)O3)=C[#6]C(=O)O[#6]/C=C/C=C/[#6@H]2C"
        },
        "16.e.2": {
            "tier": 2,
            "name": "Macrocyclic Lactone - Spiramycin (Key Feature)",
            "description": "Identifies the sugar-substituted 16-membered macrocyclic lactone core of Spiramycins.",
            "smarts": "[#6]1OC(=O)C[#6](O)[#6]C=CC=C[#6@H](C)[#6@H](O[#6]2O[#6](C)[#6@H](O[#6]3O[#6](C)C[#6@H](N(C)C)[#6@H]3O)[#6@@H]2O)[#6@H]1C"
        },
        "16.f.1": {
            "tier": 2,
            "name": "Polychlorinated Cyclodiene (Key Feature)",
            "description": "Identifies the polychlorinated, bridged norbornene scaffold characteristic of cyclodiene insecticides.",
            "smarts": "ClC1=C(Cl)C2(Cl)C3C=CC(C3)C12Cl"
        }
    },
    "Category 17: Non-aromatic Heterocyclic Compounds": {
        "17.a.1": {
            "tier": 1,
            "name": "Barbiturate (Core)",
            "description": "Identifies the core barbituric acid (pyrimidinetrione) ring.",
            "smarts": "C1(=O)NC(=O)NC(=O)C1"
        },
        "17.a.2": {
            "tier": 2,
            "name": "5,5-Disubstituted Barbiturate (Key Feature)",
            "description": "Identifies the key feature of 5,5-dialkyl substitution on the barbiturate ring.",
            "smarts": "C1(C(=O)NC(=O)NC1=O)(C)C"
        },
        "17.b.1": {
            "tier": 2,
            "name": "Cyclic Thiourea (Key Feature)",
            "description": "Identifies the cyclic thiourea moiety, a key feature of thyroid toxicants like ETU.",
            "smarts": "[#7;R][#6;R](=[#16;R0])[#7;R]"
        },
        "17.c.1": {
            "tier": 2,
            "name": "Oxazolidinedione (Key Feature)",
            "description": "Identifies the oxazolidine-2,4-dione scaffold of anticonvulsants like dimethadione.",
            "smarts": "C1C(=O)N(C(=O)O1)"
        },
        "17.c.2": {
            "tier": 2,
            "name": "Hydantoin (Key Feature)",
            "description": "Identifies the imidazolidine-2,4-dione (hydantoin) scaffold.",
            "smarts": "C1C(=O)NC(=O)N1"
        },
        "17.d.1": {
            "tier": 1,
            "name": "Common Saturated 6-Membered Heterocycles (Core)",
            "description": "A low-specificity alert for common 6-membered heterocycles (piperazine, dioxane, morpholine). Included for completeness.",
            "smarts": "C1COCCN1.C1COCOC1.C1CCNCC1"
        }
    },
    "Category 18: Miscellaneous Non-aromatic Cyclic Chemicals": {
        "18.a.1": {
            "tier": 2,
            "name": "Hydantoin (Scaffold)",
            "description": "Alert for the hydantoin (imidazolidine-2,4-dione) ring. Overlaps with other alerts.",
            "smarts": "C1C(=O)NC(=O)N1"
        },
        "18.b.1": {
            "tier": 2,
            "name": "Cycloheximide (Scaffold)",
            "description": "Alert for the glutarimide-substituted dimethyl-cyclohexanone scaffold of cycloheximide.",
            "smarts": "CC1CC(C)CC(C1)C(C(CC(=O)NC(=O)C2)C2)O"
        },
        "18.c.1": {
            "tier": 2,
            "name": "Tropolone (Scaffold)",
            "description": "Alert for the seven-membered tropolone ring scaffold of compounds like hinokitiol.",
            "smarts": "c1(C(C)C)c(O)c(=O)cccc1"
        },
        "18.d.1": {
            "tier": 2,
            "name": "Cyclohexanedione Oxime Ether (Scaffold)",
            "description": "Alert for the cyclohexanedione oxime ether moiety found in certain herbicides.",
            "smarts": "[#6]1(=O)C[#6](C(=NO))C[#6](=O)C1"
        },
        "18.e.1": {
            "tier": 2,
            "name": "Mutilin-like (Scaffold)",
            "description": "Alert for the tricyclic core of pleuromutilin antibiotics.",
            "smarts": "C12CC(C(C1C3C(C(C=C)C(O3)C)C)C(=O))C(C2(C)C)OC(=O)C"
        }
    },
    "Category 19: Dithiocarbamates, Sulfonates, and Perfluorinated Compounds": {
        "19.a.1": {
            "tier": 1,
            "name": "Dithiocarbamate (Core)",
            "description": "Identifies the general dithiocarbamate functional group (N-C(=S)-S).",
            "smarts": "[#7]C(=S)S"
        },
        "19.a.2": {
            "tier": 2,
            "name": "Thiuram Disulfide (Key Feature)",
            "description": "Identifies the thiuram disulfide structure (-S-CS-N(R)R), a key feature of fungicides like thiram.",
            "smarts": "C(SC(=S)N(C)C)SC(=S)N(C)C"
        },
        "19.b.1": {
            "tier": 1,
            "name": "Alkyl Sulfonate (Core)",
            "description": "Identifies monofunctional alkyl sulfonate esters, which are alkylating agents.",
            "smarts": "[#6;!a]S(=O)(=O)O[#6;!a]"
        },
        "19.b.2": {
            "tier": 2,
            "name": "Busulfan-like (Key Feature)",
            "description": "Identifies bifunctional bis-alkyl sulfonate esters, potent cross-linking agents like busulfan.",
            "smarts": "CS(=O)(=O)OCCCCOS(=O)(=O)C"
        },
        "19.c.1": {
            "tier": 1,
            "name": "Perfluorinated Chain (Core, >=C4)",
            "description": "Identifies compounds with a moderately long (>=C4) perfluorinated alkyl chain.",
            "smarts": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)"
        },
        "19.c.2a": {
            "tier": 2,
            "name": "Perfluoroalkylsulfonate (>=C6) (Key Feature)",
            "description": "Identifies the key feature of a long (>=C6) perfluorinated chain attached to a sulfonate head group.",
            "smarts": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"
        },
        "19.c.2b": {
            "tier": 2,
            "name": "Perfluoroalkylcarboxylate (>=C6) (Key Feature)",
            "description": "Identifies the key feature of a long (>=C6) perfluorinated chain attached to a carboxylate head group.",
            "smarts": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
        }
    },
    "Category 20: Miscellaneous Non-cyclic Chemicals": {
        "20.a.1": {
            "tier": 2,
            "name": "Methylazoxy methyl acetate (Scaffold)",
            "description": "Alert for the specific structure of methylazoxy methyl acetate. Intended for exact scaffold match.",
            "smarts": "C[N+]([O-])=NCC(=O)OC"
        },
        "20.b.1": {
            "tier": 2,
            "name": "Hexan-2-one (Scaffold)",
            "description": "Alert for the specific structure of the neurotoxic solvent hexan-2-one. Intended for exact scaffold match.",
            "smarts": "CCCCC(=O)C"
        },
        "20.c.1": {
            "tier": 2,
            "name": "2,5-Hexanedione (Scaffold)",
            "description": "Alert for the specific neurotoxic gamma-diketone 2,5-hexanedione. Intended for exact scaffold match.",
            "smarts": "CC(=O)CCC(=O)C"
        },
        "20.d.1": {
            "tier": 2,
            "name": "Multihalogenated Acetone (Key Feature)",
            "description": "Identifies hexahaloacetones, where a carbonyl is flanked by two trihalomethyl groups.",
            "smarts": "C([F,Cl,Br])([F,Cl,Br])([F,Cl,Br])C(=O)C([F,Cl,Br])([F,Cl,Br])([F,Cl,Br])"
        },
        "20.e.1": {
            "tier": 2,
            "name": "Meprobamate (Scaffold)",
            "description": "Alert for the specific structure of the anxiolytic drug meprobamate. Intended for exact scaffold match.",
            "smarts": "CCCC(C)(COC(=O)N)COC(=O)N"
        }
    },
    "Category 21: Reactive Carbonyls, Amides, and Related Functional Groups": {
        "21.a.1": {
            "tier": 1,
            "name": "Michael Acceptor (Core)",
            "description": "Identifies the alpha,beta-unsaturated carbonyl moiety, a reactive electrophile.",
            "smarts": "[#6]=[#6]C(=O)"
        },
        "21.b.1a": {
            "tier": 1,
            "name": "Amide (Core)",
            "description": "Identifies the core amide functional group (-CONH-).",
            "smarts": "[#7]C(=O)[#6]"
        },
        "21.b.1b": {
            "tier": 1,
            "name": "Thioamide (Core)",
            "description": "Identifies the core thioamide functional group (-CSNH-).",
            "smarts": "[#7]C(=S)[#6]"
        },
        "21.b.2a": {
            "tier": 1,
            "name": "Urea (Core)",
            "description": "Identifies the core urea functional group (-NHCONH-).",
            "smarts": "[#7]C(=O)[#7]"
        },
        "21.b.2b": {
            "tier": 1,
            "name": "Thiourea (Core)",
            "description": "Identifies the core thiourea functional group (-NHCSNH-).",
            "smarts": "[#7]C(=S)[#7]"
        },
        "21.b.3": {
            "tier": 1,
            "name": "Guanidine (Core)",
            "description": "Identifies the general guanidine functional group.",
            "smarts": "[#7]C(=[#7])[#7]"
        },
        "21.b.4": {
            "tier": 1,
            "name": "Carbonate (Core)",
            "description": "Identifies the general carbonate ester functional group.",
            "smarts": "[#8]C(=O)[#8]"
        },
        "21.b.5": {
            "tier": 1,
            "name": "Carbamate (Core)",
            "description": "Identifies the general carbamate functional group.",
            "smarts": "[#7]C(=O)O"
        },
        "21.b.6": {
            "tier": 2,
            "name": "N-Nitrosourea (Key Feature)",
            "description": "Identifies the potent N-nitrosourea alkylating agent moiety.",
            "smarts": "[#7]N(N=O)C(=O)"
        },
        "21.b.7": {
            "tier": 2,
            "name": "Methylaminocarbonyl Carbamate (Key Feature)",
            "description": "Identifies the O-(methylamino)carbonyl carbamate scaffold, which acts via AChE inhibition.",
            "smarts": "CNC(=O)OC(=N)"
        }
    },
    "Category 22: Alpha-Substituted Carboxylic Acids and Adipates": {
        "22.a.1": {
            "tier": 2,
            "name": "alpha-Haloacetic Acid (Key Feature)",
            "description": "Identifies carboxylic acids/esters substituted at the alpha-carbon with a chlorine or bromine atom.",
            "smarts": "[#6]([Cl,Br])C(=O)O"
        },
        "22.b.1": {
            "tier": 2,
            "name": "alpha-Hydroxyacetic Acid (Key Feature)",
            "description": "Identifies carboxylic acids/esters substituted at the alpha-carbon with a hydroxyl group (glycolates).",
            "smarts": "[#6](O)C(=O)O"
        },
        "22.b.2": {
            "tier": 1,
            "name": "alpha-Alkoxyacetic Acid (Core)",
            "description": "Identifies the general alpha-alkoxyacetic acid scaffold.",
            "smarts": "[#6](O[#6;!a])C(=O)O"
        },
        "22.b.3a": {
            "tier": 2,
            "name": "alpha-Alkoxyacetic Acid - Methoxy (Key Feature)",
            "description": "Identifies the key feature of a short-chain (C1) alkoxy group on an acetic acid, which confers higher potency.",
            "smarts": "[#6](O[CH3])C(=O)O"
        },
        "22.b.3b": {
            "tier": 2,
            "name": "alpha-Alkoxyacetic Acid - Ethoxy (Key Feature)",
            "description": "Identifies the key feature of a short-chain (C2) alkoxy group on an acetic acid, which confers higher potency.",
            "smarts": "[#6](OCC)C(=O)O"
        },
        "22.c.1": {
            "tier": 1,
            "name": "alpha-Alkylcarboxylic Acid (Core)",
            "description": "Identifies carboxylic acids/esters with a tertiary alpha-carbon (mono-alkyl substitution).",
            "smarts": "[#6;X3](C(=O)O)[#6]"
        },
        "22.c.2": {
            "tier": 2,
            "name": "alpha,alpha-Dialkylcarboxylic Acid (Key Feature)",
            "description": "Identifies carboxylic acids/esters with a quaternary alpha-carbon (e.g., valproic acid).",
            "smarts": "[#6;X4](C(=O)O)(C)C"
        },
        "22.d.1": {
            "tier": 2,
            "name": "Adipate Derivatives",
            "description": "Identifies the six-carbon dicarboxylic acid/ester backbone of adipates.",
            "smarts": "C(CCC(=O)O)C(=O)O"
        }
    },
    "Category 23: Small Halogenated Aliphatics and Mustards": {
        "23.a.1": {
            "tier": 1,
            "name": "Haloaliphatic (Core)",
            "description": "Identifies short-chain alkanes, alkenes, or ethers with at least one Cl or Br substituent.",
            "smarts": "[#6;!a;!R][Cl,Br]"
        },
        "23.a.2": {
            "tier": 2,
            "name": "Polyhalogenated - Trihalomethyl (Key Feature)",
            "description": "Identifies the key feature of a trihalomethyl group (-CX3).",
            "smarts": "C([F,Cl,Br])([F,Cl,Br])[F,Cl,Br]"
        },
        "23.c.1": {
            "tier": 2,
            "name": "Haloacetonitrile (Key Feature)",
            "description": "Identifies acetonitriles halogenated on the alpha-carbon.",
            "smarts": "[#6]([Cl,Br,F,I])C#N"
        },
        "23.d.1": {
            "tier": 2,
            "name": "Nitrogen Mustard (Key Feature)",
            "description": "Identifies the bis(2-haloethyl)amine core of nitrogen mustards.",
            "smarts": "[#7](CC[Cl,Br])(CC[Cl,Br])"
        },
        "23.d.2": {
            "tier": 2,
            "name": "Sulfur Mustard (Key Feature)",
            "description": "Identifies the bis(2-haloethyl)sulfide core of sulfur mustards.",
            "smarts": "[#16;X2](CC[Cl,Br])CC[Cl,Br]"
        }
    },
    "Category 24: Poly-functional Short Aliphatic Chains and Chelators": {
        "24.a.1": {
            "tier": 1,
            "name": "1,2-Disubstituted Ethane (Core)",
            "description": "Identifies an ethane chain substituted on each carbon with a polar heteroatom group (O, N, or S), like ethylene glycol.",
            "smarts": "[O,N,S]CC[O,N,S]"
        },
        "24.b.1": {
            "tier": 1,
            "name": "Glycol Ether (Core)",
            "description": "Identifies the general glycol ether moiety (C-O-C-C-O).",
            "smarts": "[#6]O[#6][#6]O"
        },
        "24.b.2": {
            "tier": 2,
            "name": "Glycol Ether with Primary Alcohol (Key Feature)",
            "description": "Identifies glycol ethers that are primary alcohols, which can be metabolized to toxic alkoxyacetic acids.",
            "smarts": "[#6]O[#6][CH2][OH1]"
        },
        "24.b.3": {
            "tier": 1,
            "name": "Polyether (Core)",
            "description": "Identifies the repeating oxyethylene unit of poly(ethylene glycols), requiring at least three units.",
            "smarts": "OCCOCCO"
        },
        "24.c.1": {
            "tier": 2,
            "name": "Chelator - Aminopolycarboxylic Acid (Key Feature)",
            "description": "Identifies the aminopolycarboxylic acid scaffold of chelators like EDTA.",
            "smarts": "N(CC(=O)O)CC(=O)O"
        },
        "24.c.2": {
            "tier": 2,
            "name": "Chelator - Vicinal Dithiol (Key Feature)",
            "description": "Identifies the vicinal dithiol (-SH) scaffold of chelators like BAL and DMSA.",
            "smarts": "[#16;H1]C(C[#16;H1])"
        },
        "24.c.3": {
            "tier": 2,
            "name": "Chelator - Aminothiol (Key Feature)",
            "description": "Identifies the aminothiol scaffold of chelators like penicillamine.",
            "smarts": "[#16;H1]C(C)(C)C(N)C(=O)O"
        }
    },
    "Category 25: Short-Chain Alcohols and Nitriles": {
        "25.a.1": {
            "tier": 1,
            "name": "Short-Chain Primary Alcohol (Core)",
            "description": "Identifies short-chain (C1-C4) non-branched primary alcohols.",
            "smarts": "[CH2;D2;!R][OH1]"
        },
        "25.a.2": {
            "tier": 2,
            "name": "beta-Alkyl Substituted Alcohol (Key Feature)",
            "description": "Identifies the key feature of a primary alcohol where the beta-carbon is branched (tertiary or quaternary).",
            "smarts": "[#6;X3,X4]C[OH1]"
        },
        "25.b.1": {
            "tier": 1,
            "name": "Aliphatic Nitrile (Core)",
            "description": "Identifies the general aliphatic nitrile functional group.",
            "smarts": "[#6;!a]C#N"
        },
        "25.b.2": {
            "tier": 2,
            "name": "Vinyl Nitrile (Key Feature)",
            "description": "Identifies the acrylonitrile scaffold (C=C-C#N), a feature associated with higher potency.",
            "smarts": "[#6]=[#6]C#N"
        }
    }
}

flat_alert_lookup = {
    alert_id: {**info, 'smarts_pattern': Chem.MolFromSmarts(info['smarts']), 'category': alert_id.split('.')[0], 'subcategory': alert_id.split('.')[1]}
    for category, alerts in dart_alerts_dictionary.items() for alert_id, info in alerts.items() if Chem.MolFromSmarts(info['smarts'])
}


# In[3]:


# --- Alert SUB-MODULE 1: AMES MUTAGENICITY ---
def calculate_mutagenicity_score(target_smiles: str, surrogate_smiles: str) -> float:
    target_alerts = _get_ames_alerts(target_smiles)
    surrogate_alerts = _get_ames_alerts(surrogate_smiles)
    target_set, surrogate_set = set(target_alerts), set(surrogate_alerts)
    denominator = max(len(target_set), len(surrogate_set))
    if denominator == 0: return 1.0
    if not target_alerts or not surrogate_alerts: return 0.0
    return len(target_set.intersection(surrogate_set)) / denominator

def _get_ames_alerts(smiles: str) -> list[str]:
    df = _run_toxtree_module(smiles, "toxtree.plugins.ames.AmesMutagenicityRules")
    if df is None:
        return ["Error: Toxtree execution failed"]
        
    # Your existing mapping:
    ames_alert_lookup = {
        'SA1_Ames': 'Acyl halides',
        'SA2_Ames': 'Alkyl ester of sulphonic acid / phosphonic acid',
        'SA3_Ames': 'N-methylol derivatives',
        'SA4_Ames': 'Monohaloalkene',
        'SA5_Ames': 'S or N mustard',
        'SA6_Ames': 'propiolactones / propiosultones',
        'SA7_Ames': 'epoxides and aziridines',
        'SA8_Ames': 'Aliphatic halogens',
        'SA9_Ames': 'alkyl nitrite',
        'SA10_Ames': 'alpha, beta unsaturated carbonyls',
        'SA11_Ames': 'Simple aldehyde',
        'SA12_Ames': 'Quinones1 or Quinones2',
        'SA13_Ames': 'Hydrazine not N-N=[O,N]',
        'SA14_Ames': 'Aliphatic azo or azoxy',
        'SA15_Ames': 'Isocyanate and isothiocyanate groups',
        'SA16_Ames': 'Alkyl carbamate and thiocarbamate',
        'SA18_Ames': 'Polycyclic Aromatic Hydrocarbons',
        'SA19_Ames': 'Heterocyclic Polycyclic Aromatic Hydrocarbons',
        'SA21_Ames': 'Alkyl and aryl N-nitroso groups',
        'SA22_Ames': 'Azide of triazene',
        'SA23_Ames': 'Aliphatic N-nitro',
        'SA24_Ames': 'alpha, beta unsaturated alkoxy',
        'SA25_Ames': 'Aromatic nitroso group',
        'SA26_Ames': 'Aromatic ring N-oxide',
        'SA27_Ames': 'Nitro aromatic',
        'SA28_Ames': 'Primary aromatic amine, hydroxyl amine and its derived esters (with restrictions)',
        'SA28bis_Ames': 'Aromatic mono- and dialkylamine',
        'SA28ter_Ames': 'Aromatic N-acyl amine',
        'SA29_Ames': 'Aromatic diazo',
        'SA30_Ames': 'Coumarins and Furocoumarins',
        'SA37_Ames': 'Pyrrolizidine Alkaloids',
        'SA38_Ames': 'Alkenylbenzenes',
        'SA39_Ames': 'Steroidal estrogens',
        'SA57_Ames': 'DNA Intercalating Agents with a basic side chain',
        'SA58_Ames': 'Haloalkene cysteine S-conjugates',
        'SA59_Ames': 'Xanthones, Thioxanthones, Acridones',
        'SA60_Ames': 'Flavonoids',
        'SA61_Ames': 'Alkyl hydroperoxides',
        'SA62_Ames': 'N-acyloxy-N-alkoxybenzamides',
        'SA63_Ames': 'N-aryl-N-acetoxyacetamides',
        'SA64_Ames': 'Hydroxamic acid derivatives',
        'SA65_Ames': 'Halofuranones',
        'SA66_Ames': 'Anthrones',
        'SA67_Ames': 'Triphenylimidazole and related',
        'SA68_Ames': '9,10 - dihydrophenanthrenes',
        'SA69_Ames': 'Fluorinated quinolines'
    }

    try:
        mut_col = 'Structural Alert for S. typhimurium  mutagenicity'
        if str(df.iloc[0].get(mut_col)).strip().upper() == "YES":
            triggered = [
                name for col, name in ames_alert_lookup.items()
                if col in df.columns and str(df.iloc[0].get(col)).strip().upper() == "YES"
            ]
            return triggered if triggered else ["Alert identified (unspecified)"]
        return []
    except Exception:
        return ["Error parsing results"]



# In[4]:


# --- Alert SUB-MODULE 2: DART ALERTS ---
def calculate_dart_similarity_score(target_smiles: str, surrogate_smiles: str) -> float:
    target_alerts = _screen_for_dart_alerts(target_smiles)
    surrogate_alerts = _screen_for_dart_alerts(surrogate_smiles)
    dissimilarity = _calculate_structural_dissimilarity(target_alerts, surrogate_alerts)
    return 1.0 - dissimilarity

def _screen_for_dart_alerts(smiles: str) -> set:
    triggered = set()
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        for alert_id, info in flat_alert_lookup.items():
            if mol.HasSubstructMatch(info["smarts_pattern"]):
                triggered.add(alert_id)
    return triggered

def _calculate_structural_dissimilarity(target_alerts: set, surrogate_alerts: set) -> float:
    all_alerts = target_alerts.union(surrogate_alerts)
    if not all_alerts: return 0.0
    total_dissimilarity = 0
    for alert_id in all_alerts:
        in_target, in_surrogate = alert_id in target_alerts, alert_id in surrogate_alerts
        if in_target == in_surrogate: min_penalty = 0.0
        elif in_target: min_penalty = min([_calculate_penalty(alert_id, s_id) for s_id in surrogate_alerts] or [1.0])
        else: min_penalty = min([_calculate_penalty(alert_id, t_id) for t_id in target_alerts] or [1.0])
        total_dissimilarity += min_penalty
    return total_dissimilarity / len(all_alerts)

def _calculate_penalty(a1_id: str, a2_id: str) -> float:
    a1_info, a2_info = flat_alert_lookup[a1_id], flat_alert_lookup[a2_id]
    if a1_info['subcategory'] == a2_info['subcategory']: return 0.25
    if a1_info['category'] == a2_info['category']: return 0.5
    return 1.0


# In[5]:


# --- Alert SUB-MODULE 3: CRAMER PATH SIMILARITY ---
def calculate_cramer_path_score(target_smiles: str, surrogate_smiles: str) -> float:
    target_path, _ = _get_cramer_decision_path(target_smiles)
    surrogate_path, _ = _get_cramer_decision_path(surrogate_smiles)
    score, _ = _calculate_path_divergence_score(target_path, surrogate_path)
    return score

def _get_cramer_decision_path(smiles: str) -> tuple[str, str]:
    df = _run_toxtree_module(smiles, "toxtree.tree.cramer3.RevisedCramerDecisionTree")
    if df is None:
        return "Error", "Error"
    path = df.iloc[0].get("toxtree.tree.cramer3.CDTResult", "Path not found")
    klass = df.iloc[0].get("RevisedCDT", "Class not found")
    return str(path), str(klass)


def _calculate_path_divergence_score(path1: str, path2: str) -> tuple:
    if path1 == "Error" or path2 == "Error" or path1 == path2: return 1.0, "Identical Path"
    answers1, answers2 = path1.split(','), path2.split(',')
    for a1, a2 in zip(answers1, answers2):
        if a1 != a2:
            q_num_match = re.match(r'(\d+)', a1)
            if q_num_match:
                q_num = int(q_num_match.group(1))
                divergence_point = f"Q{q_num}"
                if q_num <= 10: return 0.0, divergence_point
                if q_num <= 20: return 0.5, divergence_point
                return 0.8, divergence_point
    return 0.9, "Path lengths differ"

# --- Generic Toxtree Helper ---
from pathlib import Path
import os, subprocess, tempfile, pandas as pd

# Detect Java for tools
JAVA_BIN_DEFAULT = os.environ.get("JAVA_BIN") or shutil.which("java")

# Allow a *separate* Java for Toxtree (falls back to default if not set)
JAVA_BIN_TOXTREE = (
    os.environ.get("JAVA_BIN_TOXTREE")
    or ("/usr/lib/jvm/java-11-openjdk-amd64/bin/java" if Path("/usr/lib/jvm/java-11-openjdk-amd64/bin/java").exists() else None)
    or JAVA_BIN_DEFAULT
)

def _run_toxtree_module(smiles: str, module_klass: str) -> pd.DataFrame | None:
    """
    Run a Toxtree module (headless) on one SMILES.
    CWD = TX_DIR so Toxtree can see ext/index.properties and all plugin JARs.
    """
    from pathlib import Path
    import tempfile, subprocess, pandas as pd

    # TX_DIR and TX_JAR should already be resolved as:
    # TX_DIR = /app/tools/toxtree/Toxtree-v3.1.0.1851/Toxtree
    # TX_JAR = /app/tools/toxtree/.../Toxtree-3.1.0.1851.jar

    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        in_csv  = tmp / "input.csv"
        out_csv = tmp / "output.csv"
        pd.DataFrame([{"SMILES": smiles}]).to_csv(in_csv, index=False)

        cmd = [
            JAVA_BIN_TOXTREE,                   # /opt/java/temurin-11/bin/java
            f"-Xmx{TT_HEAP_MB}m",
            "-Djava.awt.headless=true",
            "-jar", str(TX_JAR),
            "-n",
            "-i", str(in_csv),                  # absolute path
            "-o", str(out_csv),                 # absolute path
            "-m", module_klass,                 # e.g. toxtree.plugins.ames.AmesMutagenicityRules
        ]

        print(f"[TOXTREE] CMD: {' '.join(cmd)}  CWD={TX_DIR}")
        res = subprocess.run(cmd, cwd=str(TX_DIR), capture_output=True, text=True)

        if res.returncode != 0:
            print("[TOXTREE] returncode:", res.returncode)
            print("[TOXTREE] stdout:\n", res.stdout)
            print("[TOXTREE] stderr:\n", res.stderr)
            return None

        if not out_csv.exists():
            print(f"-> Toxtree did not produce output.csv at expected path: {out_csv}")
            print("[TOXTREE] stdout:\n", res.stdout)
            print("[TOXTREE] stderr:\n", res.stderr)
            return None

        try:
            return pd.read_csv(out_csv)
        except Exception as e:
            print("[TOXTREE] Failed to read output.csv:", e)
            print("[TOXTREE] stdout:\n", res.stdout)
            print("[TOXTREE] stderr:\n", res.stderr)
            return None
# In[6]:


# --- Alter Module Functions ---

def get_mutagenicity_results(target_smiles: str, surrogate_smiles: str) -> tuple:
    """Helper to get both the score and the alert lists for mutagenicity."""
    target_alerts = _get_ames_alerts(target_smiles) # Assumes this helper exists
    surrogate_alerts = _get_ames_alerts(surrogate_smiles)
    target_set, surrogate_set = set(target_alerts), set(surrogate_alerts)
    denominator = max(len(target_set), len(surrogate_set))
    if denominator == 0: score = 1.0
    elif not target_alerts or not surrogate_alerts: score = 0.0
    else: score = len(target_set.intersection(surrogate_set)) / denominator
    return score, target_alerts, surrogate_alerts

def get_dart_results(target_smiles: str, surrogate_smiles: str) -> tuple:
    """Helper to get the DART score and alert lists."""
    target_alerts_ids = _screen_for_dart_alerts(target_smiles) # Assumes this helper exists
    surrogate_alerts_ids = _screen_for_dart_alerts(surrogate_smiles)
    dissimilarity = _calculate_structural_dissimilarity(target_alerts_ids, surrogate_alerts_ids) # Assumes this helper exists
    score = 1.0 - dissimilarity
    target_names = [flat_alert_lookup[aid]['name'] for aid in target_alerts_ids]
    surrogate_names = [flat_alert_lookup[aid]['name'] for aid in surrogate_alerts_ids]
    return score, target_names, surrogate_names

def get_cramer_results(target_smiles: str, surrogate_smiles: str) -> tuple:
    """Helper to get the Cramer score and path details."""
    target_path, target_class = _get_cramer_decision_path(target_smiles) # Assumes this helper exists
    surrogate_path, surrogate_class = _get_cramer_decision_path(surrogate_smiles)
    score, divergence = _calculate_path_divergence_score(target_path, surrogate_path) # Assumes this helper exists
    return score, f"Class: {target_class}, Path: {target_path}", f"Class: {surrogate_class}, Path: {surrogate_path}", divergence

# --- 3. FINAL INTEGRATED ANALYSIS FUNCTION with Detailed Reporting ---
def run_structural_alert_analysis(target_smiles: str, surrogate_smiles: str) -> tuple:
    """
    Runs all three structural alert modules and returns the final score and
    all underlying detailed results as a tuple.
    """
    # --- Gather all results ---
    mutagenicity_score, target_ames, surrogate_ames = get_mutagenicity_results(target_smiles, surrogate_smiles)
    dart_score, target_dart, surrogate_dart = get_dart_results(target_smiles, surrogate_smiles)
    cramer_score, target_cramer, surrogate_cramer, cramer_divergence = get_cramer_results(target_smiles, surrogate_smiles)
    
    # --- Apply the 50/30/20 weighting ---
    final_score = (mutagenicity_score * 0.5) + (dart_score * 0.3) + (cramer_score * 0.2)

    # --- UPDATED: Return the individual values as a tuple ---
    return (
        final_score,
        mutagenicity_score,
        dart_score,
        cramer_score,
        target_ames,
        surrogate_ames,
        target_dart,
        surrogate_dart,
        cramer_divergence,
        target_cramer,
        surrogate_cramer
    )


# In[ ]:





# In[7]:


# --- Phys Chem Mod Functions ---
def is_voc(smiles_string):
    """Checks if a compound is a VOC based on carbon count (<= 12 carbons)."""
    if not isinstance(smiles_string, str): return False
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol: return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        return carbon_count > 0 and carbon_count <= 12
    except:
        return False

# --- PubChem API Helper Functions ---
def get_pubchem_by_name(compound_name: str) -> dict:
    """Retrieves properties from PubChem by name."""
    if not compound_name or pd.isna(compound_name):
        return {'molecular_weight': None, 'xlogp': None}
    encoded_name = quote(compound_name)
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
    mw_url = f"{base_url}/{encoded_name}/property/MolecularWeight/txt"
    xlogp_url = f"{base_url}/{encoded_name}/property/XLogP/txt"
    return _fetch_pubchem_props(mw_url, xlogp_url, compound_name)

def get_pubchem_by_smiles(smiles: str) -> dict:
    """Retrieves properties from PubChem by SMILES."""
    if not smiles or pd.isna(smiles):
        return {'molecular_weight': None, 'xlogp': None}
    encoded_smiles = quote(smiles)
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"
    mw_url = f"{base_url}/{encoded_smiles}/property/MolecularWeight/txt"
    xlogp_url = f"{base_url}/{encoded_smiles}/property/XLogP/txt"
    return _fetch_pubchem_props(mw_url, xlogp_url, smiles)

def _fetch_pubchem_props(mw_url, xlogp_url, identifier):
    """Generic fetcher for PubChem properties with a 5-second timeout."""
    results = {'molecular_weight': None, 'xlogp': None}
    headers = {'User-Agent': 'Python-Toxicity-Tool'}
    try:
        mw_response = requests.get(mw_url, timeout=5, headers=headers)
        if mw_response.status_code == 200:
            results['molecular_weight'] = float(mw_response.text.strip())
        
        time.sleep(0.2) # To respect PubChem's rate limits (max 5 requests/sec)

        xlogp_response = requests.get(xlogp_url, timeout=5, headers=headers)
        if xlogp_response.status_code == 200:
            results['xlogp'] = float(xlogp_response.text.strip())
    except (requests.exceptions.RequestException, ValueError):
        # Fail silently on network errors or conversion errors
        pass
    return results

# --- Main Comparison Function ---
def run_physicochemical_analysis(target_name, target_smiles, surrogate_name, surrogate_smiles):
    """
    Performs a 1-to-1 physicochemical comparison and returns the score
    and the properties of the target and surrogate.
    """
    print(f"\n--- Running Physicochemical Analysis ---")
    
    target_mol = Chem.MolFromSmiles(target_smiles)
    surr_mol = Chem.MolFromSmiles(surrogate_smiles)
    
    if not target_mol or not surr_mol:
        print("❌ Invalid SMILES provided. Cannot perform comparison.")
        # Return default values on failure
        return 0.0, {}, {}

    target_props = {
        'MW': Descriptors.MolWt(target_mol),
        'logP': Crippen.MolLogP(target_mol),
        'Charge': rdmolops.GetFormalCharge(target_mol),
        'is_VOC': is_voc(target_smiles)
    }
    surrogate_props = {
        'MW': Descriptors.MolWt(surr_mol),
        'logP': Crippen.MolLogP(surr_mol),
        'Charge': rdmolops.GetFormalCharge(surr_mol),
        'is_VOC': is_voc(surrogate_smiles)
    }

    # API Layer
    target_api = get_pubchem_by_name(target_name) or get_pubchem_by_smiles(target_smiles)
    surrogate_api = get_pubchem_by_name(surrogate_name) or get_pubchem_by_smiles(surrogate_smiles)
    
    if target_api.get('molecular_weight'): target_props['MW'] = target_api['molecular_weight']
    if target_api.get('xlogp'): target_props['logP'] = target_api['xlogp']
    if surrogate_api.get('molecular_weight'): surrogate_props['MW'] = surrogate_api['molecular_weight']
    if surrogate_api.get('xlogp'): surrogate_props['logP'] = surrogate_api['xlogp']

    # Scoring
    matches = 0
    if abs(target_props['MW'] - surrogate_props['MW']) <= 0.20 * target_props['MW']: matches += 1
    if abs(target_props['logP'] - surrogate_props['logP']) <= 1.0: matches += 1
    if (target_props['Charge'] == 0) == (surrogate_props['Charge'] == 0): matches += 1
    if target_props['is_VOC'] == surrogate_props['is_VOC']: matches += 1
        
    score_map = {4: 1.0, 3: 0.6, 2: 0.33}
    final_score = score_map.get(matches, 0.0)

    # --- ADDED: The missing return statement ---
    return final_score, target_props, surrogate_props


# In[ ]:





# In[8]:


# --- Metabolite Mod Functions ---
def _get_sygma_metabolites(smiles: str, score_cutoff: float) -> dict:
    """Internal function to run SyGMa and get high-plausibility metabolites."""
    metabolites = {}
    try:
        scenario = sygma.Scenario([
            [sygma.ruleset['phase1'], 1],
            [sygma.ruleset['phase2'], 1]
        ])
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return {}
        metabolic_tree = scenario.run(mol)
        metabolic_tree.calc_scores()
        metabolite_list = metabolic_tree.to_smiles()
        for metabolite_smi, score in metabolite_list[1:]:
            if score >= score_cutoff:
                metabolites[metabolite_smi] = score
        return metabolites
    except Exception:
        return {}

def _calculate_set_tanimoto(set1_smiles: list, set2_smiles: list) -> float:
    """Internal function to calculate unweighted set-based Tanimoto similarity."""
    if not set1_smiles or not set2_smiles: return 0.0
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=2048) for s in set1_smiles]
    fps2 = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=2048) for s in set2_smiles]
    sum_best_matches = 0
    for fp1 in fps1:
        similarities = BulkTanimotoSimilarity(fp1, fps2)
        sum_best_matches += max(similarities) if similarities else 0
    for fp2 in fps2:
        similarities = BulkTanimotoSimilarity(fp2, fps1)
        sum_best_matches += max(similarities) if similarities else 0
    denominator = len(fps1) + len(fps2)
    score = sum_best_matches / denominator if denominator > 0 else 0.0
    return score / 2.0

def _assign_score_value(raw_score: float) -> float:
    """Applies the established scoring rubric to the raw similarity score."""
    if raw_score >= 0.40:
        return 1.0
    elif raw_score >= 0.30:
        return 0.8
    elif raw_score >= 0.20:
        return 0.5
    else:
        return 0.0

# --- Main Production Function ---
#def run_metabolic_similarity_analysis(target_smiles: str, surrogate_smiles: str, cutoff: float = 0.01) -> pd.DataFrame:
#    """
#    Performs a complete metabolic similarity analysis between a target and
#    surrogate, returning a structured DataFrame with the results.
#    """
#    # 1. Get high-plausibility metabolites for both compounds
#    target_metabolites = _get_sygma_metabolites(target_smiles, score_cutoff=cutoff)
#    surrogate_metabolites = _get_sygma_metabolites(surrogate_smiles, score_cutoff=cutoff)
#    
#    target_smiles_list = list(target_metabolites.keys())
#    surrogate_smiles_list = list(surrogate_metabolites.keys())
#
#    # 2. Calculate the raw similarity score
#    raw_score = _calculate_set_tanimoto(target_smiles_list, surrogate_smiles_list)
#    
#    # 3. Assign the final binned score
#    assigned_value = _assign_score_value(raw_score)
#
#    return (
#        assigned_value,
#        len(target_smiles_list),
#        len(surrogate_smiles_list),
#        target_smiles_list,
#        surrogate_smiles_list
#    )


# In[ ]:





# In[9]:


# --- Reactive Metabolite Mod Functions---
# Define the reactive metabolite pathways from Lester et al. Table S1A
REACTIVE_PATHWAY_KEYWORDS = {
    "epoxidation",
    "glutathionation",
    "N-hydroxylation",
    "gamma-diketone",
    "diphenoquinone",
    "o-quinone",
    "o-quinone-methide",
    "o-quinoneimine",
    "p-quinone_methide",
    "p-quinoneimine"
}

def standardize_smiles(smiles: str) -> str:
    """Returns a standardized, neutralized, canonical tautomer SMILES."""
    if not isinstance(smiles, str): return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        clean_mol = rdMolStandardize.Cleanup(mol)
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
        uncharger = rdMolStandardize.Uncharger()
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        te = rdMolStandardize.TautomerEnumerator()
        canonical_mol = te.Canonicalize(uncharged_parent_clean_mol)
        return Chem.MolToSmiles(canonical_mol)
    except:
        return None

def run_biotransformer_and_get_reactions(smiles: str) -> set[str]:
    """
    Run BioTransformer (full distro required: JAR + config.json + KBs).
    CWD = BT_DIR so the app finds config.json, KBs, etc.
    """
    from pathlib import Path
    import tempfile, subprocess, csv
    reactions = set()

    standardized = standardize_smiles(smiles)
    if not standardized:
        return reactions

    m = Chem.MolFromSmiles(standardized)
    if not m:
        return reactions

    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        in_sdf  = tmp / "input.sdf"
        out_csv = tmp / "output.csv"

        w = Chem.SDWriter(str(in_sdf))
        w.write(m)
        w.close()

        cmd = [
            JAVA_BIN_BT,                     # e.g. /usr/bin/java (Java 17)
            f"-Xmx{BT_HEAP_MB}m",
            "-jar", str(BT_JAR),
            "-k", "pred",
            "-b", "allHuman",
            "-isdf", str(in_sdf),           # absolute
            "-ocsv", str(out_csv),          # absolute
            "-s", "2",
            "-cm", "3",
        ]
        print(f"[BT] CMD: {' '.join(cmd)}  CWD={BT_DIR}")
        res = subprocess.run(cmd, cwd=str(BT_DIR), capture_output=True, text=True)

        if res.returncode != 0:
            print("[BT] returncode:", res.returncode)
            print("[BT] stdout:\n", res.stdout)
            print("[BT] stderr:\n", res.stderr)
            return reactions

        if not out_csv.exists():
            print(f"-> BioTransformer did not produce output CSV at: {out_csv}")
            print("[BT] stdout:\n", res.stdout)
            print("[BT] stderr:\n", res.stderr)
            return reactions

        try:
            with out_csv.open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if "Reaction" in row and row["Reaction"]:
                        reactions.add(row["Reaction"])
        except Exception as e:
            print("[BT] Failed reading output CSV:", e)

    return reactions

def _extract_alerts(
    reactions: Iterable,
    reactive_keywords: Union[Iterable[str], Dict[str, Iterable[str]]]
) -> Set[str]:
    alerts: Set[str] = set()
    if not reactions:
        return alerts

    # Normalize reactive keywords
    if isinstance(reactive_keywords, dict):
        # canonical -> patterns
        canon_to_patterns = {str(k).strip().lower(): [str(p).lower() for p in v] for k, v in reactive_keywords.items()}
        for rxn in reactions:
            text = str(rxn or "").lower()
            if not text:
                continue
            for canon, pats in canon_to_patterns.items():
                if any(p in text for p in pats):
                    alerts.add(canon)
    else:
        # simple list/set of tokens; use each token as the canonical name
        tokens = [str(k).strip().lower() for k in reactive_keywords]
        for rxn in reactions:
            text = str(rxn or "").lower()
            if not text:
                continue
            for tok in tokens:
                if tok and tok in text:
                    alerts.add(tok)
    return alerts


def run_reactive_metabolite_analysis(target_smiles: str, surrogate_smiles: str):
    """
    Screens for reactive-metabolite pathways and returns:
      (score:int, target_alerts:List[str], surrogate_alerts:List[str])

    score = 1 iff the sets of canonical alert categories are identical (both empty counts as a match).
    """
    target_reactions = run_biotransformer_and_get_reactions(target_smiles)
    surrogate_reactions = run_biotransformer_and_get_reactions(surrogate_smiles)

    # Use your global REACTIVE_PATHWAY_KEYWORDS (either Iterable[str] or Dict[str, Iterable[str]])
    t_alerts = _extract_alerts(target_reactions, REACTIVE_PATHWAY_KEYWORDS)
    s_alerts = _extract_alerts(surrogate_reactions, REACTIVE_PATHWAY_KEYWORDS)

    # Binary match
    score = 1 if t_alerts == s_alerts else 0

    # Always return sorted lists (stable, easy to serialize)
    return score, sorted(t_alerts), sorted(s_alerts)

def run_metabolic_similarity_analysis(target_smiles: str, surrogate_smiles: str, cutoff: float = 0.01) -> tuple:
    """
    Performs a complete metabolic similarity analysis and returns the score
    and detailed results as a tuple.
    """
    # 1. Get high-plausibility metabolites for both compounds
    target_metabolites = _get_sygma_metabolites(target_smiles, score_cutoff=cutoff)
    surrogate_metabolites = _get_sygma_metabolites(surrogate_smiles, score_cutoff=cutoff)
    
    target_smiles_list = list(target_metabolites.keys())
    surrogate_smiles_list = list(surrogate_metabolites.keys())

    #New code to handle empty lists (no metabolites predicted)
    coverage_flag = "ok"

    if len(target_smiles_list) == 0 and len(surrogate_smiles_list) == 0:
        # Parent-as-metabolite on both sides
        target_smiles_list    = [target_smiles]      # weight 1.0 within its (singleton) set
        surrogate_smiles_list = [surrogate_smiles]   # ditto
        coverage_flag = "parent_fallback_both_empty"
    
    elif len(target_smiles_list) == 0:
        target_smiles_list = [target_smiles]
        coverage_flag = "parent_fallback_one_empty_target"
    
    elif len(surrogate_smiles_list) == 0:
        surrogate_smiles_list = [surrogate_smiles]
        coverage_flag = "parent_fallback_one_empty_surrogate"

    # 2. Calculate the raw similarity score
    raw_score = _calculate_set_tanimoto(target_smiles_list, surrogate_smiles_list)
    
    # 3. Assign the final binned score
    assigned_value = _assign_score_value(raw_score)

    # --- UPDATED: Return the individual values as a tuple ---
    return (
        assigned_value,
        len(target_smiles_list),
        len(surrogate_smiles_list),
        target_smiles_list,
        surrogate_smiles_list
    )


# In[10]:


#new version 2 metabolite prediction 

from typing import List, Tuple, Dict
import numpy as np
import sygma

from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rdMD, rdFMCS

def _get_sygma_metabolites(smiles: str, score_cutoff: float) -> dict:
    """Internal function to run SyGMa and get high-plausibility metabolites."""
    metabolites = {}
    try:
        scenario = sygma.Scenario([
            [sygma.ruleset['phase1'], 1],
            [sygma.ruleset['phase2'], 1]
        ])
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return {}
        metabolic_tree = scenario.run(mol)
        metabolic_tree.calc_scores()
        metabolite_list = metabolic_tree.to_smiles()
        for metabolite_smi, score in metabolite_list[1:]:
            if score >= score_cutoff:
                metabolites[metabolite_smi] = score
        return metabolites
    except Exception:
        return {}

# ---------- Standardization (safe fallbacks if unavailable) ----------
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _taut_enum = rdMolStandardize.TautomerEnumerator(
        rdMolStandardize.TautomerEnumeratorParams()
    )
    def _standardize(m: Chem.Mol) -> Chem.Mol:
        if m is None: return None
        lfs = rdMolStandardize.LargestFragmentChooser()
        m = lfs.choose(m)
        m = rdMolStandardize.Reionize(m)
        m = rdMolStandardize.Uncharger()(m)
        m = _taut_enum.Canonicalize(m)
        Chem.SanitizeMol(m)
        return m
except Exception:
    def _standardize(m: Chem.Mol) -> Chem.Mol:
        return m

def _mol_from_smiles(smi: str) -> Chem.Mol:
    return _standardize(Chem.MolFromSmiles(smi))

# ---------- Phase-II “aglycone” stripping (simple heuristics) ----------
SMARTS_CONJUGATES = [
    # O-glucuronide (rough; most common core)
    "[OX2;H0]-C1OC(O)C(O)C(O)C1=O",
    # sulfate (OSO3-)
    "OS(=O)(=O)[O-]",
    # glutathione (very rough motif)
    "NCC(=O)N[C@@H](CS)C(=O)O"
]
_CONJ_PATS = [Chem.MolFromSmarts(s) for s in SMARTS_CONJUGATES if Chem.MolFromSmarts(s)]

def _strip_conjugates(m: Chem.Mol) -> Chem.Mol:
    # For production, you may want a cleaner “cleave at linker” approach
    if m is None: return None
    to_delete = set()
    for pat in _CONJ_PATS:
        for match in m.GetSubstructMatches(pat):
            to_delete.update(match)
    if not to_delete:
        return m
    em = Chem.EditableMol(Chem.Mol(m))
    for idx in sorted(to_delete, reverse=True):
        try:
            em.RemoveAtom(idx)
        except Exception:
            pass
    out = em.GetMol()
    Chem.SanitizeMol(out, catchErrors=True)
    return out

# ---------- Fingerprints ----------
def _fp_ecfp_counts(m: Chem.Mol, radius=2, useFeatures=False):
    # Sparse count vector (collision-free-ish)
    return rdMD.GetMorganFingerprint(m, radius, useFeatures=useFeatures)

def _fp_ap_bitvect(m: Chem.Mol, nBits=4096):
    return rdMD.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nBits)

# Delta (metabolite – parent) counts: only “gains” to reflect added/modified local envs

def _delta_counts_dict(parent: Chem.Mol, metabolite: Chem.Mol, radius=2, useFeatures=False) -> dict[int,int]:
    fp_p = rdMD.GetMorganFingerprint(parent, radius, useFeatures=useFeatures)
    fp_m = rdMD.GetMorganFingerprint(metabolite, radius, useFeatures=useFeatures)
    p = fp_p.GetNonzeroElements()  # dict: {hash:int count}
    m = fp_m.GetNonzeroElements()
    # keep only “gains” to represent what metabolism added/changed
    d = {}
    for k, vm in m.items():
        vp = p.get(k, 0)
        if vm > vp:
            d[k] = vm - vp
    return d

# ---------- Pairwise similarity matrices ----------
def _sim_mat_ecfp(molsA: List[Chem.Mol], molsB: List[Chem.Mol], radius=2, useFeatures=False):
    fA = [_fp_ecfp_counts(m, radius, useFeatures) for m in molsA]
    fB = [_fp_ecfp_counts(m, radius, useFeatures) for m in molsB]
    S = np.zeros((len(fA), len(fB)), dtype=float)
    for i, fa in enumerate(fA):
        for j, fb in enumerate(fB):
            S[i, j] = DataStructs.TanimotoSimilarity(fa, fb)
    return S

def _sim_mat_ap(molsA: List[Chem.Mol], molsB: List[Chem.Mol], nBits=4096):
    fA = [_fp_ap_bitvect(m, nBits) for m in molsA]
    fB = [_fp_ap_bitvect(m, nBits) for m in molsB]
    S = np.zeros((len(fA), len(fB)), dtype=float)
    for i, fa in enumerate(fA):
        S[i, :] = DataStructs.BulkTanimotoSimilarity(fa, fB)
    return S

#Replaced to account for empty lists
#def _tani_counts_dict(d1: dict, d2: dict) -> float:
#    if not d1 and not d2:
#        return 0.0
#    keys = d1.keys() | d2.keys()
#    inter = 0
#   union = 0
#    for k in keys:
#        a = d1.get(k, 0)
#        b = d2.get(k, 0)
#        inter += a if a < b else b
#        union += a if a > b else b
#    return (inter / union) if union else 0.0

def _sim_mat_delta(parentA: Chem.Mol, metsA: list[Chem.Mol],
                   parentB: Chem.Mol, metsB: list[Chem.Mol],
                   radius=2, useFeatures=False):
    dA = [_delta_counts_dict(parentA, m, radius, useFeatures) for m in metsA]
    dB = [_delta_counts_dict(parentB, m, radius, useFeatures) for m in metsB]
    S = np.zeros((len(dA), len(dB)), dtype=float)
    for i, da in enumerate(dA):
        for j, db in enumerate(dB):
            S[i, j] = _tani_counts_dict(da, db)
    return S

def _mcs_tani_atoms(a: Chem.Mol, b: Chem.Mol, timeout=1):
    res = rdFMCS.FindMCS([a, b],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=timeout
    )
    if res.canceled or res.numAtoms == 0:
        return 0.0
    c = res.numAtoms
    return c / (a.GetNumAtoms() + b.GetNumAtoms() - c)

def _sim_mat_mcs(molsA: List[Chem.Mol], molsB: List[Chem.Mol], timeout=1):
    S = np.zeros((len(molsA), len(molsB)), dtype=float)
    for i, a in enumerate(molsA):
        for j, b in enumerate(molsB):
            S[i, j] = _mcs_tani_atoms(a, b, timeout=timeout)
    return S

# ---------- Set aggregators ----------
def _normalize_weights(plaus_scores: List[float], temp: float = 1.5) -> np.ndarray:
    x = np.array(plaus_scores, dtype=float)
    if x.size == 0:
        return x
    if temp != 1.0:
        x = x**temp
    s = x.sum()
    if s <= 0:
        return np.ones_like(x) / len(x)
    return x / s

def _chamfer_symmetric(S: np.ndarray, wA: np.ndarray, wB: np.ndarray) -> float:
    # S in [0,1], wA and wB sum to 1 each. Returns [0,1].
    if S.size == 0:
        return 0.0
    a_part = float((wA * S.max(axis=1)).sum()) if S.shape[0] else 0.0
    b_part = float((wB * S.max(axis=0)).sum()) if S.shape[1] else 0.0
    return 0.5 * (a_part + b_part)

def _assignment_score(S: np.ndarray, wA: np.ndarray, wB: np.ndarray) -> float:
    # Optional: exact 1:1 matching with SciPy, greedy otherwise.
    try:
        from scipy.optimize import linear_sum_assignment
        cost = 1.0 - S
        r, c = linear_sum_assignment(cost)
        if len(r) == 0:
            return 0.0
        pair_w = (wA[r] + wB[c]) / 2.0
        return float((pair_w * S[r, c]).sum())
    except Exception:
        # Greedy fallback
        usedA, usedB = set(), set()
        total = 0.0
        while True:
            best, bi, bj = -1.0, -1, -1
            for i in range(S.shape[0]):
                if i in usedA: continue
                for j in range(S.shape[1]):
                    if j in usedB: continue
                    if S[i, j] > best:
                        best, bi, bj = S[i, j], i, j
            if best < 0:
                break
            wt = 0.5 * (wA[bi] + wB[bj])
            total += wt * best
            usedA.add(bi); usedB.add(bj)
        return float(total)

# ---------- Top-level: compute fused metabolite-set similarity ----------
def compare_metabolite_sets(
    parentA_smiles: str,
    metsA: List[Tuple[str, float]],
    parentB_smiles: str,
    metsB: List[Tuple[str, float]],
    *,
    weight_temp: float = 1.75,
    aggregator: str = "chamfer",   # "chamfer" or "assignment"
    include_aglycone: bool = True,
    use_mcs: bool = True
) -> Dict[str, float]:
    """
    Returns a dict of component scores (0..1) and 'fused'.
    """
    # standardize parents
    pA = _mol_from_smiles(parentA_smiles)
    pB = _mol_from_smiles(parentB_smiles)
    if pA is None or pB is None:
        return {"ecfp": 0, "fcfp": 0, "ap": 0, "delta": 0, "mcs": 0, "fused": 0}

    # standardize metabolites, keep plausibility
    MA, wA = [], []
    for smi, sc in metsA:
        m = _mol_from_smiles(smi)
        if m is not None:
            MA.append(m); wA.append(float(sc))
    MB, wB = [], []
    for smi, sc in metsB:
        m = _mol_from_smiles(smi)
        if m is not None:
            MB.append(m); wB.append(float(sc))
    if len(MA) == 0 or len(MB) == 0:
        return {"ecfp": 0, "fcfp": 0, "ap": 0, "delta": 0, "mcs": 0, "fused": 0}

    wA = _normalize_weights(wA, temp=weight_temp)
    wB = _normalize_weights(wB, temp=weight_temp)

    # Make aglycone sets (optional)
    if include_aglycone:
        MA_ag = [_strip_conjugates(m) for m in MA]
        MB_ag = [_strip_conjugates(m) for m in MB]
    else:
        MA_ag, MB_ag = MA, MB

    # Pairwise matrices (full)
    S_ecfp  = _sim_mat_ecfp(MA, MB, radius=2, useFeatures=False)
    S_fcfp  = _sim_mat_ecfp(MA, MB, radius=2, useFeatures=True)
    S_ap    = _sim_mat_ap(MA, MB, nBits=4096)
    S_delta = _sim_mat_delta(pA, MA, pB, MB, radius=2, useFeatures=False)
    S_mcs   = _sim_mat_mcs(MA, MB, timeout=1) if use_mcs else np.zeros_like(S_ap)

    # Pairwise matrices (aglycone)
    if include_aglycone:
        S_ecfp_ag  = _sim_mat_ecfp(MA_ag, MB_ag, radius=2, useFeatures=False)
        S_fcfp_ag  = _sim_mat_ecfp(MA_ag, MB_ag, radius=2, useFeatures=True)
        S_ap_ag    = _sim_mat_ap(MA_ag, MB_ag, nBits=4096)
        S_delta_ag = _sim_mat_delta(pA, MA_ag, pB, MB_ag, radius=2, useFeatures=False)
        S_mcs_ag   = _sim_mat_mcs(MA_ag, MB_ag, timeout=1) if use_mcs else np.zeros_like(S_ap_ag)

    agg = _assignment_score if aggregator == "assignment" else _chamfer_symmetric

    # Aggregate (take max(full, aglycone) per component)
    comps = {}
    for key, S_full, S_ag in [
        ("ecfp",  S_ecfp,  S_ecfp_ag if include_aglycone else None),
        ("fcfp",  S_fcfp,  S_fcfp_ag if include_aglycone else None),
        ("ap",    S_ap,    S_ap_ag if include_aglycone else None),
        ("delta", S_delta, S_delta_ag if include_aglycone else None),
        ("mcs",   S_mcs,   S_mcs_ag if include_aglycone else None),
    ]:
        val_full = agg(S_full, wA, wB)
        if S_ag is not None:
            val_ag = agg(S_ag, wA, wB)
            comps[key] = max(val_full, val_ag)
        else:
            comps[key] = val_full

    # Simple fusion (tune on your calibration set)
    fused = (
        0.45*comps["ecfp"] +
        0.20*comps["fcfp"] +
        0.20*comps["ap"]   +
        0.10*comps["delta"]+
        0.05*comps["mcs"]
    )
    comps["fused"] = float(fused)
    return comps

# ---- Handle empty lists (no metabolites)
def _tani_counts_dict(d1: dict, d2: dict) -> float:
    # NEW: treat “both zero” (no change vs no change) as identical
    if not d1 and not d2:
        return 1.0
    keys = d1.keys() | d2.keys()
    inter = union = 0
    for k in keys:
        a = d1.get(k, 0)
        b = d2.get(k, 0)
        inter += a if a < b else b
        union  += a if a > b else b
    return (inter / union) if union else 0.0



def run_metabolic_similarity_analysis_v2(
    target_smiles: str,
    surrogate_smiles: str,
    cutoff: float = 0.01,
    *,
    weight_temp: float = 1.75,
    aggregator: str = "chamfer"
):
    # 1) Get high-plausibility metabolites (your existing function returns dict {smi: score})
    target_mets = _get_sygma_metabolites(target_smiles, score_cutoff=cutoff)
    sur_mets    = _get_sygma_metabolites(surrogate_smiles, score_cutoff=cutoff)

    # --- coverage-aware parent fallback (keep dict shape!) ---
    coverage_flag = "ok"
    if len(target_mets) == 0 and len(sur_mets) == 0:
        target_mets = {target_smiles: 1.0}
        sur_mets    = {surrogate_smiles: 1.0}
        coverage_flag = "parent_fallback_both_empty"
    elif len(target_mets) == 0:
        target_mets = {target_smiles: 1.0}
        coverage_flag = "parent_fallback_one_empty_target"
    elif len(sur_mets) == 0:
        sur_mets = {surrogate_smiles: 1.0}
        coverage_flag = "parent_fallback_one_empty_surrogate"

    # 2) Convert to [(smi, score), ...] for the comparator
    target_list    = list(target_mets.items())
    surrogate_list = list(sur_mets.items())

    # 3) Compare sets (fused and components in [0,1])
    comps = compare_metabolite_sets(
        parentA_smiles=target_smiles,
        metsA=target_list,
        parentB_smiles=surrogate_smiles,
        metsB=surrogate_list,
        weight_temp=weight_temp,
        aggregator=aggregator,
        include_aglycone=True,
        use_mcs=True
    )

    # 4) Bin the fused score
    fused = float(comps.get("fused", 0.0))
    assigned = _assign_score_value_v2(fused)

    # 5) Extract component scores (defensive .get() in case a key is missing)
    ecfp  = float(comps.get("ecfp",  0.0))
    fcfp  = float(comps.get("fcfp",  0.0))
    ap    = float(comps.get("ap",    0.0))
    delta = float(comps.get("delta", 0.0))
    mcs   = float(comps.get("mcs",   0.0))

    # Return tuple extended with component scores + fused (minimal API change)
    return (
        assigned,                      # 1
        len(target_list),              # 2
        len(surrogate_list),           # 3
        [s for s,_ in target_list],    # 4
        [s for s,_ in surrogate_list], # 5
        ecfp,                          # 6
        fcfp,                          # 7
        ap,                            # 8
        delta,                         # 9
        mcs,                           #10
        fused                          #11
        # (optionally: coverage_flag if you decide to expose it later)
    )


def _assign_score_value_v2(fused: float) -> float:
    # Provisional mapping 
    if fused >= 0.80: return 1.0
    if fused >= 0.60: return 0.8
    if fused >= 0.40: return 0.5
    return 0.0


# In[11]:


# --- Simple Tanimoto Similarity Mod Functions ---
def calculate_tanimoto_similarity(smiles1: str, smiles2: str) -> float:
    """
    Calculates the structural similarity between two SMILES strings using the
    Tanimoto coefficient based on Morgan fingerprints.

    Args:
        smiles1 (str): The SMILES string of the first molecule.
        smiles2 (str): The SMILES string of the second molecule.

    Returns:
        float: The Tanimoto similarity score (between 0.0 and 1.0),
               or 0.0 if either SMILES is invalid.
    """
    try:
        # Convert SMILES to RDKit molecule objects
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        if mol1 is None or mol2 is None:
            print("   -> Warning: Could not parse one or both SMILES strings.")
            return 0.0

        # Generate Morgan fingerprints for each molecule
        # Radius 2 is a standard setting for general-purpose similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

        # Calculate and return the Tanimoto similarity
        return TanimotoSimilarity(fp1, fp2)

    except:
        return 0.0


# In[13]:


def run_full_read_across_assessment(target_name: str, target_smiles: str, surrogate_name: str, surrogate_smiles: str) -> pd.DataFrame:
    """
    Executes all analysis modules and compiles a comprehensive, report-ready DataFrame.
    """
    # --- Run all individual modules to get scores and detailed results ---
    pchem_score, target_pchem_props, surrogate_pchem_props = run_physicochemical_analysis(target_name, target_smiles, surrogate_name, surrogate_smiles)
    (metabolic_score, target_met_count, surrogate_met_count, target_met_list, surrogate_met_list, ms_ecfp, ms_fcfp, ms_ap, ms_delta, ms_mcs, ms_fused) = run_metabolic_similarity_analysis_v2(target_smiles, surrogate_smiles)
    struct_score, mut_score, dart_score, cramer_score, target_ames, surrogate_ames, target_dart, surrogate_dart, cramer_divergence, target_cramer, surrogate_cramer= run_structural_alert_analysis(target_smiles, surrogate_smiles)
    reactive_metabolite_match, target_alerts, surrogate_alerts = run_reactive_metabolite_analysis(target_smiles, surrogate_smiles)
    tanimoto_score = calculate_tanimoto_similarity(target_smiles, surrogate_smiles)

    # --- Calculate Final Total Score ---
    total_score = pchem_score + metabolic_score + struct_score

    # --- Build the report row by row ---
    report = []
    report.append({'Parameter': 'Target Name', 'Target Result': target_name, 'Surrogate Result': ''})
    report.append({'Parameter': 'Target SMILES', 'Target Result': target_smiles, 'Surrogate Result': ''})
    report.append({'Parameter': 'Surrogate Name', 'Target Result': '', 'Surrogate Result': surrogate_name})
    report.append({'Parameter': 'Surrogate SMILES', 'Target Result': '', 'Surrogate Result': surrogate_smiles})

    report.append({'Parameter': '---', 'Target Result': '', 'Surrogate Result': ''})
    report.append({'Parameter': 'Generic Tanimoto Similarity', 'Target Result': tanimoto_score, 'Surrogate Result': tanimoto_score})
    report.append({'Parameter': '---', 'Target Result': '', 'Surrogate Result': ''})

    # Physicochemical Section
    report.append({'Parameter': 'P-Chem Module Score', 'Target Result': pchem_score, 'Surrogate Result': pchem_score})
    for prop, val in target_pchem_props.items(): report.append({'Parameter': f'  - {prop}', 'Target Result': val, 'Surrogate Result': surrogate_pchem_props.get(prop, 'N/A')})
    
    # Metabolic Similarity Section
    report.append({'Parameter': 'Metabolic Similarity Module Score', 'Target Result': metabolic_score, 'Surrogate Result': metabolic_score})
    report.append({'Parameter': '  - High-Plausibility Metabolite Count', 'Target Result': target_met_count, 'Surrogate Result': surrogate_met_count})

    # NEW: show metabolite SMILES lists (wrapped in Excel)
    report.append({'Parameter': '  - Target Metabolites (SMILES)', 'Target Result': '; '.join(target_met_list) if target_met_list else 'None', 'Surrogate Result': ''})
    report.append({'Parameter': '  - Surrogate Metabolites (SMILES)', 'Target Result': '', 'Surrogate Result': '; '.join(surrogate_met_list) if surrogate_met_list else 'None'})

    
    # NEW: component scores
    report.append({'Parameter': '  - ECFP',  'Target Result': ms_ecfp,  'Surrogate Result': ms_ecfp})
    report.append({'Parameter': '  - FCFP',  'Target Result': ms_fcfp,  'Surrogate Result': ms_fcfp})
    report.append({'Parameter': '  - AP',    'Target Result': ms_ap,    'Surrogate Result': ms_ap})
    report.append({'Parameter': '  - Delta', 'Target Result': ms_delta, 'Surrogate Result': ms_delta})
    report.append({'Parameter': '  - MCS',   'Target Result': ms_mcs,   'Surrogate Result': ms_mcs})
    report.append({'Parameter': '  - Fused', 'Target Result': ms_fused, 'Surrogate Result': ms_fused})
    
    # Structural Alerts Section
    report.append({'Parameter': 'Structural Alert Module Score', 'Target Result': struct_score, 'Surrogate Result': struct_score})
    report.append({'Parameter': '  - Mutagenicity Similarity Score (50%)', 'Target Result': mut_score, 'Surrogate Result': mut_score})
    report.append({'Parameter': '    - Target Ames Alerts', 'Target Result': ', '.join(target_ames) or 'None', 'Surrogate Result': ''})
    report.append({'Parameter': '    - Surrogate Ames Alerts', 'Target Result': '', 'Surrogate Result': ', '.join(surrogate_ames) or 'None'})
    report.append({'Parameter': '  - DART Similarity Score (30%)', 'Target Result': dart_score, 'Surrogate Result': dart_score})
    report.append({'Parameter': '    - Target DART Alerts', 'Target Result': ', '.join(target_dart) or 'None', 'Surrogate Result': ''})
    report.append({'Parameter': '    - Surrogate DART Alerts', 'Target Result': '', 'Surrogate Result': ', '.join(surrogate_dart) or 'None'})
    report.append({'Parameter': '  - Cramer Path Similarity Score (20%)', 'Target Result': cramer_score, 'Surrogate Result': cramer_score})
    report.append({'Parameter': '    - Cramer Path Divergence', 'Target Result': cramer_divergence, 'Surrogate Result': cramer_divergence})

    # NEW: explicit Cramer class + divergence point
    report.append({'Parameter': '    - Target Cramer Class', 'Target Result': target_cramer, 'Surrogate Result': ''})
    report.append({'Parameter': '    - Surrogate Cramer Class', 'Target Result': '', 'Surrogate Result': surrogate_cramer})
    report.append({'Parameter': '    - Cramer Path Divergence Point', 'Target Result': cramer_divergence, 'Surrogate Result': cramer_divergence})

    # Reactive Metabolite Section
    report.append({'Parameter': 'Reactive Metabolite Match (Binary)', 'Target Result': reactive_metabolite_match, 'Surrogate Result': reactive_metabolite_match})
    report.append({'Parameter': '  - Target Reactive Metabolite Alerts', 'Target Result': ', '.join(target_alerts) if target_alerts else 'None', 'Surrogate Result': ''})
    report.append({'Parameter': '  - Surrogate Reactive Metabolite Alerts', 'Target Result': '', 'Surrogate Result': ', '.join(surrogate_alerts) if surrogate_alerts else 'None'})

    report.append({'Parameter': '---', 'Target Result': '', 'Surrogate Result': ''})
    report.append({'Parameter': 'TOTAL SCORE (Sum of Modules)', 'Target Result': total_score, 'Surrogate Result': total_score})

    return pd.DataFrame(report).set_index('Parameter')


# In[14]:


import pandas as pd
import re
import xlsxwriter
from typing import List, Tuple, Optional

# ---------- Helpers to parse your vertical DF into sections ----------



_TOP_LEVEL_ANCHORS = {
    "P-Chem Module Score": "Phys Chem",
    "Metabolic Similarity Module Score": "Metabolism",
    "Structural Alert Module Score": "Structural Alerts",
    "Reactive Metabolite Match (Binary)": "Reactive Metabolites",
    "TOTAL SCORE (Sum of Modules)": "Total",
}

_IDENTITY_KEYS = [
    "Target Name",
    "Target SMILES",
    "Surrogate Name",
    "Surrogate SMILES",
]

def _is_voc_label(label: str) -> bool:
    s = str(label or "")
    s = re.sub(r"^[\-\s]+", "", s).strip().lower()  # strip leading bullets/dashes
    return s in {"is_voc", "is voc", "voc", "voc (0/1)"}

def _to_01(v):
    # Normalize a value to 0/1 for VOC display; leave empty as ""
    if v is None or (isinstance(v, str) and v.strip() == ""):
        return ""
    # Booleans
    if isinstance(v, bool):
        return 1 if v else 0
    # Numbers
    if isinstance(v, (int, float)):
        # clamp any numeric to {0,1}
        return 1 if float(v) >= 0.5 else 0
    # Strings
    s = str(v).strip().lower()
    if s in {"1", "true", "yes", "y", "t"}:
        return 1
    if s in {"0", "false", "no", "n", "f", "none"}:
        return 0
    # fallback: not a recognizable token → leave original (writer will print text)
    return v

def _is_break(param: str) -> bool:
    return param.strip() == "---"

def _looks_like_subitem(param: str) -> bool:
    # Your DF uses leading spaces + dashes for detail lines (e.g., "  - ...", "    - ...")
    return bool(re.match(r"^\s+-", param))

def _df_to_sections(df_vertical: pd.DataFrame) -> list:
    """
    Convert your vertical DF (index='Parameter', columns=['Target Result','Surrogate Result'])
    into sections suitable for writing. Returns a list of:
      [ (section_title, [(label, target_value, surrogate_value), ...]), ... ]
    """
    # Ensure 'Parameter' is the index
    if df_vertical.index.name != "Parameter":
        if "Parameter" in df_vertical.columns:
            df_vertical = df_vertical.set_index("Parameter", drop=True)
        else:
            raise ValueError("Expected a 'Parameter' column or index in the DataFrame.")

    # Find the value columns (case-insensitive)
    colmap = {c.lower(): c for c in df_vertical.columns}
    if "target result".lower() not in colmap or "surrogate result".lower() not in colmap:
        raise ValueError("Expected columns 'Target Result' and 'Surrogate Result' in the DataFrame.")
    tcol = colmap["target result"]
    scol = colmap["surrogate result"]

    sections = []

    # Identity section
    ident_rows = []
    for k in _IDENTITY_KEYS:
        if k in df_vertical.index:
            row = df_vertical.loc[k]
            ident_rows.append((k, row.get(tcol, ""), row.get(scol, "")))
    if ident_rows:
        sections.append(("Identity", ident_rows))

    # Generic Tanimoto (if present)
    if "Generic Tanimoto Similarity" in df_vertical.index:
        row = df_vertical.loc["Generic Tanimoto Similarity"]
        sections.append((
            "Generic Tanimoto",
            [("Tanimoto (parent)", row.get(tcol, ""), row.get(scol, ""))]
        ))

    # Module sections (walk rows in order)
    current = None
    bucket = []

    def flush():
        nonlocal current, bucket, sections
        if current and bucket:
            sections.append((current, bucket))
        current, bucket = None, []

    for param, row in df_vertical.iterrows():  # <-- iterrows fixes the unpack error
        if param in _IDENTITY_KEYS or param == "Generic Tanimoto Similarity" or _is_break(param):
            continue

        if param in _TOP_LEVEL_ANCHORS:
            flush()
            current = _TOP_LEVEL_ANCHORS[param]
            # top-line “Module Score” → show as Final Score
            bucket.append(("Final Score", row.get(tcol, ""), row.get(scol, "")))
            continue

        if param.startswith("TOTAL SCORE"):
            flush()
            current = "Total"
            bucket.append(("Total Score", row.get(tcol, ""), row.get(scol, "")))
            continue

        if current:
            label = param.strip()
            bucket.append((label, row.get(tcol, ""), row.get(scol, "")))

    flush()
    return sections

# ---------- Excel writer ----------

def _write_vertical_sheet(
    wb: xlsxwriter.Workbook,
    sheet_name: str,
    sections: list,
    col_widths=(26, 48, 48)
):
    ws = wb.add_worksheet(sheet_name[:31])

    # Formats
    f_hdr = wb.add_format({"bold": True, "font_size": 12, "bottom": 1})
    f_sec = wb.add_format({"bold": True, "bg_color": "#EFEFEF", "border": 1})
    f_lab = wb.add_format({"bold": True})
    f_txt = wb.add_format({"text_wrap": True})
    f_num = wb.add_format({"num_format": "0.000"})

    # Columns
    ws.set_column(0, 0, col_widths[0])
    ws.set_column(1, 1, col_widths[1])
    ws.set_column(2, 2, col_widths[2])

    # Header
    ws.write(0, 0, "", f_hdr)
    ws.write(0, 1, "Target", f_hdr)
    ws.write(0, 2, "Surrogate", f_hdr)

    r = 1
    for sec_title, rows in sections:
        ws.merge_range(r, 0, r, 2, sec_title, f_sec)
        r += 1
        for label, tval, sval in rows:
            pretty_label = label.strip()
            ws.write(r, 0, pretty_label, f_lab)

            # Writers
            def write_val(c, v):
                if v is None or v == "":
                    ws.write(r, c, "", f_txt)
                else:
                    try:
                        vv = float(v)
                        ws.write_number(r, c, vv, f_num)
                    except Exception:
                        ws.write(r, c, str(v), f_txt)

            # Special case: is_VOC → 0/1
            if _is_voc_label(pretty_label):
                tv = _to_01(tval)
                sv = _to_01(sval)
                # write as number if int-like, else as text (fallback)
                if isinstance(tv, (int, float)):
                    ws.write_number(r, 1, tv, wb.add_format({"num_format": "0"}))
                else:
                    ws.write(r, 1, str(tv), f_txt)
                if isinstance(sv, (int, float)):
                    ws.write_number(r, 2, sv, wb.add_format({"num_format": "0"}))
                else:
                    ws.write(r, 2, str(sv), f_txt)
            else:
                write_val(1, tval)
                write_val(2, sval)

            r += 1
        r += 1  # blank line between sections

def create_excel_report(
    pairs: List[Tuple[str, str, str, str]],
    out_path: str,
    *,
    sheet_name_fmt: Optional[str] = "{i:02d} - {tname} vs {sname}"
):
    """
    Run your full assessment for each (target_name, target_smiles, surrogate_name, surrogate_smiles),
    then write a vertically formatted, side-by-side Excel report.

    Parameters
    ----------
    pairs : list of tuples
        [(target_name, target_smiles, surrogate_name, surrogate_smiles), ...]
    out_path : str
        Path to the .xlsx file to write.
    sheet_name_fmt : str
        Optional format for sheet names (31 char limit applies after formatting).
        Tokens: {i}, {tname}, {sname}
    """
    wb = xlsxwriter.Workbook(out_path)

    for i, (tname, tsmi, sname, ssmi) in enumerate(pairs, start=1):
        # Run your existing end-to-end function (returns vertical DF)
        df_vertical = run_full_read_across_assessment(tname, tsmi, sname, ssmi)

        # Convert to sections
        sections = _df_to_sections(df_vertical)

        # Build a safe sheet name
        safe_t = (tname or "Target")[:15]
        safe_s = (sname or "Surrogate")[:15]
        sheet_name = sheet_name_fmt.format(i=i, tname=safe_t, sname=safe_s)

        # Write sheet
        _write_vertical_sheet(wb, sheet_name, sections)

    wb.close()

