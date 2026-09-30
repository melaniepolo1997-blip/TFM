# %%
#======================
#Cargar librerías 
#======================

import pandas as pd
import requests
import time
import re
from Bio import Entrez
import numpy as np
from pathlib import Path


Entrez.email = "melaniepolo1997@gmail.com"

#%%
Base_DIR = Path(__file__).resolve().parent

Data_DIR = Base_DIR / "data"
Results_DIR = Base_DIR / "results"

Results_DIR.mkdir(exist_ok=True)

# %%
ruta = Base_DIR/"Organismo_referencia_positivos y negativos.xlsx"

organismos_referencia = pd.read_excel(ruta)

print("Archivo cargado correctamente")
print("Filas:", len(organismos_referencia))
print("Columnas:", organismos_referencia.columns.tolist())

print("\nPrimeros 10 organismos:")
print(organismos_referencia.head(10))
# %%
print("Total de organismos:", len(organismos_referencia))

print("\nValores nulos:")
print(organismos_referencia.isna().sum())

print("\nLabels:")
print(organismos_referencia["label"].value_counts(dropna=False))

print("\nOrganismos duplicados:")
duplicados = organismos_referencia[
    organismos_referencia["Organism"].duplicated(keep=False)
]

print(duplicados)
# %%

#===================
#Búsqeuda Assembly 
#===================
def buscar_assemblies(organismo, retmax=50):
    
    try:
        handle = Entrez.esearch(
            db="assembly",
            term=f'"{organismo}"[Organism]',
            retmax=retmax
        )
        
        resultado = Entrez.read(handle)
        handle.close()
        
        return resultado["IdList"]
    
    except Exception as e:
        print(f"Error buscando {organismo}: {e}")
        return []
    
# %%
#================================
#Obtención de datos en Assembly 
#================================
def obtener_datos_assembly(assembly_id):
    
    try:
        handle = Entrez.esummary(
            db="assembly",
            id=assembly_id,
            report="full"
        )
        
        resultado = Entrez.read(handle)
        handle.close()
        
        return resultado["DocumentSummarySet"]["DocumentSummary"][0]
    
    except Exception as e:
        print(f"Error con Assembly {assembly_id}: {e}")
        return None


def extraer_datos_assembly(assembly_id):
    
    datos = obtener_datos_assembly(assembly_id)
    
    if datos is None:
        return None
    
    return {
        "Assembly_ID": assembly_id,
        "Assembly_accession": datos.get("AssemblyAccession", ""),
        "Organism_NCBI": datos.get("Organism", ""),
        "Species_NCBI": datos.get("SpeciesName", ""),
        "Taxid": datos.get("Taxid", ""),
        "Assembly_status": datos.get("AssemblyStatus", ""),
        "RefSeq": datos.get("Synonym", {}).get("RefSeq", ""),
        "GenBank": datos.get("Synonym", {}).get("Genbank", ""),
        "Contig_N50": datos.get("ContigN50", ""),
        "Scaffold_N50": datos.get("ScaffoldN50", ""),
        "Ftp_RefSeq": datos.get("FtpPath_RefSeq", ""),
        "Ftp_GenBank": datos.get("FtpPath_GenBank", "")
    }
# %%
resultados_todos = []

for i, organismo in enumerate(organismos_referencia["Organism"], start=1):
    
    print(f"[{i}/{len(organismos_referencia)}] {organismo}")
    
    ids = buscar_assemblies(organismo, retmax=50)
    
    if len(ids) == 0:
        
        resultados_todos.append({
            "Organism_reference": organismo,
            "Assembly_ID": None,
            "Assembly_accession": None,
            "Organism_NCBI": None,
            "Species_NCBI": None,
            "Taxid": None,
            "Assembly_status": None,
            "RefSeq": None,
            "GenBank": None,
            "Contig_N50": None,
            "Scaffold_N50": None,
            "Ftp_RefSeq": None,
            "Ftp_GenBank": None
        })
    
    else:
        
        for assembly_id in ids:
            
            datos = extraer_datos_assembly(assembly_id)
            
            if datos is not None:
                
                resultados_todos.append({
                    "Organism_reference": organismo,
                    **datos
                })
            
            time.sleep(0.35)
    time.sleep(0.35)
    
# %%
tabla_assemblies = pd.DataFrame(resultados_todos)

print("Filas:", len(tabla_assemblies))
print("Organismos:", tabla_assemblies["Organism_reference"].nunique())

tabla_assemblies.head()

# %%
resumen_assemblies = (
    tabla_assemblies
    .groupby("Organism_reference", dropna=False)
    .agg(
        N_Assembly=("Assembly_ID", lambda x: x.notna().sum())
    )
    .reset_index()
)

print(resumen_assemblies.to_string(index=False))
# %%
sin_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 0
]

print("ORGANISMOS SIN ASSEMBLY:")
print(len(sin_assembly))

print(sin_assembly.to_string(index=False))
# %%
un_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 1
]

print("Organismos con exactamente 1 Assembly:")
print(len(un_assembly))

print(un_assembly.to_string(index=False))

# %%
multiples_assemblies = resumen_assemblies[
    resumen_assemblies["N_Assembly"] > 1
]

print("Organismos con múltiples Assembly:")
print(len(multiples_assemblies))

print(multiples_assemblies.to_string(index=False))

# %%
def normalizar_nombre(nombre):
    
    if pd.isna(nombre):
        return ""
    
    nombre = str(nombre).lower()
    
    # eliminar información entre paréntesis
    nombre = re.sub(r"\([^)]*\)", "", nombre)
    
    # normalizar espacios
    nombre = re.sub(r"\s+", " ", nombre).strip()
    
    return nombre
# %%
tabla_assemblies["Organism_ref_norm"] = (
    tabla_assemblies["Organism_reference"]
    .apply(normalizar_nombre)
)

tabla_assemblies["Organism_ncbi_norm"] = (
    tabla_assemblies["Organism_NCBI"]
    .apply(normalizar_nombre)
)
# %%
def clasificar_match(row):
    
    if pd.isna(row["Assembly_ID"]):
        return "SIN_ASSEMBLY"
    
    ref = row["Organism_ref_norm"]
    ncbi = row["Organism_ncbi_norm"]
    
#Coincidencia exacta
    if ref == ncbi:
        return "VALIDO"
    
#El nombre de referencia aparece dentro del nombre NCBI
    if ref in ncbi:
        return "VALIDO"
    
#El nombre NCBI aparece dentro de la referencia
    if ncbi in ref and ncbi != "":
        return "REVISAR"
    
    return "REVISAR"


tabla_assemblies["Match_status"] = (
    tabla_assemblies.apply(clasificar_match, axis=1)
)
# %%
tabla_assemblies.to_excel(
    "tabla_assemblies_todos_organismos.xlsx",
    index=False
)

print("Archivo guardado.")
print("Total organismos:", organismos_referencia["Organism"].nunique())
print(
    resumen_assemblies["N_Assembly"]
    .value_counts()
    .sort_index()
)
print(
    resumen_assemblies["N_Assembly"]
    .value_counts()
    .sort_index()
)
# %%
pd.set_option("display.max_rows", 200)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)

print(
    tabla_assemblies[
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank"
        ]
    ].to_string(index=False)
)
# %%
casos_multiples = (
    tabla_assemblies[
        tabla_assemblies["Organism_reference"].isin(
            multiples_assemblies["Organism_reference"]
        )
    ]
    [
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank"
        ]
    ]
)

print(casos_multiples.to_string(index=False))

# %%
duplicados_referencia = (
    organismos_referencia[
        organismos_referencia["Organism"].duplicated(keep=False)
    ]
    .sort_values("Organism")
)

print(duplicados_referencia.to_string(index=False))
print("Organismos totales:", len(organismos_referencia))
print("Organismos únicos:", organismos_referencia["Organism"].nunique())
# %%
print("Filas:", len(tabla_assemblies))
print("Columnas:")
print(tabla_assemblies.columns.tolist())

print("\nOrganismos de referencia:")
print(tabla_assemblies["Organism_reference"].nunique())

# %%
resumen_assemblies = (
    tabla_assemblies
    .groupby("Organism_reference", dropna=False)
    .agg(
        N_Assembly=("Assembly_ID", lambda x: x.notna().sum())
    )
    .reset_index()
)

print(resumen_assemblies.to_string(index=False))
# %%
sin_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 0
]

un_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 1
]

multiples_assemblies = resumen_assemblies[
    resumen_assemblies["N_Assembly"] > 1
]

print("SIN ASSEMBLY:", len(sin_assembly))
print("1 ASSEMBLY:", len(un_assembly))
print("MÚLTIPLES ASSEMBLY:", len(multiples_assemblies))

#HASTA AQUI BÚSQUEDA DE ASSEMBLY
# %%
import unicodedata

def normalizar_nombre(nombre):
    if pd.isna(nombre):
        return ""
    
    nombre = str(nombre)
    
#Normalizar caracteres Unicode
    nombre = unicodedata.normalize("NFKC", nombre)
    
#Unificar diferentes tipos de guion
    nombre = re.sub(r"[‐-‒–—−]", "-", nombre)
    
#Eliminar espacios repetidos
    nombre = re.sub(r"\s+", " ", nombre).strip()
    
#Mayúsculas para comparar
    return nombre.upper()


organismos_referencia["Organism_normalizado"] = (
    organismos_referencia["Organism"]
    .apply(normalizar_nombre)
)

duplicados_normalizados = (
    organismos_referencia[
        organismos_referencia["Organism_normalizado"].duplicated(keep=False)
    ]
    .sort_values("Organism_normalizado")
)

print(duplicados_normalizados[
    ["Organism", "label", "Reference"]
].to_string(index=False))

# %%
import pandas as pd

organismos_referencia_limpios = (
    organismos_referencia
    .groupby(["Organism_normalizado", "label"], as_index=False)
    .agg(
        Organism=("Organism", "first"),
        Reference=(
            "Reference",
            lambda x: "; ".join(
                x.dropna().astype(str).unique()
            )
        )
    )
)

#Reordenar columnas
organismos_referencia_limpios = organismos_referencia_limpios[
    ["Organism", "label", "Reference"]
]

print("Organismos antes:", len(organismos_referencia))
print("Organismos después:", len(organismos_referencia_limpios))

print("\nResultado:")
print(organismos_referencia_limpios.to_string(index=False))
# %%
print(
    organismos_referencia_limpios[
        organismos_referencia_limpios["Organism"].duplicated(keep=False)
    ]
)
# %%
print("Organismos referencia limpios:", len(organismos_referencia_limpios))
print("Organismos con Assembly buscado:", tabla_assemblies["Organism_reference"].nunique())
# %%
organismos_sin_assembly = organismos_referencia_limpios[
    ~organismos_referencia_limpios["Organism"].isin(
        tabla_assemblies["Organism_reference"].unique()
    )
]

print("Organismos sin coincidencia en tabla_assemblies:")
print(len(organismos_sin_assembly))

print(
    organismos_sin_assembly[
        ["Organism", "label", "Reference"]
    ].to_string(index=False)
)
# %%
resumen_assembly = (
    tabla_assemblies
    .groupby("Organism_reference")
    .agg(
        N_Assembly=("Assembly_accession", "count"),
        Organismos_NCBI=(
            "Organism_NCBI",
            lambda x: " | ".join(x.dropna().astype(str).unique())
        ),
        Estados=(
            "Assembly_status",
            lambda x: " | ".join(x.dropna().astype(str).unique())
        )
    )
    .reset_index()
)

print(resumen_assembly.to_string(index=False))
# %%
un_assembly = resumen_assembly[
    resumen_assembly["N_Assembly"] == 1
]

print("Organismos con un único Assembly:", len(un_assembly))

print(
    un_assembly.to_string(index=False)
)
# %%
multiples_assembly = resumen_assembly[
    resumen_assembly["N_Assembly"] > 1
]

print("Organismos con múltiples Assembly:", len(multiples_assembly))

print(
    multiples_assembly.to_string(index=False)
)
# %%
un_assembly = resumen_assembly[
    resumen_assembly["N_Assembly"] == 1
].copy()

print("Total:", len(un_assembly))
print(un_assembly.to_string(index=False))

# %%
casos_un_assembly = tabla_assemblies[
    tabla_assemblies["Organism_reference"].isin(
        un_assembly["Organism_reference"]
    )
].copy()

#Eliminar duplicados
casos_un_assembly = casos_un_assembly.drop_duplicates(
    subset=["Organism_reference", "Assembly_accession"]
)

print("Filas:", len(casos_un_assembly))

print(
    casos_un_assembly[
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank"
        ]
    ].to_string(index=False)
)

#%%
#Eliminar duplicados de Assembly

tabla_assemblies_unica = (
    tabla_assemblies
    .drop_duplicates(
        subset=["Organism_reference", "Assembly_accession"]
    )
    .copy()
)

print("Registros únicos de Assembly:", len(tabla_assemblies_unica))
#%%
organismos_genomas = tabla_assemblies_unica.merge(
    organismos_referencia_limpios[
        ["Organism", "label", "Reference"]
    ],
    left_on="Organism_reference",
    right_on="Organism",
    how="inner"
)

print("Registros después del cruce:", len(organismos_genomas))
#%%
organismos_genomas["Coincidencia_nombre"] = (
    organismos_genomas["Organism_reference"]
    .str.lower()
    .str.replace("‐", "-", regex=False)
    .apply(
        lambda x: x.split(" ")[0] if isinstance(x, str) else ""
    )
)
#%%
multiples_detalle = tabla_assemblies_unica[
    tabla_assemblies_unica["Organism_reference"].isin(
        multiples_assembly["Organism_reference"]
    )
].copy()

print("Registros:", len(multiples_detalle))

print(
    multiples_detalle[
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank"
        ]
    ].to_string(index=False)
)

#%% 
#Tabla de candidatos

#Eliminar duplicados de la búsqueda
tabla_assemblies_unica = (
    tabla_assemblies
    .drop_duplicates(
        subset=["Organism_reference", "Assembly_accession"]
    )
    .copy()
)

#Unir con la información de Label
tabla_candidatos = tabla_assemblies_unica.merge(
    organismos_referencia_limpios[
        ["Organism", "label", "Reference"]
    ],
    left_on="Organism_reference",
    right_on="Organism",
    how="inner"
)

print("Registros Assembly:", len(tabla_assemblies_unica))
print("Registros después del cruce:", len(tabla_candidatos))
print(
    "Organismos únicos:",
    tabla_candidatos["Organism_reference"].nunique()
)
# %% 
#Contar ASSEMBLY por organismo

conteo_assembly = (
    tabla_candidatos
    .groupby("Organism_reference")["Assembly_accession"]
    .nunique()
    .reset_index(name="N_Assembly")
)

print(conteo_assembly["N_Assembly"].value_counts().sort_index())
# %% 
#Separar grupos

grupo_unico = conteo_assembly[
    conteo_assembly["N_Assembly"] == 1
]

grupo_multiple = conteo_assembly[
    conteo_assembly["N_Assembly"] > 1
]

print("Organismos con 1 Assembly:", len(grupo_unico))
print("Organismos con múltiples Assembly:", len(grupo_multiple))
# %% 
#Orgnaismos con múltiples ASSEMBLY

lista_multiples = (
    tabla_candidatos[
        tabla_candidatos["Organism_reference"].isin(
            grupo_multiple["Organism_reference"]
        )
    ]
    .sort_values(
        ["Organism_reference", "Assembly_accession"]
    )
    .copy()
)

print(
    lista_multiples[
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank",
            "Match_status"
        ]
    ].to_string(index=False)
)
#%% 
#Selección atomática del ASSEMBLY
# Prioridad del nivel de ensamblaje
prioridad_status = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1
}

seleccion = tabla_candidatos.copy()

# 1. PRIORIDAD DEL NIVEL DE ASSEMBLY

seleccion["Prioridad_status"] = (
    seleccion["Assembly_status"]
    .map(prioridad_status)
    .fillna(0)
)

# 2. PRIORIDAD REFSEQ
#    GCF = RefSeq
#    GCA = GenBank

seleccion["Es_RefSeq"] = (
    seleccion["Assembly_accession"]
    .astype(str)
    .str.startswith("GCF_")
    .astype(int)
)

# 3. VERSIÓN DEL ASSEMBLY

seleccion["Version_Assembly"] = (
    seleccion["Assembly_accession"]
    .astype(str)
    .str.extract(r"\.(\d+)$")[0]
    .fillna("0")
    .astype(int)
)

# 4. SELECCIÓN DEL MEJOR ASSEMBLY

seleccion = seleccion.sort_values(
    [
        "Organism_reference",
        "Prioridad_status",
        "Es_RefSeq",
        "Version_Assembly"
    ],
    ascending=[
        True,
        False,
        False,
        False
    ]
)

tabla_assembly_final = (
    seleccion
    .drop_duplicates(
        subset=["Organism_reference"],
        keep="first"
    )
    .copy()
)

print("Organismos totales:", len(tabla_assembly_final))

print(
    tabla_assembly_final[
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "RefSeq",
            "GenBank",
            "label",
            "Match_status"
        ]
    ].to_string(index=False)
)
#%% 
#Clasificar selección

tabla_assembly_final["Tipo_seleccion"] = "Revisar"

#Negativos
tabla_assembly_final.loc[
    tabla_assembly_final["label"] == 0,
    "Tipo_seleccion"
] = "Negativo - representante"

#Positivos
tabla_assembly_final.loc[
    tabla_assembly_final["label"] == 1,
    "Tipo_seleccion"
] = "Positivo - Assembly seleccionado"

print(
    tabla_assembly_final[
        [
            "Organism_reference",
            "Assembly_accession",
            "Assembly_status",
            "label",
            "Tipo_seleccion"
        ]
    ].to_string(index=False)
)

print("\nResumen:")
print(
    tabla_assembly_final["Tipo_seleccion"]
    .value_counts()
)
#%% 
#Resumen ASSEMBLY

print(
    tabla_assembly_final[
        "Assembly_status"
    ].value_counts()
)

print("\nRefSeq vs GenBank:")

print(
    tabla_assembly_final[
        "Assembly_accession"
    ]
    .astype(str)
    .str.startswith("GCF_")
    .map({
        True: "RefSeq",
        False: "GenBank"
    })
    .value_counts()
)

print("\nLabels:")

print(
    tabla_assembly_final["label"]
    .value_counts(dropna=False)
)
#%% 
#QC ASSEMBLY final

#Resumen general
print("TOTAL ORGANISMOS:", len(tabla_assembly_final))
print("\nLABEL:")
print(tabla_assembly_final["label"].value_counts())

print("\nASSEMBLY STATUS:")
print(tabla_assembly_final["Assembly_status"].value_counts())

print("\nREFSEQ / GENBANK:")
print(
    tabla_assembly_final["Assembly_accession"]
    .astype(str)
    .str.startswith("GCF_")
    .map({
        True: "RefSeq",
        False: "GenBank"
    })
    .value_counts()
)

#Positivos y negativos por calidad
print("\nSTATUS POR LABEL:")
print(
    pd.crosstab(
        tabla_assembly_final["label"],
        tabla_assembly_final["Assembly_status"]
    )
)

#Organismos que quedaron en Contig
print("\nORGANISMOS EN CONTIG:")
print(
    tabla_assembly_final[
        tabla_assembly_final["Assembly_status"] == "Contig"
    ][
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Taxid",
            "label",
            "Match_status"
        ]
    ].to_string(index=False)
)

#Organismos que quedaron en Scaffold
print("\nORGANISMOS EN SCAFFOLD:")
print(
    tabla_assembly_final[
        tabla_assembly_final["Assembly_status"] == "Scaffold"
    ][
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Taxid",
            "label",
            "Match_status"
        ]
    ].to_string(index=False)
)
#%% MARCAR ASSEMBLY PARA REVISION

tabla_assembly_final["Revision_manual"] = "NO"

organismos_revision = [
    "Cycloclasticus pugetii PS‐1",
    "Pseudomonas putida CSV86",
    "Burkholderia sp. Ch1‐1",
    "Burkholderia sp. RP007",
    "Pseudomonas saccharophilia"
]

tabla_assembly_final.loc[
    tabla_assembly_final["Organism_reference"].isin(
        organismos_revision
    ),
    "Revision_manual"
] = "SI"

print(
    tabla_assembly_final[
        tabla_assembly_final["Revision_manual"] == "SI"
    ][
        [
            "Organism_reference",
            "Assembly_accession",
            "Organism_NCBI",
            "Species_NCBI",
            "Taxid",
            "Assembly_status",
            "label",
            "Match_status"
        ]
    ].to_string(index=False)
)
# %% EXCLUIR ORGANISMO AMBIGUO

tabla_assembly_final = tabla_assembly_final[
    tabla_assembly_final["Organism_reference"] !=
    "Pseudomonas putida CSV86"
].copy()

print("Total organismos:", len(tabla_assembly_final))
print("\nLabels:")
print(tabla_assembly_final["label"].value_counts())
# %% TABLA DEFINITIVA DE ORGANISMOS

dataset_referencia = tabla_assembly_final[
    [
        "Organism_reference",
        "label",
        "Reference",
        "Assembly_ID",
        "Assembly_accession",
        "Organism_NCBI",
        "Species_NCBI",
        "Taxid",
        "Assembly_status",
        "RefSeq",
        "GenBank",
        "Contig_N50",
        "Scaffold_N50",
        "Ftp_RefSeq",
        "Ftp_GenBank"
    ]
].copy()

dataset_referencia = dataset_referencia.rename(
    columns={
        "Organism_reference": "Organism"
    }
)

print(dataset_referencia.shape)
print(dataset_referencia.head())
# %% GUARDAR DATASET MAESTRO

ruta_salida = "C:/Users/USUARIO/Desktop/TFM VARIOS/Dataset_referencia_144_organismos.xlsx"

dataset_referencia.to_excel(
    ruta_salida,
    index=False
)

print("Archivo guardado correctamente")
print(ruta_salida)
# %% VERIFICACIÓN FINAL

print("Total:", len(dataset_referencia))

print("\nPositivos y negativos:")
print(dataset_referencia["label"].value_counts())

print("\nEstados de ensamblaje:")
print(dataset_referencia["Assembly_status"].value_counts())

print("\nOrganismos únicos:")
print(dataset_referencia["Organism"].nunique())

print("\nAssembly únicos:")
print(dataset_referencia["Assembly_accession"].nunique())
# %% DETECTAR ASSEMBLY DUPLICADOS

assembly_duplicados = (
    dataset_referencia
    .groupby("Assembly_accession")
    .agg(
        N_organismos=("Organism", "nunique"),
        Organismos=("Organism", lambda x: " | ".join(x)),
        Labels=("label", lambda x: " | ".join(x.astype(str)))
    )
    .reset_index()
)

assembly_duplicados = assembly_duplicados[
    assembly_duplicados["N_organismos"] > 1
].copy()

print(
    assembly_duplicados.to_string(index=False)
)

print(
    "\nAssembly compartidos:",
    len(assembly_duplicados)
)
# %% COMPROBAR ASSEMBLY VACÍOS

print("Assembly vacíos (NaN):",
      dataset_referencia["Assembly_accession"].isna().sum())

print("Assembly vacíos (''):",
      dataset_referencia["Assembly_accession"].eq("").sum())

print("\nOrganismos sin Assembly:")

sin_assembly = dataset_referencia[
    dataset_referencia["Assembly_accession"].isna() |
    dataset_referencia["Assembly_accession"].eq("")
][
    [
        "Organism",
        "label",
        "Assembly_accession",
        "Organism_NCBI",
        "Taxid"
    ]
]

print(sin_assembly.to_string(index=False))
# %% BUSCAR ASSEMBLIES DE LOS 71 ORGANISMOS FALTANTES

def buscar_assemblies(organismo, retmax=20):
    try:
        handle = Entrez.esearch(
            db="assembly",
            term=f'"{organismo}"[Organism]',
            retmax=retmax
        )
        resultado = Entrez.read(handle)
        handle.close()
        return resultado["IdList"]

    except Exception as e:
        print(f"Error buscando {organismo}: {e}")
        return []


def obtener_datos_assembly(assembly_id):

    try:
        handle = Entrez.esummary(
            db="assembly",
            id=assembly_id,
            report="full"
        )

        resultado = Entrez.read(handle)
        handle.close()

        return resultado["DocumentSummarySet"]["DocumentSummary"][0]

    except Exception as e:
        print(f"Error con Assembly {assembly_id}: {e}")
        return None


def extraer_datos_assembly(assembly_id):

    datos = obtener_datos_assembly(assembly_id)

    if datos is None:
        return None

    return {
        "Assembly_ID": assembly_id,
        "Assembly_accession": datos.get("AssemblyAccession", ""),
        "Assembly_name": datos.get("AssemblyName", ""),
        "Organism_NCBI": datos.get("Organism", ""),
        "Species_NCBI": datos.get("SpeciesName", ""),
        "Taxid": datos.get("Taxid", ""),
        "Assembly_status": datos.get("AssemblyStatus", ""),
        "RefSeq": datos.get("Synonym", {}).get("RefSeq", ""),
        "GenBank": datos.get("Synonym", {}).get("Genbank", ""),
        "Contig_N50": datos.get("ContigN50", ""),
        "Scaffold_N50": datos.get("ScaffoldN50", ""),
        "Ftp_RefSeq": datos.get("FtpPath_RefSeq", ""),
        "Ftp_GenBank": datos.get("FtpPath_GenBank", "")
    }


# Organismos que todavía no tienen Assembly
organismos_faltantes = (
    dataset_referencia[
        dataset_referencia["Assembly_accession"].isna()
    ]["Organism"]
    .dropna()
    .unique()
    .tolist()
)

print("Organismos a buscar:", len(organismos_faltantes))


resultados_faltantes = []

for i, organismo in enumerate(organismos_faltantes, start=1):

    print(f"[{i}/{len(organismos_faltantes)}] {organismo}")

    ids = buscar_assemblies(organismo, retmax=20)

    if not ids:

        resultados_faltantes.append({
            "Organism_original": organismo,
            "Assembly_ID": "",
            "Assembly_accession": "",
            "Assembly_name": "",
            "Organism_NCBI": "",
            "Species_NCBI": "",
            "Taxid": "",
            "Assembly_status": "",
            "RefSeq": "",
            "GenBank": "",
            "Contig_N50": "",
            "Scaffold_N50": "",
            "Ftp_RefSeq": "",
            "Ftp_GenBank": ""
        })

    else:

        for assembly_id in ids:

            datos = extraer_datos_assembly(assembly_id)

            if datos is not None:

                datos["Organism_original"] = organismo

                resultados_faltantes.append(datos)

    time.sleep(0.34)


tabla_faltantes = pd.DataFrame(resultados_faltantes)


print("\nRESULTADOS")
print("Organismos buscados:", len(organismos_faltantes))
print(
    "Organismos con algún resultado:",
    tabla_faltantes.loc[
        tabla_faltantes["Assembly_accession"] != "",
        "Organism_original"
    ].nunique()
)
print(
    "Organismos sin ningún resultado:",
    tabla_faltantes.loc[
        tabla_faltantes["Assembly_accession"] == "",
        "Organism_original"
    ].nunique()
)

print("\nTABLA DE RESULTADOS:")
print(tabla_faltantes.to_string(index=False))
# %% CONTAR POSITIVOS Y NEGATIVOS CON ASSEMBLY

con_assembly = dataset_referencia[
    dataset_referencia["Assembly_accession"].notna() &
    dataset_referencia["Assembly_accession"].ne("")
].copy()

print("ORGANISMOS CON ASSEMBLY")
print("-----------------------")
print("Total:", len(con_assembly))

print("\nPor etiqueta:")
print(con_assembly["label"].value_counts().sort_index())

print("\nDetalle:")
print(
    con_assembly
    .groupby("label")["Organism"]
    .count()
    .rename({0: "Negativos", 1: "Positivos"})
)

# Lista de organismos
print("\nPOSITIVOS CON ASSEMBLY:")
print(
    con_assembly[
        con_assembly["label"] == 1
    ]["Organism"].to_string(index=False)
)

print("\nNEGATIVOS CON ASSEMBLY:")
print(
    con_assembly[
        con_assembly["label"] == 0
    ]["Organism"].to_string(index=False)
)
# %% DATASET CANDIDATOS PARA ENTRENAMIENTO

dataset_candidatos = dataset_referencia[
    dataset_referencia["Assembly_accession"].notna() &
    dataset_referencia["Assembly_accession"].ne("")
].copy()

print("Total candidatos:", len(dataset_candidatos))

print("\nPor etiqueta:")
print(dataset_candidatos["label"].value_counts())

print("\nAssembly únicos:",
      dataset_candidatos["Assembly_accession"].nunique())
# %% ELIMINAR DUPLICADO PYR-1

dataset_candidatos = dataset_candidatos[
    dataset_candidatos["Organism"] != "Mycobacterium sp. PYR-1"
].copy()

print("Total después de eliminar duplicado:", len(dataset_candidatos))

print("\nPor etiqueta:")
print(dataset_candidatos["label"].value_counts())

print("\nAssembly únicos:",
      dataset_candidatos["Assembly_accession"].nunique())
# %% REVISAR NEGATIVOS

negativos = dataset_candidatos[
    dataset_candidatos["label"] == 0
].copy()

print(
    negativos[
        [
            "Organism",
            "Assembly_accession",
            "Organism_NCBI",
            "Taxid",
            "Assembly_status"
        ]
    ].to_string(index=False)
)
# %% REVISAR POSITIVOS

positivos = dataset_candidatos[
    dataset_candidatos["label"] == 1
].copy()

print(
    positivos[
        [
            "Organism",
            "Assembly_accession",
            "Organism_NCBI",
            "Taxid",
            "Assembly_status"
        ]
    ].to_string(index=False)
)
# %% CREAR DATASET DE ENTRENAMIENTO 29 + 29

# Copia de los organismos que tienen Assembly
dataset_entrenamiento = dataset_candidatos.copy()

# --------------------------------------------------
# 1. Eliminar duplicado de PYR-1
# --------------------------------------------------

dataset_entrenamiento = dataset_entrenamiento[
    dataset_entrenamiento["Organism"] != "Mycobacterium sp. PYR-1"
].copy()


# --------------------------------------------------
# 2. Separar positivos y negativos
# --------------------------------------------------

positivos_train = dataset_entrenamiento[
    dataset_entrenamiento["label"] == 1
].copy()

negativos_disponibles = dataset_entrenamiento[
    dataset_entrenamiento["label"] == 0
].copy()


# --------------------------------------------------
# 3. Seleccionar 29 negativos
# --------------------------------------------------

negativos_train = (
    negativos_disponibles
    .sample(
        n=29,
        random_state=42
    )
    .copy()
)


# --------------------------------------------------
# 4. Unir positivos + negativos
# --------------------------------------------------

dataset_train = pd.concat(
    [
        positivos_train,
        negativos_train
    ],
    ignore_index=True
)


# --------------------------------------------------
# 5. Ordenar
# --------------------------------------------------

dataset_train = dataset_train.sort_values(
    by=["label", "Organism"]
).reset_index(drop=True)


# --------------------------------------------------
# 6. Comprobar
# --------------------------------------------------

print("================================")
print("DATASET DE ENTRENAMIENTO")
print("================================")

print("\nTotal organismos:",
      len(dataset_train))

print("\nDistribución:")
print(
    dataset_train["label"]
    .value_counts()
    .sort_index()
)

print("\nAssembly únicos:",
      dataset_train["Assembly_accession"].nunique())

print("\nOrganismos únicos:",
      dataset_train["Organism"].nunique())

print("\nPositivos:")
print(
    dataset_train[
        dataset_train["label"] == 1
    ]["Organism"].to_string(index=False)
)

print("\nNegativos:")
print(
    dataset_train[
        dataset_train["label"] == 0
    ]["Organism"].to_string(index=False)
)
# %% IDENTIFICAR DUPLICADO EN DATASET TRAIN

duplicados_train = (
    dataset_train
    .groupby("Assembly_accession")
    .agg(
        N_organismos=("Organism", "nunique"),
        Organismos=("Organism", lambda x: " | ".join(x)),
        Labels=("label", lambda x: " | ".join(x.astype(str)))
    )
    .reset_index()
)

print(
    duplicados_train[
        duplicados_train["N_organismos"] > 1
    ].to_string(index=False)
)
# %% ELIMINAR ASSEMBLY DUPLICADO

# Identificar el Assembly que está repetido
conteo = dataset_train["Assembly_accession"].value_counts()

assemblies_duplicados = conteo[
    conteo > 1
].index.tolist()

print("Assembly duplicados:")
print(assemblies_duplicados)


# Mostrar las filas duplicadas
print("\nFilas que vamos a revisar:")

print(
    dataset_train[
        dataset_train["Assembly_accession"].isin(assemblies_duplicados)
    ][
        [
            "Organism",
            "label",
            "Assembly_accession"
        ]
    ].to_string(index=False)
)
      # %% ELIMINAR UNA COPIA DEL ASSEMBLY DUPLICADO

for assembly in assemblies_duplicados:

    indices = dataset_train[
        dataset_train["Assembly_accession"] == assembly
    ].index.tolist()

    # Conservamos la primera aparición
    # y eliminamos las siguientes
    if len(indices) > 1:
        dataset_train = dataset_train.drop(
            indices[1:]
        )


dataset_train = dataset_train.reset_index(drop=True)


# Comprobación final
print("================================")
print("DATASET TRAIN CORREGIDO")
print("================================")

print("Total organismos:",
      len(dataset_train))

print("\nDistribución:")
print(
    dataset_train["label"]
    .value_counts()
    .sort_index()
)

print("\nAssembly únicos:",
      dataset_train["Assembly_accession"].nunique())

print("\nOrganismos únicos:",
      dataset_train["Organism"].nunique())
# %% CARGAR MATRIZ X DE 107 KO

ruta_X = r"C:\Users\USUARIO\Desktop\TFM VARIOS\X_ML_107_KO.csv"

X_107 = pd.read_csv(ruta_X)

print("Dimensiones de X_107:")
print(X_107.shape)

print("\nColumnas:")
print(X_107.columns.tolist())

print("\nPrimeras filas:")
print(X_107.head().to_string())
# %% REVISAR IDENTIFICADORES DE X_107

print("\nTipos de datos:")
print(X_107.dtypes.head(15))

print("\nPrimera columna:")
print(X_107.iloc[:, 0].name)

print("\nValores de la primera columna:")
print(
    X_107.iloc[:, 0]
    .head(20)
    .to_string(index=False)
)
# %% GUARDAR DATASET DE ENTRENAMIENTO

ruta_train = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\dataset_entrenamiento_58.csv"
)

dataset_train.to_csv(
    ruta_train,
    index=False
)

print("Dataset guardado en:")
print(ruta_train)

print("\nDimensiones:", dataset_train.shape)
print(
    dataset_train["label"]
    .value_counts()
    .sort_index()
)

# %% DESCARGAR LOS 58 GENOMAS DE ENTRENAMIENTO

import os
import subprocess
import pandas as pd

# Carpeta principal
carpeta_base = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\genomas_entrenamiento"
)

os.makedirs(carpeta_base, exist_ok=True)

# Assembly únicos del dataset de entrenamiento
assemblies = (
    dataset_train[["Organism", "Assembly_accession", "label"]]
    .dropna(subset=["Assembly_accession"])
    .drop_duplicates(subset=["Assembly_accession"])
    .copy()
)

print("Número de Assembly a descargar:", len(assemblies))
print()

resultados_descarga = []

for i, fila in enumerate(assemblies.itertuples(index=False), start=1):

    organismo = fila.Organism
    assembly = fila.Assembly_accession
    label = fila.label

    print("=" * 70)
    print(f"[{i}/{len(assemblies)}] {organismo}")
    print(f"Assembly: {assembly}")
    print(f"Label: {label}")

    # Carpeta individual
    carpeta_assembly = os.path.join(
        carpeta_base,
        assembly
    )

    os.makedirs(carpeta_assembly, exist_ok=True)

    # Comprobar si ya existe el protein.faa
    ruta_proteinas = os.path.join(
        carpeta_assembly,
        "ncbi_dataset",
        "data",
        assembly,
        "protein.faa"
    )

    if os.path.exists(ruta_proteinas):
        print("✓ Ya descargado. Se omite.")

        resultados_descarga.append({
            "Organism": organismo,
            "Assembly_accession": assembly,
            "label": label,
            "Estado": "Ya descargado"
        })

        continue

    # Archivo ZIP temporal
    zip_path = os.path.join(
        carpeta_assembly,
        f"{assembly}.zip"
    )

    comando = [
        "datasets",
        "download",
        "genome",
        "accession",
        assembly,
        "--include",
        "genome,protein,gff3",
        "--filename",
        zip_path
    ]

    try:

        resultado = subprocess.run(
            comando,
            capture_output=True,
            text=True
        )

        if resultado.returncode != 0:

            print("✗ ERROR en la descarga")
            print(resultado.stderr)

            resultados_descarga.append({
                "Organism": organismo,
                "Assembly_accession": assembly,
                "label": label,
                "Estado": "Error descarga"
            })

            continue

        print("✓ Descarga completada")

        # Extraer ZIP
        resultado_extract = subprocess.run(
            [
                "tar",
                "-xf",
                zip_path,
                "-C",
                carpeta_assembly
            ],
            capture_output=True,
            text=True
        )

        if resultado_extract.returncode != 0:

            print("✗ Error extrayendo ZIP")
            print(resultado_extract.stderr)

            resultados_descarga.append({
                "Organism": organismo,
                "Assembly_accession": assembly,
                "label": label,
                "Estado": "Error extracción"
            })

            continue

        print("✓ Archivos extraídos")

        # Eliminar ZIP para ahorrar espacio
        if os.path.exists(zip_path):
            os.remove(zip_path)

        # Verificar proteína
        if os.path.exists(ruta_proteinas):

            print("✓ protein.faa encontrado")

            resultados_descarga.append({
                "Organism": organismo,
                "Assembly_accession": assembly,
                "label": label,
                "Estado": "OK"
            })

        else:

            print("⚠ protein.faa NO encontrado")

            resultados_descarga.append({
                "Organism": organismo,
                "Assembly_accession": assembly,
                "label": label,
                "Estado": "Sin protein.faa"
            })

    except Exception as e:

        print("✗ Error inesperado:", e)

        resultados_descarga.append({
            "Organism": organismo,
            "Assembly_accession": assembly,
            "label": label,
            "Estado": "Error inesperado"
        })


# Convertir resultados a DataFrame
resultados_descarga = pd.DataFrame(resultados_descarga)

print("\n" + "=" * 70)
print("RESUMEN DE DESCARGAS")
print("=" * 70)

print(
    resultados_descarga["Estado"]
    .value_counts()
)

print("\nTotal procesados:", len(resultados_descarga))
# %%
print(tabla_assembly_final.shape)
print(tabla_assembly_final.columns.tolist())
print(tabla_assembly_final[["Organism", "Assembly_accession"]].head(10))
# %%
print(dataset_train[[
    "Organism",
    "Assembly_accession",
    "label"
]].to_string(index=False))
# %%
print(resultados_descarga["Estado"].value_counts(dropna=False))
print(
    resultados_descarga[
        resultados_descarga["Estado"] != "OK"
    ][["Organism", "Assembly_accession", "label", "Estado"]]
    .to_string(index=False)
)
# %%
from pathlib import Path

carpeta = Path(
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\genomas_entrenamiento"
)

protein_files = list(carpeta.rglob("protein.faa"))

print("Protein.faa encontrados:", len(protein_files))

for p in protein_files:
    print(p)

# %%
from pathlib import Path

protein_files = list(carpeta.rglob("protein.faa"))

accesiones = []

for p in protein_files:
    coincidencias = re.findall(r"GCA_\d+\.\d+|GCF_\d+\.\d+", str(p))
    if coincidencias:
        accesiones.append(coincidencias[-1])

print("Protein.faa:", len(protein_files))
print("Accesiones totales:", len(accesiones))
print("Accesiones únicas:", len(set(accesiones)))

print("\nDuplicados:")
from collections import Counter

for acc, n in Counter(accesiones).items():
    if n > 1:
        print(acc, "→", n)
        
# %%
protein_files = list(carpeta.rglob("protein.faa"))

registros_proteinas = []

for p in protein_files:
    
    # Buscar accession en la ruta
    coincidencias = re.findall(r"GCA_\d+\.\d+|GCF_\d+\.\d+", str(p))
    
    if not coincidencias:
        continue
    
    accession = coincidencias[-1]
    
    registros_proteinas.append({
        "Assembly_accession": accession,
        "protein_faa": str(p)
    })

tabla_proteinas = pd.DataFrame(registros_proteinas)

# Eliminar el duplicado de Alcanivorax
tabla_proteinas = tabla_proteinas.drop_duplicates(
    subset="Assembly_accession"
).reset_index(drop=True)

# Incorporar organismo y label
tabla_proteinas = tabla_proteinas.merge(
    dataset_train[
        ["Organism", "Assembly_accession", "label"]
    ],
    on="Assembly_accession",
    how="left"
)

# Ordenar columnas
tabla_proteinas = tabla_proteinas[
    ["Organism", "Assembly_accession", "label", "protein_faa"]
]

print("Genomas únicos:", len(tabla_proteinas))
print("\nDistribución de labels:")
print(tabla_proteinas["label"].value_counts(dropna=False))

print("\nPrimeros registros:")
print(tabla_proteinas.head(10).to_string(index=False))
# %%
print("Organismos sin correspondencia:",
      tabla_proteinas["Organism"].isna().sum())

print("Labels sin correspondencia:",
      tabla_proteinas["label"].isna().sum())

print("Accesiones únicas:",
      tabla_proteinas["Assembly_accession"].nunique())

# %%
from pathlib import Path
import shutil

carpeta_kofam = carpeta / "proteinas_kofam"
carpeta_kofam.mkdir(exist_ok=True)

for _, fila in tabla_proteinas.iterrows():

    accession = fila["Assembly_accession"]
    origen = Path(fila["protein_faa"])

    destino = carpeta_kofam / f"{accession}.faa"

    if not destino.exists():
        shutil.copy2(origen, destino)

print("Archivos preparados:",
      len(list(carpeta_kofam.glob("*.faa"))))

print("Carpeta:")
print(carpeta_kofam)

# %%
from pathlib import Path
import pandas as pd
import shutil
import re

# ============================================================
# CONFIGURACIÓN
# ============================================================

carpeta_base = Path(
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\genomas_entrenamiento"
)

carpeta_kofam = carpeta_base / "proteinas_kofam"
carpeta_kofam.mkdir(exist_ok=True)

# ============================================================
# LOCALIZAR TODOS LOS protein.faa
# ============================================================

protein_files = list(carpeta_base.rglob("protein.faa"))

print("Protein.faa encontrados:", len(protein_files))

# ============================================================
# CONSTRUIR TABLA DE ARCHIVOS
# ============================================================

registros = []

for p in protein_files:

    coincidencias = re.findall(
        r"GCA_\d+\.\d+|GCF_\d+\.\d+",
        str(p)
    )

    if not coincidencias:
        continue

    accession = coincidencias[-1]

    registros.append({
        "Assembly_accession": accession,
        "protein_original": str(p)
    })

tabla_archivos = pd.DataFrame(registros)

# ============================================================
# ELIMINAR DUPLICADOS
# ============================================================

tabla_archivos = (
    tabla_archivos
    .drop_duplicates(subset="Assembly_accession")
    .reset_index(drop=True)
)

# ============================================================
# UNIR CON DATASET DE ENTRENAMIENTO
# ============================================================

tabla_archivos = tabla_archivos.merge(
    dataset_train[
        ["Organism", "Assembly_accession", "label"]
    ],
    on="Assembly_accession",
    how="left"
)

# ============================================================
# COPIAR A CARPETA KOFAM
# ============================================================

rutas_kofam = []

for _, fila in tabla_archivos.iterrows():

    accession = fila["Assembly_accession"]
    origen = Path(fila["protein_original"])

    destino = carpeta_kofam / f"{accession}.faa"

    if not destino.exists():
        shutil.copy2(origen, destino)

    rutas_kofam.append(str(destino))

tabla_archivos["protein_kofam"] = rutas_kofam

# ============================================================
# COMPROBACIONES
# ============================================================

print("\n========================================")
print("RESUMEN")
print("========================================")

print("Genomas únicos:", len(tabla_archivos))
print(
    "Positivos:",
    (tabla_archivos["label"] == 1).sum()
)
print(
    "Negativos:",
    (tabla_archivos["label"] == 0).sum()
)

print(
    "Sin organismo:",
    tabla_archivos["Organism"].isna().sum()
)

print(
    "Sin label:",
    tabla_archivos["label"].isna().sum()
)

print(
    "Archivos Kofam:",
    len(list(carpeta_kofam.glob("*.faa")))
)

print("\nPrimeros registros:")
print(
    tabla_archivos[
        ["Organism", "Assembly_accession", "label", "protein_kofam"]
    ].head(10).to_string(index=False)
)

# %%
# ============================================================
# LISTA MAESTRA PARA KOFAMKOALA
# ============================================================

tabla_kofam = tabla_archivos[
    ["Organism", "Assembly_accession", "label", "protein_kofam"]
].copy()

tabla_kofam = tabla_kofam.sort_values(
    "Assembly_accession"
).reset_index(drop=True)

print("TOTAL:", len(tabla_kofam))
print("\nDistribución:")
print(tabla_kofam["label"].value_counts())

print("\nLISTA COMPLETA:")
print(
    tabla_kofam.to_string(index=False)
)
# %%
archivo_lista_kofam = (
    carpeta_base / "lista_53_genomas_KofamKOALA.csv"
)

tabla_kofam.to_csv(
    archivo_lista_kofam,
    index=False,
    encoding="utf-8-sig"
)

print("\nGuardado en:")
print(archivo_lista_kofam)

ruta = r"C:\Users\USUARIO\Desktop\TFM VARIOS\X_ML_107_KO.csv"

X = pd.read_csv(ruta)

print(X.shape)
print(X.columns.tolist())

ruta_ko = r"C:\Users\USUARIO\Desktop\TFM VARIOS\genomas_entrenamiento\KO_objetivo_107.txt"

with open(ruta_ko, "w") as f:
    for ko in X.columns:
        f.write(ko + "\n")

print("Archivo creado:")
print(ruta_ko)

print("Número de KO:", len(X.columns))
# %%
import numpy as np

ruta = r"C:\Users\USUARIO\Desktop\TFM VARIOS\MATRIZ_KOFAMSCAN_53x78.csv"

df = pd.read_csv(ruta)

print("Dimensiones:", df.shape)

print("\nDistribución de etiquetas:")
print(df["Label"].value_counts())

# Separar variable objetivo
y = df["Label"].astype(int)

# Metadatos que NO entran al modelo
columnas_metadata = [
    "Organism",
    "Genome_ID",
    "Label"
]

# Matriz de variables predictoras
X = df.drop(columns=columnas_metadata)

print("\nDimensiones de X:", X.shape)
print("Dimensiones de y:", y.shape)

print("\nTipos de datos de X:")
print(X.dtypes.value_counts())

print("\nValores únicos de X:")
print(np.unique(X.values))
# %%
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import make_scorer, precision_score, recall_score, f1_score

# --------------------------------------------------
# VALIDACIÓN CRUZADA
# --------------------------------------------------

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

# --------------------------------------------------
# MODELO RANDOM FOREST
# --------------------------------------------------

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

# --------------------------------------------------
# MÉTRICAS
# --------------------------------------------------

scoring = {
    "accuracy": "accuracy",
    "precision": make_scorer(precision_score, zero_division=0),
    "recall": make_scorer(recall_score, zero_division=0),
    "f1": make_scorer(f1_score, zero_division=0),
    "roc_auc": "roc_auc"
}

# --------------------------------------------------
# VALIDACIÓN
# --------------------------------------------------

resultados_rf = cross_validate(
    rf,
    X,
    y,
    cv=cv,
    scoring=scoring,
    return_train_score=False,
    n_jobs=-1
)

# --------------------------------------------------
# RESULTADOS POR FOLD
# --------------------------------------------------

print("RESULTADOS RANDOM FOREST")
print("-" * 40)

for metrica in scoring:
    valores = resultados_rf[f"test_{metrica}"]
    
    print(f"\n{metrica.upper()}:")
    print("Por fold:", np.round(valores, 3))
    print("Media:", round(valores.mean(), 3))
    print("Desviación estándar:", round(valores.std(), 3))

# %%
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import (
    make_scorer,
    precision_score,
    recall_score,
    f1_score
)

# --------------------------------------------------
# VALIDACIÓN CRUZADA
# --------------------------------------------------

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

# --------------------------------------------------
# MODELO SVM
# --------------------------------------------------

svm = Pipeline([
    ("scaler", StandardScaler()),
    ("classifier", SVC(
        kernel="rbf",
        probability=True,
        class_weight="balanced",
        random_state=42
    ))
])

# --------------------------------------------------
# MÉTRICAS
# --------------------------------------------------

scoring = {
    "accuracy": "accuracy",
    "precision": make_scorer(
        precision_score,
        zero_division=0
    ),
    "recall": make_scorer(
        recall_score,
        zero_division=0
    ),
    "f1": make_scorer(
        f1_score,
        zero_division=0
    ),
    "roc_auc": "roc_auc"
}

# --------------------------------------------------
# VALIDACIÓN
# --------------------------------------------------

resultados_svm = cross_validate(
    svm,
    X,
    y,
    cv=cv,
    scoring=scoring,
    return_train_score=False,
    n_jobs=-1
)

# --------------------------------------------------
# RESULTADOS
# --------------------------------------------------

print("RESULTADOS SVM")
print("-" * 40)

for metrica in scoring:

    valores = resultados_svm[f"test_{metrica}"]

    print(f"\n{metrica.upper()}:")
    print("Por fold:", np.round(valores, 3))
    print("Media:", round(valores.mean(), 3))
    print("Desviación estándar:", round(valores.std(), 3))
    
# %%
# ============================================================
# RANDOM FOREST FINAL - IMPORTANCIA DE VARIABLES
# ============================================================

rf_final = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

# Entrenamiento con los 53 organismos
rf_final.fit(X, y)

# Importancia de cada KO
importancias_rf = pd.DataFrame({
    "KO": X.columns,
    "Importance": rf_final.feature_importances_
})

# Ordenar de mayor a menor
importancias_rf = importancias_rf.sort_values(
    by="Importance",
    ascending=False
).reset_index(drop=True)

print("TOP 20 KO MÁS IMPORTANTES")
print("-" * 40)
print(importancias_rf.head(20).to_string(index=False))

# %%
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

# Predicciones out-of-fold
pred_rf = cross_val_predict(
    rf,
    X,
    y,
    cv=cv,
    method="predict",
    n_jobs=-1
)

# Matriz de confusión
cm_rf = confusion_matrix(y, pred_rf)

print("MATRIZ DE CONFUSIÓN - RANDOM FOREST")
print("-" * 40)
print(cm_rf)

print("\nInterpretación:")
print("Verdaderos negativos (TN):", cm_rf[0, 0])
print("Falsos positivos (FP):", cm_rf[0, 1])
print("Falsos negativos (FN):", cm_rf[1, 0])
print("Verdaderos positivos (TP):", cm_rf[1, 1])

# Visualización
disp = ConfusionMatrixDisplay(
    confusion_matrix=cm_rf,
    display_labels=["Negativo", "Positivo"]
)

disp.plot()
plt.title("Matriz de confusión - Random Forest")
plt.show()

# %%
# ============================================================
# MODELO FINAL RANDOM FOREST
# ============================================================

rf_final = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

# Entrenar con todos los organismos de referencia
rf_final.fit(X, y)

print("Modelo Random Forest entrenado correctamente.")
print("Número de organismos:", X.shape[0])
print("Número de variables KO:", X.shape[1])

# %%
from pathlib import Path

carpeta = Path(r"C:\Users\USUARIO\Desktop\TFM VARIOS")

print("Archivos relacionados con matrices/candidatos:\n")

for archivo in sorted(carpeta.rglob("*")):
    if archivo.is_file():
        nombre = archivo.name.lower()
        
        if any(palabra in nombre for palabra in [
            "candidato",
            "matrix",
            "matriz",
            "koala",
            "kofam",
            "genoma"
        ]):
            print(archivo)
            
# %%
import pandas as pd
from pathlib import Path

archivos = [
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv",
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\MATRIZ_REFERENCIAS_KO.csv",
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\organismos_referencia_candidatos.csv",
    r"C:\Users\USUARIO\Desktop\TFM VARIOS\matriz_resumen_funcional_141.csv",
]

for ruta_archivo in archivos:
    
    print("\n" + "=" * 70)
    print(Path(ruta_archivo).name)
    print("=" * 70)
    
    tabla = pd.read_csv(ruta_archivo)
    
    print("Dimensiones:", tabla.shape)
    print("Primeras columnas:")
    print(tabla.columns[:10].tolist())
    print("Últimas columnas:")
    print(tabla.columns[-10:].tolist())
    
    print("\nPrimeras 3 filas:")
    print(tabla.head(3))
    
# %%
# ============================================================
# PREDICCIÓN SOBRE LOS 141 CANDIDATOS
# ============================================================
# ------------------------------------------------------------
# 1. Cargar matriz de candidatos
# ------------------------------------------------------------

ruta_candidatos = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv"
)

candidatos = pd.read_csv(ruta_candidatos)

print("Matriz de candidatos:", candidatos.shape)

# ------------------------------------------------------------
# 2. Obtener los 78 KO utilizados por el modelo
# ------------------------------------------------------------

KO_modelo = X.columns.tolist()

print("KO utilizados por el modelo:", len(KO_modelo))

# ------------------------------------------------------------
# 3. Comprobar que los 78 KO existen en candidatos
# ------------------------------------------------------------

KO_faltantes = [
    ko for ko in KO_modelo
    if ko not in candidatos.columns
]

print("\nKO faltantes en candidatos:", len(KO_faltantes))

if len(KO_faltantes) > 0:
    print(KO_faltantes)
    raise ValueError(
        "Faltan KO necesarios para realizar la predicción."
    )

# ------------------------------------------------------------
# 4. Construir X de candidatos
# ------------------------------------------------------------

X_candidatos = candidatos[KO_modelo].copy()

print("\nDimensiones X_candidatos:", X_candidatos.shape)

# ------------------------------------------------------------
# 5. Comprobar valores
# ------------------------------------------------------------

print("\nValores únicos:")
print(np.unique(X_candidatos.values))

# ------------------------------------------------------------
# 6. Predicción
# ------------------------------------------------------------

predicciones = rf_final.predict(X_candidatos)

# Probabilidad de pertenecer a clase positiva
probabilidades = rf_final.predict_proba(X_candidatos)[:, 1]

# ------------------------------------------------------------
# 7. Crear tabla de resultados
# ------------------------------------------------------------

resultados_candidatos = candidatos[
    ["Organism_code"]
].copy()

resultados_candidatos["Prediction"] = predicciones
resultados_candidatos["Probability_positive"] = probabilidades

# Porcentaje
resultados_candidatos["Probability_positive_%"] = (
    resultados_candidatos["Probability_positive"] * 100
)

# ------------------------------------------------------------
# 8. Ordenar por probabilidad
# ------------------------------------------------------------

resultados_candidatos = resultados_candidatos.sort_values(
    by="Probability_positive",
    ascending=False
).reset_index(drop=True)

# ------------------------------------------------------------
# 9. Ranking
# ------------------------------------------------------------

resultados_candidatos.insert(
    0,
    "Rank",
    range(1, len(resultados_candidatos) + 1)
)

# ------------------------------------------------------------
# 10. Mostrar resultados
# ------------------------------------------------------------

print("\nRESULTADOS DE LOS CANDIDATOS")
print("=" * 60)

print(
    resultados_candidatos.head(20).to_string(index=False)
)

# ------------------------------------------------------------
# 11. Distribución de predicciones
# ------------------------------------------------------------

print("\nDistribución de predicciones:")
print(
    resultados_candidatos["Prediction"].value_counts()
)

# ------------------------------------------------------------
# 12. Guardar resultados
# ------------------------------------------------------------

ruta_salida = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\predicciones_141_candidatos_RF.csv"
)

resultados_candidatos.to_csv(
    ruta_salida,
    index=False
)

print("\nArchivo generado:")
print(ruta_salida)

# %%
# ============================================================
# TABLA FINAL DE CANDIDATOS + PREDICCIÓN RF
# ============================================================

ruta_resumen = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\matriz_resumen_funcional_141.csv"
)

# Cargar resumen funcional
resumen = pd.read_csv(ruta_resumen)

print("Resumen funcional:", resumen.shape)

# ------------------------------------------------------------
# Combinar predicciones con resumen funcional
# ------------------------------------------------------------

tabla_final_candidatos = resultados_candidatos.merge(
    resumen,
    on="Organism_code",
    how="left"
)

# ------------------------------------------------------------
# Ordenar por probabilidad
# ------------------------------------------------------------

tabla_final_candidatos = tabla_final_candidatos.sort_values(
    by="Probability_positive",
    ascending=False
).reset_index(drop=True)

# Actualizar ranking
tabla_final_candidatos["Rank"] = range(
    1,
    len(tabla_final_candidatos) + 1
)

# ------------------------------------------------------------
# Mostrar primeros 30
# ------------------------------------------------------------

print("\nTABLA FINAL DE CANDIDATOS")
print("=" * 80)

print(
    tabla_final_candidatos.head(30).to_string(index=False)
)

# ------------------------------------------------------------
# Guardar
# ------------------------------------------------------------

ruta_salida_final = (
    r"C:\Users\USUARIO\Desktop\TFM VARIOS"
    r"\ranking_candidatos_RF_final.csv"
)

tabla_final_candidatos.to_csv(
    ruta_salida_final,
    index=False
)

print("\nArchivo generado:")
print(ruta_salida_final)

# %%
ruta_ranking = r"C:\Users\USUARIO\Desktop\TFM VARIOS\ranking_candidatos_RF_final.csv"

ranking = pd.read_csv(ruta_ranking)

print("=" * 80)
print("ANÁLISIS PARA PRIORIZACIÓN DE CANDIDATOS")
print("=" * 80)

print("\nDimensiones:")
print(ranking.shape)

print("\nColumnas:")
print(ranking.columns.tolist())


# ------------------------------------------------------------
# 1. DISTRIBUCIÓN DE PROBABILIDAD DEL RANDOM FOREST
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("1. PROBABILITY_POSITIVE")
print("=" * 80)

print(
    ranking["Probability_positive"]
    .describe()
    .to_string()
)


# ------------------------------------------------------------
# 2. CATEGORÍAS DE PROBABILIDAD
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("2. DISTRIBUCIÓN DE PROBABILIDAD")
print("=" * 80)

categorias_prob = pd.cut(
    ranking["Probability_positive"],
    bins=[-0.001, 0.50, 0.80, 0.95, 1.001],
    labels=[
        "<50%",
        "50–80%",
        "80–95%",
        ">95%"
    ]
)

print(categorias_prob.value_counts().sort_index())


# ------------------------------------------------------------
# 3. DISTRIBUCIÓN DE KO
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("3. KO_PRESENTES")
print("=" * 80)

print(
    ranking["KO_presentes"]
    .describe()
    .to_string()
)


# ------------------------------------------------------------
# 4. DISTRIBUCIÓN DE MÓDULOS
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("4. MODULES_PRESENTES")
print("=" * 80)

print(
    ranking["Modules_presentes"]
    .describe()
    .to_string()
)


# ------------------------------------------------------------
# 5. DISTRIBUCIÓN DE RUTAS
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("5. PATHWAYS_PRESENTES")
print("=" * 80)

print(
    ranking["Pathways_presentes"]
    .describe()
    .to_string()
)


# ------------------------------------------------------------
# 6. TABLA RESUMEN DE LOS CANDIDATOS
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("6. CORRELACIÓN ENTRE VARIABLES")
print("=" * 80)

variables = [
    "Probability_positive",
    "KO_presentes",
    "Modules_presentes",
    "Pathways_presentes",
    "Genes_objetivo"
]

print(
    ranking[variables]
    .corr()
    .round(3)
    .to_string()
)


# ------------------------------------------------------------
# 7. CANDIDATOS CON MAYOR COBERTURA FUNCIONAL
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("7. MAYOR COBERTURA FUNCIONAL")
print("=" * 80)

top_funcional = ranking.sort_values(
    ["Pathways_presentes", "Modules_presentes", "KO_presentes"],
    ascending=False
).head(20)

print(
    top_funcional[
        [
            "Rank",
            "Organism_code",
            "Probability_positive",
            "KO_presentes",
            "Modules_presentes",
            "Pathways_presentes"
        ]
    ].to_string(index=False)
)
