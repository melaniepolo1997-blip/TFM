# %%
#======================
#Cargar librerías 
#======================

import pandas as pd
import time
import re
from Bio import Entrez
from pathlib import Path
import unicodedata

Entrez.email = "melaniepolo1997@gmail.com"

#%%
Base_DIR = Path(__file__).resolve().parent

Data_DIR = Base_DIR / "data"
Results_DIR = Base_DIR / "results"

Results_DIR.mkdir(exist_ok=True)

# %%

#==============
#PRIMERA PARTE
#==============

#=================================
#Búsqueda/selección de assemblies
#=================================

ruta = Data_DIR/"Organismo_referencia_positivos y negativos.xlsx"

organismos_referencia = pd.read_excel(ruta)

print("Archivo cargado correctamente")
print("Filas:", len(organismos_referencia))
print("Columnas:", organismos_referencia.columns.tolist())

print("\nPrimeros 10 organismos:")
print(organismos_referencia.head(10))

# %%
#==========================================
#Revisar valores nulos, duplicados y labels
#==========================================

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
#%%

#=======================
#Extraer la información
#=======================

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
#Creamos un dataframe 
tabla_assemblies = pd.DataFrame(resultados_todos)

print("Filas:", len(tabla_assemblies))
print("Organismos:", tabla_assemblies["Organism_reference"].nunique())

tabla_assemblies.head()   
#%%
#Revisamos los resultados del Assembly
resumen_assemblies = (
    tabla_assemblies
    .groupby("Organism_reference", dropna=False)
    .agg(
        N_Assembly=("Assembly_ID", lambda x: x.notna().sum())
    )
    .reset_index()
)

print(resumen_assemblies.to_string(index=False))

#Organismos con 1 Assembly
un_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 1
]

print("Organismos con 1 Assembly:")
print(len(un_assembly))

print(un_assembly.to_string(index=False))

#Organismos sin Assembly
sin_assembly = resumen_assemblies[
    resumen_assemblies["N_Assembly"] == 0
]

print("ORGANISMOS SIN ASSEMBLY:")
print(len(sin_assembly))

print(sin_assembly.to_string(index=False))

#Organismos con múltiples Aseembly 
multiples_assemblies = resumen_assemblies[
    resumen_assemblies["N_Assembly"] > 1
]

print("Organismos con múltiples Assembly:")
print(len(multiples_assemblies))

print(multiples_assemblies.to_string(index=False))

#%%
#===================
#Normalizar nombres 
#===================

def normalizar_nombre(nombre):
    
    if pd.isna(nombre):
        return ""
    
    nombre = str(nombre).lower()
    
    # eliminar información entre paréntesis
    nombre = re.sub(r"\([^)]*\)", "", nombre)
    
    # normalizar espacios
    nombre = re.sub(r"\s+", " ", nombre).strip()
    
    return nombre

tabla_assemblies["Organism_ref_norm"] = (
    tabla_assemblies["Organism_reference"]
    .apply(normalizar_nombre)
)

tabla_assemblies["Organism_ncbi_norm"] = (
    tabla_assemblies["Organism_NCBI"]
    .apply(normalizar_nombre)
)

#%%

#=============================
#Clasificar las coincidencias
#=============================

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
tabla_assemblies.to_excel(
    Results_DIR / "tabla_assemblies_todos_organismos.xlsx",
    index=False
)

print("Total organismos:", organismos_referencia["Organism"].nunique())
print(
    resumen_assemblies["N_Assembly"]
    .value_counts()
    .sort_index()
)

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

#%%

#Casos con múltiples Aseembly 

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

#%%

#==============
#SEGUNDA PARTE
#==============

#====================================
#Descarga y preparación de proteínas
#====================================

#Normalizar nombres 

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

#Agrupamos

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

print(organismos_referencia_limpios.to_string(index=False))

#Organismos sin Assembly 
organismos_sin_assembly = organismos_referencia_limpios[
    ~organismos_referencia_limpios["Organism"].isin(
        tabla_assemblies["Organism_reference"].unique()
    )
]

print(
    organismos_sin_assembly[
        ["Organism", "label", "Reference"]
    ].to_string(index=False)
)

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

#Organismos con un único Assembly 
un_assembly = resumen_assembly[
    resumen_assembly["N_Assembly"] == 1
]

print("Organismos con un único Assembly:", len(un_assembly))

print(
    un_assembly.to_string(index=False)
)
#Organismos con múltiples Assembly 
multiples_assembly = resumen_assembly[
    resumen_assembly["N_Assembly"] > 1
]

print("Organismos con múltiples Assembly:", len(multiples_assembly))

print(
    multiples_assembly.to_string(index=False)
)

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

tabla_assemblies_unica = (
    tabla_assemblies
    .drop_duplicates(
        subset=["Organism_reference", "Assembly_accession"]
    )
    .copy()
)

print("Registros únicos de Assembly:", len(tabla_assemblies_unica))

#%%
#Creamos la tabla de organismos

organismos_genomas = tabla_assemblies_unica.merge(
    organismos_referencia_limpios[
        ["Organism", "label", "Reference"]
    ],
    left_on="Organism_reference",
    right_on="Organism",
    how="inner"
)

organismos_genomas["Coincidencia_nombre"] = (
    organismos_genomas["Organism_reference"]
    .str.lower()
    .str.replace("‐", "-", regex=False)
    .apply(
        lambda x: x.split(" ")[0] if isinstance(x, str) else ""
    )
)

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

print(
    "Organismos únicos:",
    tabla_candidatos["Organism_reference"].nunique()
)

#Contar assembly por organismo

conteo_assembly = (
    tabla_candidatos
    .groupby("Organism_reference")["Assembly_accession"]
    .nunique()
    .reset_index(name="N_Assembly")
)

print(conteo_assembly["N_Assembly"].value_counts().sort_index())

#Separar grupos

grupo_unico = conteo_assembly[
    conteo_assembly["N_Assembly"] == 1
]

grupo_multiple = conteo_assembly[
    conteo_assembly["N_Assembly"] > 1
]

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

#%% 
#Selección automática del ASSEMBLY
# Prioridad del nivel de ensamblaje
prioridad_status = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1
}

seleccion = tabla_candidatos.copy()

#Prioridad en el nivel de Assembly

seleccion["Prioridad_status"] = (
    seleccion["Assembly_status"]
    .map(prioridad_status)
    .fillna(0)
)

#Prioridad de Refseq: GCF = RefSeq y GCA = GenBank

seleccion["Es_RefSeq"] = (
    seleccion["Assembly_accession"]
    .astype(str)
    .str.startswith("GCF_")
    .astype(int)
)

#Versión de assembly

seleccion["Version_Assembly"] = (
    seleccion["Assembly_accession"]
    .astype(str)
    .str.extract(r"\.(\d+)$")[0]
    .fillna("0")
    .astype(int)
)

#Seleccionar el mejor Assembly

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
    tabla_assembly_final["Tipo_seleccion"]
    .value_counts()
)

#%% 

#Marcar assembly para revisión

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

#Excluir organismos ambiguo 
tabla_assembly_final = tabla_assembly_final[
    tabla_assembly_final["Organism_reference"] !=
    "Pseudomonas putida CSV86"
].copy()

#Tabla final de organismos 
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

#Guardar
ruta_salida = Results_DIR / "Dataset_referencia_144_organismos.xlsx"

dataset_referencia.to_excel(
    ruta_salida,
    index=False
)

#%%

#Detectar Assembly duplicados 

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

#Detectar Assembly vacíos 

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

#Buscar Assemblies de los organismos faltantes 
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

#Extraer información

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

#Guardar

tabla_faltantes = pd.DataFrame(resultados_faltantes)

print(tabla_faltantes.to_string(index=False))

#%%

#======================================
#Dataset candidatos para entrenamiento
#======================================

dataset_candidatos = dataset_referencia[
    dataset_referencia["Assembly_accession"].notna() &
    dataset_referencia["Assembly_accession"].ne("")
].copy()

print("\nAssembly únicos:",
      dataset_candidatos["Assembly_accession"].nunique())

#Eliminar el duplicado 

dataset_candidatos = dataset_candidatos[
    dataset_candidatos["Organism"] != "Mycobacterium sp. PYR-1"
].copy()

print("\nAssembly únicos:",
      dataset_candidatos["Assembly_accession"].nunique())

#Revisar organismos de referencia positivos y negativos 

negativos = dataset_candidatos[
    dataset_candidatos["label"] == 0
].copy()

positivos = dataset_candidatos[
    dataset_candidatos["label"] == 1
].copy()

#%%
#===============================
#Crear dataset de entrenamiento 
#===============================

#Copia de los organismos que tienen Assembly
dataset_entrenamiento = dataset_candidatos.copy()

dataset_entrenamiento = dataset_entrenamiento[
    dataset_entrenamiento["Organism"] != "Mycobacterium sp. PYR-1"
].copy()

#Separar los positivos y negativos

positivos_train = dataset_entrenamiento[
    dataset_entrenamiento["label"] == 1
].copy()

negativos_disponibles = dataset_entrenamiento[
    dataset_entrenamiento["label"] == 0
].copy()

#Seleccionar 29 organismos negativos 
negativos_train = (
    negativos_disponibles
    .sample(
        n=29,
        random_state=42
    )
    .copy()
)

#Unir organismos positivos y negativos
dataset_train = pd.concat(
    [
        positivos_train,
        negativos_train
    ],
    ignore_index=True
)
dataset_train = dataset_train.sort_values(
    by=["label", "Organism"]
).reset_index(drop=True)

dataset_train.to_csv(
    Results_DIR / "Dataset_train.csv",
    index=False,
    encoding="utf-8-sig"
)

#Identificar duplicados en el set de entrenamiento 
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
#Eliminar duplicados 
# Identificar el Assembly repetido
conteo = dataset_train["Assembly_accession"].value_counts()

assemblies_duplicados = conteo[
    conteo > 1
].index.tolist()

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

#Guardar dataset de entrenamiento 

ruta_train = (
    Data_DIR,
    Results_DIR/"dataset_entrenamiento_58.csv"
)

dataset_train.to_csv(
    ruta_train,
    index=False)

