#%%

import os
import subprocess
import pandas as pd
from pathlib import Path
import requests as re
import shutil
from collections import Counter

#%%
Base_DIR = Path(__file__).resolve().parent

Data_DIR = Base_DIR / "data"
Results_DIR = Base_DIR / "results"

Results_DIR.mkdir(exist_ok=True)

#%%

#======================================
#Decargar los genomas de entrenamiento
#======================================


# Carpeta principal
carpeta_base = (
    Data_DIR / "genomas_entrenamiento"
)

dataset_train = pd.read_csv(
    Results_DIR / "Dataset_train.csv"
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

#Carpeta individual
    carpeta_assembly = os.path.join(
        carpeta_base,
        assembly
    )

    os.makedirs(carpeta_assembly, exist_ok=True)

#Comprobar si existe el protein.faa
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

#Archivo ZIP temporal
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

#Verificar proteína
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


#Convertir resultados a DataFrame
resultados_descarga = pd.DataFrame(resultados_descarga)

#%%
carpeta = Path(
    Data_DIR / "genomas_entrenamiento"
)

protein_files = list(carpeta.rglob("protein.faa"))

print("Protein.faa encontrados:", len(protein_files))

for p in protein_files:
    print(p)


accesiones = []

for p in protein_files:
    coincidencias = re.findall(r"GCA_\d+\.\d+|GCF_\d+\.\d+", str(p))
    if coincidencias:
        accesiones.append(coincidencias[-1])

for acc, n in Counter(accesiones).items():
    if n > 1:
        print(acc, "→", n)

protein_files = list(carpeta.rglob("protein.faa"))

registros_proteinas = []

for p in protein_files:
    
#Buscar accession en la ruta
    coincidencias = re.findall(r"GCA_\d+\.\d+|GCF_\d+\.\d+", str(p))
    
    if not coincidencias:
        continue
    
    accession = coincidencias[-1]
    
    registros_proteinas.append({
        "Assembly_accession": accession,
        "protein_faa": str(p)
    })

tabla_proteinas = pd.DataFrame(registros_proteinas)

#Eliminar duplicado
tabla_proteinas = tabla_proteinas.drop_duplicates(
    subset="Assembly_accession"
).reset_index(drop=True)

#Incorporar organismo y label
tabla_proteinas = tabla_proteinas.merge(
    dataset_train[
        ["Organism", "Assembly_accession", "label"]
    ],
    on="Assembly_accession",
    how="left"
)

#Ordenar columnas
tabla_proteinas = tabla_proteinas[
    ["Organism", "Assembly_accession", "label", "protein_faa"]
]

#%%
#======
#Kofam
#======

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

