import csv
import glob
import os
from pathlib import Path

#Rutas

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent

LISTA_GENOMAS = (
    REPO_DIR
    / "data"
    / "lista_53_genomas_KofamKOALA.csv"
)

KO_FILE = SCRIPT_DIR / "KO_objetivo_107.txt"

RESULTS_DIR = SCRIPT_DIR / "resultados_53"

OUTPUT_FILE = (
    REPO_DIR
    / "results"
    / "MATRIZ_KOFAMSCAN_53x107.csv"
)


#Directorio de salida

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)


#Leer lista de genomas

with open(
    LISTA_GENOMAS,
    "r",
    encoding="utf-8-sig",
    newline=""
) as f:

    reader = csv.DictReader(f)
    lista_53 = list(reader)

print("Genomas en lista:", len(lista_53))

#Leer lista de KO objetivo
with open(
    KO_FILE,
    "r",
    encoding="utf-8"
) as f:

    KO_objetivo = [
        line.strip()
        for line in f
        if line.strip()
    ]

KO_set = set(KO_objetivo)

print("KO objetivo:", len(KO_objetivo))

#Resultados Kofamscan

resultados = sorted(
    glob.glob(
        os.path.join(
            RESULTS_DIR,
            "*_kofam107.tsv"
        )
    )
)

print(
    "Resultados KofamScan encontrados:",
    len(resultados)
)

#Crear diccionario de accessión

resultado_por_accession = {}

for archivo in resultados:

    nombre = os.path.basename(archivo)

    accession = nombre.replace(
        "_kofam107.tsv",
        ""
    )

    resultado_por_accession[accession] = archivo


Construir Matriz

filas = []
advertencias = []


for registro in lista_53:

    organism = registro["Organism"].strip()
    accession = registro["Assembly_accession"].strip()
    label = registro["label"].strip()

#Verificar los resultados
  
    if accession not in resultado_por_accession:

        advertencias.append(
            f"Sin resultado KofamScan: {accession}"
        )

        continue

    archivo_kofam = (
        resultado_por_accession[accession]
    )

Inicializar los KO en 0
  
    ko_presentes = {
        ko: 0
        for ko in KO_objetivo
    }

#Leero resultados de Kofam

    with open(
        archivo_kofam,
        "r",
        encoding="utf-8"
    ) as f:

        for linea in f:

            if linea.startswith("#"):
                continue

            linea = linea.strip()

            if not linea:
                continue

            campos = linea.split()

#hits aceptados
          
            if campos[0] == "*":

                if len(campos) >= 3:

                    ko = campos[2]

                    if ko in KO_set:
                        ko_presentes[ko] = 1

#Construir fila 

    fila = {
        "Organism": organism,
        "Genome_ID": accession,
        "Label": label
    }

    fila.update(ko_presentes)

    filas.append(fila)


#Guardar la matriz

fieldnames = [
    "Organism",
    "Genome_ID",
    "Label"
] + KO_objetivo


with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8-sig",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(filas)

#Comprobar

print()
print("Filas generadas:", len(filas))
print("Columnas generadas:", len(fieldnames))

#Número de KO por organismo

print()
print("KO aceptados por organismo:")

for fila in filas:

    n_ko = sum(
        int(fila[ko])
        for ko in KO_objetivo
    )

    print(
        f"{fila['Genome_ID']}: "
        f"{n_ko} KO | "
        f"{fila['Label']} | "
        f"{fila['Organism']}"
    )


#Advertencias 

if advertencias:

    print()
    print("ADVERTENCIAS:")

    for advertencia in advertencias:
        print(advertencia)

#Archivo final

print()
print("Archivo generado:")
print(OUTPUT_FILE)
