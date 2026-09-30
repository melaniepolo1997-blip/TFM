#%%
#======================
#Cargar librerías 
#======================
from pathlib import Path
import pandas as pd
import shutil
import re
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import make_scorer, precision_score, recall_score, f1_score
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt


#%%
Base_DIR = Path(__file__).resolve().parent

Data_DIR = Base_DIR / "data"
Results_DIR = Base_DIR / "results"
carpeta_base = Data_DIR / "genomas_entrenamiento"
dataset_train = pd.read_csv(
    Results_DIR / "Dataset_train.csv"
)
carpeta = Path(
    Data_DIR / "genomas_entrenamiento"
)
carpeta_kofam = carpeta / "proteinas_kofam"

Results_DIR.mkdir(exist_ok=True)

#%%

#==============
#Entrenamiento
#==============

#Localizar archivos protein.faa
protein_files = list(carpeta_base.rglob("protein.faa"))

#Construir la trabla de archivos

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

#Eliminar duplicados
tabla_archivos = (
    tabla_archivos
    .drop_duplicates(subset="Assembly_accession")
    .reset_index(drop=True)
)

#Unir con el dataset de entrenamiento 

tabla_archivos = tabla_archivos.merge(
    dataset_train[
        ["Organism", "Assembly_accession", "label"]
    ],
    on="Assembly_accession",
    how="left"
)

#Copiar en la carpeta Kofam
rutas_kofam = []

for _, fila in tabla_archivos.iterrows():

    accession = fila["Assembly_accession"]
    origen = Path(fila["protein_original"])

    destino = carpeta_kofam / f"{accession}.faa"

    if not destino.exists():
        shutil.copy2(origen, destino)

    rutas_kofam.append(str(destino))

tabla_archivos["protein_kofam"] = rutas_kofam

#%%
#Configuración 
#Localizar los protein.faa
protein_files = list(carpeta_base.rglob("protein.faa"))

#Construir la tabla de archivos 
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

#Eliminar duplicados
tabla_archivos = (
    tabla_archivos
    .drop_duplicates(subset="Assembly_accession")
    .reset_index(drop=True)
)

#Unir con el dataset de entrenamiento 
tabla_archivos = tabla_archivos.merge(
    dataset_train[
        ["Organism", "Assembly_accession", "label"]
    ],
    on="Assembly_accession",
    how="left"
)

#Copiar en la carpeta kofam 

rutas_kofam = []

for _, fila in tabla_archivos.iterrows():

    accession = fila["Assembly_accession"]
    origen = Path(fila["protein_original"])

    destino = carpeta_kofam / f"{accession}.faa"

    if not destino.exists():
        shutil.copy2(origen, destino)

    rutas_kofam.append(str(destino))

tabla_archivos["protein_kofam"] = rutas_kofam

#Lista para el Kofamkoala 

tabla_kofam = tabla_archivos[
    ["Organism", "Assembly_accession", "label", "protein_kofam"]
].copy()

tabla_kofam = tabla_kofam.sort_values(
    "Assembly_accession"
).reset_index(drop=True)

print(
    tabla_kofam.to_string(index=False)
)

archivo_lista_kofam = (
    carpeta_base / "lista_53_genomas_KofamKOALA.csv"
)

tabla_kofam.to_csv(
    archivo_lista_kofam,
    index=False,
    encoding="utf-8-sig"
)

ruta = Results_DIR / "X_ML_107_KO.csv"

X = pd.read_csv(ruta)

ruta_ko = Results_DIR /"KO_objetivo_107.txt"

with open(ruta_ko, "w") as f:
    for ko in X.columns:
        f.write(ko + "\n")
        
matriz_kofamscan = pd.read_csv(
    Results_DIR / "MATRIZ_KOFAMSCAN_53x78.csv"
)
df = pd.read_csv(matriz_kofamscan)

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

#%%
#Validación cruzada 

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

#===================================
#Entrenamiento modelo Random forest
#===================================
rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

#Métricas y validación 

scoring = {
    "accuracy": "accuracy",
    "precision": make_scorer(precision_score, zero_division=0),
    "recall": make_scorer(recall_score, zero_division=0),
    "f1": make_scorer(f1_score, zero_division=0),
    "roc_auc": "roc_auc"
}

resultados_rf = cross_validate(
    rf,
    X,
    y,
    cv=cv,
    scoring=scoring,
    return_train_score=False,
    n_jobs=-1
)

for metrica in scoring:
    valores = resultados_rf[f"test_{metrica}"]
    
    print(f"\n{metrica.upper()}:")
    print("Por fold:", np.round(valores, 3))
    print("Media:", round(valores.mean(), 3))
    print("Desviación estándar:", round(valores.std(), 3))
    
#===================================
#Entrenamiento modelo SVM
#===================================

svm = Pipeline([
    ("scaler", StandardScaler()),
    ("classifier", SVC(
        kernel="rbf",
        probability=True,
        class_weight="balanced",
        random_state=42
    ))
])

#Métricas
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

#Validación

resultados_svm = cross_validate(
    svm,
    X,
    y,
    cv=cv,
    scoring=scoring,
    return_train_score=False,
    n_jobs=-1
)
for metrica in scoring:
    valores = resultados_svm[f"test_{metrica}"]
    print(f"\n{metrica.upper()}:")
    print("Por fold:", np.round(valores, 3))
    print("Media:", round(valores.mean(), 3))
    print("Desviación estándar:", round(valores.std(), 3))
    
#%%
#Random Forest-Importancia variables 
rf_final = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

#Entrenamiento con los 53 organismos
rf_final.fit(X, y)

#Importancia de cada KO
importancias_rf = pd.DataFrame({
    "KO": X.columns,
    "Importance": rf_final.feature_importances_
})

#Ordenar de mayor a menor
importancias_rf = importancias_rf.sort_values(
    by="Importance",
    ascending=False
).reset_index(drop=True)

print(importancias_rf.head(20).to_string(index=False))
 
#Predicciones out-of-fold
pred_rf = cross_val_predict(
    rf,
    X,
    y,
    cv=cv,
    method="predict",
    n_jobs=-1
)

#Matriz de confusión
cm_rf = confusion_matrix(y, pred_rf)
print(cm_rf)

print("\nInterpretación:")
print("Verdaderos negativos (TN):", cm_rf[0, 0])
print("Falsos positivos (FP):", cm_rf[0, 1])
print("Falsos negativos (FN):", cm_rf[1, 0])
print("Verdaderos positivos (TP):", cm_rf[1, 1])

#Visualización
disp = ConfusionMatrixDisplay(
    confusion_matrix=cm_rf,
    display_labels=["Negativo", "Positivo"]
)

disp.plot()
plt.title("Matriz de confusión - Random Forest")
plt.show()

#Modeli final Random Forest

rf_final = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1
)

# Entrenar con todos los organismos de referencia
rf_final.fit(X, y)
#%%

#============================================
#Entrenamiento con los organismos candidatos 
#============================================

ruta_candidatos = (
    Results_DIR/"MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv"
)

candidatos = pd.read_csv(ruta_candidatos)

KO_modelo = X.columns.tolist()

#Comprobar que todos los candidatos tengan KO

KO_faltantes = [
    ko for ko in KO_modelo
    if ko not in candidatos.columns
]

if len(KO_faltantes) > 0:
    print(KO_faltantes)
    raise ValueError(
        "Faltan KO necesarios para realizar la predicción."
    )
    
#Contruir matriz X de candidatos 
X_candidatos = candidatos[KO_modelo].copy()

#%%

#===========
#Predicción
#===========

predicciones = rf_final.predict(X_candidatos)

#Probabilidad de pertenecer a clase positiva
probabilidades = rf_final.predict_proba(X_candidatos)[:, 1]

#Crear tabla de resultados 

resultados_candidatos = candidatos[
    ["Organism_code"]
].copy()

resultados_candidatos["Prediction"] = predicciones
resultados_candidatos["Probability_positive"] = probabilidades

# Porcentaje
resultados_candidatos["Probability_positive_%"] = (
    resultados_candidatos["Probability_positive"] * 100
)

#Ordenar por probabilidad
resultados_candidatos = resultados_candidatos.sort_values(
    by="Probability_positive",
    ascending=False
).reset_index(drop=True)

#Ranking 
resultados_candidatos.insert(
    0,
    "Rank",
    range(1, len(resultados_candidatos) + 1)
)

#Resultados 

print("Resultados Obtenidos:")
print(
    resultados_candidatos.head(20).to_string(index=False)
)

#Guardar resultados 
resultados_candidatos.to_csv(
    Results_DIR,
    index=False
)

#%%

#======================================================
#Contrucción de tabla final de candidatos y modelo RF
#======================================================


#Combinar predicciones con resumen funcional
tabla_final_candidatos = resultados_candidatos.merge(
    resumen,
    on="Organism_code",
    how="left"
)

#Ordenar por probabilidad
tabla_final_candidatos = tabla_final_candidatos.sort_values(
    by="Probability_positive",
    ascending=False
).reset_index(drop=True)

#Actualizar ranking
tabla_final_candidatos["Rank"] = range(
    1,
    len(tabla_final_candidatos) + 1
)

print(
    tabla_final_candidatos.head(30).to_string(index=False)
)

#Guardar

ruta_salida_final = (
    Results_DIR/"ranking_candidatos_RF_final.csv"
)

tabla_final_candidatos.to_csv(
    ruta_salida_final,
    index=False
)
print(ruta_salida_final)


#Distribución de probabilidad de RF

ruta_ranking = Results_DIR/"ranking_candidatos_RF_final.csv"
ranking = pd.read_csv(ruta_ranking)


print(
    ranking["Probability_positive"]
    .describe()
    .to_string()
)

#Categorías de probabilidad 

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

#Distribución KO

print(
    ranking["Modules_presentes"]
    .describe()
    .to_string()
)

#Distribución de rutas

print(
    ranking["Pathways_presentes"]
    .describe()
    .to_string()
)

#Tabla resumen 

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

#Organismos con mayor cobertura funcional 
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