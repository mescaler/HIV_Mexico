import pandas as pd
import os
from glob import glob
from datetime import datetime

# =========================
# Configuración
# =========================

# Ruta de las NUEVAS REDES
input_dir = "/data/shaday_data/GIT/movilidad/output/redes_nuevas/"
movilidad_files = sorted(glob(os.path.join(input_dir, "ageb_od_cvegeo_*_FullDay.csv")))

# Archivo de mapeo (se mantiene igual)
agebs_path = "/data/shaday_data/GIT/movilidad/agebs_cdmx_edomex.csv"

# Archivo de salida (nombre actualizado para no sobrescribir el anterior)
output_file = "/data/shaday_data/GIT/movilidad/output/fuerza_cdmx_edomex_nuevas_redes.csv"

# =========================
# Cargar agebs con su estado
# =========================
print(f"--- Cargando mapeo de estados desde {agebs_path} ---")
df_agebs = pd.read_csv(agebs_path, dtype={"CVEGEO": str})
# Limpieza de espacios en blanco por si acaso
df_agebs["CVEGEO"] = df_agebs["CVEGEO"].str.strip()
cvegeo_to_state = dict(zip(df_agebs["CVEGEO"], df_agebs["state"]))

# =========================
# Procesamiento
# =========================
resultados = []
print(f"--- Procesando {len(movilidad_files)} archivos ---")

for mov_file in movilidad_files:
    basename = os.path.basename(mov_file)
    
    # 1. Extraer fecha (formato nuevo: ageb_od_cvegeo_2020_08_21_FullDay.csv)
    try:
        date_str = basename.replace("ageb_od_cvegeo_", "").replace("_FullDay.csv", "")
        date = datetime.strptime(date_str, "%Y_%m_%d").date()
    except ValueError:
        print(f"[!] Saltando archivo con nombre inesperado: {basename}")
        continue

    # 2. Leer archivo
    df = pd.read_csv(mov_file, dtype={"source": str, "target": str})
    
    # Limpieza de códigos GEO
    df["source"] = df["source"].str.strip()
    df["target"] = df["target"].str.strip()

    # 3. Detectar columna de peso dinámicamente
    weight_col = None
    for cand in ["w", "dwell_10", "weight", "viajes", "count"]:
        if cand in df.columns:
            weight_col = cand
            break
    
    if weight_col is None:
        print(f"[!] Saltando {basename}: No se encontró columna de peso (w, dwell_10, etc)")
        continue

    # 4. Mapear estados
    df["state_source"] = df["source"].map(cvegeo_to_state)
    df["state_target"] = df["target"].map(cvegeo_to_state)

    # Filtrar solo filas donde ambos (origen y destino) tienen estado identificado
    df_valid = df.dropna(subset=["state_source", "state_target"])
    
    # Suma total de flujo en el día (solo entre los estados de interés mapeados)
    total_w = df_valid[weight_col].sum()

    # 5. Calcular flujos específicos
    # CDMX -> EDOMEX
    cdmx_a_edomex = df_valid.query("state_source == 'Mexico City' and state_target == 'State of Mexico'")[weight_col].sum()
    
    # EDOMEX -> CDMX
    edomex_a_cdmx = df_valid.query("state_source == 'State of Mexico' and state_target == 'Mexico City'")[weight_col].sum()

    # 6. Normalización general (Respecto al flujo total válido del día)
    cdmx_a_edomex_norm = cdmx_a_edomex / total_w if total_w > 0 else 0
    edomex_a_cdmx_norm = edomex_a_cdmx / total_w if total_w > 0 else 0

    # 7. Normalización bilateral (Solo considerando la interacción entre estas dos entidades)
    # Es decir: De todo el tráfico que hubo ENTRE Cdmx y Edomex, ¿cuánto fue en cada dirección?
    total_bilateral = cdmx_a_edomex + edomex_a_cdmx
    cdmx_a_edomex_bilateral_norm = cdmx_a_edomex / total_bilateral if total_bilateral > 0 else 0
    edomex_a_cdmx_bilateral_norm = edomex_a_cdmx / total_bilateral if total_bilateral > 0 else 0

    resultados.append({
        "fecha": date,
        "CDMX_a_EDOMEX": cdmx_a_edomex,
        "EDOMEX_a_CDMX": edomex_a_cdmx,
        "CDMX_a_EDOMEX_norm_total": cdmx_a_edomex_norm,
        "EDOMEX_a_CDMX_norm_total": edomex_a_cdmx_norm,
        "CDMX_a_EDOMEX_pct_total": cdmx_a_edomex_norm * 100,
        "EDOMEX_a_CDMX_pct_total": edomex_a_cdmx_norm * 100,
        "CDMX_a_EDOMEX_norm_bilateral": cdmx_a_edomex_bilateral_norm,
        "EDOMEX_a_CDMX_norm_bilateral": edomex_a_cdmx_bilateral_norm
    })

# === Guardar resultados ===
if resultados:
    df_resultados = pd.DataFrame(resultados).sort_values("fecha")
    df_resultados.to_csv(output_file, index=False)
    print(f"[✔] Procesamiento terminado. Archivo guardado: {output_file}")
else:
    print("[!] No se generaron resultados. Verifica los archivos de entrada.")
