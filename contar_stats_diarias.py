#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging, os, sys, time, datetime as dt
import pandas as pd, duckdb
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- CONFIGURACIÓN BASADA EN TU SCRIPT ANTERIOR ---
DATASET_ROOT = "/mnt/mobile_data/movement_parquet"
LOCATION_FILE = "/home/ShadayGuerrero/CRASH/redes/ageb_cdmx_edomex_urbanos.csv"
OUTPUT_FILE = "estadisticas_diarias_movilidad_2020_2021.csv"
MAX_WORKERS = 10       # Puedes subir esto ya que la tarea requiere menos RAM que las redes
DUCK_THREADS = 4

logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s")
logger = logging.getLogger(__name__)

def load_agebs():
    """Carga los AGEBs para asegurar que filtramos la misma zona que las redes"""
    if not os.path.exists(LOCATION_FILE):
        logger.error(f"No se encuentra: {LOCATION_FILE}")
        sys.exit(1)
    loc_df = pd.read_csv(LOCATION_FILE)
    cvegeo_list = sorted(loc_df["CVEGEO"].astype(str).unique())
    return "'" + "','".join(cvegeo_list) + "'"

def count_day_stats(cur_date, cvegeo_in):
    date_str = cur_date.strftime("%Y-%m-%d")
    year, month, day = cur_date.strftime("%Y %m %d").split()
    
    # Ruta al archivo parquet
    parquet_glob = f"{DATASET_ROOT}/year={year}/month={month}/day={day}/*"
    
    # Verificación rápida de existencia
    path_check = f"{DATASET_ROOT}/year={year}/month={month}/day={day}"
    if not os.path.exists(path_check):
        return {"date": date_str, "devices": 0, "records": 0, "status": "No Data"}

    try:
        con = duckdb.connect()
        con.execute(f"PRAGMA threads={DUCK_THREADS};")
        
        # Query para contar Pings (records) y Dispositivos (devices) en una sola pasada
        # Filtramos por cvegeo IN (...) para que coincida con tu estudio
        query = f"""
            SELECT 
                COUNT(*) as num_records,
                COUNT(DISTINCT caid) as num_devices
            FROM read_parquet('{parquet_glob}')
            WHERE cvegeo IN ({cvegeo_in})
        """
        
        df = con.execute(query).df()
        
        recs = int(df["num_records"][0]) if not df.empty else 0
        devs = int(df["num_devices"][0]) if not df.empty else 0
        
        return {"date": date_str, "devices": devs, "records": recs, "status": "OK"}

    except Exception as e:
        return {"date": date_str, "devices": 0, "records": 0, "status": f"Error: {e}"}

def main():
    logger.info("--- Iniciando conteo de estadísticas ---")
    cvegeo_in = load_agebs()
    
    # Rango de fechas: Todo 2020 y hasta Marzo 2021 (según tu paper)
    d_start = dt.date(2020, 1, 1)
    d_end = dt.date(2021, 12, 31) 
    
    dates = [d_start + dt.timedelta(days=i) for i in range((d_end - d_start).days + 1)]
    results = []
    
    t0 = time.perf_counter()
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(count_day_stats, d, cvegeo_in): d for d in dates}
        
        for i, f in enumerate(as_completed(futures)):
            res = f.result()
            results.append(res)
            if i % 10 == 0:
                logger.info(f"Procesado {res['date']}: {res['devices']} devs, {res['records']} pings")

    # Guardar CSV final
    df_res = pd.DataFrame(results)
    df_res = df_res.sort_values("date")
    df_res.to_csv(OUTPUT_FILE, index=False)
    
    logger.info(f"Guardado en: {OUTPUT_FILE}")
    logger.info("Tiempo total: %.2f s" % (time.perf_counter()-t0))

if __name__ == "__main__":
    main()
