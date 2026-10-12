#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
import sys
import time
import datetime as dt
import pandas as pd
import duckdb
from concurrent.futures import ProcessPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(asctime)s:%(message)s")
logger = logging.getLogger(__name__)

# --- Configuración fija ---
DATASET_ROOT = "/data/covid/aws/movement_parquet"
LOCATION_FILE = "/data/shaday_data/GIT/movilidad/database/AGEB_CDMX_EDOMEX_INEGI_2020_COMPLETO.csv"
OUTPUT_DIR = "/data/shaday_data/GIT/movilidad/output/redes_nuevas"
DWELL_TIME = 0   # segundos
SUFFIX = "ageb"
MAX_WORKERS = 4    # ajusta al número de núcleos de tu servidor

# --- Fechas ---
D_START = dt.date(2020, 1, 1)
D_END   = dt.date(2021, 12, 31)


def load_agebs():
    loc_df = pd.read_csv(LOCATION_FILE)
    if "CVEGEO" not in loc_df.columns:
        print("El CSV de 'location' debe contener la columna 'CVEGEO'.")
        sys.exit(1)
    cvegeo_list = sorted(loc_df["CVEGEO"].astype(str).unique())
    if not cvegeo_list:
        print("Lista de CVEGEO vacía.")
        sys.exit(1)
    return "'" + "','".join(cvegeo_list) + "'"


def process_day(cur_date, cvegeo_in):
    """Procesa un día completo y genera CSV"""
    date_str = cur_date.strftime("%Y-%m-%d")
    year, month, day = cur_date.strftime("%Y %m %d").split()

    parquet_glob = f"{DATASET_ROOT}/year={year}/month={month}/day={day}/*"
    tbl_day = f"movement_{year}{month}{day}"
    tbl_tmp1 = f"_tmp1_{year}{month}{day}"

    try:
        con = duckdb.connect()
        con.execute("PRAGMA threads=16;")  # DuckDB usa 16 hilos internos

        # 1) Carga parquet y filtra CVEGEOs
        con.execute(f"""
            CREATE OR REPLACE TABLE {tbl_day} AS
            SELECT *
            FROM read_parquet('{parquet_glob}')
            WHERE cvegeo IN ({cvegeo_in});
        """)

        # 2) Conteo de dispositivos del día
        dfc = con.execute(f"SELECT COUNT(DISTINCT caid) AS total_caid FROM {tbl_day}").df()
        total_caid = int(dfc["total_caid"][0]) if not dfc.empty else 0
        if total_caid == 0:
            return f"[{date_str}] sin dispositivos, se omite."

        # 3) Precálculo de transiciones
        con.execute(f"""
            CREATE OR REPLACE TABLE {tbl_tmp1} AS
            WITH base AS (
                SELECT
                    caid,
                    cvegeo AS source_cvegeo,
                    utc_timestamp,
                    LAG(cvegeo, -1) OVER (PARTITION BY caid ORDER BY utc_timestamp) AS target_cvegeo,
                    LAG(utc_timestamp, -1) OVER (PARTITION BY caid ORDER BY utc_timestamp) AS target_time
                FROM {tbl_day}
            )
            SELECT
                caid,
                source_cvegeo,
                target_cvegeo,
                utc_timestamp,
                target_time,
                (target_time - utc_timestamp) AS dwell_time
            FROM base;
        """)

        # 4) Full day: dwell_time aplicado
        df_od = con.execute(f"""
            SELECT
                source_cvegeo AS source,
                target_cvegeo AS target,
                (COUNT(*) * 1000000.0) / {total_caid} AS w
            FROM {tbl_tmp1}
            WHERE target_cvegeo IS NOT NULL
              AND source_cvegeo <> target_cvegeo
              AND dwell_time > {DWELL_TIME}
            GROUP BY source_cvegeo, target_cvegeo;
        """).df()

        # 5) Exportar CSV
        fname = "_".join([
            SUFFIX + "_od_cvegeo",
            date_str.replace("-", "_"),
            "FullDay"
        ]) + ".csv"
        fpath = os.path.join(OUTPUT_DIR, fname)
        df_od.to_csv(fpath, index=False)

        return f"[{date_str}] listo ({len(df_od)} aristas, {total_caid} dispositivos)."

    except Exception as e:
        return f"[{date_str}] ERROR: {e}"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    cvegeo_in = load_agebs()

    # lista de fechas
    dates = [D_START + dt.timedelta(days=i) for i in range((D_END - D_START).days + 1)]

    t0 = time.perf_counter()
    results = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_day, d, cvegeo_in): d for d in dates}
        for f in as_completed(futures):
            results.append(f.result())

    for r in sorted(results):
        print(r)

    elapsed = time.perf_counter() - t0
    logger.info("FINISHED in %.2f seconds." % elapsed)


if __name__ == "__main__":
    main()

