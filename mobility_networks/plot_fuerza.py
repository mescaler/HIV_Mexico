import pandas as pd
import matplotlib.pyplot as plt

# === Ruta del archivo de entrada ===
csv_path = "/data/shaday_data/GIT/movilidad/output/fuerza_cdmx_edomex_nuevas_redes.csv"

# === Cargar datos ===
df = pd.read_csv(csv_path, parse_dates=["fecha"])
df = df.sort_values("fecha")

# === Gráfica: fuerza absoluta ===
plt.figure(figsize=(10, 5))
plt.plot(df["fecha"], df["CDMX_a_EDOMEX"], label="CDMX → EDOMEX")
plt.plot(df["fecha"], df["EDOMEX_a_CDMX"], label="EDOMEX → CDMX")
plt.title("Fuerza absoluta de movilidad")
plt.ylabel("Movilidad (w)")
plt.xlabel("Fecha")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("/data/shaday_data/GIT/movilidad/output/plot_fuerza_absoluta.svg", format="svg")
plt.close()

# === Gráfica: normalización respecto al total de red ===
plt.figure(figsize=(10, 5))
plt.plot(df["fecha"], df["CDMX_a_EDOMEX_pct_total"], label="CDMX → EDOMEX (%)")
plt.plot(df["fecha"], df["EDOMEX_a_CDMX_pct_total"], label="EDOMEX → CDMX (%)")
plt.title("Movilidad normalizada por total de la red")
plt.ylabel("Porcentaje (%)")
plt.xlabel("Fecha")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("/data/shaday_data/GIT/movilidad/output/plot_normalizado_total.svg", format="svg")
plt.close()

# === Gráfica: normalización bilateral ===
plt.figure(figsize=(10, 5))
plt.plot(df["fecha"], df["CDMX_a_EDOMEX_norm_bilateral"], label="CDMX → EDOMEX")
plt.plot(df["fecha"], df["EDOMEX_a_CDMX_norm_bilateral"], label="EDOMEX → CDMX")
plt.title("Distribución dentro de movilidad bilateral")
plt.ylabel("Proporción (0–1)")
plt.xlabel("Fecha")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("/data/shaday_data/GIT/movilidad/output/plot_normalizado_bilateral.svg", format="svg")
plt.close()

print("[✔] Gráficas SVG actualizadas sin marcadores.")
