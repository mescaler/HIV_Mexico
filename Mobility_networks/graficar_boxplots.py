import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates

# Configuración de estilo similar al paper
sns.set_theme(style="ticks")
plt.rcParams['font.family'] = 'sans-serif'

# 1. Cargar datos
input_file = "estadisticas_diarias_movilidad_2020_2021.csv"
df = pd.read_csv(input_file)

# Filtrar días sin datos o errores
df = df[df["status"] == "OK"].copy()
df["date"] = pd.to_datetime(df["date"])

# Crear columna de Año-Mes para el eje X
df["year_month"] = df["date"].dt.to_period("M").astype(str)
df["year"] = df["date"].dt.year.astype(str) # Para el color (hue)

# Convertir a millones para que se vea igual a la figura
df["devices_m"] = df["devices"] / 1_000_000
df["records_m"] = df["records"] / 1_000_000

# === GRAFICAR ===
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: Dispositivos (Devices)
sns.boxplot(
    data=df, x="year_month", y="devices_m", hue="year",
    ax=axes[0], palette=["#8fceb6", "#f2a87d"],  # Colores verde/naranja similares a tu imagen
    dodge=False, flierprops={"marker": "o", "markersize": 3}
)
axes[0].set_title("a) Devices observed", loc="left", fontsize=12, fontweight="bold")
axes[0].set_ylabel("Devices (Million)", fontsize=11)
axes[0].set_xlabel("Date", fontsize=11)
axes[0].tick_params(axis='x', rotation=45)
axes[0].legend_.remove() # Quitamos leyenda interna para limpiar

# Panel B: Registros (Records/Pings)
sns.boxplot(
    data=df, x="year_month", y="records_m", hue="year",
    ax=axes[1], palette=["#8fceb6", "#f2a87d"],
    dodge=False, flierprops={"marker": "o", "markersize": 3}
)
axes[1].set_title("b) Data records observed daily", loc="left", fontsize=12, fontweight="bold")
axes[1].set_ylabel("Records (Million)", fontsize=11)
axes[1].set_xlabel("Date", fontsize=11)
axes[1].tick_params(axis='x', rotation=45)
axes[1].legend(title="Year", loc="upper right")

plt.tight_layout()
plt.savefig("figura_resumen_dispositivos.svg", dpi=300)
print("Gráfica guardada como figura_resumen_dispositivos.png")
plt.show()
