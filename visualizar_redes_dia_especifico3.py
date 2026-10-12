import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from shapely.geometry import LineString
import numpy as np
import sys
import os

# ================= CONFIGURACIÓN =================

# 1. FECHA Y ARCHIVOS
FECHA_OBJETIVO = "2020_01_14" 
redes_dir = "/data/shaday_data/GIT/movilidad/output/redes_nuevas/"
archivo_red = os.path.join(redes_dir, f"ageb_od_cvegeo_{FECHA_OBJETIVO}_FullDay.csv")
shapefile_path = "/data/shaday_data/GIT/inegi-duckdb/data/MG_2020_Integrado/conjunto_de_datos/00a.shp"

# 2. CONFIGURACIÓN VISUAL
FILTRO_ENTIDADES = ["09", "15"]  
MIN_PESO_LINEA = 2              # <--- ESTE ES EL FILTRO QUE SE MARCARÁ
ANCHO_MAPA = 24                  

# Estilos
COLOR_FONDO_MAPA = "black"       
COLOR_RELLENO_AGEB = "#101010"   
COLOR_BORDE_AGEB = "#333333"     
ANCHO_BORDE_AGEB = 0.3
COLOR_MAP_ARCOS = "plasma"       # Paleta: Morado -> Rosa -> Amarillo
CURVATURA = 0.2                  

# ================= FUNCIONES =================

def crear_curva(p1, p2, curvatura=0.2, resolucion=20):
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    vx, vy = x2 - x1, y2 - y1
    dist = np.sqrt(vx**2 + vy**2)
    if dist == 0: return LineString([p1, p2])
    nx, ny = -vy, vx
    offset_x = (nx / dist) * (dist * curvatura)
    offset_y = (ny / dist) * (dist * curvatura)
    cx, cy = mx + offset_x, my + offset_y
    puntos = []
    for t in np.linspace(0, 1, resolucion):
        px = (1-t)**2 * x1 + 2*(1-t)*t * cx + t**2 * x2
        py = (1-t)**2 * y1 + 2*(1-t)*t * cy + t**2 * y2
        puntos.append((px, py))
    return LineString(puntos)

# ================= PROCESAMIENTO =================

print(f"--- 1. Cargando Mapa Base ---")
try:
    gdf = gpd.read_file(shapefile_path)
except Exception as e:
    sys.exit(f"[ERROR] {e}")

if "CVEGEO" not in gdf.columns and "CVE_AGEB" in gdf.columns:
    gdf["CVEGEO"] = gdf["CVE_ENT"] + gdf["CVE_MUN"] + gdf["CVE_LOC"] + gdf["CVE_AGEB"]

gdf = gdf[gdf["CVE_ENT"].isin(FILTRO_ENTIDADES)].copy()
print(f"   Zona filtrada: {len(gdf)} AGEBs.")
gdf["centroid"] = gdf.geometry.centroid
geo_dict = gdf.set_index("CVEGEO")["centroid"].to_dict()

print(f"\n--- 2. Cargando Datos de Movilidad ({FECHA_OBJETIVO}) ---")
if not os.path.exists(archivo_red):
    sys.exit(f"[ERROR] No existe: {archivo_red}")

df_edges = pd.read_csv(archivo_red, dtype=str)
col_peso = "w" if "w" in df_edges.columns else "weight"
if col_peso not in df_edges.columns:
    df_edges["w"] = 1
    col_peso = "w"
else:
    df_edges[col_peso] = pd.to_numeric(df_edges[col_peso])

# === FILTRO DE CLÚSTERES (OPCIONAL) ===
# Si quieres solo clústeres, descomenta y adapta esto:
# cluster_map_path = "/data/shaday_data/GIT/movilidad/database/cieni_geo_data_7074_FINAL_COMPLETO.csv"
# df_c = pd.read_csv(cluster_map_path, dtype=str)
# agebs_cluster = set(df_c["CVEGEO_fixed"].str.strip().unique())
# df_edges = df_edges[df_edges["source"].isin(agebs_cluster) & df_edges["target"].isin(agebs_cluster)]

# ================= GENERACIÓN DE ARCOS =================
print("\n--- 3. Generando Arcos Curvos ---")

valid_cvegeos = set(gdf["CVEGEO"].unique())
df_plot = df_edges[
    (df_edges[col_peso] > MIN_PESO_LINEA) &
    (df_edges["source"].isin(valid_cvegeos)) &
    (df_edges["target"].isin(valid_cvegeos))
].copy()

arcos_geom = []
arcos_w = []

for _, row in df_plot.iterrows():
    pt_src = geo_dict.get(row["source"])
    pt_tgt = geo_dict.get(row["target"])
    if pt_src and pt_tgt:
        arcos_geom.append(crear_curva(pt_src, pt_tgt, curvatura=CURVATURA))
        arcos_w.append(row[col_peso])

gdf_arcos = gpd.GeoDataFrame({"weight": arcos_w}, geometry=arcos_geom, crs=gdf.crs)
print(f"   Se dibujarán {len(gdf_arcos)} arcos.")

# ================= GRAFICAR =================
print("\n--- 4. Renderizando Mapa Oscuro con Leyenda ---")

fig, ax = plt.subplots(figsize=(ANCHO_MAPA, ANCHO_MAPA), facecolor='black')
ax.set_facecolor('black')

# CAPA 1: Polígonos de AGEBs
gdf.plot(ax=ax, color=COLOR_RELLENO_AGEB, edgecolor=COLOR_BORDE_AGEB, linewidth=ANCHO_BORDE_AGEB, zorder=1)

# CAPA 2: Arcos de Movilidad
if not gdf_arcos.empty:
    gdf_arcos = gdf_arcos.sort_values("weight")
    
    # Normalización para colores y grosores
    vmin, vmax = gdf_arcos["weight"].min(), gdf_arcos["weight"].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Grosor dinámico
    linewidths = (gdf_arcos["weight"] / vmax) * 3.5 + 0.3
    
    # Dibujar arcos
    plot_arcos = gdf_arcos.plot(ax=ax, 
                   column="weight",
                   cmap=COLOR_MAP_ARCOS, 
                   norm=norm,
                   linewidth=linewidths, 
                   alpha=0.75,
                   zorder=2)

    # === AGREGAR BARRA DE COLOR (LEYENDA) ===
    # Creamos un "ScalarMappable" falso solo para la barra de color
    sm = cm.ScalarMappable(cmap=COLOR_MAP_ARCOS, norm=norm)
    sm._A = [] # Array vacío necesario para matplotlib
    
    # Posición de la barra (ajustar shrink para tamaño)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, fraction=0.03)
    
    # Estilo de la barra para fondo oscuro
    cbar.set_label('Intensidad de Flujo (Viajes)', color='white', fontsize=14, labelpad=10)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white', fontsize=12)
    
    # Marcar el mínimo en la barra (opcional)
    # cbar.ax.text(1.3, 0, f"> {MIN_PESO_LINEA}", transform=cbar.ax.transAxes, color='gray', va='center')

# Estética
ax.set_axis_off()
plt.title(f"Intensidad de Movilidad ZMVM: {FECHA_OBJETIVO}", fontsize=24, color="white", weight="bold")

# === LEYENDA DEL FILTRO ===
# Texto explícito en la parte inferior indicando el umbral
plt.figtext(0.5, 0.08, 
            f"Filtro aplicado: Solo se muestran conexiones con > {MIN_PESO_LINEA} viajes\nCurvatura de arcos: {CURVATURA}", 
            ha="center", fontsize=14, color="#AAAAAA")

# Guardar
out_file = f"mapa_arcos_legend_{FECHA_OBJETIVO}.svg"
plt.savefig(out_file, dpi=300, bbox_inches="tight", facecolor='black')
print(f"✅ Mapa guardado: {out_file}")
