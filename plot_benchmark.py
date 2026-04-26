import matplotlib.pyplot as plt
import numpy as np

# Datos del benchmark
modules = ['vcf-stats\n(Diagnóstico)', 'vcf-filter\n(Filtrado)', 'pop-structure\n(PCA/Kinship)', 'snp-tree\n(Filogenia)']
biocenicana_times = [2.37, 1.85, 11.04, 3.07]
ngsep_times = [17.02, 121.97, 25.04, 10.50]

x = np.arange(len(modules))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

# Usar colores vibrantes (verde para BioCenicana, gris oscuro para NGSEP)
rects1 = ax.bar(x - width/2, biocenicana_times, width, label='BioCenicana', color='#2ecc71', edgecolor='black', linewidth=1.5)
rects2 = ax.bar(x + width/2, ngsep_times, width, label='NGSEP 5.1.0', color='#7f8c8d', edgecolor='black', linewidth=1.5)

# Añadir etiquetas, título y leyenda
ax.set_ylabel('Tiempo de Ejecución (Segundos)', fontsize=12, fontweight='bold')
ax.set_title('Benchmark de Rendimiento: BioCenicana vs NGSEP\nGenoma Poliploide (10x, ~50K SNPs)', fontsize=14, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(modules, fontsize=11, fontweight='bold')
ax.legend(fontsize=12)

# Añadir escala logarítmica si es necesario, o anotaciones directas
ax.bar_label(rects1, padding=3, fmt='%.2fs', fontsize=10, fontweight='bold')
ax.bar_label(rects2, padding=3, fmt='%.2fs', fontsize=10, fontweight='bold')

# Mejorar el diseño del grid
ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Escala logarítmica opcional en Y para hacer visible el 121 vs 1.85,
# pero en gráfico de barras normal se ve más dramático.
# plt.yscale('symlog')

fig.tight_layout()

# Guardar la imagen
plt.savefig('docs/assets/benchmark_plot.png', dpi=300, bbox_inches='tight')
print("Gráfico generado con éxito en docs/assets/benchmark_plot.png")
