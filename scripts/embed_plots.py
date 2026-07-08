#!/usr/bin/env python3
"""
embed_plots.py
Reemplaza los <iframe src="plots/xxx.html"> en el visor_sintenia.html
por el contenido HTML de los archivos correspondientes, incrustado directamente
dentro de un <div>. Esto elimina la dependencia de un servidor HTTP y permite
abrir el visor haciendo doble-clic (protocolo file://).
"""

import sys
import re
import os

def extract_body_content(html_content):
    """Extrae solo el contenido entre <body>...</body>; si no hay body, devuelve todo."""
    m = re.search(r'<body[^>]*>(.*?)</body>', html_content, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1)
    return html_content

def extract_head_scripts(html_content):
    """Extrae los bloques <script> y <style> del <head> de un HTML de Plotly."""
    scripts = []
    for m in re.finditer(r'<script[^>]*>.*?</script>', html_content, re.DOTALL | re.IGNORECASE):
        scripts.append(m.group(0))
    for m in re.finditer(r'<style[^>]*>.*?</style>', html_content, re.DOTALL | re.IGNORECASE):
        scripts.append(m.group(0))
    return '\n'.join(scripts)

def embed_plots(visor_path):
    visor_dir = os.path.dirname(os.path.abspath(visor_path))
    plots_dir = os.path.join(visor_dir, 'plots')

    with open(visor_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Idempotency guard: if the plots were already inlined in a previous run,
    # do NOT embed again (that would inject duplicate <script> blocks and
    # duplicate Plotly.newPlot calls into the same page).
    if 'id="embed-iframe-' in content:
        print('  [SKIP] El visor ya está auto-contenido (gráficos incrustados). Nada que hacer.')
        return

    # Map of iframe IDs to their plot filenames
    iframe_map = {
        'iframe-orientation': 'blocks_orientation.html',
        'iframe-dotplot':     'dotplot_synteny_sv.html',
        'iframe-invsize':     'inversions_sizes.html',
        'iframe-sankey':      'gene_sv_snp_network.html',
        'iframe-table':       'sv_gene_table.html',
    }

    injected_scripts = []

    for iframe_id, plot_file in iframe_map.items():
        plot_path = os.path.join(plots_dir, plot_file)
        if not os.path.exists(plot_path):
            print(f'  [WARN] Plot not found, skipping: {plot_path}')
            continue

        with open(plot_path, 'r', encoding='utf-8') as f:
            plot_html = f.read()

        # Extract <script> and <style> blocks to inject in head later
        head_content = extract_head_scripts(plot_html)
        if head_content:
            injected_scripts.append(f'\n<!-- Embedded scripts from {plot_file} -->\n{head_content}')

        # Extract body content
        body = extract_body_content(plot_html)

        # Build inline div replacement
        inline_div = (
            f'<div id="embed-{iframe_id}" '
            f'style="width:100%; border-radius:12px; background:#ffffff; '
            f'box-shadow:var(--shadow-lg); overflow:auto;">\n'
            f'{body}\n'
            f'</div>'
        )

        # Replace <iframe id="iframe-xxx" ...></iframe> with the inline div
        pattern = rf'<iframe\s+id="{re.escape(iframe_id)}"[^>]*>\s*</iframe>'
        inline_div_final = inline_div  # capture for lambda
        new_content, n = re.subn(pattern, lambda m: inline_div_final, content, flags=re.DOTALL)
        if n == 0:
            # Try with src="" (lazy-loaded iframes have empty src initially)
            pattern2 = rf'<iframe\s+id="{re.escape(iframe_id)}"[^>]*/?>'
            new_content, n = re.subn(pattern2, lambda m: inline_div_final, content, flags=re.DOTALL)
        if n > 0:
            content = new_content
            print(f'  [OK] Embedded {plot_file} -> #{iframe_id}')
        else:
            print(f'  [WARN] iframe #{iframe_id} not found in visor HTML')

    # Inject collected scripts just before </head>
    if injected_scripts:
        all_scripts = '\n'.join(injected_scripts)
        content = content.replace('</head>', f'{all_scripts}\n</head>', 1)

    # Write back
    with open(visor_path, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f'\nVisor auto-contenido escrito en: {visor_path}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Uso: python3 embed_plots.py <ruta/visor_sintenia.html>')
        sys.exit(1)
    embed_plots(sys.argv[1])
