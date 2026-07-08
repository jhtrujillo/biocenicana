#!/usr/bin/env python3
"""
inline_plotly.py
Hace que el visor sea 100% autónomo: reemplaza las etiquetas
<script src="https://cdn.plot.ly/plotly-X.Y.Z.min.js"> por la librería
Plotly incrustada localmente. Así el visor carga con doble-clic (file://)
sin servidor y sin conexión a internet.

La librería se toma de scripts/vendor/plotly-<version>.min.js.

Uso: python3 inline_plotly.py <ruta/visor_sintenia.html>
"""

import sys
import re
import os

VENDOR_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vendor')

# Detecta cualquier <script ... src="...cdn.plot.ly/plotly-X.Y.Z.min.js" ...></script>
CDN_RE = re.compile(
    r'<script\b[^>]*\bsrc="https?://cdn\.plot\.ly/(plotly-[\d.]+\.min\.js)"[^>]*>\s*</script>',
    re.IGNORECASE,
)


def inline_plotly(visor_path):
    with open(visor_path, 'r', encoding='utf-8') as f:
        content = f.read()

    if 'PLOTLY-INLINED' in content:
        print('  [SKIP] Plotly ya está incrustado en este visor. Nada que hacer.')
        return

    matches = CDN_RE.findall(content)
    if not matches:
        print('  [WARN] No se encontraron referencias a Plotly en el CDN. Nada que reemplazar.')
        return

    lib_name = matches[0]  # p.ej. plotly-3.6.0.min.js
    lib_path = os.path.join(VENDOR_DIR, lib_name)
    if not os.path.exists(lib_path):
        print(f'  [ERROR] No se encontró la librería local: {lib_path}')
        print(f'          Descárgala con: curl -sL -o {lib_path} '
              f'https://cdn.plot.ly/{lib_name}')
        sys.exit(1)

    with open(lib_path, 'r', encoding='utf-8') as f:
        plotly_src = f.read()

    # Bloque incrustado (una sola vez). El marcador PLOTLY-INLINED sirve de guarda.
    inline_block = (
        '<script>/* PLOTLY-INLINED: ' + lib_name + ' */\n'
        + plotly_src +
        '\n</script>'
    )

    # Quita TODAS las referencias al CDN (había varias por la incrustación de gráficos).
    n_removed = len(CDN_RE.findall(content))
    content = CDN_RE.sub('', content)

    # Inserta la librería una sola vez, justo después de abrir <head>, para que
    # window.Plotly esté disponible antes de cualquier llamada a Plotly.newPlot.
    def insert_after_head(m):
        return m.group(0) + '\n' + inline_block + '\n'

    new_content, n = re.subn(r'<head[^>]*>', insert_after_head, content, count=1,
                             flags=re.IGNORECASE)
    if n == 0:
        # No hay <head>: insertar al inicio del documento.
        new_content = inline_block + '\n' + content

    with open(visor_path, 'w', encoding='utf-8') as f:
        f.write(new_content)

    print(f'  [OK] {n_removed} referencia(s) al CDN eliminadas; '
          f'{lib_name} incrustado en línea.')
    print(f'\nVisor 100% autónomo escrito en: {visor_path}')
    print('Ahora se puede abrir con doble-clic (file://) sin servidor ni internet.')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Uso: python3 inline_plotly.py <ruta/visor_sintenia.html>')
        sys.exit(1)
    inline_plotly(sys.argv[1])
