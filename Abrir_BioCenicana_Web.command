#!/bin/bash
cd "$(dirname "$0")"
echo "================================================="
echo "   Iniciando BioCenicana Web Control Center      "
echo "================================================="
echo "1. Arrancando servidor en el puerto 9090..."
# Iniciar el servidor en segundo plano
java -jar target/biocenicana-1.0.jar web &
SERVER_PID=$!

# Esperar un par de segundos a que el servidor suba
sleep 2

echo "2. Abriendo el navegador en http://localhost:9090"
open http://localhost:9090

echo "-------------------------------------------------"
echo "Presiona Ctrl+C para detener el servidor."
echo "-------------------------------------------------"

# Mantener el script vivo para que el servidor no muera
wait $SERVER_PID
