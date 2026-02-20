#!/usr/bin/env bash
##############################################################################
# Lancia Basket e poi Umbrella in sequenza
##############################################################################

set -e

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== BASKET: inizio $(date) ==="
Rscript "$PROJECT_DIR/R/Basket/run_basket_parallel.R" 2>&1 | tee "$PROJECT_DIR/basket_log.txt"
echo "=== BASKET: fine $(date) ==="

echo ""

echo "=== UMBRELLA: inizio $(date) ==="
Rscript "$PROJECT_DIR/R/Umbrella/run_umbrella_parallel.R" 2>&1 | tee "$PROJECT_DIR/umbrella_log.txt"
echo "=== UMBRELLA: fine $(date) ==="

echo ""
echo "=== TUTTE LE SIMULAZIONI COMPLETATE $(date) ==="
