#!/usr/bin/env python3
import os
import sys
import csv

def compute_balanced_accuracy(results_file):
    TP_1, FN_1 = 0, 0
    TP_2, FN_2 = 0, 0
    
    with open(results_file, 'r') as f:
        lines = f.readlines()

    # Ignoramos la primera línea si es cabecera (comprobar según formato exacto)
    # Asumiendo que 'results.sexcheck' tiene encabezado en la primera línea
    for line in lines[1:]:
        cols = line.strip().split('\t')
        if len(cols) < 6:
            continue

        # Extraemos PEDSEX (4ª columna) y SNPSEX (5ª columna)
        pedsex = cols[3]
        snpsex = cols[4]

        if pedsex == '1':
            if snpsex == '1':
                TP_1 += 1
            else:
                FN_1 += 1
        elif pedsex == '2':
            if snpsex == '2':
                TP_2 += 1
            else:
                FN_2 += 1

    # Calculamos TPR para cada clase
    tpr_1 = TP_1 / (TP_1 + FN_1) if (TP_1 + FN_1) > 0 else 0
    tpr_2 = TP_2 / (TP_2 + FN_2) if (TP_2 + FN_2) > 0 else 0

    # Balanced Accuracy = (TPR_1 + TPR_2) / 2
    ba = (tpr_1 + tpr_2) / 2
    return ba

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze-results.py <output_directory>")
        sys.exit(1)

    output_dir = sys.argv[1]
    results = []

    # Obtenemos la lista de items en el directorio
    # Solo nos interesan los subdirectorios (ej. "1", "2", "3", etc.)
    for item in sorted([dir for dir in os.listdir(output_dir) if dir.isdigit()], key=lambda x: int(x)):
        print(item)
        path_item = os.path.join(output_dir, item)
        if os.path.isdir(path_item):
            results_file = os.path.join(path_item, "results.sexcheck")
            if os.path.isfile(results_file):
                ba = compute_balanced_accuracy(results_file)
                results.append((item, ba))

    csv_path = os.path.join(output_dir, "results.csv")
    with open(csv_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["sample", "balanced_acc"])
        for (test_name, ba_value) in results:
            writer.writerow([test_name, f"{ba_value:.3f}"])

    print(f"Analysis finished. Results CSV: {csv_path}")

if __name__ == "__main__":
    main()
