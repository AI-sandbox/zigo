#!/usr/bin/env python3
import sys

def compute_balanced_accuracy(sexcheck_file):
    TP_1 = 0
    FN_1 = 0
    TP_2 = 0
    FN_2 = 0

    with open(sexcheck_file, 'r') as f:
        lines = f.readlines()

    # Se asume que la primera línea es la cabecera.
    for line in lines[1:]:
        cols = line.strip().split()
        if len(cols) < 6:
            continue
        # Asumiendo:
        # cols[0] = FID, cols[1] = IID, cols[2] = PEDSEX, cols[3] = SNPSEX, cols[4] = STATUS, cols[5] = F
        pedsex = cols[2]
        snpsex = cols[3]

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

    # Calcular TPR para cada clase
    tpr_1 = TP_1 / (TP_1 + FN_1) if (TP_1 + FN_1) > 0 else 0
    tpr_2 = TP_2 / (TP_2 + FN_2) if (TP_2 + FN_2) > 0 else 0

    # Balanced Accuracy = (TPR_1 + TPR_2) / 2
    balanced_accuracy = (tpr_1 + tpr_2) / 2
    return balanced_accuracy

def main():
    if len(sys.argv) < 2:
        print("Uso: python compute_balanced_accuracy.py <archivo.sexcheck>")
        sys.exit(1)

    sexcheck_file = sys.argv[1]
    ba = compute_balanced_accuracy(sexcheck_file)
    print(f"Balanced Accuracy: {ba:.3f}")

if __name__ == "__main__":
    main()

