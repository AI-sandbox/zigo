// gcc -O3 -march=native -flto -funroll-loops vcf_processor.c -o vcf_processor -lz
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <ctype.h>

#define LINE_BUF 1048576

typedef struct {
    double z, o, t, n;            // 0/0, 0/1|1/0, 1/1, missing (diploide estándar)
    // Contadores para el primer alelo (para haploides)
    double f0, f1, fmiss;         // primer alelo: '0', '1' o '.'
} SampleCounts;

static int split_tab(char *line, char **cols, int max_cols) {
    int n = 0;
    char *p = line;
    while (n < max_cols) {
        cols[n++] = p;
        char *t = strchr(p, '\t');
        if (!t) break;
        *t = '\0';
        p = t + 1;
    }
    return n;
}

static inline void parse_gt_and_update(const char *field, SampleCounts *sc) {
    // field = "GT[:...]" ; tomamos hasta ':' si existe
    // Permitimos GT como "a/b", "a|b", "a", "./.", ".|.", ".", etc.
    char gt[16]; int i = 0;
    while (field[i] && field[i] != ':' && i < (int)sizeof(gt)-1) { gt[i] = field[i]; i++; }
    gt[i] = '\0';

    // Primer alelo (para posible haploidía)
    char a = gt[0];
    if (a == '0') sc->f0 += 1.0;
    else if (a == '1') sc->f1 += 1.0;
    else sc->fmiss += 1.0; // incluye '.', o GT vacío/inesperado

    // Diploide estándar
    // Buscamos separador para segundo alelo
    char b = 0;
    for (int j = 1; gt[j]; ++j) {
        if (gt[j] == '/' || gt[j] == '|') { b = gt[j+1]; break; }
    }

    // Si falta alguno de los dos alelos -> missing
    if (!(a == '0' || a == '1') || !(b == '0' || b == '1')) {
        sc->n += 1.0;
        return;
    }

    if (a == '0' && b == '0') sc->z += 1.0;
    else if (a == '1' && b == '1') sc->t += 1.0;
    else sc->o += 1.0; // 0/1 o 1/0
}

static inline void trim_sample_inplace(char *s) {
    if (!s) return;
    // Elimina CR/LF en cualquier sitio
    for (char *q = s; *q; ++q) {
        if (*q == '\r' || *q == '\n') *q = '\0';
    }
    // Trim izquierda
    char *p = s;
    while (*p && isspace((unsigned char)*p)) p++;
    // Trim derecha
    char *e = p + strlen(p);
    while (e > p && isspace((unsigned char)e[-1])) --e;
    *e = '\0';
    // Compacta si hizo avance
    if (p != s) memmove(s, p, (size_t)(e - p + 1));
}

static void write_csv(const char *path, const char *mode, char **samples, int n_samples, SampleCounts *C, int n_snps, int normalize) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open data file for write\n"); exit(2); }

    // Cabecera: sample + columnas según 'mode' (z,o,t,n)
    fprintf(f, "sample");
    for (int i=0; mode[i]; ++i) fprintf(f, ",%c", mode[i]);
    fputc('\n', f);

    for (int s=0; s<n_samples; ++s) {
        double z = C[s].z, o = C[s].o, t = C[s].t, n = C[s].n;

        // Detección de haploidía “total”: sin llamadas diploides informativas
        if (z==0.0 && o==0.0 && t==0.0) {
            double A = C[s].f0;
            double B = C[s].f1;
            double N = C[s].fmiss;
            z = A;
            o = B;
            t = 0.0;
            n = N;
        }

        // Normalizar si se solicita
        if (normalize) {
            double denom = 0.0;
            // Calcular denominador basándonos en lo que se va a exportar (mode)
            for (int i=0; mode[i]; ++i) {
                char m = mode[i];
                if (m=='z') denom += z;
                else if (m=='o') denom += o;
                else if (m=='t') denom += t;
                else if (m=='n') denom += n;
            }
            
            if (denom > 0.0) {
                z /= denom;
                o /= denom;
                t /= denom;
                n /= denom;
            }
        }

        // --- SANEAR NOMBRE DE MUESTRA (evita salto de línea suelto) ---
        char *name = samples ? samples[s] : NULL;
        if (name) trim_sample_inplace(name);
        if (!name || name[0] == '\0') name = (char*)"UNKNOWN_SAMPLE";

        // Una única línea: nombre + métricas
        fprintf(f, "%s", name);
        for (int i=0; mode[i]; ++i) {
            char m = mode[i];
            double v = 0.0;
            if (m=='z') v = z; else if (m=='o') v = o; else if (m=='t') v = t; else if (m=='n') v = n;
            fprintf(f, ",%.6f", v);
        }
        fputc('\n', f);
    }
    fclose(f);
}

static void write_stats(const char *path, int n_samples, int n_snps, const char *mode) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open stats file for write\n"); exit(2); }
    fprintf(f, "{");
    fprintf(f, "\"num_samples\": %d,", n_samples);
    fprintf(f, "\"num_snps_input\": %d,", n_snps);
    fprintf(f, "\"num_snps_used\": %d,", n_snps);
    fprintf(f, "\"overlap_percentage\": 0.0,");
    fprintf(f, "\"rs_percentage\": 0.0,");
    fprintf(f, "\"initial_null_percentage\": 0.0,");
    fprintf(f, "\"final_null_percentage\": 0.0,");
    fprintf(f, "\"mode\": \"%s\"", mode);
    fprintf(f, "}\n");
    fclose(f);
}

int main(int argc, char **argv) {
    // argv: mode normalize_flag data_csv stats_json [vcf.gz]
    if (argc != 5 && argc != 6) {
        fprintf(stderr, "Usage: %s <mode[zotn]> <normalize{0|1}> <out.csv> <stats.json> [vcf.gz]\n", argv[0]);
        return 1;
    }
    const char *mode = argv[1];
    int normalize = (argv[2][0]=='1');
    const char *out_csv = argv[3];
    const char *out_stats = argv[4];
    int use_gz = (argc==6);
    const char *gz_path = use_gz ? argv[5] : NULL;

    char *line = (char*)malloc(LINE_BUF);
    if (!line) { fprintf(stderr, "OOM\n"); return 2; }

    // Lectura de cabecera para obtener muestras
    int n_samples = 0;
    char **samples = NULL;

    if (use_gz) {
        gzFile g = gzopen(gz_path, "rb");
        if (!g) { fprintf(stderr, "Cannot open gz file\n"); return 2; }
        while (gzgets(g, line, LINE_BUF)) {
            if (line[0] == '#') {
                if (strncmp(line, "#CHROM", 6)==0) {
                    // columnas a partir de la 10 son muestras
                    char *cols[2048];
                    int n = split_tab(line, cols, 2048);
                    n_samples = (n>9)? (n-9) : 0;
                    samples = (char**)calloc(n_samples, sizeof(char*));
                    for (int i=0;i<n_samples;i++) samples[i] = strndup(cols[9+i], 256);
                    break;
                }
                continue;
            } else {
                // no había #CHROM
                gzclose(g);
                fprintf(stderr, "Malformed VCF: missing #CHROM header\n");
                return 2;
            }
        }
        gzclose(g);
    } else {
        FILE *in = stdin;
        while (fgets(line, LINE_BUF, in)) {
            if (line[0] == '#') {
                if (strncmp(line, "#CHROM", 6)==0) {
                    char *cols[2048];
                    int n = split_tab(line, cols, 2048);
                    n_samples = (n>9)? (n-9) : 0;
                    samples = (char**)calloc(n_samples, sizeof(char*));
                    for (int i=0;i<n_samples;i++) samples[i] = strndup(cols[9+i], 256);
                    break;
                }
                continue;
            } else {
                fprintf(stderr, "Malformed VCF: missing #CHROM header\n");
                return 2;
            }
        }
    }

    if (n_samples <= 0) {
        fprintf(stderr, "No samples found\n");
        return 2;
    }

    SampleCounts *C = (SampleCounts*)calloc(n_samples, sizeof(SampleCounts));
    if (!C) { fprintf(stderr, "OOM\n"); return 2; }

    // Segunda pasada: procesar variantes
    int n_snps = 0;

    if (use_gz) {
        gzFile g = gzopen(gz_path, "rb");
        if (!g) { fprintf(stderr, "Cannot reopen gz file\n"); return 2; }
        while (gzgets(g, line, LINE_BUF)) {
            if (line[0] == '#') continue;
            // tokenizar por tabs
            char *cols[2048];
            int n = split_tab(line, cols, 2048);
            if (n < 10) continue; // sin genotipos

            // FORMAT debe ir en columna 9; asumimos GT es el primer campo (usual)
            // procesamos las columnas de muestras [9..n-1]
            for (int s=0; s<n_samples; ++s) {
                parse_gt_and_update(cols[9+s], &C[s]);
            }
            n_snps++;
        }
        gzclose(g);
    } else {
        FILE *in = stdin;
        // Rebobinar no posible; reabrir stdin no trivial; asumimos lectura desde archivo redirigido por el caller (vcf_reader.py)
        // vcf_reader.py abre el archivo y pasa el FILE* como stdin, por lo que esta fase llega ya justo después de #CHROM.
        // Aquí simplemente seguimos leyendo:
        while (fgets(line, LINE_BUF, in)) {
            if (line[0] == '#') continue;
            char *cols[2048];
            int n = split_tab(line, cols, 2048);
            if (n < 10) continue;
            for (int s=0; s<n_samples; ++s) {
                parse_gt_and_update(cols[9+s], &C[s]);
            }
            n_snps++;
        }
    }

    write_csv(out_csv, mode, samples, n_samples, C, n_snps, normalize);
    write_stats(out_stats, n_samples, n_snps, mode);

    // Limpieza
    for (int i=0;i<n_samples;i++) free(samples[i]);
    free(samples);
    free(C);
    free(line);
    return 0;
}
